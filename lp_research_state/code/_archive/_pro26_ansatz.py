"""PRO-26 Phase 2a — Rechnitzer-style ansatz for μ's UB side.

Goal: parameterize h on [0, 2] with a small ansatz (~32 params), optimize
in mpmath at 50-digit precision, evaluate sup_t M(h; t) at high precision
via FFT-accelerated correlation, BFGS to minimize.

Starting point: Together's h* with regions A_+ (h≈1, 62 cells), A_0 (h≈0,
168 cells), A_int (interior, 370 cells) — already characterized in PRO-23.

PHASE 2a strategy (per PRO26_RECHNITZER_ANALYSIS.md §4):
1. Fix region assignments A_+, A_0, A_int from Together's h*.
2. Parameterize the INTERIOR by a low-dimensional ansatz.
3. mpmath BFGS at 50 digits on the ansatz coefficients.
4. Evaluate M(t) at all grid shifts via correlation; sup is max.

PHASE 2a — this file: build the scaffolding (load + region detection +
high-precision M evaluation). Optimization itself is the next step.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from mpmath import mp, mpf

mp.dps = 60  # 60-digit working precision (60 > 50 target with safety margin)

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def load_together_h() -> np.ndarray:
    """Load Together's h* (n=600 piecewise-constant on [0, 2])."""
    with open(DATA) as f:
        data = json.load(f)
    return np.array(data["together"]["values"], dtype=np.float64)


def detect_regions(h: np.ndarray,
                   eps_low: float = 0.01,
                   eps_high: float = 0.99) -> dict[str, np.ndarray]:
    """Partition cells into A_+ (h ≈ 1), A_0 (h ≈ 0), A_int (interior).

    Matches PRO-23's characterization: at eps_low=0.01, eps_high=0.99
    gives 168 / 370 / 62 cells respectively for Together's h*.
    """
    n = len(h)
    idx = np.arange(n)
    A_0 = idx[h <= eps_low]
    A_int = idx[(h > eps_low) & (h < eps_high)]
    A_plus = idx[h >= eps_high]
    return {"A_0": A_0, "A_int": A_int, "A_plus": A_plus}


def compute_M_mp(h_mp: list[mpf], L: mpf, max_shifts: int) -> list[mpf]:
    """Evaluate M(jL) for j = 0..max_shifts-1 at mpmath precision.

    M(jL) = L * Σ_{i=0..n-1-j} h[i] * (1 - h[i+j])  (Together convention)
    """
    n = len(h_mp)
    Ms = []
    for j in range(max_shifts):
        total = mpf(0)
        for i in range(n - j):
            total += h_mp[i] * (1 - h_mp[i + j])
        Ms.append(L * total)
    return Ms


def compute_M_fast(h: np.ndarray, L: float) -> np.ndarray:
    """Fast float64 version of M at ALL 2n-1 lags. Uses np.correlate.
    BUG FIXED 2026-05-18: previously this only returned positive lags (n
    values), which is wrong for asymmetric h since M(t) ≠ M(-t). Together's
    `compute_overlap_from_f` (together_loader.py) uses full corr.max(), so
    the published UB 0.380871 IS the max over all lags.
    """
    n = len(h)
    corr = np.correlate(h, 1.0 - h, mode="full")
    return L * corr  # shape (2n-1,)


def baseline_status():
    """Print the starting point: Together's h* and its M values."""
    h = load_together_h()
    n = len(h)
    L = 2.0 / n
    print(f"Together h*: n = {n}, L = {L}, sum = {h.sum()}, ∫h = L·sum = {L*h.sum()}")

    regions = detect_regions(h)
    for name, idx in regions.items():
        print(f"  {name}: {len(idx)} cells")
    assert sum(len(idx) for idx in regions.values()) == n

    # Fast M evaluation at float64
    Ms = compute_M_fast(h, L)
    M_max = Ms.max()
    j_max = int(np.argmax(Ms))
    print(f"\nM(jL) max = {M_max:.16f} at j = {j_max} (k = {j_max*L})")
    print(f"Together reported UB = 0.380871 (6 digits)")
    print(f"Float64 precise UB   = {M_max:.16f}")

    # Top-5 active shifts
    sorted_idx = np.argsort(Ms)[::-1]
    print("\nTop-5 active shifts:")
    for r, idx in enumerate(sorted_idx[:5]):
        print(f"  rank {r+1}: j={int(idx):4d}, k={idx*L:.6f}, M={Ms[idx]:.16f}")
    # Active set sizes at multiple tolerances
    print("\nActive-set sizes:")
    for tol in [1e-15, 1e-10, 1e-7, 1e-4]:
        n_active = int((M_max - Ms <= tol).sum())
        print(f"  tol={tol:.0e}: |S| = {n_active}")

    return h, L, regions


def estimate_M_at_mp_precision(h_f64: np.ndarray, L_f64: float,
                                shifts_to_check: list[int] | None = None,
                                target_dps: int = 50) -> dict:
    """Re-evaluate M at chosen shifts in mpmath at target_dps precision.

    This is for VERIFYING the float64 value at high precision (anchoring
    the 16-digit UB at >50 digits for downstream PSLQ / closed-form work).
    """
    old_dps = mp.dps
    mp.dps = target_dps + 5  # working precision
    try:
        h_mp = [mpf(str(x)) for x in h_f64]  # exact decimal from float64
        L_mp = mpf(2) / len(h_mp)

        if shifts_to_check is None:
            # default: float64 argmax + a few around it
            Ms_f64 = compute_M_fast(h_f64, L_f64)
            j_max = int(np.argmax(Ms_f64))
            shifts_to_check = sorted(set([0, j_max - 1, j_max, j_max + 1,
                                           min(len(h_f64) - 1, j_max + 2)]))
            shifts_to_check = [j for j in shifts_to_check if 0 <= j < len(h_f64)]

        result = {}
        for j in shifts_to_check:
            n = len(h_mp)
            total = mpf(0)
            for i in range(n - j):
                total += h_mp[i] * (1 - h_mp[i + j])
            M_j = L_mp * total
            result[j] = mp.nstr(M_j, target_dps)
        return result
    finally:
        mp.dps = old_dps


def main():
    print("=" * 78)
    print("PRO-26 Phase 2a — Scaffolding: regions, baseline M, mpmath evaluation")
    print("=" * 78)
    h, L, regions = baseline_status()

    print("\n--- High-precision M values at top active shifts ---")
    print(f"(target dps = 50, working dps = 55)")
    hp_values = estimate_M_at_mp_precision(h, L, target_dps=50)
    for j, val in hp_values.items():
        print(f"  j={j:4d}: M = {val}")


if __name__ == "__main__":
    main()
