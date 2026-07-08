"""PRO-26 Phase 2a v2: Piecewise ansatz on Together's h*.

v1 (Chebyshev on full h) failed because h has bang-bang structure that smooth
basis can't represent without Gibbs oscillations.

v2 fixes the bang-bang structure by:
- HOLDING h = 0 on A_0 (168 cells where h_init ≤ 0.01)
- HOLDING h = 1 on A_plus (62 cells where h_init ≥ 0.99)
- Optimizing only the 370 interior cells via a (low-dim) ansatz

For the interior: parameterize via either
(a) the 370 cell values directly (full DOF — already tried in PRO-4 LP descent, stalled)
(b) a low-dim Chebyshev expansion on A_int's index range (smooth interior)
(c) a sinusoidal expansion (Fourier on interior)

This v2 tries (b) — Chebyshev on the interior's index range.

The key new constraint: total mass over the interior cells must equal
target_int_sum = (target_sum) - |A_plus| (the upper-active cells contribute 1
each; lower-active 0).
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def load_together_h() -> np.ndarray:
    with open(DATA) as f:
        return np.array(json.load(f)["together"]["values"], dtype=np.float64)


def split_regions(h: np.ndarray, eps_low: float = 0.01,
                  eps_high: float = 0.99) -> dict:
    n = len(h)
    idx = np.arange(n)
    A_0 = idx[h <= eps_low]
    A_int = idx[(h > eps_low) & (h < eps_high)]
    A_plus = idx[h >= eps_high]
    return {
        "A_0": A_0, "A_int": A_int, "A_plus": A_plus,
        "n": n, "n_0": len(A_0), "n_int": len(A_int), "n_plus": len(A_plus),
    }


def build_h_from_interior(a_int: np.ndarray, regions: dict, n: int) -> np.ndarray:
    """Reconstruct full-h: h=0 on A_0, h=1 on A_plus, h=a_int on A_int."""
    h = np.zeros(n)
    h[regions["A_plus"]] = 1.0
    h[regions["A_int"]] = a_int
    return h


def cheb_basis_on_indices(indices: np.ndarray, K: int) -> np.ndarray:
    """Chebyshev T_0..T_K evaluated at normalized interior indices.
    Returns shape (K+1, len(indices))."""
    i_min, i_max = indices[0], indices[-1]
    x = 2 * (indices - i_min) / (i_max - i_min) - 1   # in [-1, 1]
    B = np.zeros((K + 1, len(indices)))
    for k in range(K + 1):
        c = np.zeros(K + 1)
        c[k] = 1.0
        B[k] = cheb.chebval(x, c)
    return B


def project_interior(h_init: np.ndarray, regions: dict, K: int) -> np.ndarray:
    """Project h_init[A_int] onto Chebyshev T_0..T_K."""
    A_int = regions["A_int"]
    interior_vals = h_init[A_int]
    B = cheb_basis_on_indices(A_int, K)
    a, *_ = np.linalg.lstsq(B.T, interior_vals, rcond=None)
    return a


def interior_from_coefs(coefs: np.ndarray, A_int: np.ndarray, K: int,
                        target_sum: float) -> np.ndarray:
    """Reconstruct interior cell values from Chebyshev coefs.
    Adjusts constant term to enforce Σ(interior) = target_sum."""
    B = cheb_basis_on_indices(A_int, K)
    a_int = coefs @ B  # shape (n_int,)
    # Shift to satisfy sum constraint
    shift = (target_sum - a_int.sum()) / len(a_int)
    return a_int + shift


def compute_M_fast(h: np.ndarray, L: float) -> np.ndarray:
    """Return M at ALL 2n-1 lags (positive + negative). For asymmetric h,
    the sup_t M(h, t) may be at a negative lag — Together's convention uses
    the full correlate output (see lp_research_state/code/together_loader.py:457).
    BUG FIXED 2026-05-18: previously this only returned positive lags, which
    made BFGS find h with bad behavior at negative lags.
    """
    n = len(h)
    corr = np.correlate(h, 1.0 - h, mode="full")
    return L * corr  # shape (2n-1,), index n-1+j = lag j, index n-1-j = lag -j


def soft_max(values: np.ndarray, tau: float) -> float:
    M_max = values.max()
    return M_max + np.log(np.sum(np.exp(tau * (values - M_max)))) / tau


def objective(coefs: np.ndarray, n: int, K: int, L: float, regions: dict,
              target_int_sum: float, tau: float, box_penalty: float) -> float:
    a_int = interior_from_coefs(coefs, regions["A_int"], K, target_int_sum)
    h = build_h_from_interior(a_int, regions, n)
    # Box penalty on interior only (since A_0/A_plus are fixed in-range)
    box_viol = (np.maximum(0, -a_int).sum() + np.maximum(0, a_int - 1.0).sum())
    Ms = compute_M_fast(h, L)
    M_smax = soft_max(Ms, tau)
    return M_smax + box_penalty * box_viol


def evaluate(coefs: np.ndarray, n: int, K: int, L: float, regions: dict,
             target_int_sum: float) -> tuple[float, np.ndarray]:
    a_int = interior_from_coefs(coefs, regions["A_int"], K, target_int_sum)
    a_int = np.clip(a_int, 0.0, 1.0)
    # Re-adjust after clipping
    if abs(a_int.sum() - target_int_sum) > 1e-9 and a_int.sum() > 0:
        scale = target_int_sum / a_int.sum()
        a_int = np.clip(a_int * scale, 0.0, 1.0)
    h = build_h_from_interior(a_int, regions, n)
    Ms = compute_M_fast(h, L)
    return float(Ms.max()), h


def main():
    print("=" * 78)
    print("PRO-26 Phase 2a v2: Piecewise interior-only Chebyshev opt")
    print("=" * 78)
    h_init = load_together_h()
    n = len(h_init)
    L = 2.0 / n
    target_sum = float(h_init.sum())
    M_init = compute_M_fast(h_init, L).max()
    print(f"n = {n}, L = {L}, total ∫h = {L*target_sum} = 1")
    print(f"M(h_init) = {M_init:.10f}\n")

    regions = split_regions(h_init)
    print(f"Regions: A_0 = {regions['n_0']} cells, "
          f"A_int = {regions['n_int']} cells, "
          f"A_plus = {regions['n_plus']} cells")
    target_int_sum = target_sum - regions["n_plus"]  # interior must sum to this
    print(f"Target interior sum = {target_int_sum:.4f} "
          f"(= {target_sum} - {regions['n_plus']})\n")

    for K in [5, 10, 20, 50, 100]:
        print(f"--- K = {K} ({K+1} Chebyshev coefs on A_int) ---")
        a0 = project_interior(h_init, regions, K)
        M_proj, h_proj = evaluate(a0, n, K, L, regions, target_int_sum)
        # Reconstruction error on the interior
        a_int_proj = interior_from_coefs(a0, regions["A_int"], K, target_int_sum)
        proj_err = np.abs(a_int_proj - h_init[regions["A_int"]]).max()
        print(f"  Interior proj error = {proj_err:.3e}")
        print(f"  M(proj) = {M_proj:.10f}  Δ vs h_init = {M_proj - M_init:+.3e}")
        # Optimize
        res = minimize(objective, a0,
                       args=(n, K, L, regions, target_int_sum, 1e6, 100.0),
                       method="BFGS",
                       options={"maxiter": 200, "gtol": 1e-12})
        M_opt, h_opt = evaluate(res.x, n, K, L, regions, target_int_sum)
        improve = M_opt < M_init - 1e-10
        marker = "  ✓ IMPROVED" if improve else "  (no improvement)"
        print(f"  After BFGS ({res.nit} iters): M = {M_opt:.10f}  "
              f"Δ vs init = {M_opt - M_init:+.3e}{marker}\n")


if __name__ == "__main__":
    main()
