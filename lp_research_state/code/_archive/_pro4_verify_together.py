"""PRO-4: verify Together's reported UB from saved h*, then prep for refinement.

The objective (Haugland/Together): M(k) = ∫_0^2 h(x)(1 - h(x+k)) dx with zero extension.

For piecewise-constant h on n cells of width L = 2/n, M(k) is piecewise-LINEAR in k
(since corr(k) = ∫h(x)h(x+k)dx is piecewise-linear), so max over k ∈ [0, 2] occurs
at one of the n+1 grid points k = j*L (j = 0, ..., n).
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def load_together():
    with open(DATA) as f:
        data = json.load(f)
    tg = data["together"]
    h = np.array(tg["values"])
    n = len(h)
    L = 2.0 / n
    return h, L, n, tg["meta"]


def M_grid(h: np.ndarray, L: float) -> np.ndarray:
    """Return Together-convention M(j*L) for j = 0, 1, ..., n-1.

    Together's M (per `np.correlate(h, 1-h, mode='full')`) only sums over the
    overlap region — i.e., the integral is restricted to where both h(x) and
    h(x+jL) lie in [0, 2]. So
        M(j*L) = L * Σ_{i=0..n-1-j} h[i] * (1 - h[i+j])
               = L * (Σ_{i=0..n-1-j} h[i]) - L * Σ_{i=0..n-1-j} h[i]*h[i+j]
    For j = 0: M(0) = L * Σ h[i]*(1 - h[i]).
    """
    n = len(h)
    Ms = np.zeros(n)
    for j in range(n):
        h_lo = h[: n - j] if j > 0 else h
        h_hi = h[j:] if j > 0 else h
        Ms[j] = L * (h_lo.sum() - (h_lo * h_hi).sum())
    return Ms


def main():
    h, L, n, meta = load_together()
    print(f"n={n}, L={L}, sum_h={h.sum()}, ∫h = L*sum_h = {L*h.sum()}")
    print(f"meta claimed_bound = {meta['claimed_bound']}")

    Ms = M_grid(h, L)
    j_max = int(np.argmax(Ms))
    print(f"\nmax_j M(j*L) = {Ms[j_max]:.10f} at j = {j_max} (k = {j_max*L:.6f})")
    print(f"Reported UB = {meta['claimed_bound']}")
    diff = Ms[j_max] - meta["claimed_bound"]
    print(f"diff = {diff:+.3e}")

    # Print top-K shifts and active set at multiple tolerances
    sorted_idx = np.argsort(Ms)[::-1]
    print("\nTop-10 M(k) values:")
    for r, idx in enumerate(sorted_idx[:10]):
        print(f"  rank {r+1}: j={idx:4d}, k={idx*L:.6f}, M={Ms[idx]:.10f}, slack={Ms[j_max]-Ms[idx]:.3e}")

    for tol in [1e-12, 1e-10, 1e-9, 1e-8, 1e-6, 1e-4]:
        active = np.where(Ms[j_max] - Ms <= tol)[0]
        print(f"  active at tol={tol:.0e}: |S|={len(active)}")

    print("\nh distribution:")
    print(f"  h ≤ 0.01: {(h <= 0.01).sum()}")
    print(f"  0.01 < h < 0.99: {((h > 0.01) & (h < 0.99)).sum()}")
    print(f"  h ≥ 0.99: {(h >= 0.99).sum()}")


if __name__ == "__main__":
    main()
