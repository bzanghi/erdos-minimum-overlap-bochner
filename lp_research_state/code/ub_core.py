"""
ub_core — clean, self-contained upper-bound machinery for the Erdős minimum
overlap constant µ.

Formulation (Haugland 2016 / Together, verbatim as in `together_loader.py`):

    h : [0, 2] -> [0, 1] piecewise constant on n equal cells of width 2/n
    mass       :  (2/n) * sum(h) = 1        <=>  sum(h) = n/2
    objective  :  M(h) = sup_k  int h(x) (1 - h(x+k)) dx

with BOTH h and g := 1 - h supported on [0, 2] (zero-extended outside).  The
correlation of two step functions on a common grid is piecewise linear in k
with breakpoints on the grid, so the sup over real k is attained at an integer
cell shift and

    M(h) = (2/n) * max_j  M_j,      M_j := sum_i h_i * (1 - h)_{i+j}
                                        = np.correlate(h, 1-h, mode='full')

`np.correlate(a, v, 'full')` returns lags j = -(n-1) .. (n-1) — BOTH signs.
This matters: see the repo memory note "positive-lag correlation trap".  Every
sup here is taken over the full array, never a positive-lag slice.

For any feasible h,  µ <= M(h).  Refining the grid n -> 2n by cell doubling
preserves M exactly, so M_n is non-increasing in n along doublings and
µ = inf_n M_n.

Gradient.  M_j is exactly quadratic in h:
    M_j = sum_i h_i (1 - h_{i+j})  =  sum(h) - sum_i h_i h_{i+j}
    dM_j/dh_k = (1 - h_{k+j})·1{k+j in range} - h_{k-j}·1{k-j in range}
(the indicator terms are the ones PRO-23 dropped; see PRO33_KKT_CORRECTION.md).
Scaled objective M(h) = (2/n)·max_j M_j carries the same 2/n factor.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "overlap_profile",
    "overlap_value",
    "grad_lag",
    "project_feasible",
    "cell_double",
    "load_together",
    "ANCHOR",
]

# Canonical Together h* objective value, float64, sup over BOTH lag signs
# (argmax at s = -33 under the np.correlate(h, 1-h) convention).
ANCHOR = 0.3808703105862199


def overlap_profile(h: np.ndarray) -> np.ndarray:
    """All 2n-1 lag values, as `np.correlate(h, 1-h, 'full')`.

    CONVENTION (verified by `selftest`, not assumed): entry `m` of this array
    is the lag `j = (n-1) - m`, i.e. the axis is MIRRORED relative to the naive
    reading `m - (n-1)`.  Use `lag_of_index` / `index_of_lag` rather than doing
    the arithmetic inline — getting this backwards is the repo's documented
    "positive-lag correlation trap".
    """
    h = np.asarray(h, dtype=np.float64)
    return np.correlate(h, 1.0 - h, mode="full")


def lag_of_index(m: int, n: int) -> int:
    """Signed lag j such that profile[m] == sum_i h_i (1-h)_{i+j}."""
    return (n - 1) - m


def index_of_lag(j: int, n: int) -> int:
    """Inverse of `lag_of_index`."""
    return (n - 1) - j


def overlap_value(h: np.ndarray) -> float:
    """M(h) = (2/n) * max_j M_j.  Sup over BOTH lag signs."""
    h = np.asarray(h, dtype=np.float64)
    n = h.size
    return (2.0 / n) * float(overlap_profile(h).max())


def grad_lag(h: np.ndarray, j: int) -> np.ndarray:
    """Exact gradient of the *scaled* lag value (2/n)*M_j w.r.t. h.

    j is the signed lag in the np.correlate 'full' convention, i.e. index
    `j + (n-1)` into `overlap_profile(h)`.
    """
    h = np.asarray(h, dtype=np.float64)
    n = h.size
    g = np.zeros(n, dtype=np.float64)
    # term 1: (1 - h_{k+j}) for k+j in [0, n)
    lo, hi = max(0, -j), min(n, n - j)
    if hi > lo:
        g[lo:hi] += 1.0 - h[lo + j : hi + j]
    # term 2: -h_{k-j} for k-j in [0, n)
    lo2, hi2 = max(0, j), min(n, n + j)
    if hi2 > lo2:
        g[lo2:hi2] -= h[lo2 - j : hi2 - j]
    return (2.0 / n) * g


def project_feasible(h: np.ndarray, tol: float = 1e-14, iters: int = 200) -> np.ndarray:
    """Project onto {0 <= h <= 1, sum h = n/2} (exact; bisection on the shift)."""
    h = np.asarray(h, dtype=np.float64).copy()
    n = h.size
    target = n / 2.0

    def s(t):
        return np.clip(h + t, 0.0, 1.0).sum()

    lo, hi = -1.0, 1.0
    while s(lo) > target:
        lo *= 2.0
    while s(hi) < target:
        hi *= 2.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if s(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return np.clip(h + 0.5 * (lo + hi), 0.0, 1.0)


def cell_double(h: np.ndarray, factor: int = 2) -> np.ndarray:
    """Value-preserving grid refinement n -> factor*n.  M is exactly preserved."""
    return np.repeat(np.asarray(h, dtype=np.float64), factor)


def resample(h: np.ndarray, m: int) -> np.ndarray:
    """Cell-average re-projection of a step function from n cells onto m cells.

    Exact (no interpolation): each new cell gets the mean of the old function
    over it, computed from the cumulative integral.  Mass is preserved, so the
    result is feasible.  Unlike `cell_double` this does NOT preserve M when
    m is not a multiple of n — that is the point: a 600-cell step function is
    not representable on an 800-cell grid, so re-optimising there explores a
    genuinely different function space.
    """
    h = np.asarray(h, dtype=np.float64)
    n = h.size
    cum = np.concatenate([[0.0], np.cumsum(h) / n])  # integral over [0, k/n]
    edges = np.arange(m + 1) / m
    x = edges * n
    lo = np.floor(x).astype(int)
    lo = np.minimum(lo, n - 1)
    frac = x - lo
    vals = cum[lo] + frac * h[np.minimum(lo, n - 1)] / n
    vals[-1] = cum[-1]
    return np.clip(np.diff(vals) * m, 0.0, 1.0)


def load_together() -> np.ndarray:
    """Together's n=600 h* from the repo's canonical JSON."""
    import json
    from pathlib import Path

    p = Path(__file__).parent.parent / "data" / "together_f_star.json"
    d = json.loads(p.read_text())["together"]
    h = np.asarray(d["values"], dtype=np.float64)
    assert h.size == 600 and abs(h.sum() - 300.0) < 1e-12, "unexpected h*"
    return h


def selftest(verbose: bool = True) -> None:
    """Assert the lag convention, the gradient, and the published anchor.

    Run this before trusting any result from this module.
    """
    rng = np.random.default_rng(1)
    n = 9
    h = rng.random(n)
    g = 1.0 - h
    prof = overlap_profile(h)

    def direct(j):
        return sum(h[i] * g[i + j] for i in range(n) if 0 <= i + j < n)

    for m in range(prof.size):
        j = lag_of_index(m, n)
        assert abs(prof[m] - direct(j)) < 1e-13, f"lag convention wrong at m={m}"

    # gradient vs finite difference, on a lag with real support
    for j in (-3, 0, 2):
        gr = grad_lag(h, j)
        d = rng.standard_normal(n) * 1e-6
        f0 = (2.0 / n) * prof[index_of_lag(j, n)]
        f1 = (2.0 / n) * overlap_profile(h + d)[index_of_lag(j, n)]
        assert abs(gr @ d - (f1 - f0)) < 1e-11, f"gradient wrong at lag {j}"

    # feasibility projection
    hp = project_feasible(rng.standard_normal(50))
    assert abs(hp.sum() - 25.0) < 1e-10 and hp.min() >= -1e-15 and hp.max() <= 1 + 1e-15

    # cell doubling preserves the objective (to float64 rounding: 1 ulp)
    ht = load_together()
    assert abs(overlap_value(cell_double(ht)) - overlap_value(ht)) < 1e-15

    # published anchor, bit-for-bit
    v = overlap_value(ht)
    assert v == ANCHOR, f"anchor drift: {v!r} != {ANCHOR!r}"
    jstar = lag_of_index(int(overlap_profile(ht).argmax()), ht.size)
    assert jstar == -33, f"expected argmax lag -33, got {jstar}"

    if verbose:
        print(f"ub_core selftest OK  |  M(h*) = {v:.16f}  argmax lag = {jstar}")


if __name__ == "__main__":
    selftest()
