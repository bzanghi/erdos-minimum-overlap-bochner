"""
ub_local — local minimisation of  M(h) = (2/n) max_j sum_i h_i (1-h)_{i+j}
over  {0 <= h <= 1, sum h = n/2}.

Two stages, both validated against `ub_core` (which owns the lag convention):

  `smooth_descent`  softmax-smoothed minimax + Adam + exact feasible
                    projection.  Cheap, batched over many starts, used for
                    landscape exploration.
  `slp_polish`      trust-region sequential LP (`highs-ipm`, per PRO-34:
                    simplex hangs on these ~900-way degenerate minimax LPs).
                    Used to certify a candidate is first-order stationary.

Every `M` reported by this module is the exact objective from
`ub_core.overlap_value` — the smoothing is only ever used to produce search
directions, never to report a bound.
"""
from __future__ import annotations

import numpy as np

from ub_core import (
    overlap_profile,
    overlap_value,
    grad_lag,
    project_feasible,
    lag_of_index,
    index_of_lag,
)

__all__ = ["smooth_descent", "slp_polish", "batch_overlap_values", "random_start"]


# ---------------------------------------------------------------- batched eval

def batch_overlap_values(H: np.ndarray) -> np.ndarray:
    """M(h) for each row of H (B, n), via FFT.  Matches ub_core.overlap_value."""
    H = np.atleast_2d(np.asarray(H, dtype=np.float64))
    B, n = H.shape
    L = 1
    while L < 2 * n:
        L *= 2
    G = 1.0 - H
    # correlate(h, g) == ifft( fft(h) * conj(fft(g)) ), lags recovered cyclically
    fh = np.fft.rfft(H, L, axis=1)
    fg = np.fft.rfft(G, L, axis=1)
    c = np.fft.irfft(fh * np.conjugate(fg), L, axis=1)
    # lags 0..n-1 live at c[:, 0:n]; lags -(n-1)..-1 live at c[:, L-(n-1):L]
    m = np.maximum(c[:, :n].max(axis=1), c[:, L - (n - 1):].max(axis=1))
    return (2.0 / n) * m


def _softmax_grad(h: np.ndarray, beta: float, topk: int = 256):
    """(max value, gradient) of the softmax-smoothed max over the top-`topk` lags.

    Vectorised: the weighted sum  sum_j w_j grad M_j  collapses to two length-2n
    convolutions rather than `topk` explicit gradients.  `_softmax_grad_ref` is
    the literal loop; `selftest` asserts they agree.
    """
    n = h.size
    prof = overlap_profile(h) * (2.0 / n)
    if topk < prof.size:
        idx = np.argpartition(prof, -topk)[-topk:]
    else:
        idx = np.arange(prof.size)
    vals = prof[idx]
    vmax = vals.max()
    w = np.exp(beta * (vals - vmax))
    w /= w.sum()

    # Wv is indexed by t = j + (n-1), i.e. Wv[index_of_lag(j)] is *not* w_j —
    # profile index m already equals index_of_lag(lag), so scatter on m directly
    # after re-mirroring into the t-axis: t = (2n-2) - m.
    Wv = np.zeros(2 * n - 1)
    Wv[(2 * n - 2) - idx] = w

    g_arr = 1.0 - h
    term1 = np.convolve(Wv[::-1], g_arr)[n - 1: 2 * n - 1]
    term2 = np.convolve(Wv, h)[n - 1: 2 * n - 1]
    return vmax, (2.0 / n) * (term1 - term2)


def _softmax_grad_ref(h: np.ndarray, beta: float, topk: int = 256):
    """Literal reference implementation of `_softmax_grad` (slow)."""
    n = h.size
    prof = overlap_profile(h) * (2.0 / n)
    if topk < prof.size:
        idx = np.argpartition(prof, -topk)[-topk:]
    else:
        idx = np.arange(prof.size)
    vals = prof[idx]
    vmax = vals.max()
    w = np.exp(beta * (vals - vmax))
    w /= w.sum()
    g = np.zeros(n)
    for wi, mi in zip(w, idx):
        g += wi * grad_lag(h, lag_of_index(int(mi), n))
    return vmax, g


def smooth_descent(
    h0: np.ndarray,
    iters: int = 400,
    lr: float = 3e-3,
    beta0: float = 300.0,
    beta1: float = 20000.0,
    topk: int = 256,
    seed: int = 0,
) -> tuple[np.ndarray, float]:
    """Adam on the smoothed minimax, with beta annealed up.  Returns (h, exact M)."""
    h = project_feasible(np.asarray(h0, dtype=np.float64))
    m = np.zeros_like(h)
    v = np.zeros_like(h)
    best_h, best_M = h.copy(), overlap_value(h)
    for t in range(1, iters + 1):
        beta = beta0 * (beta1 / beta0) ** (t / iters)
        _, g = _softmax_grad(h, beta, topk)
        m = 0.9 * m + 0.1 * g
        v = 0.999 * v + 0.001 * (g * g)
        step = lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-12)
        h = project_feasible(h - step)
        cur = overlap_value(h)
        if cur < best_M:
            best_M, best_h = cur, h.copy()
    return best_h, best_M


# ------------------------------------------------------------------ SLP polish

def slp_polish(
    h0: np.ndarray,
    rounds: int = 40,
    r0: float = 1e-3,
    r_min: float = 1e-9,
    slack: float = 1e-3,
    verbose: bool = False,
):
    """Trust-region sequential LP.  Returns (h, M, certified_first_order_gain).

    The final `gain` is the LP-certified maximum first-order improvement at
    radius `r_min`: a value at machine-noise level certifies first-order
    stationarity (cf. PRO-33's +1.94e-10 for Together's h*).
    """
    from scipy.optimize import linprog

    h = project_feasible(np.asarray(h0, dtype=np.float64))
    n = h.size
    M = overlap_value(h)
    r = r0
    gain = np.nan

    for _ in range(rounds):
        if r < r_min:
            break
        prof = overlap_profile(h) * (2.0 / n)
        Mmax = prof.max()
        # A step of radius r moves any scaled M_j by at most 4r (|grad|_1 <= 4).
        keep = np.nonzero(prof >= Mmax - (4.0 * r + slack))[0]
        A = np.empty((keep.size, n))
        for row, mi in enumerate(keep):
            A[row] = grad_lag(h, lag_of_index(int(mi), n))
        # vars: [delta (n), u].  Model the linearised max as Mmax + u:
        #   M_j + grad_j . delta <= Mmax + u   <=>   grad_j . delta - u <= Mmax - M_j
        A_ub = np.hstack([A, -np.ones((keep.size, 1))])
        b_ub = Mmax - prof[keep]
        A_eq = np.hstack([np.ones((1, n)), np.zeros((1, 1))])
        b_eq = np.zeros(1)
        lo = np.maximum(-r, -h)
        hi = np.minimum(r, 1.0 - h)
        bounds = list(zip(lo, hi)) + [(-1.0, 1.0)]
        c = np.zeros(n + 1)
        c[-1] = 1.0

        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=bounds, method="highs-ipm")
        if not res.success:
            r *= 0.5
            continue

        delta, u = res.x[:n], res.x[-1]
        gain = -u  # predicted decrease in max_j M_j
        # validity of the lag restriction: dropped lags cannot overtake
        dropped = np.delete(prof, keep)
        if dropped.size and dropped.max() + 4.0 * r >= Mmax + u:
            r *= 0.5
            continue

        h_new = project_feasible(h + delta)
        M_new = overlap_value(h_new)
        if M_new < M - 1e-15:
            h, M = h_new, M_new
            r = min(r0, r * 1.6)
        else:
            r *= 0.4
        if verbose:
            print(f"  r={r:.2e}  M={M:.16f}  gain={gain:.3e}  lags={keep.size}")

    return h, M, gain


# --------------------------------------------------------------------- starts

def random_start(n: int, kind: str, rng: np.random.Generator) -> np.ndarray:
    """Structured / random initialisations spanning qualitatively distinct basins."""
    if kind == "uniform":
        h = rng.random(n)
    elif kind == "bernoulli":
        h = (rng.random(n) < 0.5).astype(float)
    elif kind == "blocks":
        k = int(rng.integers(2, 40))
        blk = (rng.random(k) < 0.5).astype(float)
        h = np.repeat(blk, int(np.ceil(n / k)))[:n].astype(float)
        h = h + rng.normal(0, 0.05, n)
    elif kind == "fractional":
        # match h*'s observed 28% lower / 62% interior / 10% upper cell mix
        u = rng.random(n)
        h = np.where(u < 0.28, 0.0, np.where(u < 0.90, rng.random(n), 1.0))
    elif kind == "smooth":
        k = int(rng.integers(1, 12))
        x = np.linspace(0, 2, n, endpoint=False)
        h = 0.5 + 0.4 * np.sin(np.pi * k * x / 2 + rng.random() * 2 * np.pi)
    elif kind == "symmetric":
        half = rng.random((n + 1) // 2)
        h = np.concatenate([half, half[::-1]])[:n]
    elif kind == "antisym":
        half = rng.random((n + 1) // 2)
        h = np.concatenate([half, 1.0 - half[::-1]])[:n]
    elif kind == "sparse":
        h = np.zeros(n)
        h[rng.permutation(n)[: n // 2]] = 1.0
    else:
        raise ValueError(kind)
    return project_feasible(h)


START_KINDS = ("uniform", "bernoulli", "blocks", "fractional", "smooth",
               "symmetric", "antisym", "sparse")
