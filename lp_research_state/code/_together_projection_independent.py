"""
Independent re-implementation of step-function -> Fourier projection.

Method: per-cell high-order Gauss-Legendre quadrature of f(x) cos(omega_k x)
and f(x) sin(omega_k x). This is a DIFFERENT algorithm from
together_diagnostic.project_step_function (which uses closed-form sin/cos
per-cell integrals). The two must agree to >= 10 significant digits --
project cross-verification rule.

Convention (per lp_research_state/code/_fourier_convention_notes.md):
    - f lives on [-2, 2], period 4
    - omega_k = pi * k / 2  for k = 1..T
    - c[k] = integral_{-2}^{2} f(x) cos(omega_k x) dx
    - d[k] = integral_{-2}^{2} f(x) sin(omega_k x) dx
    - c[0] = 0.5 (placeholder for f-hat(0) in White's half-period normalization)
    - d[0] = 0.0

Why per-cell Gauss-Legendre is a clean independent check:
  - On each cell [b_i, b_{i+1}] the integrand is wv_i * cos(omega x) (or sin),
    a single-mode trig function times a constant. A 64-point Gauss-Legendre
    rule integrates polynomials of degree up to 127 exactly. The Taylor
    expansion of cos at our scales (|omega| <= pi*T/2 ~ 6283 for T=4000,
    cell width ~ 4/n_cells ~ 0.007) means each cell has |omega * (b-a)/2|
    bounded by ~22 -- the GL error decays like (e * |omega| * h / (4n))^(2n)
    which is well below machine eps at n=64 for almost all (i, k). The few
    high-k modes where this saturates produce sub-1e-12 noise, comfortably
    below the 1e-10 cross-check threshold.
  - The numerics use only leggauss nodes/weights and numpy cos/sin. The
    closed form uses analytic sin / omega. The only shared infrastructure
    is the convention constants (period, basis, omega_k); the arithmetic
    paths are disjoint.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


def project_step_function_quad(
    breakpoints: np.ndarray,
    values: np.ndarray,
    T: int,
    n_gauss: int = 64,
):
    """Project a step function onto White's truncated Fourier basis via
    per-cell Gauss-Legendre quadrature.

    Parameters
    ----------
    breakpoints : array of shape (n_cells + 1,)
        Cell endpoints in [-2, 2], strictly increasing.
    values : array of shape (n_cells,)
        Constant value of f on each cell.
    T : int
        Number of Fourier modes to project onto (k = 1..T).
    n_gauss : int
        Order of the Gauss-Legendre quadrature per cell. Default 64
        (degree-127-exact polynomial integration; ample for our regime).

    Returns
    -------
    c, d : np.ndarray of shape (T + 1,)
        c[0] = 0.5, d[0] = 0.0 (convention placeholders).
        c[k] = integral_{-2}^{2} f(x) cos(pi k x / 2) dx,  k = 1..T
        d[k] = integral_{-2}^{2} f(x) sin(pi k x / 2) dx,  k = 1..T
    """
    breakpoints = np.asarray(breakpoints, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    n_cells = len(values)
    if len(breakpoints) != n_cells + 1:
        raise ValueError(
            f"breakpoints has length {len(breakpoints)}, expected "
            f"{n_cells + 1} (one more than values)."
        )
    if T < 1:
        raise ValueError(f"T must be >= 1; got {T}.")

    nodes, weights = leggauss(n_gauss)  # nodes, weights on [-1, 1]
    c = np.zeros(T + 1, dtype=np.float64)
    d = np.zeros(T + 1, dtype=np.float64)
    c[0] = 0.5
    d[0] = 0.0

    k_arr = np.arange(1, T + 1, dtype=np.float64)
    omega_arr = np.pi * k_arr / 2.0  # shape (T,)

    for i in range(n_cells):
        wv = values[i]
        if wv == 0.0:
            continue
        a = breakpoints[i]
        b = breakpoints[i + 1]
        half = (b - a) / 2.0
        mid = (a + b) / 2.0
        x_cell = mid + half * nodes  # shape (n_gauss,)
        # arg[k_idx, g_idx] = omega_{k+1} * x_g
        arg = np.outer(omega_arr, x_cell)  # shape (T, n_gauss)
        # Per-cell integral: half * sum_g w_g cos/sin(omega_k * x_g)
        cos_int = (np.cos(arg) @ weights) * half  # shape (T,)
        sin_int = (np.sin(arg) @ weights) * half  # shape (T,)
        c[1:] += wv * cos_int
        d[1:] += wv * sin_int

    return c, d


# --- Cross-check against the closed-form projection ---------------------------


def cross_verify():
    # Prefer absolute import (matches project convention).
    try:
        from lp_research_state.code.together_loader import (
            load_together_raw,
            to_white_convention_even,
            to_white_convention_direct,
        )
        from lp_research_state.code.together_diagnostic import (
            project_step_function,
        )
    except ImportError:
        sys.path.insert(0, str(Path(__file__).parent))
        from together_loader import (  # type: ignore
            load_together_raw,
            to_white_convention_even,
            to_white_convention_direct,
        )
        from together_diagnostic import project_step_function  # type: ignore

    bp_t, vals_t, dom, _ = load_together_raw()
    embeddings = [
        ("even", to_white_convention_even),
        ("direct", to_white_convention_direct),
    ]
    overall_ok = True
    for kind, embed in embeddings:
        wb, wv = embed(bp_t, vals_t, dom)
        c_fast, d_fast = project_step_function(wb, wv, T=4000)
        c_slow, d_slow = project_step_function_quad(
            wb, wv, T=4000, n_gauss=64
        )
        max_diff_c = float(np.max(np.abs(c_fast - c_slow)))
        max_diff_d = float(np.max(np.abs(d_fast - d_slow)))
        max_diff = max(max_diff_c, max_diff_d)
        denom = max(
            float(np.max(np.abs(c_fast))), float(np.max(np.abs(d_fast))), 1e-300
        )
        rel_diff = max_diff / denom
        print(
            f"{kind}: max abs diff = {max_diff:.3e}  "
            f"(c: {max_diff_c:.3e}, d: {max_diff_d:.3e}), "
            f"max rel diff = {rel_diff:.3e}"
        )
        if max_diff >= 1e-10:
            overall_ok = False
            print(
                f"  !! {kind} cross-verify FAILED: max_diff "
                f"{max_diff:.3e} >= 1e-10"
            )
        assert max_diff < 1e-10, (
            f"{kind} cross-verify FAILED at {max_diff:.3e}; should be < 1e-10"
        )
    if overall_ok:
        print(
            "OK: Fourier projections agree to >= 10 digits on both embeddings"
        )


if __name__ == "__main__":
    cross_verify()
