"""
Diagnostic: evaluate every constraint in our SDP at Together's f*.

Task 3 (this file, initial skeleton): project Together's f* into White's
truncated Fourier basis to produce (c, d) arrays of length T+1 = 4001.

Outputs (this task):
  lp_research_state/data/together_f_star_fourier_even.npz   — (c, d) from
                                                             even embedding
  lp_research_state/data/together_f_star_fourier_direct.npz — (c, d) from
                                                             direct embedding

Future tasks will populate:
  lp_research_state/data/together_diagnostic_results.json   — constraint slacks
  lp_research_state/data/together_gap_function.npz          — f̃(x) - f*(x)

Fourier convention (per `_fourier_convention_notes.md`, verified against
`white_full_convex.py:230` and surrounding code):
  Domain: [-2, 2], length 4. Basis: cos(π k x / 2), sin(π k x / 2), k=1..T.
  f̂(k) := (1/2) ∫_{-2}^{2} f(x) e^{-i π k x / 2} dx.
  white_full_convex.py:230 states  f̂(0) = 1/2,  f̂(k) = (c[k-1] − i·d[k-1]) / 2.

  Therefore:
      c[k-1] = ∫_{-2}^{2} f(x) cos(π k x / 2) dx
      d[k-1] = ∫_{-2}^{2} f(x) sin(π k x / 2) dx        for k = 1..T.

  For a step function with values wv_i on cells [b_i, b_{i+1}], with ω = π k / 2:
      c[k-1] = Σ_i wv_i · (sin(ω b_{i+1}) − sin(ω b_i)) / ω
      d[k-1] = Σ_i wv_i · (cos(ω b_i) − cos(ω b_{i+1})) / ω

Important caveat for c[0], d[0]:
  In White's variable layout, c[0] and d[0] are *parameters* with bounds
  [p1, p2], [q1, q2] (line 201) — they pin the lowest Fourier mode (k=1)
  to a small rectangle defining a "row" of the residual region. They are
  NOT a slot for f̂(0); f̂(0) = 1/2 is hardcoded.

  This `project_step_function` returns c, d as arrays of length T+1 indexed
  by k = 0..T, where:
      c[k] (k >= 1)  = the k-th cosine Fourier coefficient ∫ f cos(πkx/2) dx
                      (so c[1] corresponds to White's c[0] variable, etc.)
      c[0]          = f̂(0)-style placeholder = 0.5  (NOT used in the SDP
                      as a Fourier mode; provided for symmetry with the
                      math notation only).
      d[k] (k >= 1)  = analogous for sine.
      d[0]          = 0.0.

  Callers that map this projection back into the SDP's `c, d` arrays of
  length T should slice [1:T+1].
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

# Robust import: prefer absolute (works under `python -c "from
# lp_research_state.code.together_diagnostic import ..."` from repo root,
# the invocation pattern used in the Task-3 plan).
try:
    from lp_research_state.code.together_loader import (
        load_together_raw,
        to_white_convention_even,
        to_white_convention_direct,
    )
except ImportError:  # pragma: no cover - fallback for direct script runs
    sys.path.insert(0, str(Path(__file__).parent))
    from together_loader import (  # type: ignore
        load_together_raw,
        to_white_convention_even,
        to_white_convention_direct,
    )


DATA_DIR = Path(__file__).parent.parent / "data"


def project_step_function(
    breakpoints: np.ndarray, values: np.ndarray, T: int
):
    """Project a step function (in White's [-2, 2] domain) onto the first
    T cosine and sine Fourier modes (k = 1..T) of period 4.

    Closed-form per-cell integrals are used — NO numerical quadrature.

    Parameters
    ----------
    breakpoints : np.ndarray of shape (n_cells + 1,)
        Cell endpoints in [-2, 2], strictly increasing.
    values : np.ndarray of shape (n_cells,)
        Constant value of f on each cell.
    T : int
        Number of Fourier modes to project onto (k = 1..T).

    Returns
    -------
    c : np.ndarray of shape (T + 1,)
        c[0] = 0.5 (White's f̂(0) convention; not a Fourier mode in the
        usual sense, see module docstring).
        c[k] for k = 1..T equals  ∫_{-2}^{2} f(x) cos(π k x / 2) dx.
    d : np.ndarray of shape (T + 1,)
        d[0] = 0.0 (sine coefficient at k = 0 is identically zero).
        d[k] for k = 1..T equals  ∫_{-2}^{2} f(x) sin(π k x / 2) dx.

    Notes
    -----
    For a step function with values wv_i on cells [b_i, b_{i+1}]:
        c[k] = Σ_i wv_i · [sin(ω b_{i+1}) − sin(ω b_i)] / ω
        d[k] = Σ_i wv_i · [cos(ω b_i) − cos(ω b_{i+1})] / ω
    with ω = π k / 2.  We use vectorized numpy ops over cells.
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

    c = np.zeros(T + 1, dtype=np.float64)
    d = np.zeros(T + 1, dtype=np.float64)
    # f̂(0) = 1/2 in White's convention (see line 230 of white_full_convex.py).
    # The d[0] slot has no Fourier meaning at k = 0 (sin(0) = 0).
    c[0] = 0.5
    d[0] = 0.0

    b_left = breakpoints[:-1]   # shape (n_cells,)
    b_right = breakpoints[1:]   # shape (n_cells,)

    k_arr = np.arange(1, T + 1, dtype=np.float64)   # (T,)
    omega = np.pi * k_arr / 2.0                     # (T,)

    # Compute sin(ω b) and cos(ω b) for each k, each breakpoint endpoint.
    # Shape: (T, n_cells).
    # outer-products: omega[:, None] * b_left[None, :].
    arg_left = np.outer(omega, b_left)
    arg_right = np.outer(omega, b_right)
    sin_left = np.sin(arg_left)
    sin_right = np.sin(arg_right)
    cos_left = np.cos(arg_left)
    cos_right = np.cos(arg_right)

    # c[k] = Σ_i wv_i (sin(ω b_{i+1}) − sin(ω b_i)) / ω
    #      = (1/ω) Σ_i wv_i (sin_right_i − sin_left_i)
    delta_sin = sin_right - sin_left           # (T, n_cells)
    delta_cos_neg = cos_left - cos_right       # (T, n_cells)
    c[1:] = (delta_sin @ values) / omega
    d[1:] = (delta_cos_neg @ values) / omega

    return c, d


# --- Unit tests -----------------------------------------------------------


def _test_projection_constant():
    """Sanity check: f = 1/4 on [-2, 2] gives c[1:] = d[1:] = 0 exactly.

    For a single cell with value 1/4 on [-2, 2], ω = πk/2, and
    sin(πk) - sin(-πk) = 0 for all integer k. Same for the cosine
    formula via cos(-πk) - cos(πk) = 0. So all c[k>=1], d[k>=1] vanish.
    """
    breakpoints = np.array([-2.0, 2.0])
    values = np.array([0.25])
    T = 50
    c, d = project_step_function(breakpoints, values, T)
    assert abs(c[0] - 0.5) < 1e-15, f"c[0]={c[0]}"
    assert abs(d[0]) < 1e-15, f"d[0]={d[0]}"
    max_c = float(np.max(np.abs(c[1:])))
    max_d = float(np.max(np.abs(d[1:])))
    assert max_c < 1e-12, f"max|c[1:]|={max_c:.3e} (expected ~0)"
    assert max_d < 1e-12, f"max|d[1:]|={max_d:.3e} (expected ~0)"
    print("[OK] projection on constant f = 1/4: c[1:]=0, d[1:]=0")


def _test_projection_single_cell():
    """Single-cell step: f(x) = 1 on [0, 1], 0 elsewhere on [-2, 2].

    Closed-form Fourier coefficients (ω = πk/2):
        c[k] = ∫_0^1 cos(ωx) dx = sin(ω)/ω = (2/(πk)) · sin(πk/2)
        d[k] = ∫_0^1 sin(ωx) dx = (1 - cos(ω))/ω
                                = (2/(πk)) · (1 - cos(πk/2))
    """
    breakpoints = np.array([-2.0, 0.0, 1.0, 2.0])
    values = np.array([0.0, 1.0, 0.0])
    T = 20
    c, d = project_step_function(breakpoints, values, T)
    for k in [1, 2, 3, 5, 10]:
        expected_c = (2.0 / (np.pi * k)) * np.sin(np.pi * k / 2)
        expected_d = (2.0 / (np.pi * k)) * (1.0 - np.cos(np.pi * k / 2))
        err_c = abs(c[k] - expected_c)
        err_d = abs(d[k] - expected_d)
        assert err_c < 1e-12, (
            f"c[{k}]: got {c[k]:.12e}, expected {expected_c:.12e}, "
            f"err={err_c:.3e}"
        )
        assert err_d < 1e-12, (
            f"d[{k}]: got {d[k]:.12e}, expected {expected_d:.12e}, "
            f"err={err_d:.3e}"
        )
    print("[OK] projection on single-cell step matches closed form")


def _test_projection_even_symmetry():
    """Symmetry: a function symmetric about x = 0 (f(-x) = f(x)) must
    have d[k] = 0 for all k >= 1.

    Use a small step function symmetric about 0 to verify, before we run
    on Together's full f_even.
    """
    breakpoints = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    values = np.array([0.1, 0.4, 0.4, 0.1])
    T = 30
    c, d = project_step_function(breakpoints, values, T)
    max_d = float(np.max(np.abs(d[1:])))
    assert max_d < 1e-13, (
        f"symmetric step gave max|d[1:]|={max_d:.3e} (expected ~0)"
    )
    print("[OK] projection on symmetric step gives d[1:] = 0")


def truncation_tail_bound(breakpoints, values, T):
    """Rigorous upper bound on ||f - f_T||_{L^2([-2,2])} and ||f - f_T||_{L^1([-2,2])}
    for a step function with breakpoints in [-2, 2], where f_T is the
    truncation of f to its first T cosine + sine Fourier modes (k = 1..T)
    in White's convention (basis cos(πkx/2), sin(πkx/2), domain length 4).

    Derivation
    ----------
    Let c_k = ∫_{-2}^{2} f(x) cos(πkx/2) dx, d_k = ∫_{-2}^{2} f(x) sin(πkx/2) dx.
    For a step function with breakpoints b_0 < b_1 < ... < b_M and values
    v_0, ..., v_{M-1} (extended by 0 outside [b_0, b_M] ⊆ [-2, 2]),
    integration by parts gives (with ω = πk/2):
        c_k = (1/ω) · Σ_i (v_i - v_{i-1}) · sin(ω b_i)
    where the sum runs over ALL jumps including the boundary "jumps" from
    0 → v_0 at b_0 and v_{M-1} → 0 at b_M (since f = 0 outside support).
    Hence |c_k| ≤ (1/ω) · V = (2V/(π k)) where
        V := Σ_i |v_i - v_{i-1}|    (total variation, boundary-extended)
    and the same bound holds for |d_k|.

    By Parseval — with the orthonormal basis
        {1/√4} ∪ {cos(πkx/2)/√2, sin(πkx/2)/√2 : k ≥ 1}   on [-2, 2] —
    the inner products are √2 · c_k / 2 and √2 · d_k / 2 ... equivalently
        ||f||_{L^2}^2 = (∫f)^2 / 4 + Σ_{k≥1} (c_k^2 + d_k^2) / 2
    so that
        ||f - f_T||_{L^2}^2 = Σ_{k>T} (c_k^2 + d_k^2) / 2
                            ≤ Σ_{k>T} (2V/(πk))^2
                            = (4 V^2 / π^2) · Σ_{k>T} 1/k^2.
    Using the standard tail bound Σ_{k>T} 1/k^2 ≤ ∫_T^∞ dx/x^2 = 1/T:
        ||f - f_T||_{L^2}^2 ≤ (4 V^2)/(π^2 T)
        ||f - f_T||_{L^2}   ≤ (2 V / π) / √T.
    For L^1, Cauchy–Schwarz on a domain of length L = 4 gives:
        ||f - f_T||_{L^1} ≤ √L · ||f - f_T||_{L^2} = 2 · ||f - f_T||_{L^2}.

    Parameters
    ----------
    breakpoints : np.ndarray, shape (n_cells + 1,)
        Cell endpoints in [-2, 2].
    values : np.ndarray, shape (n_cells,)
        Constant value of f on each cell.
    T : int
        Truncation level (keep modes k = 1..T).

    Returns
    -------
    dict with keys:
        'V'             : float, total variation (boundary-extended).
        'L2_bound'      : float, rigorous bound on ||f - f_T||_{L^2}.
        'L1_bound'      : float, rigorous bound on ||f - f_T||_{L^1}.
        'T'             : int, the truncation level.
        'domain_length' : float, b_M - b_0 (should be 4.0 on [-2, 2]).
    """
    breakpoints = np.asarray(breakpoints, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    # Boundary-extended values: prepend & append 0 to capture support-edge jumps.
    v_extended = np.concatenate([[0.0], values, [0.0]])
    V = float(np.sum(np.abs(np.diff(v_extended))))
    domain_len = float(breakpoints[-1] - breakpoints[0])
    L2_bound = (2.0 * V / np.pi) / np.sqrt(float(T))
    L1_bound = np.sqrt(domain_len) * L2_bound
    return {
        "V": V,
        "L2_bound": L2_bound,
        "L1_bound": L1_bound,
        "T": int(T),
        "domain_length": domain_len,
    }


def _test_total_variation():
    """For f(x) = 1 on [0, 1], 0 elsewhere on [-2, 2]:
    Including boundary jumps (f = 0 outside the support), the total
    variation is exactly 2 — one unit jump up at x = 0, one unit jump
    down at x = 1.
    """
    bp = np.array([-2.0, 0.0, 1.0, 2.0])
    vals = np.array([0.0, 1.0, 0.0])
    tb = truncation_tail_bound(bp, vals, T=100)
    assert abs(tb["V"] - 2.0) < 1e-12, f"V should be 2, got {tb['V']}"
    print("[OK] total variation calculation correct (V = 2 for [0,1] indicator)")


def truncation_tail_exact(breakpoints, values, T):
    """Exact L² truncation tail error for a step function on [-2, 2], using Parseval.

    For f = a_0 + Σ_k [c_k cos(πkx/2) + d_k sin(πkx/2)] with the convention
    c_k = ∫ f cos(πkx/2) dx (per _fourier_convention_notes.md), Parseval gives:
        ||f||² = (∫f)² / 4 + Σ_k (c_k² + d_k²) / 2

    Hence ||f - f_T||² = ||f||² - (∫f)²/4 - Σ_{k=1}^{T} (c_k² + d_k²)/2

    This is RIGOROUS (an exact equality), strictly tighter than the V-based bound
    in `truncation_tail_bound(...)`, suitable as the trust threshold for downstream
    tasks comparing constraint slacks.

    Returns
    -------
    dict with keys 'L2_sq_exact', 'L2_exact', 'L1_bound', 'T', plus
    'f_norm_sq', 'f_integral_sq_over_4', 'truncated_energy' for sanity-checking.
    """
    breakpoints = np.asarray(breakpoints, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    # Closed-form ||f||² and ∫f for a step function:
    widths = np.diff(breakpoints)
    f_norm_sq = float(np.sum(values**2 * widths))
    f_int = float(np.sum(values * widths))

    # Use the existing project_step_function for (c, d); guard T=0 since the
    # projection helper requires T >= 1 (no modes to integrate against).
    if T == 0:
        truncated_energy = 0.0
    else:
        c, d = project_step_function(breakpoints, values, T)
        # c[0] is a placeholder; the real projection coefs are c[1..T], d[1..T]
        truncated_energy = float(0.5 * np.sum(c[1:] ** 2 + d[1:] ** 2))
    constant_energy = (f_int**2) / 4.0

    L2_sq_exact = f_norm_sq - constant_energy - truncated_energy
    # Should be nonneg up to floating-point noise; clamp to 0 if very slightly negative
    if -1e-14 < L2_sq_exact < 0:
        L2_sq_exact = 0.0
    if L2_sq_exact < 0:
        raise ValueError(
            f"L²_sq_exact negative ({L2_sq_exact:.3e}); Parseval convention "
            f"or projection has a bug"
        )

    L2_exact = float(np.sqrt(L2_sq_exact))
    domain_len = float(breakpoints[-1] - breakpoints[0])
    L1_bound = float(np.sqrt(domain_len) * L2_exact)
    return {
        "L2_sq_exact": L2_sq_exact,
        "L2_exact": L2_exact,
        "L1_bound": L1_bound,
        "T": int(T),
        "f_norm_sq": f_norm_sq,
        "f_integral_sq_over_4": constant_energy,
        "truncated_energy": truncated_energy,
    }


def _test_truncation_tail_exact_smooth():
    """For f(x) = 1/4 (constant on [-2, 2]), ||f||² = 1, (∫f)²/4 = 1, so tail = 0 for any T ≥ 1."""
    bp = np.array([-2.0, 2.0])
    vals = np.array([0.25])
    r = truncation_tail_exact(bp, vals, T=10)
    assert r["L2_sq_exact"] < 1e-20, (
        f"constant function should have zero tail, got {r['L2_sq_exact']:.3e}"
    )
    print(f"[OK] exact tail = 0 for constant function (got {r['L2_sq_exact']:.3e})")


def _test_truncation_tail_exact_one_mode():
    """For f representing a single Fourier mode (approximately), truncating at T=0 should leave that mode's energy.

    Use f(x) = 0.25 + 0.5 cos(πx/2) on [-2, 2] (approximated as a step function with many cells).
    Then c[1] ≈ 1 (coefficient of cos(πx/2) is 0.5, the integral ∫ f cos(πx/2) dx = 0.5 · 2 = 1).
    Truncating at T=0 should leave tail² ≈ (1² + 0²)/2 = 0.5.

    ||f||² closed form: ∫(0.25 + 0.5 cos(πx/2))² dx
                  = ∫(0.0625 + 0.25 cos(πx/2) + 0.25 cos²(πx/2)) dx
                  = 0.25 + 0 + 0.5 = 0.75   (using ∫_{-2}^{2} cos²(πx/2) dx = 2)
    ∫f = 0.25·4 + 0.5·0 = 1, so (∫f)²/4 = 0.25
    expected L²_sq tail = 0.75 - 0.25 = 0.5
    """
    n = 10000
    x = np.linspace(-2, 2, n + 1)
    f_vals = 0.25 + 0.5 * np.cos(np.pi * (x[:-1] + x[1:]) / 4)  # midpoint sample
    r = truncation_tail_exact(x, f_vals, T=0)
    assert abs(r["L2_sq_exact"] - 0.5) < 1e-3, (
        f"expected L²_sq ≈ 0.5, got {r['L2_sq_exact']:.6f}"
    )
    print(f"[OK] T=0 tail ≈ 0.5: got {r['L2_sq_exact']:.6f}")


def _run_all_tests():
    _test_projection_constant()
    _test_projection_single_cell()
    _test_projection_even_symmetry()
    _test_total_variation()
    _test_truncation_tail_exact_smooth()
    _test_truncation_tail_exact_one_mode()
    print("[ALL] projection tests passed")


# --- Together f* projection ----------------------------------------------


def project_together_f_star(T: int = 4000, kind: str = "even"):
    """Project Together's h* (via the chosen White embedding) onto T
    Fourier modes and save the result as an .npz file.

    Parameters
    ----------
    T : int
        Number of modes (k = 1..T). Default 4000 matches the Phase-5 SDP.
    kind : str
        "even"   — symmetric reflection f(x) = h(|x|)/2 on [-2, 2].
        "direct" — asymmetric f(x) = h(x) on [0, 2], zero on [-2, 0].

    Returns
    -------
    c, d : np.ndarray, each of shape (T+1,)

    Side effect: writes lp_research_state/data/together_f_star_fourier_<kind>.npz
    with keys c, d, T, kind.
    """
    bp_t, vals_t, dom, _ = load_together_raw()
    if kind == "even":
        wb, wv = to_white_convention_even(bp_t, vals_t, dom)
    elif kind == "direct":
        wb, wv = to_white_convention_direct(bp_t, vals_t, dom)
    else:
        raise ValueError(f"unknown kind {kind!r}; use 'even' or 'direct'")
    c, d = project_step_function(wb, wv, T)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out_path = DATA_DIR / f"together_f_star_fourier_{kind}.npz"
    np.savez(out_path, c=c, d=d, T=np.int64(T), kind=kind)
    return c, d


if __name__ == "__main__":
    _run_all_tests()
    for kind in ("even", "direct"):
        c, d = project_together_f_star(T=4000, kind=kind)
        print(
            f"{kind}: c[0]={c[0]:.4f} d[0]={d[0]:.4f} "
            f"max|c[1:]|={np.max(np.abs(c[1:])):.4e} "
            f"max|d[1:]|={np.max(np.abs(d[1:])):.4e}"
        )
