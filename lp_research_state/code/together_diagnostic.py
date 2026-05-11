"""
Diagnostic: evaluate every constraint in our SDP at Together's f*.

Task 3 (this file, initial skeleton): project Together's f* into White's
truncated Fourier basis to produce (c, d) arrays of length T+1 = 4001.

Task 6 (added below): with (c, d) PINNED to the projection, minimize Ω over
the SDP's remaining variables.  This yields the SDP's certified upper bound
on sup_t ∫ f f(·+t) dx for THAT specific f.

Results (Task 6, row 4, bochner_n=30, N=10000, T=4000, R=10):
  - row 4 box check: NEITHER embedding lands in row 4's residual box
    [p1,p2]×[q1,q2] = [0.3875,0.3875]×[-0.02,0.02].
    Together's projected f̂(1) is essentially zero (c1 ≈ -2.3e-4) for both,
    because its mass is concentrated near the middle of the integer-shifted
    overlap region, not the row's residual corner.  We override the row
    bounds (p1=p2=c1, q1=q2=d1) so (c[0], d[0]) can be pinned consistently.

  - even embedding:   status=optimal,    Ω = 0.459311  (solve 72.7 s)
  - direct embedding: status=INFEASIBLE                (solve 24.2 s)

  Comparison numbers:
    * White Phase-5 SDP optimum (row 4)        : Ω ≈ 0.380128
    * True autocorrelation of f_even           : 0.387337
    * True autocorrelation of f_direct         : 0.774675
    * Ω_SDP(f_even) at pinned (c, d)           : 0.459311  (this task)
    * Ω_SDP(f_direct) at pinned (c, d)         : infeasible

  Why direct is infeasible:  the projection violates the SDP's pre-Bochner
  trigonometric envelope.  |d̂(1)| = 0.807 > 2/π ≈ 0.637 (white_full_convex
  line 199), and Σ(c² + d²) = 1.05 > 0.5 (line 200).  These are Parseval-side
  necessary conditions for 0 ≤ f ≤ 1 with ∫f = 1 on [-2, 2] — the direct
  embedding, which concentrates all mass on [0, 2] with f = 1 there,
  trivially fails them.  The SDP is *correctly* rejecting f_direct as an
  invalid candidate for the symmetric f* class.

  Why Ω(f_even) = 0.459 ≫ 0.387 (autocorr):  the SDP's cell-tested
  autocorrelation upper bound on f_even's Ω is much *looser* than f_even's
  true autocorrelation.  Bochner-PSD at level 30 doesn't tighten that looser
  envelope back to the true value.  Concretely, the cell-kernel bounds in
  white_full_convex.py:176–190 over-estimate sup ∫f(x)f(x+t) dt for a step
  function whose Fourier tail is heavy (Together's f* has bumpy support).
  This says the SDP relaxation is structurally non-tight on Together's f*,
  which is consistent with — and helps explain — the gap between our
  rigorous SDP bound 0.380128 and Together's heuristic upper-bound proxy.

Outputs of this module:
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


# --- Task 6: Evaluate Ω(f*) by pinning (c, d) in White's SDP ---------------
#
# Conceptual note (from Task 2):
#   White's `Ω` is an SDP variable bounded from above by:
#     (i) pointwise upper bound:   Ω ≥ w_j, Ω ≥ v_j  for all j (line 141)
#     (ii) cos/sin cell-kernel autocorrelation upper bounds for m=1..2R (lines 176–190)
#     (iii) the Bochner-PSD relaxation of f ≥ 0, 1 − f ≥ 0 (lines 228–258)
#
#   When (c, d) are PINNED to Together's projected f*, the SDP minimizes Ω over
#   (w, v, eps, dlt, …) subject to all those constraints. The resulting Ω is the
#   SDP's certified upper bound on sup_t ∫ f f(·+t) dx FOR THAT specific f.
#
#   If Ω_at_f* > White's SDP-optimal Ω (≈ 0.380128 at Phase 5), Together's f is a
#   strictly worse Ω-point than White's SDP optimum — the diagnostic conclusion.
#
# Row 4 (per CLAUDE.md): (h1,h2,p1,p2,q1,q2) = (0.004, 0.004, 0.3875, 0.3875,
# −0.02, +0.02). So p is point-pinned to 0.3875 and q ranges in [−0.02, +0.02].


def check_row_box_for_projection(c_proj, d_proj, row: str = "row4"):
    """Check whether Together's projected first Fourier mode f̂(1) lies in row's box.

    SDP's c[0] = proj_c[1], SDP's d[0] = proj_d[1]. Row 4 has
    p1 = p2 = 0.3875, q1 = −0.02, q2 = +0.02.

    Returns dict with row bounds, the projection's (c1_proj, d1_proj),
    and `in_box` flag.
    """
    rows = {"row4": (0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02)}
    if row not in rows:
        raise ValueError(f"unknown row {row!r}; known: {list(rows)}")
    h1, h2, p1, p2, q1, q2 = rows[row]
    c1 = float(c_proj[1])  # SDP's c[0] (k=1 cosine Fourier coef)
    d1 = float(d_proj[1])  # SDP's d[0] (k=1 sine   Fourier coef)
    in_box = (p1 - 1e-9 <= c1 <= p2 + 1e-9) and (q1 - 1e-9 <= d1 <= q2 + 1e-9)
    return {
        "row": row,
        "h1": h1, "h2": h2,
        "p1": p1, "p2": p2, "q1": q1, "q2": q2,
        "c1_proj": c1, "d1_proj": d1,
        "in_box": bool(in_box),
    }


def evaluate_omega_at_f_star(
    c_proj,
    d_proj,
    row: str = "row4",
    bochner_n: int = 30,
    N: int = 10000,
    T: int = 4000,
    R: int = 10,
    verbose: bool = False,
    override_row_bounds: bool = True,
):
    """Pin (c, d) in White's SDP to Together's projected f* and minimize Ω.

    SDP variable mapping (verified against white_full_convex.py:135–136, 201):
        SDP's c[k]  ↔  Fourier mode k+1  ↔  proj's c[k+1]   for k = 0..T-1
        SDP's d[k]  ↔  proj's d[k+1]                          for k = 0..T-1

    The row's residual-region pins SDP's c[0] ∈ [p1, p2], d[0] ∈ [q1, q2]
    (line 201). If `override_row_bounds=True` and Together's f̂(1) doesn't sit
    inside [p1,p2] × [q1,q2], we force p1=p2=proj_c[1], q1=q2=proj_d[1] so the
    SDP can pin (c[0], d[0]) consistently with the rest of the projection.

    Args
    ----
    c_proj, d_proj : np.ndarray of length T+1, in this module's convention.
    row : str — currently only "row4" defined.
    bochner_n : int — Bochner PSD level. Default 30 matches Phase 5.
    N, T, R : SDP discretization params (defaults match the CLAUDE.md headline).
    verbose : pass to CLARABEL.
    override_row_bounds : see note above.

    Returns dict with status, `Omega_at_f_star`, and config.
    """
    import cvxpy as cp  # local import; this module is also imported in lightweight contexts

    # Robust import of build_problem.
    try:
        from lp_research_state.code.white_full_convex import build_problem
    except ImportError:  # pragma: no cover - fallback for direct script runs
        sys.path.insert(0, str(Path(__file__).parent))
        from white_full_convex import build_problem  # type: ignore

    rows = {"row4": (0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02)}
    if row not in rows:
        raise ValueError(f"unknown row {row!r}; known: {list(rows)}")
    h1, h2, p1, p2, q1, q2 = rows[row]

    c1 = float(c_proj[1])
    d1 = float(d_proj[1])
    if override_row_bounds:
        p1 = p2 = c1
        q1 = q2 = d1

    Omega, w, v, c_var, d_var, eps, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, q1, q2,
        bochner_n=bochner_n,
    )

    # Pin remaining c[1..T-1], d[1..T-1] of the SDP. (SDP's c[0], d[0] are pinned
    # via p1=p2, q1=q2 above.) SDP index k corresponds to proj index k+1.
    pin = []
    for k in range(1, T):
        pin.append(c_var[k] == float(c_proj[k + 1]))
        pin.append(d_var[k] == float(d_proj[k + 1]))

    prob = cp.Problem(cp.Minimize(Omega), cons + pin)
    prob.solve(solver=cp.CLARABEL, verbose=verbose)

    return {
        "status": prob.status,
        "Omega_at_f_star": float(Omega.value) if Omega.value is not None else None,
        "row": row,
        "row_bounds_overridden": bool(override_row_bounds),
        "c1_pinned": c1,
        "d1_pinned": d1,
        "bochner_n": bochner_n,
        "N": N, "T": T, "R": R,
    }


# =========================================================================
# Tasks 7-9: per-constraint-family slacks at Together's f*
# =========================================================================
#
# Convention reminder (CRITICAL — verified against bochner.py and
# poly_moment.py):
#   - This module's projection returns c, d as length T+1 arrays
#     indexed by k = 0..T, where c[0]=0.5 (placeholder for f̂(0))
#     and c[k] = ∫_{-2}^{2} f cos(πkx/2) dx for k = 1..T.
#   - The SDP's c[k] (k = 0..T-1) corresponds to proj_c[k+1].
#   - White's Fourier coeffs: f̂(0) = 1/2, f̂(k) = (proj_c[k] - i proj_d[k]) / 2
#     for k = 1..T (this module's projection indexing directly).
#   - poly_moment.py's m_k_SDP formula equals (1/2^k) · ∫_{-2}^{2} x^k f(x) dx
#     because it's written on the [-1, 1] basis with c, d as ∫_{-2}^{2}
#     f cos(πkx/2) dx (i.e. SDP's c, d):
#         m_k_SDP = ∫_{-1}^{1} x^k · [1/2 + Σ_j (c_j cos(πjx) + d_j sin(πjx))] dx
#         Plugging c_j = 2·(g's [-1,1] cos coef) where g(y)=f(2y) gives
#         m_k_SDP = 2 · ∫_{-1}^{1} y^k g(y) dy = (1/2^k) · ∫_{-2}^{2} x^k f(x) dx
#     (derivation in the task plan notes; verified numerically below).


# --- Task 7: Bochner-PSD slack at f* ------------------------------------


def bochner_matrix_at_f_star(c_proj, d_proj, n):
    """Construct M_n(f*) = [f̂(j-k)]_{j,k=0..n} (Hermitian (n+1)×(n+1)).

    Uses White's convention:
        f̂(0) = 1/2,  f̂(k) = (c_proj[k] - i d_proj[k]) / 2  for k ≥ 1.

    Parameters
    ----------
    c_proj, d_proj : array-like of length ≥ n+1
        Projection-indexed Fourier coefficients (c[k] is the k-th cosine
        Fourier integral; this module's `project_step_function` output).
    n : int
        Bochner level. Returned matrix has shape (n+1, n+1).

    Returns
    -------
    M : np.ndarray of shape (n+1, n+1), complex
        Hermitian Toeplitz matrix.
    """
    c_proj = np.asarray(c_proj, dtype=np.float64)
    d_proj = np.asarray(d_proj, dtype=np.float64)
    f_hat = np.empty(n + 1, dtype=complex)
    f_hat[0] = 0.5
    for k in range(1, n + 1):
        f_hat[k] = (float(c_proj[k]) - 1j * float(d_proj[k])) / 2.0
    M = np.empty((n + 1, n + 1), dtype=complex)
    for j in range(n + 1):
        for k in range(n + 1):
            diff = j - k
            M[j, k] = f_hat[diff] if diff >= 0 else np.conj(f_hat[-diff])
    return M


def _bochner_xcheck(c_proj, d_proj, n=4):
    """Cross-check bochner_matrix_at_f_star against bochner.py's encoded
    construction by evaluating the same real-form (2n+2)×(2n+2) matrix and
    asserting bit-equality.

    bochner.py builds Re_M, Im_M (n+1)x(n+1) via:
        Re M[j,k] = 1/2 if j==k, else sign * c_{|j-k|-1+1}/2 = sign * c_{|j-k|}/2
                    (SDP indexing: SDP's c[m-1] is the m-th Fourier coef
                    = proj's c_proj[m])
        Im M[j,k] = 0 if j==k, else
                     -sign * d_{|j-k|}/2  if j>k (ell>0),
                     +sign * d_{|j-k|}/2  if j<k (ell<0)
    Our M is then Re_M + i Im_M (taking sign=+1 for f≥0).
    """
    n_ = n
    c_proj = np.asarray(c_proj, dtype=np.float64)
    d_proj = np.asarray(d_proj, dtype=np.float64)
    Re_M = np.zeros((n_ + 1, n_ + 1))
    Im_M = np.zeros((n_ + 1, n_ + 1))
    sign = +1
    for j in range(n_ + 1):
        for k in range(n_ + 1):
            ell = j - k
            if ell == 0:
                Re_M[j, k] = 0.5
                Im_M[j, k] = 0.0
            else:
                aell = abs(ell)
                # bochner.py uses SDP c-indexing: SDP's c[aell-1] is the
                # aell-th Fourier mode, which in our proj-indexing is
                # c_proj[aell].
                Re_M[j, k] = sign * 0.5 * c_proj[aell]
                if ell > 0:
                    Im_M[j, k] = -sign * 0.5 * d_proj[aell]
                else:
                    Im_M[j, k] = +sign * 0.5 * d_proj[aell]
    M_from_bochner_py = Re_M + 1j * Im_M
    M_ours = bochner_matrix_at_f_star(c_proj, d_proj, n_)
    diff = np.max(np.abs(M_from_bochner_py - M_ours))
    return float(diff)


def _test_bochner_xcheck():
    rng = np.random.default_rng(0)
    c = np.concatenate([[0.5], rng.normal(size=20) * 0.1])
    d = np.concatenate([[0.0], rng.normal(size=20) * 0.1])
    for n_ in (1, 2, 4, 8):
        diff = _bochner_xcheck(c, d, n=n_)
        assert diff < 1e-14, f"bochner xcheck failed at n={n_}: diff={diff}"
    print("[OK] bochner_matrix_at_f_star matches bochner.py construction to ~1e-14")


def bochner_diagnostic(c_proj, d_proj, n=30):
    """Compute λ_min, λ_max of M_n(f) and M_n(1-f) at the given f*.

    For 1-f: (1-f)̂(0) = 1 - 1/2 = 1/2; (1-f)̂(k) = -f̂(k) for k ≥ 1.
    Equivalently: c_{1-f, k} = -c_{f, k}, d_{1-f, k} = -d_{f, k} for k ≥ 1,
    with c_{1-f}[0] kept at 0.5 to encode (1-f)̂(0) = 1/2.

    Returns dict with eigenvalue extrema for both matrices.
    """
    c_proj = np.asarray(c_proj, dtype=np.float64)
    d_proj = np.asarray(d_proj, dtype=np.float64)
    M_f = bochner_matrix_at_f_star(c_proj, d_proj, n)
    c_1mf = -c_proj.copy()
    c_1mf[0] = 0.5
    d_1mf = -d_proj.copy()
    d_1mf[0] = 0.0
    M_1mf = bochner_matrix_at_f_star(c_1mf, d_1mf, n)
    eigs_f = np.linalg.eigvalsh(M_f)
    eigs_1mf = np.linalg.eigvalsh(M_1mf)
    return {
        "lambda_min_M_n(f)": float(eigs_f.min()),
        "lambda_max_M_n(f)": float(eigs_f.max()),
        "lambda_min_M_n(1-f)": float(eigs_1mf.min()),
        "lambda_max_M_n(1-f)": float(eigs_1mf.max()),
        "trace_M_n(f)": float(np.real(np.trace(M_f))),
        "trace_M_n(1-f)": float(np.real(np.trace(M_1mf))),
        "n": int(n),
    }


# --- Task 8: Polynomial-moment slack at f* ------------------------------


def poly_moments_direct(breakpoints, values, k_list):
    """Compute the EXACT physical moments m_k^phys = ∫_{-2}^{2} x^k f(x) dx
    for a step function f on [-2, 2].

    For a single cell [a, b] with constant value v:
        ∫_a^b x^k dx = (b^{k+1} - a^{k+1}) / (k+1)

    Parameters
    ----------
    breakpoints : array of length n_cells + 1
    values      : array of length n_cells
    k_list      : iterable of int (moments to compute)

    Returns
    -------
    dict {k: m_k^phys}
    """
    bp = np.asarray(breakpoints, dtype=np.float64)
    vv = np.asarray(values, dtype=np.float64)
    out = {}
    for k in k_list:
        # antiderivative at each breakpoint, then diff per cell
        anti = bp ** (k + 1) / (k + 1)
        m = float(np.sum(vv * (anti[1:] - anti[:-1])))
        out[int(k)] = m
    return out


def poly_moment_diagnostic(breakpoints, values, c_proj, d_proj, k_max=14, T=4000):
    """For each even k in {2, 4, ..., k_max}, compute:
      - m_k^phys = ∫_{-2}^{2} x^k f(x) dx  (closed-form, exact)
      - m_k^SDP_via_cd = poly_moment.py's expression evaluated at SDP's c, d
        (which equals the projection's c_proj[1:T+1], d_proj[1:T+1]).
      - tail_bound from poly_moment.even_moment_tail_bound(k, T).
      - relation check: m_k^SDP_via_cd should equal m_k^phys / 2^k
        (derivation in module banner).
      - slack from the SDP constraint  m_k^SDP_via_cd ≥ -tail_bound
        (positive ⇒ constraint satisfied with slack; negative ⇒ violated).

    Returns
    -------
    dict { k: {m_phys, m_SDP_via_cd, m_SDP_expected_from_phys, tail_bound,
               slack, slack_per_tail} }
    """
    # Lazy import to keep this module light if poly_moment unavailable.
    try:
        from lp_research_state.code.poly_moment import (
            fourier_coeffs_of_xk,
            even_moment_tail_bound,
        )
    except ImportError:  # pragma: no cover
        sys.path.insert(0, str(Path(__file__).parent))
        from poly_moment import (  # type: ignore
            fourier_coeffs_of_xk,
            even_moment_tail_bound,
        )

    if k_max % 2 != 0:
        k_max -= 1

    c_proj = np.asarray(c_proj, dtype=np.float64)
    d_proj = np.asarray(d_proj, dtype=np.float64)
    # SDP's c-indexed Fourier coeffs are proj[1..T] (length T):
    c_sdp = c_proj[1:T + 1]
    d_sdp = d_proj[1:T + 1]
    assert len(c_sdp) == T and len(d_sdp) == T, (
        f"projection length mismatch: c[1:T+1] has length {len(c_sdp)}, expected T={T}"
    )

    k_vals = list(range(2, k_max + 1, 2))
    alpha0, alpha, beta = fourier_coeffs_of_xk(k_max, T)
    phys = poly_moments_direct(breakpoints, values, k_list=k_vals)

    out = {}
    for k in k_vals:
        m_SDP = 0.5 * alpha0[k] + alpha[k, :] @ c_sdp + beta[k, :] @ d_sdp
        m_SDP = float(m_SDP)
        m_phys = float(phys[k])
        m_SDP_expected = m_phys / (2 ** k)
        tail = float(even_moment_tail_bound(k, T))
        slack = m_SDP - (-tail)  # m_SDP >= -tail; slack ≥ 0 means satisfied
        out[int(k)] = {
            "m_phys": m_phys,
            "m_SDP_via_cd": m_SDP,
            "m_SDP_expected_from_phys": m_SDP_expected,
            "convention_residual": float(m_SDP - m_SDP_expected),
            "tail_bound": tail,
            "slack": slack,
            "slack_per_tail": slack / tail if tail > 0 else float("inf"),
        }
    return out


# --- Task 9: Hankel-PSD slack at f* -------------------------------------


def hankel_diagnostic(breakpoints, values, c_proj, d_proj, n=6, T=4000):
    """At f*: build (a) the Hankel matrix H_n[i,j] = m_{i+j}^SDP_via_cd,
    (b) the interval-positivity matrix A_n[i,j] = m_{i+j} - m_{i+j+2},
    using the SAME SDP formula as hankel_probe.py / poly_moment.py.

    Reports eigenvalue extrema and the moments used. Also reports the
    physically-anchored moments m_k_phys = ∫_{-2}^{2} x^k f(x) dx so the
    reader can see both (the SDP's m_k corresponds to (1/2^k) · m_k_phys
    when conventions are coherent).

    Parameters
    ----------
    breakpoints : step function breakpoints (length n_cells + 1)
    values      : step function values     (length n_cells)
    c_proj, d_proj : projection arrays of length ≥ T + 1
    n           : Hankel level (matrix size n+1)
    T           : LP truncation level (default 4000, the SDP's level)

    Returns dict with lambda_min/max of H_n and A_n, plus the moments.
    """
    try:
        from lp_research_state.code.hankel_probe import (
            moments_from_cd, build_hankel, build_interval_pos,
            estimate_tail_bound,
        )
    except ImportError:  # pragma: no cover
        sys.path.insert(0, str(Path(__file__).parent))
        from hankel_probe import (  # type: ignore
            moments_from_cd, build_hankel, build_interval_pos,
            estimate_tail_bound,
        )

    c_proj = np.asarray(c_proj, dtype=np.float64)
    d_proj = np.asarray(d_proj, dtype=np.float64)
    c_sdp = c_proj[1:T + 1]
    d_sdp = d_proj[1:T + 1]
    assert len(c_sdp) == T, f"projection length mismatch: got {len(c_sdp)}, expected T={T}"

    k_max = 2 * n + 2  # need up to m_{2n+2} for A_n
    m_trunc = moments_from_cd(c_sdp, d_sdp, k_max=k_max)
    m_trunc = np.asarray(m_trunc, dtype=np.float64)
    H = build_hankel(m_trunc, n)
    A = build_interval_pos(m_trunc, n)
    eigs_H = np.linalg.eigvalsh(H)
    eigs_A = np.linalg.eigvalsh(A)

    # Physical moments for context.
    phys = poly_moments_direct(breakpoints, values, k_list=list(range(0, k_max + 1)))
    moments_table = {}
    for k in range(k_max + 1):
        tail = float(estimate_tail_bound(k, T)) if k > 0 else 0.0
        moments_table[int(k)] = {
            "m_SDP_via_cd": float(m_trunc[k]),
            "m_phys_over_2k": float(phys[k]) / (2 ** k),
            "tail_bound": tail,
        }

    return {
        "lambda_min_H_n": float(eigs_H.min()),
        "lambda_max_H_n": float(eigs_H.max()),
        "H_eigs": [float(x) for x in eigs_H],
        "lambda_min_A_n": float(eigs_A.min()),
        "lambda_max_A_n": float(eigs_A.max()),
        "A_eigs": [float(x) for x in eigs_A],
        "n": int(n),
        "moments": moments_table,
    }


if __name__ == "__main__":
    _run_all_tests()
    _test_bochner_xcheck()
    for kind in ("even", "direct"):
        c, d = project_together_f_star(T=4000, kind=kind)
        print(
            f"{kind}: c[0]={c[0]:.4f} d[0]={d[0]:.4f} "
            f"max|c[1:]|={np.max(np.abs(c[1:])):.4e} "
            f"max|d[1:]|={np.max(np.abs(d[1:])):.4e}"
        )
