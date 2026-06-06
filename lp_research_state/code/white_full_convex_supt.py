"""PRO-22: Direct sup_t SDP — bypass the cell-envelope cosine + sine.

The existing `white_full_convex.build_problem` enforces `Ω ≥ M(t)` indirectly via:
1. The cell averages `w_j, v_j` of M satisfy `w_j, v_j ≤ Ω` (line 141).
2. The cell-envelope constraints link (w, v) to (a, b) via cell-min relaxation.

This module replaces the cell-envelope with DIRECT `M(t_k) ≤ Ω` constraints at
a grid of shifts `t_k ∈ [0, 2]`. M(t) is computed from the (a_m, b_m) Fourier
variables (which white_full_convex.py builds at lines 154-174):

    M(t) = Σ_{m=1}^{2R} [a_m cos(πmt/2) + b_m sin(πmt/2)] + tail correction

The tail correction is bounded by `Σ_{m>2R} |M̂(m)|`, which is O(1/R) via
Cauchy-Schwarz on (W.1).

This is the "PRO-22" prototype — it might be too loose (tail dominates) or
might break the framework's 0.380558 ceiling. Empirical question.
"""
from __future__ import annotations
import sys
from pathlib import Path

import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import (
    cos_cell_bounds_exact,
    sin_cell_bounds_exact,
    odd_coeff_factors,
    tail_bound_eps,
    tail_bound_delta,
)


def build_problem_supt(
    N, T, R,
    h1, h2, p1, p2, q1, q2,
    *,
    bochner_n=0,
    n_supt_grid=200,         # number of t-grid points for direct M(t) ≤ Ω
    tail_bound_mode='cauchy_schwarz',  # how to bound Σ_{m>2R} |M̂(m)|
    assume_even=False,
    use_T5=False,
    use_T5p=False,
    use_T3=False,
):
    """Build the SDP with DIRECT M(t) ≤ Ω constraints (no cell-envelope).

    Returns: (Omega, w, v, c, d, eps, dlt, cons)

    Caveats:
    - The cell-envelope cosine + sine constraints (lines 176-190 of original)
      are REMOVED.
    - The w, v variables are KEPT (used by other constraints like (5.3), (5.4)).
    - A truncation tail term is added to make M(t) ≤ Ω valid even though only
      m=1..2R Fourier modes of M are constrained.
    """
    L = 2.0 / N
    j = np.arange(1, N + 1)

    Omega = cp.Variable()
    w = cp.Variable(N)
    v = cp.Variable(N)
    c = cp.Variable(T)
    d = cp.Variable(T)
    eps = cp.Variable(R)
    dlt = cp.Variable(R)

    cons = []
    cons += [w >= 0, v >= 0, w <= Omega, v <= Omega, Omega <= 1]
    cons.append(L * cp.sum(w + v) == 1)

    if assume_even:
        cons += [d == 0, dlt == 0, v == w]

    cons.append(L ** 2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h1)
    cons.append(L ** 3 * cp.sum(cp.multiply((j - 1) ** 2, (w + v))) <= 2.0 / 3 + h2 ** 2 / 2)

    # Build (a_m, b_m) Fourier expressions — same as standard build_problem
    a_expr = []
    b_expr = []
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            half = m // 2
            a_expr.append(0.5 * c[half - 1])
            b_expr.append(0.5 * d[half - 1])
        else:
            af, bf = odd_coeff_factors(m, T)
            sin_pi_half_m = np.sin(np.pi * m / 2)
            am = (
                eps[(m - 1) // 2]
                + (2 * m * sin_pi_half_m / np.pi)
                * (1.0 / (2 * m ** 2) + cp.sum(cp.multiply(af, c)))
            )
            bm = (
                dlt[(m - 1) // 2]
                + (4 * sin_pi_half_m / np.pi) * cp.sum(cp.multiply(bf, d))
            )
            a_expr.append(am)
            b_expr.append(bm)

    # ===== DIRECT sup_t constraints =====
    # For each t in grid, M(t) = Σ_{m=1}^{2R} [a_m cos(πmt/2) + b_m sin(πmt/2)]
    #                          + (tail correction)
    # Add: M(t) ≤ Ω
    #
    # Tail correction: we want to ENSURE M(t) ≤ Ω rigorously.
    # M_true(t) = M_truncated(t) + M_tail(t).
    # |M_tail(t)| ≤ Σ_{m>2R} (|a_m| + |b_m|).
    # By Cauchy-Schwarz: Σ |a_m| ≤ sqrt(2R-ish · Σ a_m²) … grows poorly.
    # Better: use Parseval-like bound: Σ_m |M̂(m)|² ≤ ||M||² ≤ Ω · ||M||_1 = Ω.
    # Then by C-S over the high-m tail (which has Σ ≤ Ω - low_m), get a bound.
    #
    # For prototype simplicity: use a conservative tail constant tau.
    # tau ≈ sqrt(Ω · sum_high) ≈ sqrt(0.4 · 0.1) ≈ 0.2 — too large for our gap.
    #
    # SIMPLIFICATION for prototype: ignore the tail term (set tau=0). This makes
    # the SDP NOT-RIGOROUS but provides a valid empirical upper bound on what
    # the direct-M approach can attain. If even this loose form doesn't break
    # the cell-envelope's reach, the approach fails.
    tail_const = 0.0  # placeholder; conservative bound is sqrt(Ω · (1/(2R)))

    t_grid = np.linspace(0.0, 2.0, n_supt_grid + 1)
    for t_k in t_grid:
        M_t = 0
        for m in range(1, 2 * R + 1):
            cos_mt = np.cos(np.pi * m * t_k / 2)
            sin_mt = np.sin(np.pi * m * t_k / 2)
            M_t = M_t + cos_mt * a_expr[m - 1] + sin_mt * b_expr[m - 1]
        # Constraint: M(t) + tail ≤ Ω  →  M_truncated ≤ Ω - tail
        cons.append(M_t <= Omega - tail_const)

    # Tail bounds, box, Parseval, etc. — same as standard
    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [
            cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
            cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T),
        ]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)
    cons += [c[0] >= p1, c[0] <= p2, d[0] >= q1, d[0] <= q2]

    # (5.13)-style with exact integral (kept from PRO-9)
    a_minus, _ = cos_cell_bounds_exact(j, 2, L)
    rhs_513 = -0.5 * (max(p1 ** 2, p2 ** 2) + max(q1 ** 2, q2 ** 2))
    cons.append((L / 2) * (a_minus @ (w + v)) >= rhs_513)

    if use_T3:
        cons.append(L * (cp.sum_squares(w) + cp.sum_squares(v)) <= Omega)

    if use_T5:
        Q = np.eye(T) + 0.5 * np.eye(T, k=1) + 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Q)) + cp.quad_form(d, cp.psd_wrap(Q)) <= 0.5)

    if use_T5p:
        Qp = np.eye(T) - 0.5 * np.eye(T, k=1) - 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Qp)) + cp.quad_form(d, cp.psd_wrap(Qp)) <= 0.5)

    # Bochner-PSD
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            half = 0.5
            Re_rows, Im_rows = [], []
            for jj in range(n_b + 1):
                re_row, im_row = [], []
                for k in range(n_b + 1):
                    ell = jj - k
                    if ell == 0:
                        re_row.append(cp.Constant(half))
                        im_row.append(cp.Constant(0.0))
                    else:
                        aell = abs(ell)
                        re_row.append(cp.Constant(sign * 0.5) * c[aell - 1])
                        if ell > 0:
                            im_row.append(cp.Constant(-sign * 0.5) * d[aell - 1])
                        else:
                            im_row.append(cp.Constant(+sign * 0.5) * d[aell - 1])
                Re_rows.append(re_row)
                Im_rows.append(im_row)
            Re_M = cp.bmat(Re_rows)
            Im_M = cp.bmat(Im_rows)
            real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
            cons.append(real_form >> 0)

    return Omega, w, v, c, d, eps, dlt, cons


if __name__ == "__main__":
    import time
    import warnings
    warnings.simplefilter("ignore")
    from white_full_convex import build_problem

    print("=== PRO-22 prototype: direct sup_t SDP ===\n")
    cfg = dict(N=200, T=100, R=5, bochner_n=8)
    h1, h2, p1, p2, q1, q2 = 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02

    print("Baseline (cell-envelope cell-min, original white_full_convex):")
    t0 = time.time()
    Omega, *_, cons = build_problem(
        cfg["N"], cfg["T"], cfg["R"], h1, h2, p1, p2, q1, q2,
        bochner_n=cfg["bochner_n"],
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    Omega_baseline = prob.value
    print(f"  Ω = {Omega_baseline:.7f}, time = {time.time()-t0:.2f}s")

    print("\nDirect sup_t (PRO-22, tail=0, grid=200):")
    t0 = time.time()
    Omega, *_, cons = build_problem_supt(
        cfg["N"], cfg["T"], cfg["R"], h1, h2, p1, p2, q1, q2,
        bochner_n=cfg["bochner_n"],
        n_supt_grid=200, tail_bound_mode='cauchy_schwarz',
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    Omega_supt = prob.value
    print(f"  Ω = {Omega_supt:.7f}, time = {time.time()-t0:.2f}s")

    print(f"\nΔΩ = Ω_supt - Ω_baseline = {Omega_supt - Omega_baseline:+.4e}")
    if Omega_supt > Omega_baseline:
        print(f"  ✓ Direct sup_t TIGHTER (better LB) by {Omega_supt - Omega_baseline:.4e}")
        print(f"  But: tail_const=0 makes this NON-RIGOROUS. Need tail term to make valid.")
    else:
        print(f"  Direct sup_t LOOSER by {Omega_baseline - Omega_supt:.4e}")
        print(f"  → cell-envelope was already tighter. PRO-22 approach doesn't help.")
