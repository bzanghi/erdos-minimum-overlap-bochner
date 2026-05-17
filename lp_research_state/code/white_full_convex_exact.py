"""ERD-9: Exact-integral cell-envelope variant of White's SDP.

The original `white_full_convex.build_problem` enforces the cosine cell-envelope
constraint (W.1) via a CELL-MIN relaxation:

    (L/2) · α_m^-(j) · (w + v)  +  2(a_m² + b_m²)  −  c_m · a_m  ≤  0

where α_m^-(j) = min over cell j of cos(πmx/2). This is an UNDER-estimate of
the cell integral, making the inequality LOOSER than the analytical form.

This module builds the EXACT version, replacing the cell-min by the cell
integral:

    (1/2) · I_m(j) · (w + v)   where   I_m(j) = ∫_{(j-1)L}^{jL} cos(πmx/2) dx
                                              = (2/(πm))·[sin(πmjL/2) − sin(πm(j-1)L/2)]

Similarly for the sine cell-envelope: the two-sided inequalities (line 189-190
of `white_full_convex.py`) collapse to a single equality at the exact integral.

The Step E saturation theorem (LEVER_I_PRIME_FINAL.md) predicts that
`Ω_exact ≤ C_explicit` for any SDP instance, with C_explicit = LB + ResidualGain.
This module is the empirical falsification test of that theorem.

Usage:
    from white_full_convex_exact import build_problem_exact
    Omega, *_, cons = build_problem_exact(N, T, R, h1, h2, p1, p2, q1, q2, ...)
    prob = cvxpy.Problem(cvxpy.Minimize(Omega), cons)
    prob.solve()
    Omega_exact = prob.value
"""
from __future__ import annotations

import sys
from pathlib import Path

import cvxpy as cp
import numpy as np

# Re-use everything from the standard module except the relaxed cell-envelope
sys.path.insert(0, str(Path(__file__).resolve().parent))
from white_full_convex import (  # noqa: E402
    cos_cell_bounds_exact,
    sin_cell_bounds_exact,
    odd_coeff_factors,
    tail_bound_eps,
    tail_bound_delta,
)


def cos_cell_integral(j_arr, m, L):
    """∫_{(j-1)L}^{jL} cos(πmx/2) dx  =  (2/(πm)) · [sin(πmjL/2) − sin(πm(j-1)L/2)].

    For positive-side cells. Negative-side integrals (over [-jL, -(j-1)L])
    are EQUAL to positive-side because cos is even.
    """
    x_lo = (j_arr - 1) * L
    x_hi = j_arr * L
    return (2.0 / (np.pi * m)) * (
        np.sin(np.pi * m * x_hi / 2) - np.sin(np.pi * m * x_lo / 2)
    )


def sin_cell_integral(j_arr, m, L):
    """∫_{(j-1)L}^{jL} sin(πmx/2) dx  =  -(2/(πm)) · [cos(πmjL/2) − cos(πm(j-1)L/2)].

    For positive-side cells. Negative-side integrals (over [-jL, -(j-1)L])
    are NEGATIVE of positive-side because sin is odd.
    """
    x_lo = (j_arr - 1) * L
    x_hi = j_arr * L
    return -(2.0 / (np.pi * m)) * (
        np.cos(np.pi * m * x_hi / 2) - np.cos(np.pi * m * x_lo / 2)
    )


def build_problem_exact(
    N, T, R,
    h1, h2, p1, p2, q1, q2,
    *,
    bochner_n=0,
    assume_even=False,
    use_T5=False,
    use_T5p=False,
    use_T3=False,
    cos_bnds=cos_cell_bounds_exact,
    sin_bnds=sin_cell_bounds_exact,
):
    """Build the SDP with EXACT cell-envelope integrals (no cell-min relaxation).

    All other constraints (Bochner-PSD, eps/dlt tail, Parseval, box, etc.)
    are identical to `white_full_convex.build_problem`.

    Returns: (Omega, w, v, c, d, eps, dlt, cons)
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

    # Fourier expression for (a_m, b_m) as in standard build_problem
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

    # ===== KEY DIFFERENCE: COSINE CELL-ENVELOPE with EXACT integral =====
    # Original (relaxed): (L/2) · α_m^-(j) · (w + v) + 2(a_m² + b_m²) − c_m·a_m ≤ 0
    # Exact:              (1/2) · I_m(j)   · (w + v) + 2(a_m² + b_m²) − c_m·a_m ≤ 0
    for m in range(1, 2 * R + 1):
        am = a_expr[m - 1]; bm = b_expr[m - 1]
        I_m = cos_cell_integral(j, m, L)
        lhs = (1.0 / 2) * (I_m @ (w + v))
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
        cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)

    # ===== KEY DIFFERENCE: SINE CELL-ENVELOPE as EQUALITY at exact integral =====
    # Original (relaxed, two-sided):
    #   (L/2) · (b_minus·w − b_plus·v) ≤ rhs
    #   (L/2) · (b_plus ·w − b_minus·v) ≥ rhs
    # Exact (collapsed to equality):
    #   (1/2) · J_m(j) · (w − v) = -(8/(πm))·sin(πm/2)·b_m
    # (since negative-side ∫ sin = -J_m by oddness of sin)
    for m in range(1, 2 * R + 1):
        bm = b_expr[m - 1]
        J_m = sin_cell_integral(j, m, L)
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs = -(8.0 / (m * np.pi)) * sin_pi_half_m * bm
        cons.append((1.0 / 2) * (J_m @ (w - v)) == rhs)

    # ===== Standard tail bounds, box constraints, etc. (unchanged) =====
    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [
            cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
            cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T),
        ]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)
    cons += [c[0] >= p1, c[0] <= p2, d[0] >= q1, d[0] <= q2]

    # (5.13): use EXACT integral for the m=2 cos cell-envelope reformulation
    # Original code uses a_plus_2 (cell-max); we use the exact integral I_2(j).
    I_2 = cos_cell_integral(j, 2, L)
    rhs_513 = -0.5 * (max(p1 ** 2, p2 ** 2) + max(q1 ** 2, q2 ** 2))
    cons.append((1.0 / 2) * (I_2 @ (w + v)) >= rhs_513)

    if use_T3:
        cons.append(L * (cp.sum_squares(w) + cp.sum_squares(v)) <= Omega)

    if use_T5:
        Q = np.eye(T) + 0.5 * np.eye(T, k=1) + 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Q)) + cp.quad_form(d, cp.psd_wrap(Q)) <= 0.5)

    if use_T5p:
        Qp = np.eye(T) - 0.5 * np.eye(T, k=1) - 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Qp)) + cp.quad_form(d, cp.psd_wrap(Qp)) <= 0.5)

    # Bochner-PSD (identical to standard)
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            Re_rows, Im_rows = [], []
            half = 0.5
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
    # Smoke test: solve at small N and compare to the relaxed version
    import time
    import warnings
    warnings.simplefilter("ignore")
    from white_full_convex import build_problem

    print("=== ERD-9 smoke test: exact vs relaxed cell-envelope ===")
    cfg = dict(N=200, T=100, R=5, bochner_n=8)
    h1, h2, p1, p2, q1, q2 = 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02

    print("\nRelaxed (cell-min) baseline:")
    t0 = time.time()
    Omega, *_ , cons = build_problem(
        cfg["N"], cfg["T"], cfg["R"], h1, h2, p1, p2, q1, q2,
        bochner_n=cfg["bochner_n"],
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    Omega_relaxed = prob.value
    print(f"  Omega_relaxed = {Omega_relaxed:.7f}, status = {prob.status}, time = {time.time()-t0:.2f}s")

    print("\nExact-integral version:")
    t0 = time.time()
    Omega, *_, cons = build_problem_exact(
        cfg["N"], cfg["T"], cfg["R"], h1, h2, p1, p2, q1, q2,
        bochner_n=cfg["bochner_n"],
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    Omega_exact = prob.value
    print(f"  Omega_exact   = {Omega_exact:.7f}, status = {prob.status}, time = {time.time()-t0:.2f}s")

    print(f"\nΔΩ = Omega_exact - Omega_relaxed = {Omega_exact - Omega_relaxed:+.4e}")
    print(f"(positive means exact version proves a stronger LB, consistent with Step E theorem)")
