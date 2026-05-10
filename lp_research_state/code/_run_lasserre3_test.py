"""Standalone test runner for Lasserre level-3 at row 4.

Compares:
  (B) Bochner-on-f  alone, level n_b=20
  (C) Bochner + Lasserre-2  (T_max=10, T_loc=10)
  (D) Bochner + Lasserre-2 + Lasserre-3 trilinear lift (T3_lift=1, M3_max=2 or 3)

Adds level-3 by directly attaching it to the cvxpy constraint list AFTER
white_full_convex.build_problem returns it.

NOTE: build_problem already adds the level-2 lift if lasserre_T_max>0, but
it doesn't expose M_top. To get M_top for level-3, we instead build the
program manually with both levels added.
"""
from __future__ import annotations
import sys, os, json, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp

from white_full_convex import (
    cos_cell_bounds_exact, sin_cell_bounds_exact,
    odd_coeff_factors, tail_bound_eps, tail_bound_delta,
)
import lasserre as L2
import lasserre3_toeplitz as L3


def build_program(N, T, R, h, p, q1, q2, *,
                   bochner_n=0, lasserre_T_max=0, lasserre_T_loc=0,
                   l3_T3_lift=0, l3_M3_max=0):
    """Hand-build the white program + Bochner + Lasserre-2 + Lasserre-3."""
    L = 2.0 / N
    j_arr = np.arange(1, N + 1)

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
    cons.append(L ** 2 * cp.sum(cp.multiply(j_arr, w) - cp.multiply(j_arr - 1, v)) >= h)
    cons.append(L ** 3 * cp.sum(cp.multiply((j_arr - 1) ** 2, (w + v))) <= 2.0 / 3 + h ** 2 / 2)

    a_expr, b_expr = [], []
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            half = m // 2
            a_expr.append(0.5 * c[half - 1])
            b_expr.append(0.5 * d[half - 1])
        else:
            af, bf = odd_coeff_factors(m, T)
            sin_pi_half_m = np.sin(np.pi * m / 2)
            am = (eps[(m - 1) // 2]
                  + (2 * m * sin_pi_half_m / np.pi) * (1.0 / (2 * m ** 2) + cp.sum(cp.multiply(af, c))))
            bm = (dlt[(m - 1) // 2]
                  + (4 * sin_pi_half_m / np.pi) * cp.sum(cp.multiply(bf, d)))
            a_expr.append(am)
            b_expr.append(bm)

    for m in range(1, 2 * R + 1):
        am = a_expr[m - 1]; bm = b_expr[m - 1]
        a_minus, _ = cos_bnds_(j_arr, m, L)
        lhs = (L / 2) * (a_minus @ (w + v))
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
        cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)

    for m in range(1, 2 * R + 1):
        bm = b_expr[m - 1]
        b_minus, b_plus = sin_bnds_(j_arr, m, L)
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs = -(8.0 / (m * np.pi)) * sin_pi_half_m * bm
        cons.append((L / 2) * (b_minus @ w - b_plus @ v) <= rhs)
        cons.append((L / 2) * (b_plus @ w - b_minus @ v) >= rhs)

    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
                 cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T)]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)
    cons += [c[0] >= p, c[0] <= p, d[0] >= q1, d[0] <= q2]

    _, a_plus_2 = cos_bnds_(j_arr, 2, L)
    rhs_513 = -0.5 * (p ** 2 + max(q1 ** 2, q2 ** 2))
    cons.append((L / 2) * (a_plus_2 @ (w + v)) >= rhs_513)

    # Bochner-on-f and Bochner-on-(1-f).
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            half = 0.5
            Re_rows, Im_rows = [], []
            for jj in range(n_b + 1):
                re_row, im_row = [], []
                for kk in range(n_b + 1):
                    ell = jj - kk
                    if ell == 0:
                        re_row.append(cp.Constant(half)); im_row.append(cp.Constant(0.0))
                    else:
                        aell = abs(ell)
                        re_row.append(cp.Constant(sign * 0.5) * c[aell - 1])
                        if ell > 0:
                            im_row.append(cp.Constant(-sign * 0.5) * d[aell - 1])
                        else:
                            im_row.append(cp.Constant(+sign * 0.5) * d[aell - 1])
                Re_rows.append(re_row); Im_rows.append(im_row)
            Re_M = cp.bmat(Re_rows); Im_M = cp.bmat(Im_rows)
            real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
            cons.append(real_form >> 0)

    M_top = None
    # Lasserre-2 (lift)
    if lasserre_T_max > 0:
        T_loc = lasserre_T_loc if lasserre_T_loc > 0 else lasserre_T_max
        M_top = L2.add_lasserre2_constraint(cons, c, d, T_max=lasserre_T_max, T_loc=T_loc)

    # Lasserre-3 (trilinear/quadrilinear lift)
    if l3_T3_lift > 0:
        if M_top is None:
            raise ValueError("Lasserre-3 requires lasserre_T_max > 0.")
        T_loc3 = lasserre_T_loc if lasserre_T_loc > 0 else lasserre_T_max
        L3.add_lasserre3_toeplitz_constraint(
            cons, c, d, M_top, T_max=lasserre_T_max,
            T_loc3=T_loc3, T3_lift=l3_T3_lift,
            M3_max=l3_M3_max if l3_M3_max > 0 else None,
        )

    return Omega, cons


def cos_bnds_(j_arr, m, L):
    return cos_cell_bounds_exact(j_arr, m, L)


def sin_bnds_(j_arr, m, L):
    return sin_cell_bounds_exact(j_arr, m, L)


def run(label, **kwargs):
    t0 = time.time()
    Omega, cons = build_program(**kwargs)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    dt = time.time() - t0
    print(f"  [{label:50s}] Ω* = {prob.value:.7f}  ({prob.status}, {dt:.1f}s)")
    return {"label": label, "value": float(prob.value) if prob.value is not None else None,
            "status": prob.status, "time_s": dt}


if __name__ == "__main__":
    H = 0.004; P = 0.3875; QM, QP = -0.02, 0.02

    # Choose problem size: start small then scale up if time allows.
    SIZE = os.environ.get("L3_SIZE", "small")  # "small" | "med" | "full"
    if SIZE == "small":
        N, T, R = 500, 50, 10
        L2_T_max = 5; L2_T_loc = 5
    elif SIZE == "med":
        N, T, R = 1000, 100, 10
        L2_T_max = 10; L2_T_loc = 10
    else:
        N, T, R = 2000, 200, 10
        L2_T_max = 10; L2_T_loc = 10
    BOCH_N = 20

    print(f"=== Lasserre-3 row 4 test (SIZE={SIZE}) ===")
    print(f"  N={N}, T={T}, R={R}; L2: T_max={L2_T_max}, T_loc={L2_T_loc}; BOCH_N={BOCH_N}")
    print(f"  h={H}, p={P}, q∈[{QM},{QP}]")
    print()

    common = dict(N=N, T=T, R=R, h=H, p=P, q1=QM, q2=QP)
    results = []

    print("(A) Baseline (no extras):")
    res_A = run("baseline", **common); results.append(res_A)

    print("(B) Bochner-on-f alone (n_b={}):".format(BOCH_N))
    res_B = run("bochner", bochner_n=BOCH_N, **common); results.append(res_B)

    print("(C) Bochner + Lasserre-2:")
    res_C = run("bochner+L2", bochner_n=BOCH_N, lasserre_T_max=L2_T_max, lasserre_T_loc=L2_T_loc, **common)
    results.append(res_C)

    print("(D1) Bochner + L2 + L3 (T3_lift=1, M3_max=2):")
    res_D1 = run("bochner+L2+L3(T3=1, M3=2)",
                 bochner_n=BOCH_N, lasserre_T_max=L2_T_max, lasserre_T_loc=L2_T_loc,
                 l3_T3_lift=1, l3_M3_max=2, **common)
    results.append(res_D1)

    print("(D2) Bochner + L2 + L3 (T3_lift=1, M3_max=3):")
    res_D2 = run("bochner+L2+L3(T3=1, M3=3)",
                 bochner_n=BOCH_N, lasserre_T_max=L2_T_max, lasserre_T_loc=L2_T_loc,
                 l3_T3_lift=1, l3_M3_max=3, **common)
    results.append(res_D2)

    print("(D3) Bochner + L2 + L3 (T3_lift=2, M3_max=3):")
    res_D3 = run("bochner+L2+L3(T3=2, M3=3)",
                 bochner_n=BOCH_N, lasserre_T_max=L2_T_max, lasserre_T_loc=L2_T_loc,
                 l3_T3_lift=2, l3_M3_max=3, **common)
    results.append(res_D3)

    print()
    print("=== Summary deltas ===")
    if res_C["value"] is not None and res_B["value"] is not None:
        print(f"  Δ(C - B): L2 over Bochner    = {res_C['value'] - res_B['value']:+.7e}")
    for res_D in [res_D1, res_D2, res_D3]:
        if res_D["value"] is not None and res_C["value"] is not None:
            print(f"  Δ({res_D['label']:30s} - C) = {res_D['value'] - res_C['value']:+.7e}")

    out = {
        "config": {"N": N, "T": T, "R": R, "L2_T_max": L2_T_max, "L2_T_loc": L2_T_loc, "BOCH_N": BOCH_N,
                   "h": H, "p": P, "q1": QM, "q2": QP, "size": SIZE},
        "results": results,
    }
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "parallel_results", "lasserre3.json")
    out_path = os.path.normpath(out_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nWrote {out_path}")
