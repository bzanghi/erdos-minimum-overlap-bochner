"""
Path B (White-style ellipse-extension) with Bochner + Lasserre level-2 augmentation.

Same structure as path_b_analytical.py, but the inner SDP includes:
  - Bochner-on-f and Bochner-on-(1-f) at level n_b
  - Lasserre level-2 (M_top + Hermitian localizing matrix) at (T_max, T_loc).

The Lasserre/Bochner constraints depend only on (c, d), so they do NOT enter
the (h, p, q) sensitivity formula — Path B applies unchanged. We only need
the dual values of the (h, p, q)-dependent constraints (5.3, 5.4, 5.12, 5.13).

USAGE: see __main__ for a single-row CLI that writes JSON to disk.
"""
from __future__ import annotations
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp
import importlib.util as _ilu
from pathlib import Path as _Path
from white_full_convex import (
    cos_cell_bounds_exact, sin_cell_bounds_exact,
    odd_coeff_factors, tail_bound_eps, tail_bound_delta, WHITE_TABLE3,
)
from path_b_analytical import find_ellipse_h_p, in_ellipse


def _load_lasserre():
    here = _Path(__file__).resolve().parent
    spec = _ilu.spec_from_file_location("lasserre", here / "lasserre.py")
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_problem_with_dual_handles_BL(
    N, T, R, h1, h2, p1, p2, q1, q2,
    cell_mode="exact", bochner_n=0,
    lasserre_T_max=0, lasserre_T_loc=0,
):
    """Bochner+Lasserre-augmented program with named cvxpy handles for (h, p, q)
    constraints, used for Path B dual extraction."""
    L = 2.0 / N
    j = np.arange(1, N + 1)

    cos_bnds = cos_cell_bounds_exact
    sin_bnds = sin_cell_bounds_exact

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

    con_53 = L**2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h1
    cons.append(con_53)
    con_54 = L**3 * cp.sum(cp.multiply((j - 1)**2, (w + v))) <= 2.0/3 + h2**2 / 2
    cons.append(con_54)

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
            am = (eps[(m - 1) // 2]
                  + (2 * m * sin_pi_half_m / np.pi)
                  * (1.0 / (2 * m**2) + cp.sum(cp.multiply(af, c))))
            bm = (dlt[(m - 1) // 2]
                  + (4 * sin_pi_half_m / np.pi) * cp.sum(cp.multiply(bf, d)))
            a_expr.append(am)
            b_expr.append(bm)

    for m in range(1, 2 * R + 1):
        am = a_expr[m - 1]; bm = b_expr[m - 1]
        a_minus, _ = cos_bnds(j, m, L)
        lhs = (L / 2) * (a_minus @ (w + v))
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
        cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)

    for m in range(1, 2 * R + 1):
        bm = b_expr[m - 1]
        b_minus, b_plus = sin_bnds(j, m, L)
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

    con_512_pL = c[0] >= p1
    con_512_pU = c[0] <= p2
    con_512_qL = d[0] >= q1
    con_512_qU = d[0] <= q2
    cons += [con_512_pL, con_512_pU, con_512_qL, con_512_qU]

    _, a_plus_2 = cos_bnds(j, 2, L)
    rhs_513 = -0.5 * (max(p1**2, p2**2) + max(q1**2, q2**2))
    con_513 = (L / 2) * (a_plus_2 @ (w + v)) >= rhs_513
    cons.append(con_513)

    # Bochner moment-matrix PSD constraints
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

    # Lasserre level-2 augmentation
    if lasserre_T_max > 0:
        las = _load_lasserre()
        T_loc = lasserre_T_loc if lasserre_T_loc > 0 else lasserre_T_max
        las.add_lasserre2_constraint(cons, c, d, T_max=lasserre_T_max, T_loc=T_loc)

    handles = {
        "con_53": con_53,
        "con_54": con_54,
        "con_512_pL": con_512_pL,
        "con_512_pU": con_512_pU,
        "con_512_qL": con_512_qL,
        "con_512_qU": con_512_qU,
        "con_513": con_513,
        "Omega": Omega, "w": w, "v": v, "c": c, "d": d, "eps": eps, "dlt": dlt,
    }
    return Omega, cons, handles


def solve_and_extract_duals_BL(N, T, R, h, p, q1, q2, bochner_n,
                               lasserre_T_max, lasserre_T_loc, solver="CLARABEL"):
    Omega, cons, H = build_problem_with_dual_handles_BL(
        N, T, R, h, h, p, p, q1, q2, bochner_n=bochner_n,
        lasserre_T_max=lasserre_T_max, lasserre_T_loc=lasserre_T_loc,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=solver, verbose=False)
    elapsed = time.time() - t0

    duals = {}
    for key in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                "con_512_qL", "con_512_qU", "con_513"):
        dv = H[key].dual_value
        duals[key] = float(dv) if dv is not None else 0.0

    return {
        "value": float(prob.value) if prob.value is not None else None,
        "status": prob.status,
        "time": elapsed,
        "duals": duals,
        "h_c": h, "p_c": p, "q1": q1, "q2": q2,
        "N": N, "T": T, "R": R,
        "bochner_n": bochner_n,
        "lasserre_T_max": lasserre_T_max,
        "lasserre_T_loc": lasserre_T_loc,
    }


def run_one_row(label, h, p, qm, qp, N, T, R, bochner_n,
                lasserre_T_max, lasserre_T_loc, target=0.379005,
                out_dir=None):
    if out_dir is None:
        out_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/lasserre2_path_b"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{label}.json")
    print(f"=== {label}: h={h:.4f} p={p:.4f} q∈[{qm:.3f},{qp:.3f}] ===", flush=True)
    center = solve_and_extract_duals_BL(
        N, T, R, h, p, qm, qp, bochner_n, lasserre_T_max, lasserre_T_loc,
    )
    print(f"  primal={center['value']:.7f}  ({center['status']}, {center['time']:.1f}s)", flush=True)
    ell = find_ellipse_h_p(center, center['duals'], qm, qp, target=target)
    print(f"  ellipse: peak (h*={ell['h_star']:.4f}, p*={ell['p_star']:.4f}), "
          f"V_max={ell['V_max']:.6f}; semi (h={ell['semi_h']:.4f}, p={ell['semi_p']:.4f})", flush=True)

    out = {
        "label": label,
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n,
                   "lasserre_T_max": lasserre_T_max, "lasserre_T_loc": lasserre_T_loc},
        "h_c": center["h_c"], "p_c": center["p_c"],
        "q1": center["q1"], "q2": center["q2"],
        "primal_value_at_center": center["value"],
        "status": center["status"],
        "time_s": center["time"],
        "duals": center["duals"],
        "ellipse": {
            "semi_h": ell["semi_h"], "semi_p": ell["semi_p"],
            "h_star": ell["h_star"], "p_star": ell["p_star"], "V_max": ell["V_max"],
            "A_h2": ell["A_h2"], "A_h1": ell["A_h1"], "A_h0": ell["A_h0"],
            "A_p2": ell["A_p2"], "A_p1": ell["A_p1"], "A_p0": ell["A_p0"],
            "const_q": ell["const_q"], "V_c": ell["V_c"], "target": ell["target"],
        },
    }
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  written to {out_file}", flush=True)
    return out


def aggregate_and_verify(in_dir=None, out_file=None, target=0.379005,
                         h_range=(0.0, 0.06), p_range=(0.35, 0.45),
                         q_range=(-0.02, 0.02), n_grid=81):
    if in_dir is None:
        in_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/lasserre2_path_b"
    if out_file is None:
        out_file = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/lasserre2_path_b_summary.json"

    rows_data = []
    for label in [f"row{i}" for i in range(1, 8)]:
        f = os.path.join(in_dir, f"{label}.json")
        if not os.path.exists(f):
            print(f"  missing {f}; skipping")
            continue
        with open(f) as fh:
            rows_data.append(json.load(fh))
    print(f"  Loaded {len(rows_data)} rows.")

    qm_test, qp_test = q_range
    h_grid = np.linspace(h_range[0], h_range[1], n_grid)
    p_grid = np.linspace(p_range[0], p_range[1], n_grid)
    min_lb = np.inf
    min_loc = None
    min_row = None
    grid_uncovered = []
    for hr in h_grid:
        for pr in p_grid:
            best_obj = -np.inf
            best_row_label = None
            for r in rows_data:
                ell = r['ellipse']
                if r['q1'] <= qm_test and r['q2'] >= qp_test:
                    val = (ell["V_c"] + ell["const_q"]
                           + ell["A_h2"]*hr**2 + ell["A_h1"]*hr + ell["A_h0"]
                           + ell["A_p2"]*pr**2 + ell["A_p1"]*pr + ell["A_p0"])
                    if val > best_obj:
                        best_obj = val
                        best_row_label = r['label']
            if best_obj < min_lb:
                min_lb = best_obj
                min_loc = (float(hr), float(pr))
                min_row = best_row_label
            if best_obj < target and len(grid_uncovered) < 20:
                grid_uncovered.append({"h": float(hr), "p": float(pr),
                                       "best_obj": float(best_obj),
                                       "best_row": best_row_label})

    print(f"  GRID MIN dual obj over (5.16) ({n_grid}x{n_grid}): "
          f"{min_lb:.7f} at (h={min_loc[0]:.4f}, p={min_loc[1]:.4f}, best row={min_row})")
    print(f"  Improvement vs White's 0.379005: {min_lb - 0.379005:+.6e}")

    out = {
        "config": {"target": target, "qm_test": qm_test, "qp_test": qp_test,
                   "n_grid": n_grid, "h_range": list(h_range), "p_range": list(p_range)},
        "rows": rows_data,
        "coverage": {
            "grid_min_obj": float(min_lb),
            "grid_min_loc_h": min_loc[0],
            "grid_min_loc_p": min_loc[1],
            "grid_min_row": min_row,
            "improvement_vs_0p379005": float(min_lb - 0.379005),
            "grid_uncovered_examples": grid_uncovered,
            "grid_uncovered_count": len(grid_uncovered),
        },
    }
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  written to {out_file}")
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    if len(sys.argv) >= 2 and sys.argv[1] == "row":
        # CLI: row <label> <N> <T> <R> <bochner_n> <T_max> <T_loc>
        label = sys.argv[2]
        N = int(sys.argv[3])
        T = int(sys.argv[4])
        R = int(sys.argv[5])
        bn = int(sys.argv[6])
        Tmx = int(sys.argv[7])
        Tlc = int(sys.argv[8])
        for (h, p, qm, qp, lbl) in WHITE_TABLE3:
            if lbl == label:
                run_one_row(label, h, p, qm, qp, N, T, R, bn, Tmx, Tlc)
                break
    elif len(sys.argv) >= 2 and sys.argv[1] == "aggregate":
        aggregate_and_verify()
    else:
        print("usage: row <label> <N> <T> <R> <bochner_n> <T_max> <T_loc>  |  aggregate")
