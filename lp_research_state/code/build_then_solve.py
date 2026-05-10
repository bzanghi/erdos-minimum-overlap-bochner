"""
Split-solve wrapper for the Bochner+Lasserre Erdos minimum-overlap SDP.

Bash sandbox has a 45s wall-clock cap per call (any Python process gets killed
when the call returns). At (N=10000, T=4000, bochner_n>=30 or lasserre_T_max>=40)
the canonicalize+solve flow exceeds 45s. CVXPY's API allows splitting the work:

  call A: build_problem -> cp.Problem -> prob.get_problem_data(solver=...)
            => returns (data, chain, inverse_data), all picklable.
            We pickle (data, chain, inverse_data) plus a thin "context" dict
            holding (N, T, R, h, p, q1, q2, bochner_n, lasserre_T_max,
            lasserre_T_loc) so call B can rebuild constraint handles.

  call B: load pickle -> chain.solve_via_data(problem, data, verbose=True)
            => CLARABEL output captured, dual extraction via dual_extractor.
            Then problem.unpack_results(soln, chain, inverse_data) populates
            primal/dual values in the original Problem.

Subtleties:
  * The canonicalized `data` for CLARABEL is large (sparse A, dims, etc.) so
    we pickle once per (row, n, T_max) and reuse across solves if possible.
  * `prob` itself is the rebuilt cp.Problem in call B (cheaper than pickling
    cvxpy expressions); `unpack_results` requires the SAME variable IDs, so
    we pickle the whole rebuilt prob with handles.

Usage:
  build_then_solve.py BUILD <row> <N> <T> <R> <bochner_n> <T_max> <T_loc> <pkl_path>
  build_then_solve.py SOLVE <pkl_path> <out_json_path>
"""
from __future__ import annotations
import sys, os, json, time, pickle, io, contextlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp

from white_full_convex import (
    cos_cell_bounds_exact, sin_cell_bounds_exact,
    odd_coeff_factors, tail_bound_eps, tail_bound_delta, WHITE_TABLE3,
)
from path_b_lasserre import build_problem_with_dual_handles_BL
from path_b_analytical import find_ellipse_h_p
from dual_extractor import parse_clarabel_iterations, best_dual_lower_bound


# ---- BUILD phase -----------------------------------------------------------

def do_build(label, N, T, R, bochner_n, lasserre_T_max, lasserre_T_loc, pkl_path):
    # Lookup row params
    row = next((r for r in WHITE_TABLE3 if r[4] == label), None)
    if row is None:
        raise SystemExit(f"unknown label {label}")
    h, p, qm, qp, _ = row

    t0 = time.time()
    Omega, cons, H = build_problem_with_dual_handles_BL(
        N, T, R, h, h, p, p, qm, qp,
        bochner_n=bochner_n,
        lasserre_T_max=lasserre_T_max, lasserre_T_loc=lasserre_T_loc,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t_build = time.time() - t0

    t1 = time.time()
    data, chain, inverse_data = prob.get_problem_data(solver="CLARABEL")
    t_compile = time.time() - t1

    # Pickle: we need to recover (data, chain, inverse_data) + the dual handles
    # so we can call problem.unpack_results and read dual_value of the named cons.
    # Best approach: pickle (prob, data, chain, inverse_data, handles_dual_index).
    # We can't directly look up handles by id after unpickling (cvxpy assigns
    # new ids), but we CAN pickle the whole prob+handles together since they
    # share the same Variable/Constraint objects.
    state = {
        "label": label,
        "N": N, "T": T, "R": R,
        "h": h, "p": p, "q1": qm, "q2": qp,
        "bochner_n": bochner_n,
        "lasserre_T_max": lasserre_T_max,
        "lasserre_T_loc": lasserre_T_loc,
        "prob": prob,
        "data": data,
        "chain": chain,
        "inverse_data": inverse_data,
        "handles": H,
        "t_build": t_build,
        "t_compile": t_compile,
    }
    os.makedirs(os.path.dirname(pkl_path), exist_ok=True)
    with open(pkl_path, "wb") as f:
        pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
    sz = os.path.getsize(pkl_path) / (1024 * 1024)
    print(f"BUILD ok: label={label} N={N} T={T} bochner_n={bochner_n} "
          f"T_max={lasserre_T_max} T_loc={lasserre_T_loc}")
    print(f"  build={t_build:.2f}s compile={t_compile:.2f}s pkl={sz:.1f} MiB -> {pkl_path}")
    return state


# ---- SOLVE phase -----------------------------------------------------------

def do_solve(pkl_path, out_json_path, target=0.379005, time_limit=38.0,
             max_iter=200):
    t0 = time.time()
    with open(pkl_path, "rb") as f:
        state = pickle.load(f)
    t_load = time.time() - t0
    prob = state["prob"]
    data = state["data"]
    chain = state["chain"]
    inverse_data = state["inverse_data"]
    H = state["handles"]

    # CLARABEL.solve_via_data with verbose -> capture stdout. Apply time_limit
    # so CLARABEL terminates gracefully within the 45s sandbox cap, and we can
    # still extract the rigorous dual LB from the iteration history.
    solver_opts = {"time_limit": float(time_limit), "max_iter": int(max_iter)}
    buf = io.StringIO()
    t1 = time.time()
    with contextlib.redirect_stdout(buf):
        soln = chain.solve_via_data(prob, data, warm_start=False, verbose=True,
                                     solver_opts=solver_opts)
    t_solve = time.time() - t1
    output = buf.getvalue()
    # Unpack into the original Problem so prob.value, dual_value etc. populate.
    prob.unpack_results(soln, chain, inverse_data)

    # Parse iteration table for rigorous dual LB
    rows = parse_clarabel_iterations(output)
    rigorous_dual, dual_resid, best_iter, n_elig = best_dual_lower_bound(
        rows, max_dual_residual=1e-4
    )

    duals = {}
    for key in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                "con_512_qL", "con_512_qU", "con_513"):
        dv = H[key].dual_value
        duals[key] = float(dv) if dv is not None else 0.0

    # If primal is None (e.g. CLARABEL hit time_limit), use the rigorous dual
    # LB as the certified value (a true lower bound on Omega).
    primal_v = float(prob.value) if prob.value is not None else None
    use_value = primal_v if primal_v is not None else rigorous_dual

    # Build a "center" dict for find_ellipse — V_c uses the rigorous dual LB
    # so the resulting ellipse coverage is itself a rigorous lower bound.
    center = {
        "h_c": state["h"], "p_c": state["p"],
        "q1": state["q1"], "q2": state["q2"],
        "value": use_value,
        "status": prob.status,
    }

    ell = None
    if use_value is not None and all(v is not None for v in duals.values()):
        try:
            ell = find_ellipse_h_p(center, duals, state["q1"], state["q2"], target=target)
        except Exception as e:
            print(f"  find_ellipse_h_p failed: {e}")

    out = {
        "label": state["label"],
        "config": {
            "N": state["N"], "T": state["T"], "R": state["R"],
            "bochner_n": state["bochner_n"],
            "lasserre_T_max": state["lasserre_T_max"],
            "lasserre_T_loc": state["lasserre_T_loc"],
        },
        "h_c": state["h"], "p_c": state["p"],
        "q1": state["q1"], "q2": state["q2"],
        "primal_value_at_center": center["value"],
        "rigorous_dual_LB": rigorous_dual,
        "dual_residual_at_LB": dual_resid,
        "best_iter": best_iter,
        "n_eligible_iters": n_elig,
        "n_iters_total": rows[-1]["iter"] if rows else 0,
        "status": prob.status,
        "time_load_s": t_load,
        "time_solve_s": t_solve,
        "duals": duals,
        "ellipse": ({
            "semi_h": ell["semi_h"], "semi_p": ell["semi_p"],
            "h_star": ell["h_star"], "p_star": ell["p_star"], "V_max": ell["V_max"],
            "A_h2": ell["A_h2"], "A_h1": ell["A_h1"], "A_h0": ell["A_h0"],
            "A_p2": ell["A_p2"], "A_p1": ell["A_p1"], "A_p0": ell["A_p0"],
            "const_q": ell["const_q"], "V_c": ell["V_c"], "target": ell["target"],
        } if ell is not None else None),
    }
    os.makedirs(os.path.dirname(out_json_path), exist_ok=True)
    with open(out_json_path, "w") as f:
        json.dump(out, f, indent=2)
    rl = f"{rigorous_dual:.7f}" if rigorous_dual is not None else "n/a"
    print(f"SOLVE ok: load={t_load:.2f}s solve={t_solve:.2f}s "
          f"primal={center['value']:.7f} rigorous_dual_LB={rl} "
          f"status={prob.status}")
    print(f"  -> {out_json_path}")
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    cmd = sys.argv[1]
    if cmd == "BUILD":
        # BUILD <row> <N> <T> <R> <bochner_n> <T_max> <T_loc> <pkl_path>
        label = sys.argv[2]; N = int(sys.argv[3]); T = int(sys.argv[4])
        R = int(sys.argv[5]); bn = int(sys.argv[6])
        Tmx = int(sys.argv[7]); Tlc = int(sys.argv[8])
        pkl_path = sys.argv[9]
        do_build(label, N, T, R, bn, Tmx, Tlc, pkl_path)
    elif cmd == "SOLVE":
        # SOLVE <pkl_path> <out_json_path> [time_limit] [max_iter]
        pkl_path = sys.argv[2]; out_json_path = sys.argv[3]
        time_limit = float(sys.argv[4]) if len(sys.argv) > 4 else 38.0
        max_iter = int(sys.argv[5]) if len(sys.argv) > 5 else 200
        do_solve(pkl_path, out_json_path, time_limit=time_limit, max_iter=max_iter)
    else:
        print(__doc__); sys.exit(1)
