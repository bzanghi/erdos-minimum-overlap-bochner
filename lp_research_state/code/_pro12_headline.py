"""PRO-12 headline: Mosek at Phase-5 production scale (N=10000, T=4000, bn=20),
row 4 (binding row). Compare to CLARABEL at the same params.
"""
from __future__ import annotations
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cvxpy as cp

from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
from mosek_runner import solve_with_mosek


HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.realpath(os.path.join(HERE, "..", "parallel_results"))
OUT_JSON = os.path.join(RESULTS_DIR, "pro12_mosek_verify.json")


def build_row4(N, T, R, bn):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bn,
    )
    return cp.Problem(cp.Minimize(Omega), cons)


def headline(timeout_sec=3600.0, also_clarabel=True):
    # Phase 5 production headline params per CLAUDE.md reproduce snippet
    N, T, R, bn = 10000, 4000, 10, 20
    print(f"=== PRO-12 HEADLINE  N={N} T={T} R={R} bn={bn}  (row 4) ===")
    print(f"timeout_sec = {timeout_sec}")

    cl_block = None
    if also_clarabel:
        print("\n[1/2] CLARABEL (baseline) ...")
        t0 = time.time()
        p1 = build_row4(N, T, R, bn)
        t_build = time.time() - t0
        print(f"  built in {t_build:.1f}s")
        try:
            cl = solve_with_dual_extraction(p1)
            print(f"  reported      = {cl['reported_value']!r}")
            print(f"  rigorous LB   = {cl['rigorous_dual_LB']!r}")
            print(f"  status        = {cl['status']}")
            print(f"  iters total   = {cl['n_iters_total']}")
            print(f"  time          = {cl['time']:.1f}s")
            cl_block = {
                "reported_value": cl["reported_value"],
                "rigorous_dual_LB": cl["rigorous_dual_LB"],
                "dual_residual_at_LB": cl["dual_residual_at_LB"],
                "best_iter": cl["best_iter"],
                "n_iters_total": cl["n_iters_total"],
                "status": cl["status"],
                "time": cl["time"],
            }
        except Exception as e:
            print(f"  CLARABEL failed: {e!r}")
            cl_block = {"error": repr(e)}

    print("\n[2/2] MOSEK ...")
    t0 = time.time()
    p2 = build_row4(N, T, R, bn)
    t_build = time.time() - t0
    print(f"  built in {t_build:.1f}s")
    mk = solve_with_mosek(
        p2,
        timeout_sec=timeout_sec,
        mosek_params={
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS":    1e-10,
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS":    1e-10,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP":  1e-12,
            "MSK_DPAR_INTPNT_CO_TOL_MU_RED":   1e-12,
            "MSK_IPAR_NUM_THREADS":            0,
            "MSK_IPAR_INTPNT_MAX_ITERATIONS":  400,
        },
        keep_log=True,
        log_max_chars=120_000,
    )
    print(f"  prob_value    = {mk['cvxpy_prob_value']!r}")
    print(f"  primal_obj    = {mk['primal_obj']!r}  (LB cert)")
    print(f"  dual_obj      = {mk['dual_obj']!r}    (UB cert)")
    print(f"  duality_gap   = {mk['duality_gap']!r}")
    print(f"  primal_viol   = {mk['primal_viol']!r}")
    print(f"  dual_viol     = {mk['dual_viol']!r}")
    print(f"  rigorous LB   = {mk['rigorous_dual_LB']!r}")
    print(f"  mosek status  = {mk['mosek_problem_status']} / {mk['mosek_solution_status']}")
    print(f"  cvxpy status  = {mk['status']}")
    print(f"  iterations    = {mk['iterations']}")
    print(f"  runtime       = {mk['runtime_sec']:.1f}s")
    if mk.get("error"):
        print(f"  ERROR         = {mk['error']!r}")

    # Comparison
    if cl_block and cl_block.get("reported_value") is not None and mk["primal_obj"] is not None:
        print("\n--- Mosek vs CLARABEL (row 4 headline) ---")
        d1 = mk["dual_obj"] - cl_block["reported_value"]
        d2 = mk["primal_obj"] - (cl_block["rigorous_dual_LB"] or cl_block["reported_value"])
        print(f"  Mosek.dual_obj - CLARABEL.reported = {d1:+.3e}  (encoding consistency)")
        print(f"  Mosek.primal_obj - CLARABEL.rigLB  = {d2:+.3e}  (LB margin upgrade)")

    out = {
        "label": "HEADLINE_phase5_row4",
        "N": N, "T": T, "R": R, "bochner_n": bn,
        "centers": {"h": 0.004, "p": 0.3875, "qm": -0.02, "qp": 0.02},
        "clarabel": cl_block,
        "mosek": {
            "primal_obj": mk["primal_obj"],
            "dual_obj":   mk["dual_obj"],
            "duality_gap": mk["duality_gap"],
            "primal_viol": mk["primal_viol"],
            "dual_viol":   mk["dual_viol"],
            "rigorous_dual_LB": mk["rigorous_dual_LB"],
            "cvxpy_prob_value": mk["cvxpy_prob_value"],
            "status": mk["status"],
            "mosek_problem_status": mk["mosek_problem_status"],
            "mosek_solution_status": mk["mosek_solution_status"],
            "iterations": mk["iterations"],
            "runtime_sec": mk["runtime_sec"],
            "error": mk.get("error"),
        },
    }
    return out


if __name__ == "__main__":
    timeout = float(os.environ.get("PRO12_TIMEOUT_SEC", "3600"))
    skip_clarabel = os.environ.get("PRO12_SKIP_CLARABEL", "0") == "1"
    res = headline(timeout_sec=timeout, also_clarabel=not skip_clarabel)
    existing = []
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as fh:
            existing = json.load(fh)
    existing = [r for r in existing if r.get("label") != "HEADLINE_phase5_row4"]
    existing.append(res)
    with open(OUT_JSON, "w") as fh:
        json.dump(existing, fh, indent=2)
    print(f"\nWrote: {OUT_JSON}")
