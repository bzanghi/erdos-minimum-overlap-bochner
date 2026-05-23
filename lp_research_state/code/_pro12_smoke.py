"""PRO-12 smoke: row 4 at SMOKE scale — three-way CLARABEL vs Mosek vs SDPA-GMP.

The SDPA-GMP reference value at (N=200, T=80, R=10, bochner_n=4) is read from
lp_research_state/parallel_results/pro11_sdpa_s_serializer.json (the SMOKE row).
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
SDPA_REF = os.path.join(RESULTS_DIR, "pro11_sdpa_s_serializer.json")
OUT_JSON = os.path.join(RESULTS_DIR, "pro12_mosek_verify.json")


def build_row4(N, T, R, bn):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bn,
    )
    return cp.Problem(cp.Minimize(Omega), cons)


def smoke():
    N, T, R, bn = 200, 80, 10, 4
    print(f"=== PRO-12 SMOKE  N={N} T={T} R={R} bn={bn}  (row 4) ===")

    # ---- CLARABEL ----
    p1 = build_row4(N, T, R, bn)
    print("\n[1/3] CLARABEL ...")
    cl = solve_with_dual_extraction(p1)
    print(f"   reported      = {cl['reported_value']!r}")
    print(f"   rigorous LB   = {cl['rigorous_dual_LB']!r}")
    print(f"   status        = {cl['status']}   time = {cl['time']:.2f}s")

    # ---- Mosek ----
    p2 = build_row4(N, T, R, bn)
    print("\n[2/3] MOSEK ...")
    mk = solve_with_mosek(
        p2,
        mosek_params={
            "MSK_DPAR_INTPNT_CO_TOL_PFEAS":    1e-11,
            "MSK_DPAR_INTPNT_CO_TOL_DFEAS":    1e-11,
            "MSK_DPAR_INTPNT_CO_TOL_REL_GAP":  1e-13,
            "MSK_DPAR_INTPNT_CO_TOL_MU_RED":   1e-13,
        },
        keep_log=False,
    )
    print(f"   prob_value    = {mk['cvxpy_prob_value']!r}")
    print(f"   primal_obj    = {mk['primal_obj']!r}  (LB cert)")
    print(f"   dual_obj      = {mk['dual_obj']!r}    (UB cert)")
    print(f"   duality_gap   = {mk['duality_gap']!r}")
    print(f"   primal_viol   = {mk['primal_viol']!r}")
    print(f"   rigorous LB   = {mk['rigorous_dual_LB']!r}")
    print(f"   mosek_status  = {mk['mosek_problem_status']} / {mk['mosek_solution_status']}")
    print(f"   status        = {mk['status']}   iters = {mk['iterations']}   time = {mk['runtime_sec']:.2f}s")

    # ---- SDPA-GMP reference ----
    print("\n[3/3] SDPA-GMP reference (from pro11_sdpa_s_serializer.json) ...")
    with open(SDPA_REF) as fh:
        sdpa_runs = json.load(fh)
    sd = next(r for r in sdpa_runs if r.get("label") == "SMOKE")
    print(f"   primal_obj    = {sd['sdpa_primal_obj']!r}")
    print(f"   dual_obj      = {sd['sdpa_dual_obj']!r}")
    print(f"   duality_gap   = {sd['sdpa_duality_gap']!r}")
    print(f"   phase         = {sd['sdpa_phase']}   digits = {sd['sdpa_precision_digits']}")

    # ---- Comparison ----
    print("\n--- Three-way digit agreement (row 4 SMOKE) ---")
    refs = {
        "CLARABEL.reported":  cl["reported_value"],
        "Mosek.primal_obj":   mk["primal_obj"],
        "Mosek.dual_obj":     mk["dual_obj"],
        "Mosek.prob_value":   mk["cvxpy_prob_value"],
        "SDPA.primal_obj":    sd["sdpa_primal_obj"],
        "SDPA.dual_obj":      sd["sdpa_dual_obj"],
    }
    for k, v in refs.items():
        print(f"   {k:24} = {v!r}")
    pairs = [
        ("CLARABEL.reported", "SDPA.primal_obj"),
        ("CLARABEL.reported", "Mosek.dual_obj"),
        ("CLARABEL.reported", "Mosek.primal_obj"),
        ("Mosek.primal_obj",  "SDPA.primal_obj"),
        ("Mosek.dual_obj",    "SDPA.dual_obj"),
        ("Mosek.primal_obj",  "Mosek.dual_obj"),
    ]
    diffs = {}
    for a, b in pairs:
        d = abs(refs[a] - refs[b])
        diffs[f"{a} - {b}"] = d
        print(f"   |{a:22} - {b:22}| = {d:.3e}")

    out = {
        "label": "SMOKE",
        "N": N, "T": T, "R": R, "bochner_n": bn,
        "clarabel": {
            "reported_value": cl["reported_value"],
            "rigorous_dual_LB": cl["rigorous_dual_LB"],
            "status": cl["status"],
            "time": cl["time"],
        },
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
        },
        "sdpa_reference": {
            "primal_obj": sd["sdpa_primal_obj"],
            "dual_obj":   sd["sdpa_dual_obj"],
            "duality_gap": sd["sdpa_duality_gap"],
            "phase": sd["sdpa_phase"],
            "precision_digits": sd["sdpa_precision_digits"],
        },
        "pairwise_diffs": diffs,
    }
    return out


if __name__ == "__main__":
    res = smoke()
    # Append to or initialize the verification JSON
    existing = []
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as fh:
            existing = json.load(fh)
    # Replace any prior SMOKE row
    existing = [r for r in existing if r.get("label") != "SMOKE"]
    existing.append(res)
    with open(OUT_JSON, "w") as fh:
        json.dump(existing, fh, indent=2)
    print(f"\nWrote: {OUT_JSON}")
