"""
PRO-38 — Independently lift GATE region R16 by adding FRESH augmented dual-feasible
centers, so OUR Phi (not White's 0.38) certifies the floor over R16.

R16 box: h in [0,0.06], p in [0.33,0.45], q in [-0.025,-0.02] (narrow strip).
Current corehalo Phi floor = 0.38011094 (clears 0.380000 indep) but 1.731e-4 SHORT
of the core headline 0.380284. Worst corner found at (h=0.00375, p=0.3915, q=-0.025),
grid_min there 0.380259 dragged below target only by the Lipschitz eps_grid=1.49e-4.

METHOD (matches _verify_cover_dualext.py extraction EXACTLY):
  build_problem_with_dual_handles(N,T,R, h_c,h_c, p_c,p_c, q1,q2, bochner_n=BN)
  + build_even_moment_nonneg_constraints(H['c'],H['d'],T,k_max=PMK)
  solve_with_dual_extraction -> primal, dual_lb, dual_resid
  duals = {7 keys}.  Conservative anchor = primal - 1e-5.
  build_problem_with_dual_handles already bakes the corrected mside coeff (4.0;
  the 2026-05-31 White email correction, hardcoded at path_b_analytical.py L101-108).

A center is admitted as a VALID global LB only if the solve CONVERGED (small dual
residual) so its duals are a genuine dual-feasible point. We report dual_resid.

Centers are SINGLE-POINT (h1=h2, p1=p2) at the worst corner (and, optionally, a
couple spread along p), each at the worst |q| of the strip (q1=q2=-0.025).

Output: parallel_results/fullspace_promote_R16.json in the same center format as
the halo file (label,h_c,p_c,q1,q2,primal,dual_lb,gap,dual_resid,status,duals,time,config).
"""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction

PR = CODE.parent / "parallel_results"
OUT = PR / "fullspace_promote_R16.json"

DUAL_KEYS = ("con_53", "con_54", "con_512_pL", "con_512_pU",
             "con_512_qL", "con_512_qU", "con_513")


def solve_center(label, h_c, p_c, q1, q2, N, T, R, BN, PMK):
    """Replicate _verify_cover_dualext.py extraction exactly. Returns center dict."""
    t0 = time.time()
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=BN)
    pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=PMK)
    cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: float(H[k].dual_value) if H[k].dual_value is not None else 0.0
             for k in DUAL_KEYS}
    primal = res["reported_value"]
    dual_lb = res["rigorous_dual_LB"]
    gap = (primal - dual_lb) if (primal is not None and dual_lb is not None) else None
    rec = {
        "label": label, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "primal": primal, "dual_lb": dual_lb, "gap": gap,
        "dual_resid": res["dual_residual_at_LB"], "status": res["status"],
        "duals": duals, "time": time.time() - t0,
        "config": {"N": N, "T": T, "R": R, "bochner_n": BN, "pm_k_max": PMK},
    }
    return rec


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--T", type=int, default=2000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--BN", type=int, default=30)
    ap.add_argument("--PMK", type=int, default=20)
    # centers as semicolon-separated label,h,p,q1,q2
    ap.add_argument("--centers", type=str, required=True)
    args = ap.parse_args()

    existing = []
    if OUT.exists():
        existing = json.load(open(OUT)).get("centers", [])
    existing_labels = {c["label"] for c in existing}

    specs = []
    for tok in args.centers.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        label, h, p, q1, q2 = tok.split(",")
        specs.append((label, float(h), float(p), float(q1), float(q2)))

    results = list(existing)
    for (label, h_c, p_c, q1, q2) in specs:
        if label in existing_labels:
            print(f"[skip] {label} already present", flush=True)
            continue
        print(f"=== solving {label}: h={h_c} p={p_c} q=[{q1},{q2}] "
              f"N={args.N} T={args.T} BN={args.BN} PMK={args.PMK} ===", flush=True)
        rec = solve_center(label, h_c, p_c, q1, q2,
                           args.N, args.T, args.R, args.BN, args.PMK)
        print(f"  primal={rec['primal']:.8f}  dual_lb={rec['dual_lb']}  "
              f"gap={rec['gap']:.2e}  dual_resid={rec['dual_resid']:.2e}  "
              f"status={rec['status']}  ({rec['time']:.0f}s)", flush=True)
        results.append(rec)
        OUT.write_text(json.dumps(
            {"region": 16, "method": "fresh_augmented_centers_at_worst_corner",
             "anchor": "primal_m1e5", "centers": results}, indent=2, default=float))
        print(f"  saved -> {OUT.name}  ({len(results)} centers total)", flush=True)


if __name__ == "__main__":
    main()
