"""
PRO-38 R17 — FRESH AUGMENTED Phi-CENTERS deliverable (task-specified mechanism).

Region 17: h_range=[0,0.06], p_range=[0.33,0.45], q_range=[0.02,0.025]
(width_class=narrow). With the 12 core centers only, ours_Phi_min over R17 is
~0.3801482 (clears 0.380000 but ~1.36e-4 short of the core headline 0.380284),
worst @ (h=0.004, p=0.392, q=0.025) wit=cde_n30_iter3.

This script places FRESH dual-feasible centers along the p-range at the worst
|q| edge (q1=q2=0.025) and re-evaluates the combined cover (12 core + fresh)
with _fullspace_eval.cover_min_over_box, exactly the task's requested
center+grid+Lipschitz mechanism. The corrected mside_sin_coeff=4.0 is hardcoded
in path_b_analytical.build_problem_with_dual_handles's 5.6/5.7 RHS (literal 4.0).

Anchor convention (matches _promote_R8_centers.py): store 'primal' := the
dual-extracted rigorous LB, so cover_min_over_box's 'primal_m1e5' mode yields
dual_lb - 1e-5 (doubly conservative). Each center solved SEQUENTIALLY (memory).
"""
from __future__ import annotations
import argparse, json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
from _fullspace_eval import (load_centers, cover_min_over_box, CORE_HEADLINE,
                             WHITE_OUTSIDE_FLOOR)

OUT_CENTERS = CODE.parent / "parallel_results" / "fullspace_promote_R17.json"
TARGET = CORE_HEADLINE  # 0.380284
H_R, P_R, Q_R = (0.0, 0.06), (0.33, 0.45), (0.02, 0.025)


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        res = solve_with_dual_extraction(prob)
    except Exception as e:
        return {"status": "solver_failed", "reported_value": None,
                "rigorous_dual_LB": None, "dual_residual_at_LB": None,
                "time": 0.0, "error": f"{type(e).__name__}: {e}"}, {}
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return res, duals


def main():
    ap = argparse.ArgumentParser()
    # SDP-scale is the limiter at the HARD q=0.02 edge: at p=0.39 the primal is
    # only 0.38009 (N=8000) but 0.38063 (N=20000/bn=40), clearing target with
    # ~3.3e-4 margin. The hard edge is q=0.02 (NOT q=0.025): primal at p=0.39
    # rises monotonically with q (0.38009->0.38020->0.38036 for q=0.02->0.0225->0.025).
    # So place strong N=20000/bn=40 centers at the q=0.02 (single-point) edge.
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=40)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--margin", type=float, default=1e-5)
    # center grid: (h_c, p_c) pairs at the HARD |q| edge q1=q2=0.02.
    ap.add_argument("--h_c", type=float, default=0.0)
    ap.add_argument("--q_c", type=float, default=0.02)
    ap.add_argument("--p_grid", type=str, default="0.36,0.39,0.42")
    args = ap.parse_args()

    config = {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bn,
              "pm_k_max": args.pm_k, "mside_sin_coeff": 4.0}
    p_list = [float(x) for x in args.p_grid.split(",")]

    fresh = []
    if OUT_CENTERS.exists():
        try:
            fresh = json.load(open(OUT_CENTERS)).get("centers", [])
        except Exception:
            fresh = []
    have = {(round(c["h_c"], 5), round(c["p_c"], 5), round(c["q1"], 5))
            for c in fresh}

    print(f"=== R17 fresh centers: config N={args.N} T={args.T} bn={args.bn} "
          f"pm_k={args.pm_k} (mside_sin_coeff=4.0) ===")
    print(f"region: h{H_R} p{P_R} q{Q_R}")
    print(f"centers: h_c={args.h_c}, q_c={args.q_c}, p in {p_list}\n", flush=True)

    for p_c in p_list:
        key = (round(args.h_c, 5), round(p_c, 5), round(args.q_c, 5))
        if key in have:
            print(f"  skip existing (h={args.h_c},p={p_c},q={args.q_c})", flush=True)
            continue
        t0 = time.time()
        res, duals = solve_center(args.h_c, p_c, args.q_c, args.q_c,
                                  args.N, args.T, args.R, args.bn, args.pm_k)
        dt = time.time() - t0
        st = res["status"]
        if res["rigorous_dual_LB"] is None:
            print(f"  (h={args.h_c},p={p_c},q={args.q_c}): {st} -> dropped "
                  f"({dt:.0f}s) err={res.get('error')}", flush=True)
            continue
        anchor = res["rigorous_dual_LB"]
        rec = {
            "label": f"R17_h{args.h_c}_c{p_c}_q{args.q_c}",
            "h_c": args.h_c, "p_c": p_c, "q1": args.q_c, "q2": args.q_c,
            "primal": anchor,  # so cover_min_over_box 'primal_m1e5' = dual_lb - 1e-5
            "reported_primal": res["reported_value"], "dual_lb": anchor,
            "dual_resid": res["dual_residual_at_LB"], "status": st,
            "duals": duals, "time": dt, "config": config,
            "file": str(OUT_CENTERS),
        }
        fresh.append(rec)
        print(f"  (h={args.h_c},p={p_c},q={args.q_c}): primal={res['reported_value']:.7f} "
              f"dualLB={anchor:.7f} resid={res['dual_residual_at_LB']:.2e} "
              f"con_513={duals['con_513']:.4f} con_512qU={duals['con_512_qU']:.4f} "
              f"({dt:.0f}s)", flush=True)
        OUT_CENTERS.write_text(json.dumps(
            {"region": 17, "target": TARGET, "config": config,
             "anchor": "dual_LB - 1e-5 (conservative)",
             "h_range": list(H_R), "p_range": list(P_R), "q_range": list(Q_R),
             "centers": fresh}, indent=2, default=float))

    # ---- combined cover over the FULL R17 box ----
    core, _ = load_centers()
    combined = core + [{"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                        "q1": c["q1"], "q2": c["q2"], "primal": c["primal"],
                        "duals": c["duals"]} for c in fresh]
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        combined, "primal_m1e5", H_R, P_R, Q_R, n_h=61, n_p=161, n_q=21)
    print(f"\n=== R17 FULL-box cover ({len(core)} core + {len(fresh)} fresh) ===")
    print(f"box h{H_R} p{P_R} q{Q_R}")
    print(f"phi_min_lb={lb:.7f} grid_min={gmin:.7f} eps_grid={eps:.2e} L_max={Lm:.3f}")
    print(f"worst @ (h={pt[0]:.5f}, p={pt[1]:.5f}, q={pt[2]:.5f}) wit={wit}")
    print(f"clears 0.380000 (indep): {lb >= WHITE_OUTSIDE_FLOOR}; "
          f"clears 0.380284 (indep): {lb >= TARGET}; shortfall={TARGET-lb:.7f}")

    out = json.load(open(OUT_CENTERS)) if OUT_CENTERS.exists() else {"centers": fresh}
    out.update({
        "region": 17, "target": TARGET, "config": config,
        "h_range": list(H_R), "p_range": list(P_R), "q_range": list(Q_R),
        "n_core": len(core), "n_fresh": len(fresh),
        "cover_phi_min": lb, "cover_grid_min": gmin, "cover_eps_grid": eps,
        "cover_L_max": Lm, "cover_worst_point": pt, "cover_witness": wit,
        "clears_380000_indep": bool(lb >= WHITE_OUTSIDE_FLOOR),
        "clears_380284_indep": bool(lb >= TARGET),
    })
    OUT_CENTERS.write_text(json.dumps(out, indent=2, default=float))
    print(f"saved -> {OUT_CENTERS}")


if __name__ == "__main__":
    main()
