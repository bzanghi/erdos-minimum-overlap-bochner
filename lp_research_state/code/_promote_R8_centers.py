"""
PRO-38 R8 — FRESH AUGMENTED Phi-CENTERS deliverable (task-specified mechanism).

Companion to _promote_R8_boxlp.py (the box-LP subdivision certificate). Here we place
fresh dual-feasible centers over the FEASIBLE part of R8 and evaluate the combined cover
(12 core + fresh) with _fullspace_eval.cover_min_over_box, exactly the task's requested
center+grid+Lipschitz mechanism, with the corrected mside_sin_coeff=4.0 (hardcoded in
path_b_analytical.build_problem_with_dual_handles 5.6/5.7 RHS).

R8 feasibility (mapped by _promote_R8_feas/_lowq): the augmented program enforces
|d|<=2/pi and sum_squares(c,d)<=0.5, so the strip q in [0.05,1.0] is FEASIBLE only for
roughly p(=c1) <~ 0.6 and q <~ q_max(p) (q_max ~0.42 at p=0, shrinking with p); the rest is
cleanly infeasible (vacuously covered). Feasible-interior primals are 0.399-0.58, tightest
at the low-q edge (min ~0.3991 at h=0,p=0.25,q=0.05). So a modest center grid suffices.

Centers are placed at single-point q (q1=q2) per task. Anchor = conservative dual_LB - 1e-5
(stored as 'primal':=dual_lb so cover_min_over_box's 'primal_m1e5' = dual_lb - 1e-5).
We evaluate cover over the FEASIBLE sub-box of R8 (p in [0,P_FEAS], q in [0.05,Q_FEAS]); the
complementary infeasible part is certified vacuous by clean-infeasible solves (box-LP file).
Each center solved SEQUENTIALLY (memory). Infeasible center solves are dropped (a center must
be feasible to yield a finite anchor).
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

OUT = CODE.parent / "parallel_results" / "fullspace_promote_R8.json"  # box-LP file; centers go to a sibling
OUT_CENTERS = CODE.parent / "parallel_results" / "fullspace_promote_R8_centers.json"
TARGET = CORE_HEADLINE  # 0.380284


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        res = solve_with_dual_extraction(prob)
        st = res["status"]
    except Exception as e:
        return {"status": "solver_failed", "reported_value": None,
                "rigorous_dual_LB": None, "dual_residual_at_LB": None,
                "time": 0.0, "error": f"{type(e).__name__}"}, {}
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return res, duals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--T", type=int, default=2000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=30)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--margin", type=float, default=1e-5)
    # feasible sub-box of R8 to certify with the center cover:
    ap.add_argument("--P_FEAS", type=float, default=0.6)
    ap.add_argument("--Q_FEAS", type=float, default=0.40)
    # center grid (feasible zone). p x q at h in {0, 0.08}. low-q edge reinforced.
    ap.add_argument("--p_grid", type=str, default="0.0,0.2,0.4,0.6")
    ap.add_argument("--q_grid", type=str, default="0.08,0.18,0.30")
    ap.add_argument("--h_grid", type=str, default="0.0,0.08")
    args = ap.parse_args()

    config = {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bn,
              "pm_k_max": args.pm_k, "mside_sin_coeff": 4.0}
    p_list = [float(x) for x in args.p_grid.split(",")]
    q_list = [float(x) for x in args.q_grid.split(",")]
    h_list = [float(x) for x in args.h_grid.split(",")]

    fresh = []
    if OUT_CENTERS.exists():
        try:
            fresh = json.load(open(OUT_CENTERS)).get("centers", [])
        except Exception:
            fresh = []
    have = {(round(c["h_c"], 4), round(c["p_c"], 4), round(c["q1"], 4))
            for c in fresh}

    print(f"=== R8 fresh centers: config N={args.N} bn={args.bn} pm_k={args.pm_k} "
          f"(mside_sin_coeff=4.0) ===")
    print(f"feasible sub-box: p in [0,{args.P_FEAS}] q in [0.05,{args.Q_FEAS}]")
    print(f"center grid: h{h_list} x p{p_list} x q{q_list}\n", flush=True)

    for h_c in h_list:
        for p_c in p_list:
            for q_c in q_list:
                key = (round(h_c, 4), round(p_c, 4), round(q_c, 4))
                if key in have:
                    print(f"  skip existing (h={h_c},p={p_c},q={q_c})", flush=True)
                    continue
                t0 = time.time()
                res, duals = solve_center(h_c, p_c, q_c, q_c, args.N, args.T, args.R,
                                          args.bn, args.pm_k)
                dt = time.time() - t0
                st = res["status"]
                if res["rigorous_dual_LB"] is None:
                    print(f"  (h={h_c},p={p_c},q={q_c}): {st} -> dropped ({dt:.0f}s)",
                          flush=True)
                    continue
                anchor = res["rigorous_dual_LB"]
                rec = {
                    "label": f"R8_h{h_c}_c{p_c}_q{q_c}",
                    "h_c": h_c, "p_c": p_c, "q1": q_c, "q2": q_c,
                    "primal": anchor,  # so cover_min_over_box 'primal_m1e5' = dual_lb-1e-5
                    "reported_primal": res["reported_value"], "dual_lb": anchor,
                    "dual_resid": res["dual_residual_at_LB"], "status": st,
                    "duals": duals, "time": dt, "config": config,
                    "file": str(OUT_CENTERS),
                }
                fresh.append(rec)
                print(f"  (h={h_c},p={p_c},q={q_c}): primal={res['reported_value']:.6f} "
                      f"dualLB={anchor:.6f} resid={res['dual_residual_at_LB']:.1e} "
                      f"con_513={duals['con_513']:.4f} ({dt:.0f}s)", flush=True)
                OUT_CENTERS.write_text(json.dumps(
                    {"region": 8, "target": TARGET, "config": config,
                     "anchor": "dual_LB - 1e-5 (conservative)", "centers": fresh},
                    indent=2, default=float))

    # ---- combined cover over the FEASIBLE sub-box of R8 ----
    core, _ = load_centers()
    combined = core + [{"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                        "q1": c["q1"], "q2": c["q2"], "primal": c["primal"],
                        "duals": c["duals"]} for c in fresh]
    FB = ((0.0, 0.08), (0.0, args.P_FEAS), (0.05, args.Q_FEAS))
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        combined, "primal_m1e5", *FB, n_h=41, n_p=161, n_q=81)
    print(f"\n=== R8 FEASIBLE sub-box cover ({len(core)} core + {len(fresh)} fresh) ===")
    print(f"sub-box {FB}")
    print(f"phi_min_lb={lb:.6f} grid_min={gmin:.6f} eps_grid={eps:.2e} L_max={Lm:.3f}")
    print(f"worst @ (h={pt[0]:.4f}, p={pt[1]:.4f}, q={pt[2]:.4f}) wit={wit}")
    print(f"clears 0.380000 (indep): {lb >= WHITE_OUTSIDE_FLOOR}; "
          f"clears 0.380284 (indep): {lb >= TARGET}")

    out = json.load(open(OUT_CENTERS)) if OUT_CENTERS.exists() else {"centers": fresh}
    out.update({
        "region": 8, "target": TARGET, "config": config,
        "feasible_subbox": {"h": list(FB[0]), "p": list(FB[1]), "q": list(FB[2])},
        "n_core": len(core), "n_fresh": len(fresh),
        "feasible_cover_phi_min": lb, "feasible_cover_grid_min": gmin,
        "feasible_cover_eps_grid": eps, "feasible_cover_L_max": Lm,
        "feasible_cover_worst_point": pt, "feasible_cover_witness": wit,
        "clears_380000_indep_feasible": bool(lb >= WHITE_OUTSIDE_FLOOR),
        "clears_380284_indep_feasible": bool(lb >= TARGET),
        "note": ("center cover over the FEASIBLE sub-box only; the complementary part "
                 "of R8 (p>~0.6 / q>q_max) is cleanly SDP-infeasible (vacuous) -- see "
                 "fullspace_promote_R8.json box-LP leaves for the rigorous full-region "
                 "certificate."),
    })
    OUT_CENTERS.write_text(json.dumps(out, indent=2, default=float))
    print(f"saved -> {OUT_CENTERS}")


if __name__ == "__main__":
    main()
