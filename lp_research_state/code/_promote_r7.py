"""
PRO-38 R7 promotion: fresh AUGMENTED dual-feasible centers for region 7.

R7 box (task): h_range=[0,0.08], p_range=[0,1], q_range=[-0.05,-0.025].
Diagnostic (existing 12 + prior 11 halo, fine subdivision): the ONLY band below
0.380284 at the strip's worst q=-0.05 is p in [0.35,0.42], grid_min ~0.38002
(>White 0.380000 but <core headline). Everywhere else clears 0.380284 independently.

We place fresh PRODUCTION-style centers in that weak band at the worst h (small,
near 0) and at the strip's worst edge q1=q2=-0.05 (single point => zero q-decay
penalty there; Phi is concave in q so q=-0.05, the largest |q|, is the worst case
over the strip). EXACT same encoding as _fullspace_stage2_halo_centers.solve_center
so duals compose with the existing 12-center cover (max_c Phi_c). Conservative
anchor = dual_LB - 1e-5 (we store primal := dual_lb so anchor_value primal_m1e5
gives dual_lb - 1e-5, matching the existing convention).
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
from _fullspace_eval import load_centers, cover_min_over_box, CORE_HEADLINE, WHITE_OUTSIDE_FLOOR

OUT = CODE.parent / "parallel_results" / "fullspace_promote_R7.json"
TARGET = CORE_HEADLINE  # 0.380284
DUAL_KEYS = ("con_53", "con_54", "con_512_pL", "con_512_pU",
             "con_512_qL", "con_512_qU", "con_513")


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in DUAL_KEYS}
    return res, duals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--T", type=int, default=2000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=30)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--margin", type=float, default=1e-5)
    # weak band centers: (h_c, p_c) ; q fixed at strip worst edge -0.05 (single point)
    ap.add_argument("--centers", type=str,
                    default="0.0:0.36,0.0:0.39,0.0:0.42")
    ap.add_argument("--q", type=float, default=-0.05)
    args = ap.parse_args()

    spec = [tuple(float(x) for x in pair.split(":")) for pair in args.centers.split(",")]
    q1 = q2 = args.q

    fresh = []
    if OUT.exists():
        try:
            fresh = json.load(open(OUT)).get("centers", [])
        except Exception:
            fresh = []
    have = {(round(c["h_c"], 5), round(c["p_c"], 5), round(c["q1"], 5), round(c["q2"], 5))
            for c in fresh}

    print(f"=== PRO-38 R7 fresh centers (N={args.N} T={args.T} bn={args.bn} pm_k={args.pm_k}) ===")
    print(f"spec={spec} q1=q2={q1}\n", flush=True)

    for (h_c, p_c) in spec:
        key = (round(h_c, 5), round(p_c, 5), round(q1, 5), round(q2, 5))
        if key in have:
            print(f"  skip existing (h={h_c}, p={p_c}, q={q1})", flush=True)
            continue
        t0 = time.time()
        res, duals = solve_center(h_c, p_c, q1, q2, args.N, args.T, args.R, args.bn, args.pm_k)
        dt = time.time() - t0
        if res["rigorous_dual_LB"] is None:
            print(f"  center (h={h_c}, p={p_c}): {res['status']} NO DUAL ({dt:.0f}s)", flush=True)
            continue
        anchor = res["rigorous_dual_LB"] - args.margin  # conservative dual_LB - margin
        fresh.append({
            "label": f"r7_h{h_c}_c{p_c}_q{q1}_{q2}_N{args.N}_bn{args.bn}",
            "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
            "primal": res["reported_value"],
            "dual_lb_raw": res["rigorous_dual_LB"],
            "anchor": anchor,
            "dual_resid": res["dual_residual_at_LB"],
            "status": res["status"], "duals": duals, "time": dt,
            "config": {"N": args.N, "T": args.T, "R": args.R, "bn": args.bn,
                       "pm_k": args.pm_k, "margin": args.margin},
        })
        print(f"  center (h={h_c}, p={p_c}): primal={res['reported_value']:.7f} "
              f"dualLB={res['rigorous_dual_LB']:.7f} anchor={anchor:.7f} "
              f"dual_resid={res['dual_residual_at_LB']:.2e} con_513={duals['con_513']:.4f} "
              f"status={res['status']} ({dt:.0f}s)", flush=True)
        OUT.write_text(json.dumps({"config": vars(args), "centers": fresh},
                                  indent=2, default=float))

    print(f"\nsaved {len(fresh)} fresh centers -> {OUT}", flush=True)


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
