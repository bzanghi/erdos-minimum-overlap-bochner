"""
STAGE 2 (halo) — fresh AUGMENTED Phi-centers for the near-core "halo" regions.

The box-LP-min subdivision (in _fullspace_stage2_solve.py) cleanly clears the wide
far regions (R1-R5, R10, and the |d1|-large parts of R6-R9 are infeasible or high).
But the NEAR-CORE halo -- c1 in [~0.31, 0.44], small |d1| in [0.02, 0.05], i.e. the
problematic parts of R7, R9, R16, R17 -- cannot be cleared by box-LP-min at feasible
compute: the range relaxation's (5.13) max(p^2),max(q^2) conservatism is ~1.5e-3 even
at production for a c1-width 0.0375 box, so a narrow range-box gives 0.3791 < TARGET
while the single interior point gives 0.3808 > TARGET.

The right tool here is White's Phi-extension from CENTERS (smooth quadratic, exact at
the center, no range-relaxation conservatism). The earlier failure was center PLACEMENT
(centers at c1 in {0,0.15,0.5} missed the worst point at c1~0.39). Here we place fresh
PRODUCTION centers spanning the halo c1-range at a q-range that COVERS each halo's d1
extent, so Phi has no fatal q-decay over the halo. Phi from each center is a GLOBALLY
valid LB; combined cover = max over (existing 12 + new). We then evaluate the combined
cover's box-min over each halo region with the SAME rigorous grid+Lipschitz method as
Stage 1 (_fullspace_eval.cover_min_over_box).

Anchors: CONSERVATIVE dual-extracted, anchor = dual_LB - margin (margin 1e-5), exactly
the existing-cover convention. Each new center is independently dual-feasible, so adding
it only raises the cover (validity automatic).

CONFIG: production N=10000, T=4000, bn=40, pm_k=20 (single-point solves, ~90 s each).
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
from _fullspace_eval import load_centers, cover_min_over_box, WHITE_TABLE2, CORE_HEADLINE

OUT_CENTERS = CODE.parent / "parallel_results" / "fullspace_stage2_halo_centers.json"
TARGET = 0.380284


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return res, duals


# Halo region indices (the ones box-LP-min can't clear cheaply) and their boxes.
# (h_range, p_range, q_range)
HALO_REGIONS = {
    7:  ((0.0, 0.08), (0.0, 1.0), (-0.05, -0.025)),
    9:  ((0.0, 0.08), (0.0, 1.0), (0.025, 0.05)),
    16: ((0.0, 0.06), (0.33, 0.45), (-0.025, -0.02)),
    17: ((0.0, 0.06), (0.33, 0.45), (0.02, 0.025)),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=40)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--margin", type=float, default=1e-5)
    # center grid: c1 values spanning the halo, q-range covering halo d1 extents.
    ap.add_argument("--c1s", type=str, default="0.34,0.39,0.44")
    ap.add_argument("--hs", type=str, default="0.0")
    ap.add_argument("--qlo", type=float, default=-0.05)
    ap.add_argument("--qhi", type=float, default=0.05)
    args = ap.parse_args()

    c1_list = [float(x) for x in args.c1s.split(",")]
    h_list = [float(x) for x in args.hs.split(",")]

    # load any previously-saved halo centers (accumulate)
    halo = []
    if OUT_CENTERS.exists():
        try:
            halo = json.load(open(OUT_CENTERS)).get("centers", [])
        except Exception:
            halo = []
    have = {(round(c["h_c"], 5), round(c["p_c"], 5), round(c["q1"], 5), round(c["q2"], 5))
            for c in halo}

    print(f"=== Stage-2 halo centers (production N={args.N} bn={args.bn}) ===")
    print(f"c1 grid={c1_list} h grid={h_list} q-range=[{args.qlo},{args.qhi}]\n", flush=True)

    for h_c in h_list:
        for p_c in c1_list:
            key = (round(h_c, 5), round(p_c, 5), round(args.qlo, 5), round(args.qhi, 5))
            if key in have:
                print(f"  skip existing center (h={h_c}, c1={p_c})", flush=True)
                continue
            t0 = time.time()
            res, duals = solve_center(h_c, p_c, args.qlo, args.qhi,
                                      args.N, args.T, args.R, args.bn, args.pm_k)
            dt = time.time() - t0
            if res["rigorous_dual_LB"] is None:
                print(f"  center (h={h_c}, c1={p_c}): {res['status']} (no dual) ({dt:.0f}s)",
                      flush=True)
                continue
            anchor = res["rigorous_dual_LB"]
            halo.append({
                "label": f"halo_h{h_c}_c{p_c}_q{args.qlo}_{args.qhi}",
                "h_c": h_c, "p_c": p_c, "q1": args.qlo, "q2": args.qhi,
                "primal": res["reported_value"], "dual_lb": anchor,
                "dual_resid": res["dual_residual_at_LB"], "status": res["status"],
                "duals": duals, "time": dt,
            })
            print(f"  center (h={h_c}, c1={p_c}): primal={res['reported_value']:.6f} "
                  f"dualLB={anchor:.6f} con_513={duals['con_513']:.4f} ({dt:.0f}s)",
                  flush=True)
            OUT_CENTERS.write_text(json.dumps(
                {"config": {"N": args.N, "T": args.T, "bn": args.bn, "pm_k": args.pm_k},
                 "centers": halo}, indent=2, default=float))

    # ---- Evaluate combined cover (existing 12 + halo) over each halo region ----
    # Adapt halo centers to load_centers() schema (keys: label,h_c,p_c,q1,q2,primal,duals).
    existing, _ = load_centers()
    combined = existing + halo
    print(f"\n=== Combined cover: {len(existing)} existing + {len(halo)} halo = "
          f"{len(combined)} centers ===\n", flush=True)

    # anchor mode: existing centers use primal-1e-5; halo centers we want dual_LB-margin.
    # cover_min_over_box uses anchor_value(c, mode); for uniformity we patch halo centers'
    # 'primal' so that primal-1e-5 == dual_LB-margin (i.e. set primal := dual_lb - margin + 1e-5).
    # Cleaner: temporarily set each halo center primal to (dual_lb) so primal-1e-5 = dual_lb-1e-5.
    for c in halo:
        c["primal"] = c["dual_lb"]  # so anchor 'primal_m1e5' = dual_lb - 1e-5 (conservative)

    results = {}
    overall_min = np.inf
    overall_region = None
    for rid, (hr, pr, qr) in HALO_REGIONS.items():
        n_h = 41
        n_p = 161 if (pr[1] - pr[0]) > 0.2 else 81
        n_q = 41
        lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
            combined, "primal_m1e5", hr, pr, qr, n_h=n_h, n_p=n_p, n_q=n_q)
        cleared = lb >= TARGET
        results[rid] = {"region": rid, "h_range": list(hr), "p_range": list(pr),
                        "q_range": list(qr), "combined_phi_min": lb,
                        "grid_min": gmin, "eps_grid": eps, "L_max": Lm,
                        "worst_point": pt, "witness": wit,
                        "cleared_target": bool(cleared)}
        flag = "OK" if cleared else "BELOW"
        print(f"[R{rid}] combined-Phi-min={lb:.6f} {flag}  worst @ "
              f"(h={pt[0]:.4f}, c1={pt[1]:.4f}, d1={pt[2]:.4f}) wit={wit}", flush=True)
        if lb < overall_min:
            overall_min = lb
            overall_region = rid

    out = {
        "method": "fresh production Phi-centers (halo) + combined-cover box-min eval",
        "target": TARGET, "anchor": "dual_LB - 1e-5 (conservative)",
        "config": {"N": args.N, "T": args.T, "bn": args.bn, "pm_k": args.pm_k},
        "n_existing": len(existing), "n_halo": len(halo),
        "centers": halo,
        "region_eval": [results[r] for r in sorted(results)],
        "halo_min": overall_min, "halo_min_region": overall_region,
        "all_halo_cleared": bool(all(results[r]["cleared_target"] for r in results)),
    }
    OUT_CENTERS.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nhalo min over regions = {overall_min:.6f} at R{overall_region}")
    print(f"all halo regions cleared >= {TARGET}: {out['all_halo_cleared']}")
    print(f"saved -> {OUT_CENTERS}")


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
