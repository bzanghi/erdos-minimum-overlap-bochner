"""
Re-solve each path-B center with poly-moment constraints AND bochner_n=30,
extract new duals, recompute ellipses, recompute envelope min.

Combines Phase 2 (cover refinement + n=30) with Phase 3 (poly-moment) into the
unified rigorous LB.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles, find_ellipse_h_p
from poly_moment import build_even_moment_nonneg_constraints

H_BOX = (0.0, 0.06); P_BOX = (0.35, 0.45); TARGET = 0.379005


def solve_with_pm(N, T, R, h_c, p_c, q1, q2, bochner_n, pm_k_max):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bochner_n,
    )
    if pm_k_max > 0:
        pm_cons, tb = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k_max)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time(); prob.solve(solver="CLARABEL", verbose=False); dt = time.time() - t0
    duals = {k: float(H[k].dual_value) if H[k].dual_value is not None else 0.0
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return {"value": float(prob.value), "status": prob.status,
            "duals": duals, "time": dt}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=30)
    ap.add_argument("--pm_k_max", type=int, default=14)
    ap.add_argument("--margin", type=float, default=1e-6)
    ap.add_argument("--out", type=str, default="lp_research_state/parallel_results/cde_phase3.json")
    args = ap.parse_args()

    # Centers to re-solve
    centers = []
    for r in range(1, 8):
        d = json.load(open(CODE.parent / "parallel_results" / "path_b" / f"row{r}.json"))
        centers.append({"label": d["label"], "h_c": d["h_c"], "p_c": d["p_c"],
                        "q1": d["q1"], "q2": d["q2"]})
    n30 = json.load(open(CODE.parent / "parallel_results" / "cde_iter_n30.json"))
    for h in n30["history"]:
        if "new_center" not in h: continue
        nc = h["new_center"]
        centers.append({"label": f"cde_n30_iter{h['iter']}",
                        "h_c": nc["h_c"], "p_c": nc["p_c"],
                        "q1": -0.02, "q2": 0.02})

    print(f"=== Phase 3: re-solve all {len(centers)} centers at "
          f"bochner_n={args.bochner_n} + poly_moment k_max={args.pm_k_max} ===\n")

    results = []
    for i, c in enumerate(centers):
        print(f"[{i+1}/{len(centers)}] {c['label']:20s} (h={c['h_c']:.4f}, p={c['p_c']:.4f}) ...", flush=True)
        try:
            r = solve_with_pm(args.N, args.T, args.R, c['h_c'], c['p_c'], c['q1'], c['q2'],
                              args.bochner_n, args.pm_k_max)
            V_c_rig = r['value'] - args.margin
            center = {**c, "value": V_c_rig}
            ell = find_ellipse_h_p(center, r['duals'], c['q1'], c['q2'], target=TARGET)
            results.append({**c, "V_c": r['value'], "V_c_rigorous": V_c_rig,
                            "status": r['status'], "duals": r['duals'],
                            "ellipse": ell, "time_s": r['time']})
            print(f"   V_c = {r['value']:.7f}  ({r['status']}, {r['time']:.1f}s)", flush=True)
        except Exception as e:
            print(f"   ERROR: {type(e).__name__}: {e}", flush=True)
            results.append({**c, "error": str(e)})

    # Compute envelope min
    h_grid = np.linspace(*H_BOX, 4001); p_grid = np.linspace(*P_BOX, 4001)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing='ij')
    env = np.full_like(HH, -np.inf); witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(results):
        if "error" in r: continue
        e = r['ellipse']
        F = (r['V_c_rigorous'] + e.get('const_q', 0) + e['A_h2']*HH*HH + e['A_h1']*HH
             + e['A_h0'] + e['A_p2']*PP*PP + e['A_p1']*PP + e['A_p0'])
        mask = F > env; env[mask] = F[mask]; witness[mask] = i
    grid_min = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    L_max = 0.0
    for r in results:
        if "error" in r: continue
        e = r['ellipse']
        def lin_max_abs(c2, c1, lo, hi):
            return max(abs(2*c2*lo + c1), abs(2*c2*hi + c1))
        L = float(np.sqrt(lin_max_abs(e['A_h2'], e['A_h1'], *H_BOX)**2 + lin_max_abs(e['A_p2'], e['A_p1'], *P_BOX)**2))
        L_max = max(L_max, L)
    cell_h = (H_BOX[1]-H_BOX[0])/4000; cell_p = (P_BOX[1]-P_BOX[0])/4000
    eps_grid = L_max * 0.5 * float(np.sqrt(cell_h**2 + cell_p**2))
    rigorous_LB = grid_min - eps_grid

    print(f"\n=== Final rigorous LB ===")
    print(f"  grid_min = {grid_min:.7f}")
    print(f"  binding  = ({float(HH[arg]):.5f}, {float(PP[arg]):.5f})  witness={results[int(witness[arg])]['label']}")
    print(f"  eps_grid = {eps_grid:.2e}   L_max = {L_max:.4f}")
    print(f"  RIGOROUS LB µ ≥ {rigorous_LB:.7f}")
    print(f"  vs White (+5.4e-4 baseline): {rigorous_LB - 0.379005:+.7f}")
    print(f"  vs published 0.379544:      {rigorous_LB - 0.379544:+.7f}")

    out_path = CODE.parent.parent / args.out
    out_path.write_text(json.dumps({
        "config": {"N": args.N, "T": args.T, "R": args.R,
                   "bochner_n": args.bochner_n, "pm_k_max": args.pm_k_max,
                   "margin": args.margin},
        "rigorous_LB": rigorous_LB, "grid_min": grid_min,
        "binding_point": [float(HH[arg]), float(PP[arg])],
        "witness": results[int(witness[arg])]['label'],
        "improvement_vs_white": rigorous_LB - 0.379005,
        "improvement_vs_prior_headline": rigorous_LB - 0.379544,
        "eps_grid": eps_grid, "L_max_grad": L_max,
        "centers": results,
    }, indent=2, default=float))
    print(f"\n→ saved {out_path}")


if __name__ == "__main__":
    main()
