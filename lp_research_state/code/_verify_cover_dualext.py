"""
Re-derive the Phase-5 headline cover with RIGOROUS dual-extracted anchors.

Motivation: phase5_N20K_bn40.json anchored each ellipse at  primal_value - 1e-6.
A fresh re-solve of the binding center gives a primal ~1e-5 BELOW the recorded
primal, so the fixed 1e-6 margin was not conservative enough to guarantee a valid
lower bound at optimal_inaccurate gaps. This script re-solves all 12 centers at
the production config (N=20000, T=4000, bochner_n=40, pm_k_max=20), captures BOTH
the fresh primal and the dual-extracted rigorous LB, and recomputes the cover
headline under conservative conventions.

Sequential (one solve at a time) to bound RAM. ~12 x 200s ~ 40 min.
"""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles, find_ellipse_h_p
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction

N, T, R, BN, PMK = 20000, 4000, 10, 40, 20
H_BOX = (0.0, 0.06); P_BOX = (0.35, 0.45); TARGET = 0.379005
OUT = CODE.parent / "parallel_results" / "phase5_N20K_bn40_dualext.json"

# center list: 7 White + 5 cde (same as path_b_with_polymoment.main)
centers = []
for r in range(1, 8):
    d = json.load(open(CODE.parent / "parallel_results" / "path_b" / f"row{r}.json"))
    centers.append({"label": d["label"], "h_c": d["h_c"], "p_c": d["p_c"],
                    "q1": d["q1"], "q2": d["q2"]})
n30 = json.load(open(CODE.parent / "parallel_results" / "cde_iter_n30.json"))
for h in n30["history"]:
    if "new_center" not in h:
        continue
    nc = h["new_center"]
    centers.append({"label": f"cde_n30_iter{h['iter']}", "h_c": nc["h_c"],
                    "p_c": nc["p_c"], "q1": -0.02, "q2": 0.02})

print(f"re-solving {len(centers)} centers @ N={N} bn={BN} pmk={PMK} with dual extraction\n", flush=True)
results = []
for i, c in enumerate(centers):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, c["h_c"], c["h_c"], c["p_c"], c["p_c"], c["q1"], c["q2"], bochner_n=BN)
    pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=PMK)
    cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: float(H[k].dual_value) if H[k].dual_value is not None else 0.0
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    primal = res["reported_value"]
    dual_lb = res["rigorous_dual_LB"]
    gap = (primal - dual_lb) if (primal is not None and dual_lb is not None) else None
    results.append({**c, "primal": primal, "dual_lb": dual_lb, "gap": gap,
                    "dual_resid": res["dual_residual_at_LB"], "status": res["status"],
                    "duals": duals, "time": res["time"]})
    print(f"[{i+1:2d}/{len(centers)}] {c['label']:16s} primal={primal:.8f} "
          f"dual_lb={dual_lb if dual_lb is None else round(dual_lb,6)} "
          f"gap={gap if gap is None else format(gap,'.2e')} ({res['time']:.0f}s)", flush=True)
    OUT.write_text(json.dumps({"config": {"N": N, "T": T, "R": R, "bochner_n": BN,
                   "pm_k_max": PMK}, "centers": results}, indent=2, default=float))


def cover_headline(anchor_fn, label):
    h_grid = np.linspace(*H_BOX, 4001); p_grid = np.linspace(*P_BOX, 4001)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    L_max = 0.0
    for c in results:
        anchor = anchor_fn(c)
        syn = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"], "value": anchor}
        e = find_ellipse_h_p(syn, c["duals"], c["q1"], c["q2"], target=TARGET)
        F = (anchor + e["const_q"] + e["A_h2"]*HH*HH + e["A_h1"]*HH + e["A_h0"]
             + e["A_p2"]*PP*PP + e["A_p1"]*PP + e["A_p0"])
        np.maximum(env, F, out=env)
        lam = lambda c2, c1, lo, hi: max(abs(2*c2*lo+c1), abs(2*c2*hi+c1))
        L_max = max(L_max, float(np.hypot(lam(e["A_h2"], e["A_h1"], *H_BOX),
                                          lam(e["A_p2"], e["A_p1"], *P_BOX))))
    gmin = float(env.min())
    cell_h = (H_BOX[1]-H_BOX[0])/4000; cell_p = (P_BOX[1]-P_BOX[0])/4000
    eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    print(f"  [{label}] grid_min={gmin:.8f}  eps_grid={eps_grid:.2e}  "
          f"rigorous_LB={gmin-eps_grid:.8f}  (vs White {gmin-eps_grid-0.379005:+.2e})")
    return gmin - eps_grid

print("\n=== recomputed cover headlines ===")
cover_headline(lambda c: c["primal"] - 1e-6, "fresh primal - 1e-6 (orig convention)")
cover_headline(lambda c: c["dual_lb"], "dual-extracted LB (rigorous, coarse 5-digit print)")
cover_headline(lambda c: c["primal"] - 1e-5, "fresh primal - 1e-5 (conservative)")
print(f"\nsaved -> {OUT.name}")
