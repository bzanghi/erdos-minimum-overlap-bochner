"""PRO-12 Phase 2 (throwaway driver): full Phase-5 ellipse-extension cover with
MOSEK as the SDP backend, compared to the CLARABEL-based headline.

This replicates path_b_with_polymoment.main()'s cover logic (re-solve each center,
recompute per-center ellipse via find_ellipse_h_p, then envelope-min over box
(5.16) with the Lipschitz eps_grid margin) but routes the per-center solve through
mosek_runner via solve_with_pm(..., solver="MOSEK").

Run modes:
  --mode lite7  : 7 White centers, N=10000 T=4000 bn=20, pm/hankel OFF.
                  De-risk config (~5s/row Mosek). Compares to CLARABEL same config.
  --mode head12 : 12 centers (7 White + 5 CDE), the 0.3803027 headline machine.
                  Default N=20000 bn=40 pm_k_max=20 hankel_n=6.

The final µ LB comes from the envelope-min cover (the SAME path_b_* machinery
that produced 0.3803027), NOT from min-over-rows of raw duals.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import find_ellipse_h_p  # noqa: E402
from path_b_with_polymoment import solve_with_pm, H_BOX, P_BOX, TARGET  # noqa: E402

WHITE = 0.379005
PRIOR_HEADLINE = 0.379544
CLARABEL_HEADLINE_LB = 0.3803027  # PRO-21 Phase 8 N=20K bn=40

RESULTS = CODE.parent / "parallel_results"           # worktree (output)
# Center-definition JSON (path_b/rowX.json, cde_iter_n30.json) may live only in
# the canonical main-repo parallel_results, not the worktree checkout.  Resolve a
# read source that actually contains them.
_CENTER_CANDIDATES = [
    RESULTS,
    Path("/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/parallel_results"),
]
CENTER_SRC = next((p for p in _CENTER_CANDIDATES
                   if (p / "path_b" / "row1.json").exists()), RESULTS)
OUT_JSON = RESULTS / "pro12_phase5_mosek.json"


def load_centers(n_centers):
    """7 White + (optionally) 5 CDE centers, exactly as path_b_with_polymoment.main."""
    centers = []
    for r in range(1, 8):
        d = json.load(open(CENTER_SRC / "path_b" / f"row{r}.json"))
        centers.append({"label": d["label"], "h_c": d["h_c"], "p_c": d["p_c"],
                        "q1": d["q1"], "q2": d["q2"]})
    if n_centers > 7:
        n30 = json.load(open(CENTER_SRC / "cde_iter_n30.json"))
        for h in n30["history"]:
            if "new_center" not in h:
                continue
            nc = h["new_center"]
            centers.append({"label": f"cde_n30_iter{h['iter']}",
                            "h_c": nc["h_c"], "p_c": nc["p_c"],
                            "q1": -0.02, "q2": 0.02})
    return centers


def cover_min(results, margin, n_grid=4001):
    """Envelope-min of best-of-rows over box (5.16) + Lipschitz eps_grid margin.
    Identical to path_b_with_polymoment.main()'s cover computation."""
    h_grid = np.linspace(*H_BOX, n_grid)
    p_grid = np.linspace(*P_BOX, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(results):
        if "error" in r:
            continue
        e = r["ellipse"]
        F = (r["V_c_rigorous"] + e.get("const_q", 0)
             + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        mask = F > env
        env[mask] = F[mask]
        witness[mask] = i
    grid_min = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)

    L_max = 0.0
    for r in results:
        if "error" in r:
            continue
        e = r["ellipse"]
        def lin_max_abs(c2, c1, lo, hi):
            return max(abs(2 * c2 * lo + c1), abs(2 * c2 * hi + c1))
        L = float(np.sqrt(lin_max_abs(e["A_h2"], e["A_h1"], *H_BOX) ** 2
                          + lin_max_abs(e["A_p2"], e["A_p1"], *P_BOX) ** 2))
        L_max = max(L_max, L)
    cell_h = (H_BOX[1] - H_BOX[0]) / (n_grid - 1)
    cell_p = (P_BOX[1] - P_BOX[0]) / (n_grid - 1)
    eps_grid = L_max * 0.5 * float(np.sqrt(cell_h ** 2 + cell_p ** 2))
    return {
        "grid_min": grid_min,
        "binding_point": [float(HH[arg]), float(PP[arg])],
        "witness": results[int(witness[arg])]["label"],
        "L_max_grad": L_max,
        "eps_grid": eps_grid,
        "rigorous_LB": grid_min - eps_grid,
        "n_grid": n_grid,
    }


def run(N, T, R, bochner_n, pm_k_max, hankel_n, use_T5p, margin, solver, n_centers,
        label):
    centers = load_centers(n_centers)
    print(f"=== {label}: {len(centers)} centers @ solver={solver}, "
          f"N={N} T={T} R={R} bn={bochner_n} pm={pm_k_max} hankel={hankel_n} ===")
    results = []
    t0 = time.time()
    for i, c in enumerate(centers):
        print(f"[{i+1}/{len(centers)}] {c['label']:18s} "
              f"(h={c['h_c']:.4f}, p={c['p_c']:.4f}) ...", flush=True)
        try:
            r = solve_with_pm(N, T, R, c["h_c"], c["p_c"], c["q1"], c["q2"],
                              bochner_n, pm_k_max, hankel_n=hankel_n,
                              use_T5p=use_T5p, solver=solver)
            V_c_rig = r["value"] - margin
            center = {**c, "value": V_c_rig}
            ell = find_ellipse_h_p(center, r["duals"], c["q1"], c["q2"], target=TARGET)
            row = {**c, "V_c": r["value"], "V_c_rigorous": V_c_rig,
                   "status": r["status"], "duals": r["duals"],
                   "ellipse": ell, "time_s": r["time"], "solver": r.get("solver")}
            for k in ("rigorous_dual_LB", "primal_obj", "duality_gap",
                      "mosek_problem_status", "mosek_solution_status",
                      "primal_viol", "dual_viol", "iterations"):
                if k in r:
                    row[k] = r[k]
            results.append(row)
            extra = ""
            if "duality_gap" in r:
                extra = (f"  gap={r['duality_gap']:.2e}  "
                         f"[{r.get('mosek_solution_status')}]")
            print(f"   V_c={r['value']:.10f}  ({r['status']}, {r['time']:.1f}s){extra}",
                  flush=True)
        except Exception as e:
            print(f"   ERROR: {type(e).__name__}: {e}", flush=True)
            results.append({**c, "error": str(e)})
    total = time.time() - t0

    cov = cover_min(results, margin)
    print(f"\n--- {label} cover (envelope-min over box 5.16) ---")
    print(f"  grid_min      = {cov['grid_min']:.10f}")
    print(f"  eps_grid      = {cov['eps_grid']:.3e}  (L_max={cov['L_max_grad']:.4f})")
    print(f"  margin        = {margin:.1e}")
    print(f"  RIGOROUS µ LB = {cov['rigorous_LB']:.10f}")
    print(f"  binding       = {cov['binding_point']}  witness={cov['witness']}")
    print(f"  vs White      : {cov['rigorous_LB'] - WHITE:+.7e}")
    print(f"  vs CLARABEL headline 0.3803027: {cov['rigorous_LB'] - CLARABEL_HEADLINE_LB:+.3e}")

    return {
        "label": label,
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n,
                   "pm_k_max": pm_k_max, "hankel_n": hankel_n,
                   "use_T5p": use_T5p, "margin": margin, "solver": solver,
                   "n_centers": len(centers)},
        "centers": results,
        "cover": cov,
        "rigorous_LB": cov["rigorous_LB"],
        "improvement_vs_white": cov["rigorous_LB"] - WHITE,
        "vs_clarabel_headline_0p3803027": cov["rigorous_LB"] - CLARABEL_HEADLINE_LB,
        "total_time_s": total,
    }


def save(rec):
    existing = []
    if OUT_JSON.exists():
        existing = json.load(open(OUT_JSON))
    existing = [r for r in existing if r.get("label") != rec["label"]]
    existing.append(rec)
    OUT_JSON.write_text(json.dumps(existing, indent=2, default=float))
    print(f"\n→ wrote {OUT_JSON}")


def main():
    import warnings
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["lite7", "head12"], default="lite7")
    ap.add_argument("--solver", default="MOSEK")
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=None)
    ap.add_argument("--pm_k_max", type=int, default=None)
    ap.add_argument("--hankel_n", type=int, default=None)
    ap.add_argument("--use_T5p", action="store_true")
    ap.add_argument("--margin", type=float, default=1e-6)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    if args.mode == "lite7":
        N = args.N or 10000
        bn = args.bochner_n if args.bochner_n is not None else 20
        pm = args.pm_k_max if args.pm_k_max is not None else 0
        hk = args.hankel_n if args.hankel_n is not None else 0
        n_centers = 7
        label = args.label or f"lite7_{args.solver}_N{N}_bn{bn}"
    else:  # head12
        N = args.N or 20000
        bn = args.bochner_n if args.bochner_n is not None else 40
        pm = args.pm_k_max if args.pm_k_max is not None else 20
        hk = args.hankel_n if args.hankel_n is not None else 6
        n_centers = 12
        label = args.label or f"head12_{args.solver}_N{N}_bn{bn}"

    rec = run(N, args.T, args.R, bn, pm, hk, args.use_T5p, args.margin,
              args.solver, n_centers, label)
    save(rec)


if __name__ == "__main__":
    main()
