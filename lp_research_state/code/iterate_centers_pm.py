"""
Phase 5: iterative center addition USING poly_moment + Hankel-PSD constraints.

Same idea as iterate_centers.py but each new center is solved with the augmented
constraint set, so binding-point search reflects the actual current cover.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_with_polymoment import solve_with_pm
from path_b_analytical import find_ellipse_h_p

H_BOX = (0.0, 0.06); P_BOX = (0.35, 0.45); TARGET = 0.379005
MARGIN = 1e-6


def envelope_min(rows, n_grid=4001, margin=MARGIN):
    h_grid = np.linspace(*H_BOX, n_grid); p_grid = np.linspace(*P_BOX, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf, dtype=float)
    witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(rows):
        e = r["ellipse"]; V_c_rig = r["V_c"] - margin
        F = (V_c_rig + e.get("const_q", 0.0)
             + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        mask = F > env; env[mask] = F[mask]; witness[mask] = i
    arg = np.unravel_index(int(env.argmin()), env.shape)

    L_max = 0.0
    for r in rows:
        e = r["ellipse"]
        def lin_max_abs(c2, c1, lo, hi):
            return max(abs(2*c2*lo + c1), abs(2*c2*hi + c1))
        L = float(np.sqrt(lin_max_abs(e["A_h2"], e["A_h1"], *H_BOX)**2
                          + lin_max_abs(e["A_p2"], e["A_p1"], *P_BOX)**2))
        L_max = max(L_max, L)
    cell = (H_BOX[1]-H_BOX[0])/(n_grid-1); cell_p = (P_BOX[1]-P_BOX[0])/(n_grid-1)
    eps_grid = L_max * 0.5 * float(np.sqrt(cell*cell + cell_p*cell_p))
    rigorous_LB = float(env.min()) - eps_grid

    return {
        "grid_min": float(env.min()),
        "rigorous_LB": rigorous_LB,
        "eps_grid": eps_grid,
        "h_min": float(HH[arg]), "p_min": float(PP[arg]),
        "witness_idx": int(witness[arg]),
        "witness_label": rows[int(witness[arg])]["label"],
    }


def load_phase_centers(json_path):
    """Load centers from a path_b_with_polymoment.py output JSON."""
    d = json.load(open(json_path))
    rows = []
    for c in d["centers"]:
        if "error" in c: continue
        rows.append({"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                     "V_c": c["V_c"], "ellipse": c["ellipse"]})
    return rows, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start_from", type=str,
                    default="lp_research_state/parallel_results/cde_phase4b.json")
    ap.add_argument("--max_iters", type=int, default=5)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=30)
    ap.add_argument("--pm_k_max", type=int, default=20)
    ap.add_argument("--hankel_n", type=int, default=6)
    ap.add_argument("--min_delta", type=float, default=5e-7)
    ap.add_argument("--out", type=str,
                    default="lp_research_state/parallel_results/cde_phase5.json")
    args = ap.parse_args()

    rows, init_data = load_phase_centers(CODE.parent.parent / args.start_from)
    print(f"=== Phase 5: cover iteration with poly_moment + Hankel-PSD ===")
    print(f"  Starting cover: {len(rows)} centers from {args.start_from}")
    print(f"  bochner_n={args.bochner_n} pm_k_max={args.pm_k_max} hankel_n={args.hankel_n}\n")

    e = envelope_min(rows)
    print(f"[0] start: rigorous_LB = {e['rigorous_LB']:.7f}  binding=({e['h_min']:.5f}, "
          f"{e['p_min']:.5f}) via {e['witness_label']}")
    history = [{"iter": 0, **e}]

    last = e["rigorous_LB"]
    for it in range(1, args.max_iters + 1):
        h_b, p_b = e["h_min"], e["p_min"]
        print(f"\n[{it}] solving at binding ({h_b:.5f}, {p_b:.5f}) ...", flush=True)
        t0 = time.time()
        r = solve_with_pm(args.N, args.T, args.R, h_b, p_b, -0.02, 0.02,
                          args.bochner_n, args.pm_k_max, hankel_n=args.hankel_n)
        dt = time.time() - t0
        V_c_rig = r["value"] - MARGIN
        center = {"h_c": h_b, "p_c": p_b, "q1": -0.02, "q2": 0.02, "value": V_c_rig}
        ell = find_ellipse_h_p(center, r["duals"], -0.02, 0.02, target=TARGET)
        rows.append({"label": f"cde_phase5_iter{it}", "h_c": h_b, "p_c": p_b,
                     "V_c": r["value"], "ellipse": ell})
        e = envelope_min(rows)
        delta = e["rigorous_LB"] - last
        print(f"  V_c = {r['value']:.7f}  ({dt:.0f}s)")
        print(f"  new rigorous_LB = {e['rigorous_LB']:.7f}  Δ = {delta:+.7f}  "
              f"binding now: ({e['h_min']:.5f}, {e['p_min']:.5f}) via {e['witness_label']}")
        history.append({"iter": it, "new_center": {"h_c": h_b, "p_c": p_b,
                        "V_c": r["value"], "duals": r["duals"], "elapsed_s": dt}, **e,
                        "delta": delta})
        if delta < args.min_delta:
            print(f"\nSaturated (Δ < {args.min_delta:.0e})"); break
        last = e["rigorous_LB"]

    out = CODE.parent.parent / args.out
    out.write_text(json.dumps({
        "config": vars(args),
        "history": history,
        "final_rigorous_LB": e["rigorous_LB"],
        "improvement_vs_white": e["rigorous_LB"] - TARGET,
        "improvement_vs_phase4b_start": history[0]["rigorous_LB"] - history[-1]["rigorous_LB"],
        "n_centers_final": len(rows),
    }, indent=2, default=float))
    print(f"\n=== Phase 5 Final ===")
    print(f"  µ ≥ {e['rigorous_LB']:.7f}")
    print(f"  vs White:  +{e['rigorous_LB'] - TARGET:.7f}")
    print(f"  Saved to {out}")


if __name__ == "__main__":
    main()
