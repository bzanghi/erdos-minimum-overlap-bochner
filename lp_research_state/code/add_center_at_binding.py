"""
Phase-1 cap experiment: solve the Bochner-augmented SDP at a NEW center placed
at the current binding point of the envelope-min cover, add its ellipse to the
existing 7-row cover, recompute envelope-min.

Hypothesis: the binding point (h, p) = (0, 0.394) currently has V ≈ 0.37955.
A new SDP center there should give V_c ≈ 0.3797–0.3799 at full scale, and its
ellipse should raise the envelope min by a measurable Δµ.
"""
from __future__ import annotations
import argparse, json, os, sys, time
from pathlib import Path

import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import solve_and_extract_duals, find_ellipse_h_p

PATH_B_DIR = CODE.parent / "parallel_results" / "path_b"
H_BOX = (0.0, 0.06)
P_BOX = (0.35, 0.45)
TARGET = 0.379005  # White's published bound


def load_existing_rows():
    rows = []
    for r in range(1, 8):
        d = json.load(open(PATH_B_DIR / f"row{r}.json"))
        rows.append({
            "label": d["label"],
            "h_c": d["h_c"], "p_c": d["p_c"],
            "V_c": d["primal_value_at_center"],
            "ellipse": d["ellipse"],
            "config": d["config"],
        })
    return rows


def envelope_min(rows, n_grid=4001):
    h_grid = np.linspace(H_BOX[0], H_BOX[1], n_grid)
    p_grid = np.linspace(P_BOX[0], P_BOX[1], n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf, dtype=float)
    witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(rows):
        e = r["ellipse"]
        F = (e["V_c"] + e.get("const_q", 0.0)
             + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        mask = F > env
        env[mask] = F[mask]
        witness[mask] = i
    arg = np.unravel_index(int(env.argmin()), env.shape)
    return {
        "grid_min": float(env.min()),
        "h_min": float(HH[arg]), "p_min": float(PP[arg]),
        "witness_row": int(witness[arg]),
        "row_label": rows[int(witness[arg])]["label"],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h_c", type=float, default=0.0)
    ap.add_argument("--p_c", type=float, default=0.394)
    ap.add_argument("--q1", type=float, default=-0.02)
    ap.add_argument("--q2", type=float, default=0.02)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=20)
    ap.add_argument("--label", type=str, default="row8_binding")
    ap.add_argument("--margin", type=float, default=1e-6)
    args = ap.parse_args()

    print("=" * 70)
    print("Constraint Discovery Engine — Phase 1 cap experiment")
    print("Add SDP center at current binding point to tighten the envelope cover")
    print("=" * 70)

    rows = load_existing_rows()
    print("\n[1] Current 7-row envelope min:")
    e0 = envelope_min(rows)
    print(f"  grid_min = {e0['grid_min']:.7f}")
    print(f"  at (h, p) = ({e0['h_min']:.5f}, {e0['p_min']:.5f})")
    print(f"  witness = {e0['row_label']}")
    print(f"  vs White's 0.379005: +{e0['grid_min'] - TARGET:.7f}")
    print(f"  vs project headline 0.379544: {e0['grid_min'] - 0.379544:+.7f}")

    print(f"\n[2] Solving SDP at new center "
          f"(h_c={args.h_c}, p_c={args.p_c}, q=[{args.q1}, {args.q2}]) "
          f"at N={args.N}, T={args.T}, bochner_n={args.bochner_n}...")
    t0 = time.time()
    res = solve_and_extract_duals(
        args.N, args.T, args.R, args.h_c, args.p_c, args.q1, args.q2,
        bochner_n=args.bochner_n,
    )
    print(f"  solved in {res['time']:.1f}s, status = {res['status']}")
    print(f"  V_c (primal) = {res['value']:.10f}")
    print(f"  dual values: {res['duals']}")

    # Apply rigour margin (subtract for safety)
    V_c_rigorous = res["value"] - args.margin
    print(f"  V_c_rigorous (margin={args.margin}) = {V_c_rigorous:.10f}")

    print(f"\n[3] Computing new ellipse around (h_c, p_c) = ({args.h_c}, {args.p_c})...")
    center = {"h_c": args.h_c, "p_c": args.p_c,
              "q1": args.q1, "q2": args.q2, "value": V_c_rigorous}
    new_ellipse = find_ellipse_h_p(center, res["duals"], args.q1, args.q2, target=TARGET)
    new_row = {
        "label": args.label,
        "h_c": args.h_c, "p_c": args.p_c,
        "V_c": V_c_rigorous,
        "ellipse": new_ellipse,
        "config": {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bochner_n},
    }
    print(f"  ellipse params: semi_h={new_ellipse.get('semi_h')}, semi_p={new_ellipse.get('semi_p')}")
    print(f"  A_h2={new_ellipse['A_h2']:.4f}  A_p2={new_ellipse['A_p2']:.4f}")

    print(f"\n[4] Envelope min with 8 rows (new center added):")
    rows_aug = rows + [new_row]
    e1 = envelope_min(rows_aug)
    print(f"  grid_min = {e1['grid_min']:.7f}")
    print(f"  at (h, p) = ({e1['h_min']:.5f}, {e1['p_min']:.5f})")
    print(f"  witness = {e1['row_label']}")
    print(f"  vs White's 0.379005: +{e1['grid_min'] - TARGET:.7f}")
    print(f"  vs project headline 0.379544: {e1['grid_min'] - 0.379544:+.7f}")

    delta = e1['grid_min'] - e0['grid_min']
    print(f"\n[5] Improvement: Δgrid_min = {delta:+.7f}")
    if delta > 1e-6:
        print(f"  ✓ Center addition produced a measurable improvement")
    elif delta > 0:
        print(f"  ~ Marginal improvement (within margin/noise)")
    else:
        print(f"  ✗ No improvement (new center dominated by existing cover at binding point)")

    out = {
        "phase": "1-cap",
        "before_grid_min": e0["grid_min"],
        "after_grid_min": e1["grid_min"],
        "delta": delta,
        "new_center": {"h_c": args.h_c, "p_c": args.p_c, "q1": args.q1, "q2": args.q2,
                       "V_c": res["value"], "V_c_rigorous": V_c_rigorous,
                       "duals": res["duals"], "ellipse": new_ellipse,
                       "config": new_row["config"], "elapsed_s": res["time"]},
        "before_binding": e0, "after_binding": e1,
    }
    out_file = CODE.parent / "parallel_results" / f"cde_phase1_{args.label}.json"
    out_file.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n→ saved to {out_file}")


if __name__ == "__main__":
    main()
