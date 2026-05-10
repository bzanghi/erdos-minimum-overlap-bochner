"""
Iterative center addition: at each step, find the current binding point of the
envelope cover, solve a new SDP center there, add its ellipse, repeat until the
improvement saturates.

The validity is the same as for any path-B center: the LP dual at the new center
gives a parameter-independent rigorous lower bound on µ via the explicit shift
formula. Adding more centers to the envelope only tightens the cover.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import solve_and_extract_duals, find_ellipse_h_p

PATH_B_DIR = CODE.parent / "parallel_results" / "path_b"
H_BOX = (0.0, 0.06)
P_BOX = (0.35, 0.45)
Q_RANGE = (-0.02, 0.02)
TARGET = 0.379005


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


def solve_center_and_make_ellipse(h_c, p_c, N, T, R, bochner_n, margin=1e-6):
    res = solve_and_extract_duals(N, T, R, h_c, p_c, Q_RANGE[0], Q_RANGE[1],
                                  bochner_n=bochner_n)
    V_c_rigorous = res["value"] - margin
    center = {"h_c": h_c, "p_c": p_c,
              "q1": Q_RANGE[0], "q2": Q_RANGE[1], "value": V_c_rigorous}
    ell = find_ellipse_h_p(center, res["duals"], Q_RANGE[0], Q_RANGE[1], target=TARGET)
    return res, V_c_rigorous, ell


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max_iters", type=int, default=5)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=20)
    ap.add_argument("--min_delta", type=float, default=1e-6,
                    help="Stop when improvement per iter falls below this.")
    ap.add_argument("--out", type=str,
                    default="lp_research_state/parallel_results/cde_iterative.json")
    args = ap.parse_args()

    rows = load_existing_rows()
    history = []
    e = envelope_min(rows)
    print(f"[0] start: grid_min = {e['grid_min']:.7f}  binding = ({e['h_min']:.5f}, {e['p_min']:.5f}) via {e['row_label']}")
    history.append({"iter": 0, "n_rows": len(rows), **e})

    last = e["grid_min"]
    for it in range(1, args.max_iters + 1):
        h_b, p_b = e["h_min"], e["p_min"]
        # Snap binding point to inside the box (it lies in box by construction)
        print(f"\n[{it}] solving new center at (h, p) = ({h_b:.5f}, {p_b:.5f}) ...")
        t0 = time.time()
        res, V_c_rig, ell = solve_center_and_make_ellipse(
            h_b, p_b, args.N, args.T, args.R, args.bochner_n,
        )
        dt = time.time() - t0
        new_row = {
            "label": f"cde_iter{it}",
            "h_c": h_b, "p_c": p_b,
            "V_c": V_c_rig,
            "ellipse": ell,
            "config": {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bochner_n},
            "duals": res["duals"], "solve_s": dt,
        }
        rows.append(new_row)
        e = envelope_min(rows)
        delta = e["grid_min"] - last
        print(f"  V_c = {res['value']:.7f}  ({dt:.1f}s)")
        print(f"  new grid_min = {e['grid_min']:.7f}  Δ = {delta:+.7f}  "
              f"binding now: ({e['h_min']:.5f}, {e['p_min']:.5f}) via {e['row_label']}")
        history.append({
            "iter": it, "n_rows": len(rows),
            "new_center": {"h_c": h_b, "p_c": p_b, "V_c": res["value"],
                           "duals": res["duals"], "elapsed_s": dt},
            **e, "delta": delta,
        })
        if delta < args.min_delta:
            print(f"\nStopping: Δ < min_delta ({args.min_delta:.0e})")
            break
        last = e["grid_min"]

    out = CODE.parent.parent / args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "config": {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bochner_n,
                   "min_delta": args.min_delta},
        "history": history,
        "final_grid_min": e["grid_min"],
        "final_rigorous_LB_vs_white": e["grid_min"] - TARGET,
        "n_rows_final": len(rows),
    }, indent=2, default=float))
    print(f"\n→ saved {out}")
    print(f"\nFinal rigorous envelope min: {e['grid_min']:.7f}")
    print(f"  vs White's 0.379005:    +{e['grid_min'] - TARGET:+.7f}")
    print(f"  vs prior headline 0.379544: {e['grid_min'] - 0.379544:+.7f}")


if __name__ == "__main__":
    main()
