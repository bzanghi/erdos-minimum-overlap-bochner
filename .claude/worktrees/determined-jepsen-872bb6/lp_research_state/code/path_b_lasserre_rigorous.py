"""
Path B RIGOROUS aggregation for the Bochner+Lasserre runs.

Mirrors path_b_rigorous.py but reads /parallel_results/lasserre2_path_b/rowX.json,
applies the same V_c <- prob.value - margin treatment + vertex-closed-form min
+ Lipschitz-corrected envelope min on box (5.16).
"""
from __future__ import annotations
import sys, os, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from white_full_convex import WHITE_TABLE3
from path_b_analytical import find_ellipse_h_p
from path_b_rigorous import (
    build_rigorous_row, per_row_min_on_box,
    envelope_grad_lipschitz_bound, envelope_min_on_box,
)

WHITE_BOUND = 0.379005


def main(args):
    out_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results"
    if not os.path.isdir(out_dir):
        out_dir = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/parallel_results"
    row_dir = os.path.join(out_dir, "lasserre2_path_b")
    out_file = os.path.join(out_dir, "lasserre2_path_b_rigorous.json")

    margin = args.margin
    print(f"Rigour margin: {margin:.1e}")
    print(f"Box (5.16): h in [0, 0.06], p in [0.35, 0.45], q in [-0.02, 0.02]\n")

    rows_data = []
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        fp = os.path.join(row_dir, f"{label}.json")
        if not os.path.exists(fp):
            print(f"  MISSING {fp}; skipping")
            continue
        with open(fp) as fh:
            row_json = json.load(fh)
        new_row = build_rigorous_row(row_json, margin=margin)
        rows_data.append(new_row)

    print("=== (a) Per-row CLOSED-FORM minimum on box (5.16) ===")
    per_row_results = []
    for r in rows_data:
        pm = per_row_min_on_box(r["ellipse"])
        per_row_results.append(pm)
        print(f"  {r['label']:>5}: f_min = {pm['f_min']:.10f}  at "
              f"(h={pm['h_min']:.3f}, p={pm['p_min']:.3f})  concave={pm['concave']}")

    min_per_row = min(p["f_min"] for p in per_row_results)
    argmin_idx = int(np.argmin([p["f_min"] for p in per_row_results]))
    print(f"\n  MIN over 7 per-row minima = {min_per_row:.10f}  (row {rows_data[argmin_idx]['label']})")
    print(f"    vs White {WHITE_BOUND}: improvement = {min_per_row - WHITE_BOUND:+.6e}")

    print("\n=== (b) Best-of-7 ENVELOPE min on box (5.16) with Lipschitz LB ===")
    env = envelope_min_on_box(rows_data, n_grid=args.n_grid)
    print(f"  fine grid              : {env['n_grid']} x {env['n_grid']}")
    print(f"  envelope grid_min      : {env['grid_min']:.10f}  at (h={env['h_min']:.5f}, p={env['p_min']:.5f})")
    print(f"  witness row at min     : {rows_data[env['witness_row']]['label']}")
    print(f"  L_max_grad (over box)  : {env['L_max_grad']:.6f}")
    print(f"  half-cell diag         : {env['cell_half_diag']:.3e}")
    print(f"  eps_grid (Lipschitz)   : {env['eps_grid']:.3e}")
    print(f"  rigorous envelope LB   : {env['rigorous_envelope_min_LB']:.10f}")
    print(f"    vs White {WHITE_BOUND}: improvement = {env['rigorous_envelope_min_LB'] - WHITE_BOUND:+.6e}")

    # alt margins
    alt_margins = {}
    for alt in (1e-7, 1e-6, 5e-5):
        rows_alt = []
        for (h, p, qm, qp, label) in WHITE_TABLE3:
            fp = os.path.join(row_dir, f"{label}.json")
            if not os.path.exists(fp): continue
            with open(fp) as fh:
                row_json = json.load(fh)
            rows_alt.append(build_rigorous_row(row_json, margin=alt))
        env_alt = envelope_min_on_box(rows_alt, n_grid=args.n_grid)
        alt_margins[f"margin_{alt:.0e}"] = {
            "margin": alt,
            "rigorous_envelope_LB": env_alt["rigorous_envelope_min_LB"],
            "improvement_vs_white": env_alt["rigorous_envelope_min_LB"] - WHITE_BOUND,
            "grid_min": env_alt["grid_min"],
            "eps_grid": env_alt["eps_grid"],
        }

    summary = {
        "config": {
            "target": WHITE_BOUND,
            "h_box": [0.0, 0.06],
            "p_box": [0.35, 0.45],
            "n_grid_envelope": args.n_grid,
            "rigour_margin": margin,
            "augmentation": "Bochner_n=20 + Lasserre_T_max=25_T_loc=10",
        },
        "alt_margins": alt_margins,
        "rows": rows_data,
        "per_row_min_on_box": [
            {"label": rows_data[i]["label"],
             "f_min": per_row_results[i]["f_min"],
             "h_min": per_row_results[i]["h_min"],
             "p_min": per_row_results[i]["p_min"],
             "concave": per_row_results[i]["concave"],
             "vertex_values": per_row_results[i]["vertex_values"],
            } for i in range(len(rows_data))
        ],
        "min_over_per_row_mins": {
            "value": min_per_row,
            "row": rows_data[argmin_idx]["label"],
            "improvement_vs_white": min_per_row - WHITE_BOUND,
        },
        "envelope_min": env,
        "envelope_min_witness_row": rows_data[env["witness_row"]]["label"],
        "improvement_envelope_vs_white": env["rigorous_envelope_min_LB"] - WHITE_BOUND,
        "headline": {
            "MIN_envelope_rigorous_LB": env["rigorous_envelope_min_LB"],
            "improvement_over_white_0p379005": env["rigorous_envelope_min_LB"] - WHITE_BOUND,
        }
    }
    with open(out_file, "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "__float__") else str(o))
    print(f"\nResults written to {out_file}")
    return summary


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_grid", type=int, default=2001)
    parser.add_argument("--margin", type=float, default=1e-6)
    args = parser.parse_args()
    main(args)
