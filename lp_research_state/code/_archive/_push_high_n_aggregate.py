"""Aggregate the (n=20, T_max=40) push_high_n results with rigorous dual LB.

For each row, V_c = rigorous_dual_LB (parsed from CLARABEL verbose output) which
is a true lower bound on the LP optimum at the row center.  Then run Path B
ellipse-extension and the closed-form vertex min + Lipschitz envelope LB.
"""
from __future__ import annotations
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from white_full_convex import WHITE_TABLE3
from path_b_analytical import find_ellipse_h_p
from path_b_rigorous import (
    per_row_min_on_box,
    envelope_grad_lipschitz_bound,
    envelope_min_on_box,
)

WHITE_BOUND = 0.379005
ROW_DIR = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/push_high_n"
OUT_FILE = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/push_high_n_summary.json"


def build_row_with_dual_LB(row_json):
    """V_c = rigorous_dual_LB (parsed CLARABEL dual_obj is a true LP LB)."""
    label = row_json["label"]
    h_c = row_json["h_c"]; p_c = row_json["p_c"]
    q1 = row_json["q1"]; q2 = row_json["q2"]
    duals = row_json["duals"]
    V_dual = row_json["rigorous_dual_LB"]
    V_primal = row_json["primal_value_at_center"]
    synthetic_center = {
        "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2, "value": V_dual,
    }
    ell = find_ellipse_h_p(synthetic_center, duals, q1, q2, target=WHITE_BOUND)
    return {
        "label": label,
        "config": row_json["config"],
        "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "V_primal": V_primal,
        "V_dual_LB": V_dual,
        "dual_residual_at_LB": row_json.get("dual_residual_at_LB"),
        "best_iter": row_json.get("best_iter"),
        "n_iters_total": row_json.get("n_iters_total"),
        "status": row_json.get("status"),
        "duals": duals,
        "ellipse": {
            "V_c": V_dual,
            "semi_h": ell["semi_h"], "semi_p": ell["semi_p"],
            "h_star": ell["h_star"], "p_star": ell["p_star"], "V_max": ell["V_max"],
            "A_h2": ell["A_h2"], "A_h1": ell["A_h1"], "A_h0": ell["A_h0"],
            "A_p2": ell["A_p2"], "A_p1": ell["A_p1"], "A_p0": ell["A_p0"],
            "const_q": ell["const_q"], "target": ell["target"],
        },
    }


def main():
    rows_data = []
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        fp = os.path.join(ROW_DIR, f"{label}.json")
        if not os.path.exists(fp):
            print(f"MISSING: {fp}")
            continue
        with open(fp) as fh:
            row_json = json.load(fh)
        rows_data.append(build_row_with_dual_LB(row_json))

    print("=== Per-row vertex-closed-form min on box (5.16) ===")
    per_row_results = []
    for r in rows_data:
        pm = per_row_min_on_box(r["ellipse"])
        per_row_results.append(pm)
        print(f"  {r['label']:>5}: V_dual_LB={r['V_dual_LB']:.7f}  "
              f"f_min_box={pm['f_min']:.10f}  at (h={pm['h_min']:.3f}, "
              f"p={pm['p_min']:.3f})  concave={pm['concave']}")

    min_per_row = min(p["f_min"] for p in per_row_results)
    argmin_idx = int(np.argmin([p["f_min"] for p in per_row_results]))
    print(f"\n  MIN over per-row mins: {min_per_row:.10f}  (row {rows_data[argmin_idx]['label']})")
    print(f"    vs White {WHITE_BOUND}: improvement = {min_per_row - WHITE_BOUND:+.6e}")

    print("\n=== Best-of-7 envelope min on box, with Lipschitz LB ===")
    env = envelope_min_on_box(rows_data, n_grid=2001)
    print(f"  grid: {env['n_grid']}x{env['n_grid']}")
    print(f"  grid_min = {env['grid_min']:.10f} at "
          f"(h={env['h_min']:.5f}, p={env['p_min']:.5f}); "
          f"witness = {rows_data[env['witness_row']]['label']}")
    print(f"  L_max_grad = {env['L_max_grad']:.6f}")
    print(f"  eps_grid (Lipschitz) = {env['eps_grid']:.3e}")
    print(f"  rigorous envelope LB = {env['rigorous_envelope_min_LB']:.10f}")
    print(f"    vs White {WHITE_BOUND}: improvement = "
          f"{env['rigorous_envelope_min_LB'] - WHITE_BOUND:+.6e}")

    summary = {
        "config": {
            "augmentation": "Bochner_n=20 + Lasserre_T_max=40, T_loc=8",
            "N": 10000, "T": 4000, "R": 10,
            "method": "split-solve via cvxpy.get_problem_data + chain.solve_via_data",
            "target": WHITE_BOUND,
            "h_box": [0.0, 0.06],
            "p_box": [0.35, 0.45],
            "V_c_choice": "rigorous_dual_LB (parsed CLARABEL dual_obj — a true LP LB)",
        },
        "rows": rows_data,
        "per_row_min_on_box": [
            {"label": rows_data[i]["label"], **per_row_results[i]}
            for i in range(len(rows_data))
        ],
        "min_over_per_row_mins": {
            "value": min_per_row,
            "row": rows_data[argmin_idx]["label"],
            "improvement_vs_white": min_per_row - WHITE_BOUND,
        },
        "envelope_min": env,
        "envelope_min_witness_row": rows_data[env["witness_row"]]["label"],
        "headline": {
            "MIN_envelope_rigorous_LB": env["rigorous_envelope_min_LB"],
            "improvement_over_white_0p379005": env["rigorous_envelope_min_LB"] - WHITE_BOUND,
            "previous_best_rigorous_LB_n20_T30": 0.3798283157774792,
            "improvement_vs_previous": (env["rigorous_envelope_min_LB"]
                                         - 0.3798283157774792),
        }
    }
    with open(OUT_FILE, "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "__float__") else str(o))
    print(f"\nWritten to {OUT_FILE}")
    return summary


if __name__ == "__main__":
    main()
