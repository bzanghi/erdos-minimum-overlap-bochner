"""
Compute the rigorous envelope LB applying the SAME margin convention used in
path_b_rigorous.json (margin=1e-6 on every V_c + Lipschitz cell-error bound).

Compares baseline (7 rows) vs augmented (7 + iterated centers) under identical
rigor convention.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

H_BOX = (0.0, 0.06)
P_BOX = (0.35, 0.45)
TARGET = 0.379005
MARGIN = 1e-6


def lipschitz_grad_bound(rows):
    """Max over rows of upper bound on |grad f_r| over the box. Linear gradient
    of A_h2 h^2 + A_h1 h + ... is 2 A_h2 h + A_h1; max over box."""
    def lin_max_abs(c2, c1, lo, hi):
        return max(abs(2 * c2 * lo + c1), abs(2 * c2 * hi + c1))
    L_max = 0.0
    for r in rows:
        e = r["ellipse"]
        gh = lin_max_abs(e["A_h2"], e["A_h1"], H_BOX[0], H_BOX[1])
        gp = lin_max_abs(e["A_p2"], e["A_p1"], P_BOX[0], P_BOX[1])
        L_max = max(L_max, float(np.sqrt(gh * gh + gp * gp)))
    return L_max


def envelope_min_rigorous(rows, n_grid=4001, margin=MARGIN):
    """Envelope min with V_c -> V_c - margin uniformly + Lipschitz cell error."""
    h_grid = np.linspace(H_BOX[0], H_BOX[1], n_grid)
    p_grid = np.linspace(P_BOX[0], P_BOX[1], n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf, dtype=float)
    witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(rows):
        e = r["ellipse"]
        V_c_rig = r["V_c"] - margin
        F = (V_c_rig + e.get("const_q", 0.0)
             + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        mask = F > env
        env[mask] = F[mask]; witness[mask] = i
    grid_min = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    cell_h = (H_BOX[1] - H_BOX[0]) / (n_grid - 1)
    cell_p = (P_BOX[1] - P_BOX[0]) / (n_grid - 1)
    half_diag = 0.5 * np.sqrt(cell_h ** 2 + cell_p ** 2)
    L_max = lipschitz_grad_bound(rows)
    eps_grid = L_max * half_diag
    rigorous_LB = grid_min - eps_grid
    return {
        "grid_min_after_margin": grid_min,
        "h_min": float(HH[arg]), "p_min": float(PP[arg]),
        "witness_idx": int(witness[arg]),
        "witness_label": rows[int(witness[arg])]["label"],
        "L_max_grad": L_max,
        "eps_grid": eps_grid,
        "rigorous_LB": rigorous_LB,
        "improvement_vs_white": rigorous_LB - TARGET,
    }


def load_existing_rows():
    base = CODE.parent / "parallel_results" / "path_b"
    rows = []
    for r in range(1, 8):
        d = json.load(open(base / f"row{r}.json"))
        rows.append({"label": d["label"], "h_c": d["h_c"], "p_c": d["p_c"],
                     "V_c": d["primal_value_at_center"], "ellipse": d["ellipse"]})
    return rows


def load_iterated_rows():
    """Pull the new centers from cde_iterative.json output."""
    iter_file = CODE.parent / "parallel_results" / "cde_iterative.json"
    if not iter_file.exists():
        return []
    data = json.load(open(iter_file))
    rows = []
    for h in data["history"]:
        if "new_center" not in h: continue
        # Reconstruct ellipse via find_ellipse_h_p
        from path_b_analytical import find_ellipse_h_p
        nc = h["new_center"]
        center = {"h_c": nc["h_c"], "p_c": nc["p_c"],
                  "q1": -0.02, "q2": 0.02, "value": nc["V_c"] - MARGIN}
        ell = find_ellipse_h_p(center, nc["duals"], -0.02, 0.02, target=TARGET)
        rows.append({"label": f"cde_iter{h['iter']}",
                     "h_c": nc["h_c"], "p_c": nc["p_c"],
                     "V_c": nc["V_c"], "ellipse": ell})
    return rows


def main():
    base = load_existing_rows()
    new = load_iterated_rows()
    print("=" * 70)
    print("Rigorous envelope LB — uniform margin + Lipschitz grid error")
    print("=" * 70)

    print(f"\n[BASELINE] 7 White rows only:")
    e0 = envelope_min_rigorous(base)
    for k, v in e0.items(): print(f"  {k}: {v}")

    print(f"\n[AUGMENTED] 7 White rows + {len(new)} CDE-discovered centers:")
    e1 = envelope_min_rigorous(base + new)
    for k, v in e1.items(): print(f"  {k}: {v}")

    print(f"\n=== Delta ===")
    print(f"  Δrigorous_LB = {e1['rigorous_LB'] - e0['rigorous_LB']:+.7f}")
    print(f"  baseline:  µ ≥ {e0['rigorous_LB']:.7f}  (+{e0['improvement_vs_white']:.7f} vs White)")
    print(f"  augmented: µ ≥ {e1['rigorous_LB']:.7f}  (+{e1['improvement_vs_white']:.7f} vs White)")
    print(f"  vs prior headline 0.379544: {e1['rigorous_LB'] - 0.379544:+.7f}")

    out = CODE.parent / "parallel_results" / "cde_rigorous.json"
    out.write_text(json.dumps({
        "convention": {"margin": MARGIN, "h_box": H_BOX, "p_box": P_BOX,
                       "target_white": TARGET, "n_grid": 4001},
        "baseline_7rows": e0,
        "augmented_with_cde": e1,
        "delta_rigorous_LB": e1["rigorous_LB"] - e0["rigorous_LB"],
        "n_cde_centers": len(new),
        "cde_centers": [{"label": r["label"], "h_c": r["h_c"], "p_c": r["p_c"],
                         "V_c": r["V_c"]} for r in new],
    }, indent=2, default=float))
    print(f"\n→ saved {out}")


if __name__ == "__main__":
    main()
