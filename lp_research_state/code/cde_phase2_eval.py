"""
Evaluate the Phase 2 cover (7 White at n=20 + 5 CDE at n=30) under uniform
margin convention. Outputs the rigorous LB matching the project's published
rigor convention.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import find_ellipse_h_p

H_BOX = (0.0, 0.06); P_BOX = (0.35, 0.45); TARGET = 0.379005; MARGIN = 1e-6


def main():
    rows = []
    for r in range(1, 8):
        d = json.load(open(CODE.parent / "parallel_results" / "path_b" / f"row{r}.json"))
        rows.append({"label": d["label"], "h_c": d["h_c"], "p_c": d["p_c"],
                     "V_c": d["primal_value_at_center"], "ellipse": d["ellipse"]})
    n30 = json.load(open(CODE.parent / "parallel_results" / "cde_iter_n30.json"))
    for h in n30["history"]:
        if "new_center" not in h: continue
        nc = h["new_center"]
        center = {"h_c": nc["h_c"], "p_c": nc["p_c"],
                  "q1": -0.02, "q2": 0.02, "value": nc["V_c"] - MARGIN}
        ell = find_ellipse_h_p(center, nc["duals"], -0.02, 0.02, target=TARGET)
        rows.append({"label": f"cde_n30_iter{h['iter']}",
                     "h_c": nc["h_c"], "p_c": nc["p_c"],
                     "V_c": nc["V_c"], "ellipse": ell})

    h_grid = np.linspace(*H_BOX, 4001); p_grid = np.linspace(*P_BOX, 4001)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf); witness = np.zeros_like(HH, dtype=int)
    for i, r in enumerate(rows):
        e = r["ellipse"]; V_c_rig = r["V_c"] - MARGIN
        F = (V_c_rig + e.get("const_q", 0) + e["A_h2"]*HH*HH + e["A_h1"]*HH
             + e["A_h0"] + e["A_p2"]*PP*PP + e["A_p1"]*PP + e["A_p0"])
        mask = F > env; env[mask] = F[mask]; witness[mask] = i

    grid_min = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    L_max = 0.0
    for r in rows:
        e = r["ellipse"]
        def lin_max_abs(c2, c1, lo, hi):
            return max(abs(2*c2*lo + c1), abs(2*c2*hi + c1))
        L = float(np.sqrt(lin_max_abs(e["A_h2"], e["A_h1"], *H_BOX)**2 + lin_max_abs(e["A_p2"], e["A_p1"], *P_BOX)**2))
        L_max = max(L_max, L)
    cell_h = (H_BOX[1]-H_BOX[0])/4000; cell_p = (P_BOX[1]-P_BOX[0])/4000
    eps_grid = L_max * 0.5 * float(np.sqrt(cell_h**2 + cell_p**2))
    rigorous_LB = grid_min - eps_grid

    print("=== Phase 2 (n=30 CDE + n=20 White) under uniform margin ===")
    print(f"  n rows:         {len(rows)} (7 White + 5 CDE@n=30)")
    print(f"  grid_min:       {grid_min:.7f}")
    print(f"  binding point:  ({float(HH[arg]):.5f}, {float(PP[arg]):.5f})")
    print(f"  witness:        {rows[int(witness[arg])]['label']}")
    print(f"  L_max:          {L_max:.4f}")
    print(f"  eps_grid:       {eps_grid:.2e}")
    print(f"  RIGOROUS LB:    µ ≥ {rigorous_LB:.7f}")
    print(f"  vs White:       +{rigorous_LB - TARGET:.7f}")
    print(f"  vs 0.379544:    {rigorous_LB - 0.379544:+.7f}")
    print(f"  vs Phase 1:     {rigorous_LB - 0.3796201:+.7f}")

    out = CODE.parent / "parallel_results" / "cde_phase2_rigorous.json"
    out.write_text(json.dumps({
        "phase": 2,
        "convention": {"margin": MARGIN, "h_box": H_BOX, "p_box": P_BOX,
                       "target_white": TARGET, "n_grid": 4001},
        "rigorous_LB": rigorous_LB, "grid_min": grid_min,
        "binding_point": [float(HH[arg]), float(PP[arg])],
        "witness": rows[int(witness[arg])]["label"],
        "L_max_grad": L_max, "eps_grid": eps_grid,
        "improvement_vs_white": rigorous_LB - TARGET,
        "improvement_vs_published": rigorous_LB - 0.379544,
        "improvement_vs_phase1": rigorous_LB - 0.3796201,
        "n_centers": len(rows),
    }, indent=2, default=float))
    print(f"\n→ saved {out}")


if __name__ == "__main__":
    main()
