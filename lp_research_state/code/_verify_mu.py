"""
Independent verification of the headline rigorous LB  mu >= 0.3803027228
(Phase 5, N=20000, T=4000, bochner_n=40, pm_k_max=20).

Three checks that need NO SDP solve (fast, pure arithmetic on saved artifact):
  (A) Re-derive each ellipse from its saved duals via find_ellipse_h_p and
      compare to the recorded ellipse  -> JSON integrity / determinism.
  (B) Recompute grid_min, eps_grid, rigorous_LB from the ellipses
      -> confirms 0.3803027228.
  (C) Feed the SAME saved duals into the INDEPENDENT Phi_row reconstruction
      (path_b_independent) and recompute the box-min over a dense grid
      -> cross-checks the sign convention / linearization, second codebase.

Run:  ../../.venv/bin/python _verify_mu.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import find_ellipse_h_p
import path_b_independent as indep

ART = CODE.parent / "parallel_results" / "phase5_N20K_bn40.json"
H_BOX = (0.0, 0.06); P_BOX = (0.35, 0.45); TARGET = 0.379005
WHITE = 0.379005

data = json.load(open(ART))
centers = [c for c in data["centers"] if "error" not in c]
print(f"loaded {len(centers)} centers from {ART.name}")
print(f"recorded: rigorous_LB={data['rigorous_LB']:.10f}  grid_min={data['grid_min']:.10f}  "
      f"eps_grid={data['eps_grid']:.3e}\n")

# ---------- (A) re-derive ellipses from saved duals ----------
print("=== (A) re-derive ellipses from saved duals ===")
max_ell_err = 0.0
for c in centers:
    synthetic = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
                 "value": c["V_c_rigorous"]}
    ell = find_ellipse_h_p(synthetic, c["duals"], c["q1"], c["q2"], target=TARGET)
    for k in ("A_h2", "A_h1", "A_h0", "A_p2", "A_p1", "A_p0", "const_q"):
        max_ell_err = max(max_ell_err, abs(ell[k] - c["ellipse"][k]))
print(f"max |re-derived - recorded| ellipse coeff = {max_ell_err:.2e}  "
      f"({'OK' if max_ell_err < 1e-12 else 'MISMATCH'})\n")

# ---------- (B) recompute grid_min / eps_grid / rigorous_LB ----------
print("=== (B) recompute grid_min, eps_grid, rigorous_LB ===")
h_grid = np.linspace(*H_BOX, 4001); p_grid = np.linspace(*P_BOX, 4001)
HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
env = np.full_like(HH, -np.inf)
for c in centers:
    e = c["ellipse"]
    F = (c["V_c_rigorous"] + e.get("const_q", 0)
         + e["A_h2"]*HH*HH + e["A_h1"]*HH + e["A_h0"]
         + e["A_p2"]*PP*PP + e["A_p1"]*PP + e["A_p0"])
    np.maximum(env, F, out=env)
grid_min = float(env.min())
L_max = 0.0
for c in centers:
    e = c["ellipse"]
    lam = lambda c2, c1, lo, hi: max(abs(2*c2*lo+c1), abs(2*c2*hi+c1))
    L = float(np.hypot(lam(e["A_h2"], e["A_h1"], *H_BOX), lam(e["A_p2"], e["A_p1"], *P_BOX)))
    L_max = max(L_max, L)
cell_h = (H_BOX[1]-H_BOX[0])/4000; cell_p = (P_BOX[1]-P_BOX[0])/4000
eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
rig = grid_min - eps_grid
print(f"grid_min   = {grid_min:.10f}   (recorded {data['grid_min']:.10f}, "
      f"d={abs(grid_min-data['grid_min']):.2e})")
print(f"eps_grid   = {eps_grid:.3e}     (recorded {data['eps_grid']:.3e})")
print(f"rigorousLB = {rig:.10f}   (recorded {data['rigorous_LB']:.10f}, "
      f"d={abs(rig-data['rigorous_LB']):.2e})")
print(f"vs White 0.379005:  {rig-WHITE:+.3e}\n")

# ---------- (C) independent Phi_row cross-check ----------
print("=== (C) independent Phi_row (path_b_independent) box-min ===")
# map saved duals (con_*) -> independent's (lam_*) and build records
recs = []
for c in centers:
    d = c["duals"]
    recs.append({
        "label": c["label"],
        "value": c["V_c_rigorous"],
        "center": {"h_c": c["h_c"], "p_c": c["p_c"], "q1_c": c["q1"], "q2_c": c["q2"]},
        "duals": {"lam_53": d["con_53"], "lam_54": d["con_54"],
                  "lam_pL": d["con_512_pL"], "lam_pU": d["con_512_pU"],
                  "lam_qL": d["con_512_qL"], "lam_qU": d["con_512_qU"],
                  "lam_513": d["con_513"]},
    })
hg = np.linspace(*H_BOX, 601); pg = np.linspace(*P_BOX, 601)
qg = np.array([-0.02, 0.0, 0.02])
bmin, bloc, brow = indep.grid_min_vectorized(recs, hg, pg, qg)
print(f"independent box-min (Phi_row) = {bmin:.10f} at "
      f"(h={bloc[0]:.4f}, p={bloc[1]:.4f}, q={bloc[2]:.3f}) witness={brow}")
print(f"vs path_b grid_min {grid_min:.10f}  -> diff {abs(bmin-grid_min):.2e}")
print("(small positive diff expected: independent grid is coarser 601 vs 4001,\n"
      " and independent uses q^2 not max(q1^2,q2^2) in (5.13))")
