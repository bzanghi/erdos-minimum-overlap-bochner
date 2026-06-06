"""Sanity-check: cross-check the closed-form envelope minimum against a fine-grid
scan, and inspect the row4-row5 KKT crossing analytically."""
import json
import numpy as np
import os, sys

OUT = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results"
if not os.path.isdir(OUT):
    OUT = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/parallel_results"

margin = 1e-6
rows = []
for i in range(1, 8):
    with open(f"{OUT}/lasserre2_path_b/row{i}.json") as f:
        d = json.load(f)
    e = d["ellipse"].copy()
    e["V_c"] = d["primal_value_at_center"] - margin
    rows.append({"label": d["label"], "ellipse": e})

def Phi(e, h, p):
    return (e["V_c"] + e["const_q"]
            + e["A_h2"]*h*h + e["A_h1"]*h + e["A_h0"]
            + e["A_p2"]*p*p + e["A_p1"]*p + e["A_p0"])

# Very fine grid for cross-check
h_grid = np.linspace(0.0, 0.06, 3001)
p_grid = np.linspace(0.35, 0.45, 4001)
H, P = np.meshgrid(h_grid, p_grid, indexing='ij')
env = -np.inf*np.ones_like(H)
witness = np.zeros_like(H, dtype=int)
for i, r in enumerate(rows):
    F = Phi(r["ellipse"], H, P)
    mask = F > env
    env[mask] = F[mask]
    witness[mask] = i

idx = np.unravel_index(np.argmin(env), env.shape)
print(f"  fine grid 3001x4001 min E = {env[idx]:.10f}")
print(f"  argmin h = {H[idx]:.6f}, p = {P[idx]:.6f}, witness = row{witness[idx]+1}")

# Analytical row4-row5 binding crossover, with h fixed at 0
e4 = rows[3]["ellipse"]
e5 = rows[4]["ellipse"]
h0 = 0.0

# Phi_4(0, p) - Phi_5(0, p) = 0
# = (V_c4-V_c5) + (A_h0_4-A_h0_5) + (A_p2_4-A_p2_5) p^2 + (A_p1_4-A_p1_5) p + (A_p0_4-A_p0_5)
# (no h2/h1 since h=0)
A2 = e4["A_p2"] - e5["A_p2"]
A1 = e4["A_p1"] - e5["A_p1"]
A0 = (e4["V_c"] + e4["const_q"] + e4["A_h0"] + e4["A_p0"]) - (e5["V_c"] + e5["const_q"] + e5["A_h0"] + e5["A_p0"])
disc = A1**2 - 4*A2*A0
if disc >= 0:
    sd = np.sqrt(disc)
    p1 = (-A1 + sd) / (2*A2)
    p2 = (-A1 - sd) / (2*A2)
    print(f"\n  Phi_4 = Phi_5 at h=0:  p in {{{p1:.10f}, {p2:.10f}}}")
    for p_ in (p1, p2):
        if 0.35 <= p_ <= 0.45:
            v4 = Phi(e4, 0.0, p_)
            v5 = Phi(e5, 0.0, p_)
            env_v = max(Phi(r["ellipse"], 0.0, p_) for r in rows)
            print(f"    p = {p_:.10f}:  Phi_4 = {v4:.10f},  Phi_5 = {v5:.10f},  envelope = {env_v:.10f}")

# Check the closed-form output explicitly
import importlib.util
candidate_paths = [
    "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code/path_b_closed_form.py",
    "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code/path_b_closed_form.py",
]
mod_path = next(p for p in candidate_paths if os.path.exists(p))
spec = importlib.util.spec_from_file_location("path_b_closed_form", mod_path)
mod = importlib.util.module_from_spec(spec)
sys.path.insert(0, os.path.dirname(mod_path))
spec.loader.exec_module(mod)

# Re-run closed-form with margin=0 (no rigour margin) to see "pure" min
print("\n--- Same analysis with margin = 0 (no rigour margin) ---")
rows_pure = []
for i in range(1, 8):
    with open(f"{OUT}/lasserre2_path_b/row{i}.json") as f:
        d = json.load(f)
    e = d["ellipse"].copy()
    e["V_c"] = d["primal_value_at_center"]
    rows_pure.append({"label": d["label"], "ellipse": e})

res = mod.closed_form_envelope_min(rows_pure, verbose=False)
print(f"  margin=0 closed-form min = {res['min_E']:.10f}")
a = res["argmin"]
print(f"    type: {a['type']}, h={a['h']:.10f}, p={a['p']:.10f}, witness={a['witness_row']}")
if 'pair' in a: print(f"    pair: {a['pair']}")

# Try margin=1e-7 (tightest defensible)
print("\n--- With margin = 1e-7 (tight, only IPM gap) ---")
rows_tight = []
for i in range(1, 8):
    with open(f"{OUT}/lasserre2_path_b/row{i}.json") as f:
        d = json.load(f)
    e = d["ellipse"].copy()
    e["V_c"] = d["primal_value_at_center"] - 1e-7
    rows_tight.append({"label": d["label"], "ellipse": e})

res = mod.closed_form_envelope_min(rows_tight, verbose=False)
print(f"  margin=1e-7 closed-form min = {res['min_E']:.10f}")
