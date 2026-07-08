"""Quick test: solve at one center with small N=600, check sign convention by perturbation."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from path_b_analytical import (
    solve_and_extract_duals, dual_objective_shift,
)

N, T, R = 800, 400, 10
bochner_n = 10  # smaller for speed
h_c, p_c = 0.015, 0.381
qm, qp = -0.02, 0.02

print("Center solve at row1 (h_c=0.015, p_c=0.381)...")
center = solve_and_extract_duals(N, T, R, h_c, p_c, qm, qp, bochner_n)
print(f"value = {center['value']:.7f}  status={center['status']}  time={center['time']:.1f}s")
print(f"duals = {center['duals']}")

# Perturbation tests: vary h, p, q individually
for label, (dh, dp, dq1, dq2) in [
    ("dh+5e-3", (5e-3, 0, 0, 0)),
    ("dh-5e-3", (-5e-3, 0, 0, 0)),
    ("dp+5e-3", (0, 5e-3, 0, 0)),
    ("dp-5e-3", (0, -5e-3, 0, 0)),
    ("dq+5e-3", (0, 0, 5e-3, 5e-3)),
]:
    pert = solve_and_extract_duals(N, T, R, h_c+dh, p_c+dp, qm+dq1, qp+dq2, bochner_n)
    pred = center['value'] + dual_objective_shift(h_c+dh, p_c+dp, qm+dq1, qp+dq2,
                                                   center, center['duals'])
    print(f"  {label}: true={pert['value']:.6f}  pred={pred:.6f}  err={pert['value']-pred:+.2e}")
