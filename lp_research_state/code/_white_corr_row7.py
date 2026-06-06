"""
_white_corr_row7.py  (scratch; provenance only — do not import)

Measure White's 2026-05-31 email correction (mside_sin_coeff 8 -> 4 on the
imag/sine cell-consistency constraint 5.6/5.7) at the row7 CENTER, with coeff 8
and coeff 4 at an IDENTICAL light solver config so the only difference is the
coefficient.

row7 center: h = 0.03, p = 0.375, q in [-0.02, 0.02].
Light config: build_problem(2000, 800, 10, 0.03, 0.03, 0.375, 0.375,
                            -0.02, 0.02, bochner_n=20, mside_sin_coeff=coeff)

Report rigorous_dual_LB (preferred) and prob.value (backup) for each coeff,
and delta = LB(4) - LB(8).
"""
from __future__ import annotations
import sys, os, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cvxpy as cp
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction

# row7 center
N, T, R = 2000, 800, 10
h = 0.03
p = 0.375
q1, q2 = -0.02, 0.02
BN = 20

CONFIG = f"build_problem({N},{T},{R},{h},{h},{p},{p},{q1},{q2},bochner_n={BN},mside_sin_coeff=COEFF)"

results = {}
for coeff in (8.0, 4.0):
    Omega, w, v, c, d, eps_v, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2,
        bochner_n=BN, mside_sin_coeff=coeff,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    results[coeff] = res
    lb = res["rigorous_dual_LB"]
    val = res["reported_value"]
    print(f"[coeff={coeff}] status={res['status']}  "
          f"rigorous_dual_LB={lb!r}  value={val!r}  "
          f"resid@LB={res['dual_residual_at_LB']!r}  "
          f"iter={res['best_iter']}/{res['n_iters_total']}  t={res['time']:.1f}s")

lb8 = results[8.0]["rigorous_dual_LB"]
lb4 = results[4.0]["rigorous_dual_LB"]
v8 = results[8.0]["reported_value"]
v4 = results[4.0]["reported_value"]
delta = lb4 - lb8

print("\n==== SUMMARY (row7 center) ====")
print(f"config         = {CONFIG}")
print(f"center h,p     = {h}, {p}   q in [{q1},{q2}]")
print(f"lb_coeff8      = {lb8!r}")
print(f"lb_coeff4      = {lb4!r}")
print(f"value_coeff8   = {v8!r}")
print(f"value_coeff4   = {v4!r}")
print(f"delta(4 - 8)   = {delta!r}")

if abs(delta) <= 1e-6:
    interp = "SLACK/neutral: constraint not binding at this center; 8 and 4 give the same LB."
elif delta > 1e-6:
    interp = f"coeff 4 gives HIGHER LB by {delta:.3e}: old 8 was conservative/looser; 4 IMPROVES."
else:
    interp = f"coeff 4 gives LOWER LB by {abs(delta):.3e}: old 8 was an OVERCLAIM by that much."
print(f"interpretation = {interp}")
