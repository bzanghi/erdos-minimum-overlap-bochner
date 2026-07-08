"""
Scratch: measure White's 2026-05-31 coeff correction (8 -> 4) on the
imag/sine cell-consistency constraint (5.6/5.7) at center cde_n30_iter3.

Center cde_n30_iter3: h=0.00392, p=0.39225, d1 range q1=-0.02, q2=0.02.
Light config so each solve is fast: N=2000, T=800, R=10, bochner_n=20.
Compare mside_sin_coeff in [8.0, 4.0] at IDENTICAL config.

delta_coeff4_minus_coeff8 = lb_coeff4 - lb_coeff8
  |delta| <= ~1e-6  -> constraint effectively SLACK/neutral here
  delta > +1e-6      -> coeff 4 gives higher LB (8 was conservative; 4 improves)
  delta < -1e-6      -> coeff 4 gives lower  LB (8 was an overclaim by |delta|)
"""
import sys, os, warnings, json
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp

H = 0.00392
P = 0.39225
Q1, Q2 = -0.02, 0.02
N, T, R, BN = 2000, 800, 10, 20
CONFIG = f"build_problem({N},{T},{R},{H},{H},{P},{P},{Q1},{Q2},bochner_n={BN},mside_sin_coeff=<coeff>)"

results = {}
for coeff in (8.0, 4.0):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, H, H, P, P, Q1, Q2,
        bochner_n=BN, mside_sin_coeff=coeff,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    results[coeff] = {
        "rigorous_dual_LB": res["rigorous_dual_LB"],
        "value": res["reported_value"],
        "status": res["status"],
        "dual_residual_at_LB": res["dual_residual_at_LB"],
        "time": res["time"],
    }
    print(f"coeff={coeff}: status={res['status']} "
          f"LB={res['rigorous_dual_LB']} value={res['reported_value']} "
          f"resid={res['dual_residual_at_LB']} t={res['time']:.1f}s")

lb8 = results[8.0]["rigorous_dual_LB"]
lb4 = results[4.0]["rigorous_dual_LB"]
delta = lb4 - lb8
print("\n=== SUMMARY (JSON) ===")
print(json.dumps({
    "center": "cde_n30_iter3",
    "h": H, "p": P, "q1": Q1, "q2": Q2,
    "config": CONFIG,
    "lb_coeff8": lb8,
    "lb_coeff4": lb4,
    "value_coeff8": results[8.0]["value"],
    "value_coeff4": results[4.0]["value"],
    "delta_coeff4_minus_coeff8": delta,
}, indent=2))
