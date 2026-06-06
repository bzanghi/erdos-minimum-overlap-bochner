"""
Scratch (PRO / White-2026-05-31 correction): measure the effect of the
sine cell-consistency constraint RHS coefficient (5.6/5.7) at the row-4
ellipse CENTER. White's email said this RHS had an 8 in the numerator and
should be a 4.

We solve the SAME light-config SDP at row-4 center for coeff in {8.0, 4.0}
and compare rigorous_dual_LB. The ABSOLUTE values matter less than the
DIFFERENCE (coeff4 - coeff8) at IDENTICAL solver config.

Rigor logic: coeff 4 < 8 shrinks |rhs|, narrowing the two-sided sine band.
  - delta ~ 0   => constraint slack/neutral at this center.
  - delta > 0   => coeff 4 gives HIGHER LB (8 was conservative/looser; 4 improves).
  - delta < 0   => coeff 4 gives LOWER LB (8 was an OVERCLAIM by |delta|).
"""
import sys, os, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cvxpy as cp
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction

# Row 4 center: (h, p, q-range) = (0.004, 0.3875, [-0.02, 0.02])
N, T, R = 2000, 800, 10
h, p, q1, q2 = 0.004, 0.3875, -0.02, 0.02
BOCHNER_N = 20
CONFIG = (f"build_problem(N={N}, T={T}, R={R}, h={h}, p={p}, "
          f"q=[{q1},{q2}], bochner_n={BOCHNER_N}, exact cells)")

results = {}
for coeff in (8.0, 4.0):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2,
        bochner_n=BOCHNER_N, mside_sin_coeff=coeff,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    results[coeff] = {
        "status": res["status"],
        "rigorous_dual_LB": res["rigorous_dual_LB"],
        "reported_value": res["reported_value"],
        "dual_residual_at_LB": res["dual_residual_at_LB"],
        "best_iter": res["best_iter"],
        "n_iters_total": res["n_iters_total"],
        "time": res["time"],
    }
    print(f"[coeff={coeff}] status={res['status']} "
          f"LB={res['rigorous_dual_LB']} value={res['reported_value']} "
          f"resid={res['dual_residual_at_LB']} t={res['time']:.1f}s")

lb8 = results[8.0]["rigorous_dual_LB"]
lb4 = results[4.0]["rigorous_dual_LB"]
v8 = results[8.0]["reported_value"]
v4 = results[4.0]["reported_value"]
delta = (lb4 - lb8) if (lb4 is not None and lb8 is not None) else None

summary = {
    "center": "row4 (h=0.004, p=0.3875, q in [-0.02,0.02])",
    "h": h, "p": p, "config": CONFIG,
    "lb_coeff8": lb8, "lb_coeff4": lb4,
    "value_coeff8": v8, "value_coeff4": v4,
    "delta_coeff4_minus_coeff8": delta,
}
print("\n=== SUMMARY (JSON) ===")
print(json.dumps(summary, indent=2))
