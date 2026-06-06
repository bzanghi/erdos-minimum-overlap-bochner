"""Push the T5p composition test to the actual Phase-5 config:
bochner_n=30 + poly_moment k_max=14 + hankel_n=4. Row 4 = binding row.
"""
import sys, time, warnings
sys.path.insert(0, '/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code')
warnings.filterwarnings('ignore')
import cvxpy as cp
import numpy as np
from white_full_convex import build_problem
from poly_moment import build_even_moment_nonneg_constraints, build_even_hankel_psd
from dual_extractor import solve_with_dual_extraction

# row 4 center (binding)
PARAMS = dict(N=3000, T=1200, R=10, h1=0.004, h2=0.004, p1=0.3875, p2=0.3875, q1=-0.02, q2=0.02)

def solve(label, use_T5p=False, bochner_n=0, pm_k=0, hk_n=0):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(use_T5p=use_T5p, bochner_n=bochner_n, **PARAMS)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(c, d, PARAMS['T'], k_max=pm_k); cons.extend(pm_cons)
    if hk_n > 0:
        hk_cons, _, _ = build_even_hankel_psd(c, d, PARAMS['T'], n_hankel=hk_n); cons.extend(hk_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time(); prob.solve(solver='CLARABEL'); dt = time.time() - t0
    print(f"  {label:50s}: {prob.value:.10f}  ({prob.status}, {dt:.1f}s)")
    return prob.value

print("=== Row 4 N=3000 T=1200, T5p composition w/ Phase-5 stack ===")
v1 = solve("bochner_n=20 + pm k=10 + hk n=4",          bochner_n=20, pm_k=10, hk_n=4)
v2 = solve("bochner_n=20 + pm k=10 + hk n=4 + T5p",    bochner_n=20, pm_k=10, hk_n=4, use_T5p=True)
v3 = solve("bochner_n=30 + pm k=14 + hk n=4",          bochner_n=30, pm_k=14, hk_n=4)
v4 = solve("bochner_n=30 + pm k=14 + hk n=4 + T5p",    bochner_n=30, pm_k=14, hk_n=4, use_T5p=True)
print(f"\nDelta T5p on (bochner_20+pm_10+hk_4): {v2-v1:+.3e}")
print(f"Delta T5p on (bochner_30+pm_14+hk_4): {v4-v3:+.3e}")
