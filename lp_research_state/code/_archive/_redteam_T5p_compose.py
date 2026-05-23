"""Did anyone test use_T5p composed with bochner_n=30 + poly_moment + Hankel?
T5p is described in white_full_convex.py as 'biggest impact' but was first-gen
before Bochner was added. Maybe it's subsumed; maybe not — let's check."""
import sys, time, warnings
sys.path.insert(0, '/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code')
warnings.filterwarnings('ignore')
import cvxpy as cp
import numpy as np
from white_full_convex import build_problem
from poly_moment import build_even_moment_nonneg_constraints, build_even_hankel_psd

PARAMS = dict(N=2000, T=1000, R=10, h1=0.0, h2=0.0, p1=0.381, p2=0.381, q1=-0.02, q2=0.02)

def solve(label, **kwargs):
    pm_k = kwargs.pop('_pm_k', 0); hk_n = kwargs.pop('_hk_n', 0)
    Omega, w, v, c, d, eps, dlt, cons = build_problem(**{**PARAMS, **kwargs})
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(c, d, PARAMS['T'], k_max=pm_k); cons.extend(pm_cons)
    if hk_n > 0:
        hk_cons, _, _ = build_even_hankel_psd(c, d, PARAMS['T'], n_hankel=hk_n); cons.extend(hk_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time(); prob.solve(solver='CLARABEL'); dt = time.time() - t0
    print(f"  {label:45s}: {prob.value:.10f}  ({prob.status}, {dt:.1f}s)")
    return prob.value

print("=== Row 6 N=2000 T=1000, T5p composition test ===")
v0 = solve("baseline", )
v_t5p = solve("+ T5p", use_T5p=True)
v_b = solve("+ bochner_n=20", bochner_n=20)
v_b_t5p = solve("+ bochner_n=20 + T5p", bochner_n=20, use_T5p=True)
v_full = solve("+ bochner_n=20 + pm k=10 + hk n=4", bochner_n=20, _pm_k=10, _hk_n=4)
v_full_t5p = solve("+ bochner_n=20 + T5p + pm k=10 + hk n=4", bochner_n=20, use_T5p=True, _pm_k=10, _hk_n=4)

print(f"\nDelta T5p alone vs baseline: {v_t5p - v0:+.3e}")
print(f"Delta T5p in presence of Bochner_n=20: {v_b_t5p - v_b:+.3e}")
print(f"Delta T5p in presence of full Phase-5 stack: {v_full_t5p - v_full:+.3e}")
