"""Red-team check: can we get a tighter CONDITIONAL bound by composing
assume_even=True with the full Phase 5 augmentation stack
(bochner_n=30 + poly_moment k_max=14 + Hankel-PSD n=4)?

If yes AND the even-f conjecture is true (open), this would be a new lever.
"""
from __future__ import annotations
import sys, time
sys.path.insert(0, '/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code')
import warnings; warnings.filterwarnings('ignore')
import cvxpy as cp
import numpy as np
from white_full_convex import build_problem
from poly_moment import build_even_moment_nonneg_constraints, build_even_hankel_psd

def solve_one(label, **kwargs):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(**kwargs)
    # Add poly-moment if requested
    pm_k = kwargs.get('_pm_k_max', 0)
    hk_n = kwargs.get('_hankel_n', 0)
    T = kwargs['T']
    if pm_k > 0:
        pm_cons, tb = build_even_moment_nonneg_constraints(c, d, T, k_max=pm_k)
        cons.extend(pm_cons)
    if hk_n > 0:
        hk_cons, m_var, tails = build_even_hankel_psd(c, d, T, n_hankel=hk_n)
        cons.extend(hk_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver='CLARABEL', verbose=False)
    dt = time.time() - t0
    print(f"  {label}: status={prob.status}  Omega*={prob.value:.10f}  t={dt:.1f}s")
    return prob.value

# Row6: h=0, p=0.381, q ∈ [-0.02, 0.02]
PARAMS = dict(N=2000, T=1000, R=10, h1=0.0, h2=0.0, p1=0.381, p2=0.381, q1=-0.02, q2=0.02)

# Strip the _pm_k_max etc from build_problem call
def strip(kwargs):
    return {k: v for k, v in kwargs.items() if not k.startswith('_')}

print("=== Row 6 at small scale ===")

# baseline (no augmentation)
print("\n--- assume_even=False, no augment (baseline) ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(**PARAMS)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  baseline_uncond: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- assume_even=True, no augment ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(assume_even=True, **PARAMS)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  even_no_augment: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- assume_even=False, bochner_n=20 (unconditional reference) ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(bochner_n=20, **PARAMS)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  uncond+bochner_n=20: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- assume_even=True, bochner_n=20 ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(assume_even=True, bochner_n=20, **PARAMS)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  even+bochner_n=20: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- assume_even=True, bochner_n=20 + poly_moment k_max=10 ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(assume_even=True, bochner_n=20, **PARAMS)
pm_cons, tb = build_even_moment_nonneg_constraints(c, d, PARAMS['T'], k_max=10)
cons.extend(pm_cons)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  even+bochner+pm k=10: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- assume_even=True, bochner_n=20 + pm k=10 + hankel n=4 ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(assume_even=True, bochner_n=20, **PARAMS)
pm_cons, tb = build_even_moment_nonneg_constraints(c, d, PARAMS['T'], k_max=10)
cons.extend(pm_cons)
hk_cons, m_var, tails = build_even_hankel_psd(c, d, PARAMS['T'], n_hankel=4)
cons.extend(hk_cons)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  even+full_aug: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")

print("\n--- For comparison, UNCONDITIONAL same stack ---")
Omega, w, v, c, d, eps, dlt, cons = build_problem(assume_even=False, bochner_n=20, **PARAMS)
pm_cons, tb = build_even_moment_nonneg_constraints(c, d, PARAMS['T'], k_max=10)
cons.extend(pm_cons)
hk_cons, m_var, tails = build_even_hankel_psd(c, d, PARAMS['T'], n_hankel=4)
cons.extend(hk_cons)
prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time(); prob.solve(solver='CLARABEL'); print(f"  uncond+full_aug: {prob.value:.10f} ({prob.status}, {time.time()-t0:.1f}s)")
