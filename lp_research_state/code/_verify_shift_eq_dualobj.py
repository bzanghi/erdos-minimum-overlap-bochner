"""VERIFY (small N): cover shift reconstruction == exact change in dual objective
-b(theta)^T z for the SAME numeric dual z; and A,c are theta-independent (Lemma 10)."""
import sys, numpy as np, cvxpy as cp
sys.path.insert(0, '.')
from path_b_analytical import build_problem_with_dual_handles, dual_objective_shift

N,T,R,bn = 200, 80, 6, 6
h_c,p_c,q1,q2 = 0.004, 0.3875, -0.02, 0.02

Omega, cons, H = build_problem_with_dual_handles(N,T,R,h_c,h_c,p_c,p_c,q1,q2,bochner_n=bn)
prob = cp.Problem(cp.Minimize(Omega), cons)
prob.solve(solver=cp.CLARABEL, verbose=False)         # populates dual_value
# get canonical data + z (re-derive z from the solved problem)
data, chain, inv = prob.get_problem_data(cp.CLARABEL)
sol = chain.solve_via_data(prob, data)
z = np.asarray(sol.z); b = np.asarray(data['b']); A = data['A']; c = np.asarray(data['c'])
neg_bz_center = float(-(b @ z))
duals = {k: float(H[k].dual_value) for k in ("con_53","con_54","con_512_pL","con_512_pU","con_512_qL","con_512_qU","con_513")}
D = c + A.T @ z
print("center: prob.value=%.10f  obj_val=%.10f  -b^Tz=%.10f  defect_inf=%.2e" %
      (float(prob.value), float(sol.obj_val), neg_bz_center, np.max(np.abs(D))))
print("duals:", {k: round(v,8) for k,v in duals.items()})

def get_bAc(h,p,qq1,qq2):
    Om2, cons2, H2 = build_problem_with_dual_handles(N,T,R,h,h,p,p,qq1,qq2,bochner_n=bn)
    prob2 = cp.Problem(cp.Minimize(Om2), cons2)
    d2,_,_ = prob2.get_problem_data(cp.CLARABEL)
    return np.asarray(d2['b']), d2['A'], np.asarray(d2['c'])

worst_shift=0.0; worst_A=0.0; worst_c=0.0
rng=np.random.default_rng(0)
for _ in range(12):
    h = rng.uniform(0.0,0.06); p = rng.uniform(0.35,0.45); q = rng.uniform(-0.02,0.02)
    b2,A2,c2 = get_bAc(h,p,q,q)
    dA = abs((A2-A)).max() if A2.shape==A.shape else 9e9
    dc = float(np.max(np.abs(c2-c)))
    worst_A=max(worst_A,float(dA)); worst_c=max(worst_c,dc)
    neg_bz_pert = float(-(b2 @ z))
    exact_shift = neg_bz_pert - neg_bz_center
    recon_shift = dual_objective_shift(h,p,q,q,{"h_c":h_c,"p_c":p_c,"q1":q1,"q2":q2}, duals)
    worst_shift=max(worst_shift, abs(exact_shift-recon_shift))

print("\nLemma-10: max|A_theta - A_center| = %.2e (MUST ~0)" % worst_A)
print("Lemma-10: max|c_theta - c_center| = %.2e (MUST ~0)" % worst_c)
print("shift recon vs exact -b(theta)^Tz change: worst|Delta| = %.3e" % worst_shift)
print("\nVERDICT shift==exact:", bool(worst_shift<1e-9),
      "| A,c theta-indep:", bool(worst_A<1e-12 and worst_c<1e-12))
