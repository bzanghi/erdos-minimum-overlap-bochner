"""R9 probe B: sweep q-range tightness & N to find a center whose Phi clears the
strip p in [0.33,0.45], q in [0.025,0.05]. Key: a center at q=[0.025,0.025] (single
pt, small |q|) has higher primal; its reconstruction at query q=0.05 incurs the
(5.13)+(5.12q) q-penalty. We test whether anchor minus that penalty still > target."""
from __future__ import annotations
import sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import cvxpy as cp
CODE = Path(__file__).resolve().parent; sys.path.insert(0,str(CODE))
from path_b_analytical import build_problem_with_dual_handles, dual_objective_shift
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
NEED=("con_53","con_54","con_512_pL","con_512_pU","con_512_qL","con_512_qU","con_513")

def solve_center(h_c,p_c,q1,q2,N,T,R,bn,pm_k):
    Omega,cons,H = build_problem_with_dual_handles(N,T,R,h_c,h_c,p_c,p_c,q1,q2,bochner_n=bn)
    if pm_k>0:
        pm_cons,_=build_even_moment_nonneg_constraints(H["c"],H["d"],T,k_max=pm_k); cons.extend(pm_cons)
    prob=cp.Problem(cp.Minimize(Omega),cons)
    t0=time.time(); res=solve_with_dual_extraction(prob); dt=time.time()-t0
    duals={k:(float(H[k].dual_value) if H[k].dual_value is not None else 0.0) for k in NEED}
    return res,duals,dt

if __name__=="__main__":
    # candidate configs: (h_c,p_c,q1,q2,N,T,bn)
    cfgs = [
        (0.004, 0.39, 0.025, 0.025, 5000, 2000, 30),  # single-pt small q
        (0.004, 0.39, 0.0375, 0.0375, 5000, 2000, 30), # single-pt mid q
        (0.004, 0.39, 0.05, 0.05, 5000, 2000, 30),     # single-pt at worst q edge
    ]
    center={"h_c":0.004,"p_c":0.39}
    for (h_c,p_c,q1,q2,N,T,bn) in cfgs:
        res,duals,dt = solve_center(h_c,p_c,q1,q2,N,T,10,bn,20)
        anchor = res['rigorous_dual_LB']-1e-5
        c = {"h_c":h_c,"p_c":p_c,"q1":q1,"q2":q2}
        # reconstruct Phi at the strip's worst query: q=0.05, p=0.39, h=0.004
        phi_worst = anchor + dual_objective_shift(0.004,0.39,0.05,0.05,c,duals)
        # also at q=0.025 (other strip edge)
        phi_lo = anchor + dual_objective_shift(0.004,0.39,0.025,0.025,c,duals)
        print(f"q=[{q1},{q2}] N={N} bn={bn}: status={res['status']} primal={res['reported_value']:.6f} "
              f"dualLB={res['rigorous_dual_LB']:.6f} resid={res['dual_residual_at_LB']:.1e} ({dt:.0f}s)")
        print(f"    anchor={anchor:.6f}  Phi@(q=0.05)={phi_worst:.6f}  Phi@(q=0.025)={phi_lo:.6f}  con_513={duals['con_513']:.4f}")
