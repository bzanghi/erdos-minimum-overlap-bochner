"""PRO-38 R9 probe: one augmented dual-feasible center at the strip worst corner.
Memory-safe config N=5000,T=2000,R=10,bn=30,pm_k=20. Corrected mside coeff=4.0 is
baked into build_problem_with_dual_handles (5.6/5.7 RHS uses 4.0). Conservative
anchor = dual_LB - 1e-5.  Solve pattern mirrors _fullspace_stage2_halo_centers.solve_center.
"""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import cvxpy as cp
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction

NEED = ("con_53","con_54","con_512_pL","con_512_pU","con_512_qL","con_512_qU","con_513")

def solve_center(h_c,p_c,q1,q2,N=5000,T=2000,R=10,bn=30,pm_k=20):
    Omega,cons,H = build_problem_with_dual_handles(N,T,R,h_c,h_c,p_c,p_c,q1,q2,bochner_n=bn)
    if pm_k>0:
        pm_cons,_ = build_even_moment_nonneg_constraints(H["c"],H["d"],T,k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0=time.time(); res = solve_with_dual_extraction(prob); dt=time.time()-t0
    duals = {k:(float(H[k].dual_value) if H[k].dual_value is not None else 0.0) for k in NEED}
    return res, duals, dt

if __name__=="__main__":
    h_c,p_c,q1,q2 = 0.004, 0.39, 0.025, 0.05
    res,duals,dt = solve_center(h_c,p_c,q1,q2)
    print(f"PROBE center h={h_c} p={p_c} q=[{q1},{q2}]")
    print(f"  status={res['status']} primal={res['reported_value']:.7f} dual_LB={res['rigorous_dual_LB']}")
    print(f"  dual_resid={res['dual_residual_at_LB']:.3e} time={dt:.0f}s")
    print(f"  anchor(dualLB-1e-5)={res['rigorous_dual_LB']-1e-5:.7f}")
    print(f"  con_513={duals['con_513']:.4f} con_53={duals['con_53']:.4f} con_54={duals['con_54']:.4f}")
