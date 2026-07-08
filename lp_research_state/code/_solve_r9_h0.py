"""Solve fresh R9 centers at h=0 (worst h-edge), c1 near worst, q across strip.
N=8000 for high anchor. Cache to _r9_probe_cache.json."""
from __future__ import annotations
import json,sys,time,warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import cvxpy as cp
CODE=Path(__file__).resolve().parent; sys.path.insert(0,str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
NEED=("con_53","con_54","con_512_pL","con_512_pU","con_512_qL","con_512_qU","con_513")
CACHE=CODE.parent/"parallel_results"/"_r9_probe_cache.json"
def solve_center(h_c,p_c,q1,q2,N,T,R,bn,pm_k):
    Omega,cons,H=build_problem_with_dual_handles(N,T,R,h_c,h_c,p_c,p_c,q1,q2,bochner_n=bn)
    if pm_k>0:
        pm,_=build_even_moment_nonneg_constraints(H["c"],H["d"],T,k_max=pm_k); cons.extend(pm)
    prob=cp.Problem(cp.Minimize(Omega),cons)
    t0=time.time(); res=solve_with_dual_extraction(prob); dt=time.time()-t0
    duals={k:(float(H[k].dual_value) if H[k].dual_value is not None else 0.0) for k in NEED}
    return res,duals,dt
if __name__=="__main__":
    import ast
    specs=ast.literal_eval(sys.argv[1]); N=int(sys.argv[2]); bn=int(sys.argv[3]); T=int(sys.argv[4]) if len(sys.argv)>4 else 2000
    cache=json.load(open(CACHE)) if CACHE.exists() else {}
    for (h_c,p_c,q1,q2) in specs:
        key=f"{h_c}_{p_c}_{q1}_{q2}_{N}_{bn}"
        if key in cache:
            print(f"  cached {key}: dualLB={cache[key]['dual_lb']:.6f}"); continue
        res,duals,dt=solve_center(h_c,p_c,q1,q2,N,T,10,bn,20)
        cache[key]={"label":f"r9_h{h_c}_c{p_c}_q{q1}_{q2}_N{N}bn{bn}","h_c":h_c,"p_c":p_c,"q1":q1,"q2":q2,
            "primal":res["reported_value"],"dual_lb":res["rigorous_dual_LB"],"dual_resid":res["dual_residual_at_LB"],
            "status":res["status"],"duals":duals,"time":dt,"config":{"N":N,"T":T,"bn":bn,"pm_k":20}}
        CACHE.write_text(json.dumps(cache,indent=2,default=float))
        print(f"  solved {key}: primal={res['reported_value']:.6f} dualLB={res['rigorous_dual_LB']:.6f} resid={res['dual_residual_at_LB']:.1e} con53={duals['con_53']:.3f} con513={duals['con_513']:.3f} ({dt:.0f}s)",flush=True)
