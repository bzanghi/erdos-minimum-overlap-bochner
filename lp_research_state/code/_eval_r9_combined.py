"""Evaluate combined cover (corehalo + candidate fresh centers) over the R9 strip.
Fresh centers are passed inline (already solved in probes). We re-solve them here
once and cache to JSON so the eval is reproducible."""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np, cvxpy as cp
CODE=Path(__file__).resolve().parent; sys.path.insert(0,str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
import _fs_recon as FS
from _fullspace_eval import cover_min_over_box, CORE_HEADLINE
NEED=("con_53","con_54","con_512_pL","con_512_pU","con_512_qL","con_512_qU","con_513")

def solve_center(h_c,p_c,q1,q2,N,T,R,bn,pm_k):
    Omega,cons,H=build_problem_with_dual_handles(N,T,R,h_c,h_c,p_c,p_c,q1,q2,bochner_n=bn)
    if pm_k>0:
        pm,_=build_even_moment_nonneg_constraints(H["c"],H["d"],T,k_max=pm_k); cons.extend(pm)
    prob=cp.Problem(cp.Minimize(Omega),cons)
    t0=time.time(); res=solve_with_dual_extraction(prob); dt=time.time()-t0
    duals={k:(float(H[k].dual_value) if H[k].dual_value is not None else 0.0) for k in NEED}
    return res,duals,dt

CACHE = CODE.parent/"parallel_results"/"_r9_probe_cache.json"

def get_centers(specs, N=5000,T=2000,bn=30):
    cache = json.load(open(CACHE)) if CACHE.exists() else {}
    out=[]
    for (h_c,p_c,q1,q2) in specs:
        key=f"{h_c}_{p_c}_{q1}_{q2}_{N}_{bn}"
        if key not in cache:
            res,duals,dt=solve_center(h_c,p_c,q1,q2,N,T,10,bn,20)
            cache[key]={"label":f"r9_h{h_c}_c{p_c}_q{q1}_{q2}_N{N}","h_c":h_c,"p_c":p_c,
                "q1":q1,"q2":q2,"primal":res["reported_value"],"dual_lb":res["rigorous_dual_LB"],
                "dual_resid":res["dual_residual_at_LB"],"status":res["status"],"duals":duals,"time":dt,"config":{"N":N,"T":T,"bn":bn,"pm_k":20}}
            CACHE.write_text(json.dumps(cache,indent=2,default=float))
            print(f"  solved {key}: dualLB={res['rigorous_dual_LB']:.6f} resid={res['dual_residual_at_LB']:.1e} ({dt:.0f}s)",flush=True)
        c=dict(cache[key]); c["primal"]=c["dual_lb"]  # anchor=dual_lb-1e-5
        out.append(c)
    return out

if __name__=="__main__":
    core,_=FS.load_core(); halo=FS.load_halo(); corehalo=core+halo
    import ast
    specs = ast.literal_eval(sys.argv[1])  # list of (h,p,q1,q2)
    N=int(sys.argv[2]) if len(sys.argv)>2 else 5000
    bn=int(sys.argv[3]) if len(sys.argv)>3 else 30
    fresh=get_centers(specs,N=N,bn=bn)
    combined=corehalo+fresh
    hr=(0.0,0.08); pr=(0.33,0.45); qr=(0.025,0.05); TARGET=CORE_HEADLINE
    lb,pt,wit,gmin,eps,Lm=cover_min_over_box(combined,"primal_m1e5",hr,pr,qr,n_h=81,n_p=121,n_q=41)
    print(f"\nSTRIP p{pr} q{qr} with {len(fresh)} fresh: floor={lb:.6f} {'CLEARS' if lb>=TARGET else 'BELOW'} "
          f"(target {TARGET}) grid_min={gmin:.6f} eps={eps:.3e} Lmax={Lm:.3f} worst@h={pt[0]:.4f} p={pt[1]:.4f} q={pt[2]:.4f} wit={wit}")
