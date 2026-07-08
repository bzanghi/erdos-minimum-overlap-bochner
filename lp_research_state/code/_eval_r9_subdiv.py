"""Evaluate R9 strip with q-subdivision to cut eps_grid. Uses cached fresh centers.
A box-min over a q-partition = min of sub-box floors (rigorous)."""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
CODE=Path(__file__).resolve().parent; sys.path.insert(0,str(CODE))
import _fs_recon as FS
from _fullspace_eval import cover_min_over_box, CORE_HEADLINE
CACHE=CODE.parent/"parallel_results"/"_r9_probe_cache.json"

def load_fresh(keys):
    cache=json.load(open(CACHE)); out=[]
    for k in keys:
        c=dict(cache[k]); c["primal"]=c["dual_lb"]; out.append(c)
    return out

if __name__=="__main__":
    core,_=FS.load_core(); halo=FS.load_halo(); corehalo=core+halo
    import ast
    keys=ast.literal_eval(sys.argv[1])
    fresh=load_fresh(keys)
    combined=corehalo+fresh
    TARGET=CORE_HEADLINE; hr=(0.0,0.08); pr=(0.33,0.45)
    # subdivide q into n_sub sub-strips
    n_sub=int(sys.argv[2]) if len(sys.argv)>2 else 2
    n_psub=int(sys.argv[3]) if len(sys.argv)>3 else 1  # also subdivide p?
    qedges=np.linspace(0.025,0.05,n_sub+1)
    pedges=np.linspace(0.33,0.45,n_psub+1)
    overall=np.inf; worst=None
    for qi in range(n_sub):
        for pi in range(n_psub):
            qsub=(float(qedges[qi]),float(qedges[qi+1]))
            prsub=(float(pedges[pi]),float(pedges[pi+1]))
            lb,pt,wit,gmin,eps,Lm=cover_min_over_box(combined,"primal_m1e5",hr,prsub,qsub,n_h=81,n_p=81,n_q=41)
            print(f"  q{qsub} p{prsub}: floor={lb:.6f} {'OK' if lb>=TARGET else 'BELOW'} grid_min={gmin:.6f} eps={eps:.3e} Lmax={Lm:.3f} worst@h={pt[0]:.4f}p={pt[1]:.4f}q={pt[2]:.4f} wit={wit}")
            if lb<overall: overall=lb; worst=(qsub,prsub,pt,wit)
    print(f"\nSTRIP min over {n_sub}x{n_psub} subdiv = {overall:.6f}  {'CLEARS '+str(TARGET) if overall>=TARGET else 'BELOW '+str(TARGET)}")
    print(f"  binding sub-box q{worst[0]} p{worst[1]} @ {worst[2]} wit={worst[3]}")
