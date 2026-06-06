"""R9 strip: subdivide h x p x q into a grid of sub-boxes; min of sub-box floors is
a rigorous floor for the strip. Smaller sub-boxes => smaller eps_grid (Lipschitz)."""
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
    import ast
    keys=ast.literal_eval(sys.argv[1]); nh,np_,nq=int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    core,_=FS.load_core(); halo=FS.load_halo()
    combined=core+halo+load_fresh(keys)
    TARGET=CORE_HEADLINE
    he=np.linspace(0,0.08,nh+1); pe=np.linspace(0.33,0.45,np_+1); qe=np.linspace(0.025,0.05,nq+1)
    overall=np.inf; worst=None; nbox=0
    for i in range(nh):
        for j in range(np_):
            for k in range(nq):
                hr=(float(he[i]),float(he[i+1])); pr=(float(pe[j]),float(pe[j+1])); qr=(float(qe[k]),float(qe[k+1]))
                lb,pt,wit,gmin,eps,Lm=cover_min_over_box(combined,"primal_m1e5",hr,pr,qr,n_h=21,n_p=21,n_q=11)
                nbox+=1
                if lb<overall: overall=lb; worst=(hr,pr,qr,pt,wit,gmin,eps,Lm)
    print(f"subdiv {nh}x{np_}x{nq} = {nbox} sub-boxes; strip floor = {overall:.6f}  "
          f"{'CLEARS '+str(TARGET) if overall>=TARGET else 'BELOW'}")
    hr,pr,qr,pt,wit,gmin,eps,Lm=worst
    print(f"  binding sub-box h{hr} p{pr} q{qr}: grid_min={gmin:.6f} eps={eps:.3e} Lmax={Lm:.3f} @h={pt[0]:.4f}p={pt[1]:.4f}q={pt[2]:.4f} wit={wit}")
