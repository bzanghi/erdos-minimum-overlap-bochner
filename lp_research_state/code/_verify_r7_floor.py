"""Rigorous certified floor for R7 over the FULL box, with active-center pruning
(legitimate: a center provably below env minus the kept Lipschitz slack can be
dropped from L_max) and adaptive subdivision until per-box eps<5e-4."""
import sys; sys.path.insert(0,'.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import phi_center_grid, anchor_value, WHITE_TABLE2

hv=M.harvest_centers(); centers=hv[0]
kept=[]
for c in centers:
    if c.get('primal') is None and c.get('dual_lb') is not None:
        c=dict(c); c['primal']=float(c['dual_lb'])+1e-5
    kept.append(c)
for c in kept: c['_a']=anchor_value(c,'primal_m1e5')

def gmax(lc,qc,lo,hi): return max(abs(lc+qc*lo),abs(lc+qc*hi))
def subbox_L(c,h0,h1,p0,p1,q0,q1):
    d=c['duals']
    return float(np.sqrt(gmax(d['con_53'],-d['con_54'],h0,h1)**2
        +gmax(d['con_512_pL']-d['con_512_pU'],-d['con_513'],p0,p1)**2
        +gmax(d['con_512_qL']-d['con_512_qU'],-d['con_513'],q0,q1)**2))

def cover_lb_pruned(h0,h1,p0,p1,q0,q1,nh,npp,nq):
    hg=np.linspace(h0,h1,nh); pg=np.linspace(p0,p1,npp)
    qg=np.unique(np.concatenate([np.linspace(q0,q1,nq),[0.0] if q0<=0<=q1 else []]))
    HH,PP=np.meshgrid(hg,pg,indexing='ij')
    half_diag=0.5*float(np.sqrt(((h1-h0)/(nh-1) if nh>1 else 0)**2
        +((p1-p0)/(npp-1) if npp>1 else 0)**2+((q1-q0)/(len(qg)-1) if len(qg)>1 else 0)**2))
    env=np.full_like(HH,-np.inf); gridmax={}; wit=np.empty(HH.shape,dtype=object)
    for q in qg:
        for c in kept:
            F=phi_center_grid(c,c['_a'],HH,PP,q)
            gridmax[c['label']]=max(gridmax.get(c['label'],-np.inf),float(F.max()))
            m=F>env; env[m]=F[m]; wit[m]=c['label']
    gme=float(env.min()); a=np.unravel_index(int(env.argmin()),env.shape)
    wp=(float(HH[a]),float(PP[a]))
    Ls={c['label']:subbox_L(c,h0,h1,p0,p1,q0,q1) for c in kept}
    Lk=max(Ls.values())
    for _ in range(40):
        thr=gme-Lk*half_diag
        kc=[lab for lab,gm in gridmax.items() if gm+Ls[lab]*half_diag>=thr]
        nL=max(Ls[lab] for lab in kc) if kc else 0.0
        if abs(nL-Lk)<1e-13: break
        Lk=nL
    return gme-Lk*half_diag, gme, Lk*half_diag, Lk, len(kc), wp, str(wit[a])

def refine(box,depth,budget):
    h0,h1,p0,p1,q0,q1=box
    lb,gme,eps,Lk,nk,wp,w=cover_lb_pruned(h0,h1,p0,p1,q0,q1,9,9,7)
    if eps<5e-4 or depth>=budget:
        return lb,(box,gme,eps,Lk,nk,wp,w)
    exts=sorted([(h1-h0,'h'),(p1-p0,'p'),(q1-q0,'q')],reverse=True)
    dim=exts[0][1]
    if dim=='h': bs=[(h0,(h0+h1)/2,p0,p1,q0,q1),((h0+h1)/2,h1,p0,p1,q0,q1)]
    elif dim=='p': bs=[(h0,h1,p0,(p0+p1)/2,q0,q1),(h0,h1,(p0+p1)/2,p1,q0,q1)]
    else: bs=[(h0,h1,p0,p1,q0,(q0+q1)/2),(h0,h1,p0,p1,(q0+q1)/2,q1)]
    best=np.inf; bi=None
    for b in bs:
        l,info=refine(b,depth+1,budget)
        if l<best: best=l; bi=info
    return best,bi

hr,pr,qr,wb=WHITE_TABLE2[6]
h0,h1=hr; q0,q1=qr
# Cover ENTIRE p in [0,1]. Coarse p-scan then recurse each.
import time, json
npbox=50; pe=np.linspace(0,1,npbox+1)
glob=np.inf; ginfo=None; t0=time.time()
percentile=[]
for i in range(npbox):
    b=(h0,h1,float(pe[i]),float(pe[i+1]),q0,q1)
    l,info=refine(b,0,12)
    percentile.append((float(pe[i]),float(pe[i+1]),l))
    if l<glob: glob=l; ginfo=info
print(f"[FULL R7 certified floor, p in [0,1], pruned+adaptive, eps<5e-4] {time.time()-t0:.1f}s")
print(f"  CERTIFIED FLOOR = {glob:.7f}")
box,gme,eps,Lk,nk,wp,w=ginfo
print(f"  binding box h=[{box[0]:.4f},{box[1]:.4f}] p=[{box[2]:.4f},{box[3]:.4f}] q=[{box[4]:.5f},{box[5]:.5f}]")
print(f"  grid_min={gme:.7f} eps={eps:.3e} L_keep={Lk:.4f} n_active={nk} worst_pt=(h={wp[0]:.4f},p={wp[1]:.4f}) wit={w}")
print(f"  vs 0.3802838: {glob-0.3802838:+.3e}  => {'CONFIRM' if glob>=0.3802838 else 'REFUTE'}")
# worst few p-bands
percentile.sort(key=lambda r:r[2])
print("  worst 6 p-bands (floor):")
for a,b,l in percentile[:6]:
    print(f"    p=[{a:.3f},{b:.3f}] floor={l:.6f}")
json.dump({"certified_floor":glob,"binding_box":box,"grid_min":gme,"eps":eps,"L_keep":Lk,
           "n_active":nk,"worst_pt":wp,"witness":w},
          open("../parallel_results/_r7_floor_pruned.json","w"),indent=2)
