"""Branch-and-bound rigorous floor for R6 feasible region.
Recursively subdivide; a sub-box is CLOSED once cover_min_over_box(subbox) >= TARGET
(it can never bring the regional min below target). Only sub-boxes still < TARGET
are split further. Reports the regional floor = min over all closed/leaf sub-boxes,
and whether every sub-box cleared TARGET. This is RIGOROUS: the regional cover-min
lower bound is the min over the partition of the per-subbox rigorous lower bounds."""
import sys, os, json; sys.path.insert(0,'.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box

TARGET = 0.3802838
DROP = set(os.environ.get('R6_DROP','').split(',')) if os.environ.get('R6_DROP') else set()

hv = M.harvest_centers()
kept=[]
for c in hv[0]:
    if c.get('primal') is None:
        if c.get('dual_lb') is not None:
            c=dict(c); c['primal']=float(c['dual_lb'])+1e-5
        else: continue
    if c['label'] in DROP: continue
    kept.append(c)
print('usable centers:', len(kept), 'dropped:', sorted(DROP), flush=True)

def cm(hr,pr,qr,nh,npp,nq):
    return cover_min_over_box(kept,'primal_m1e5',hr,pr,qr,n_h=nh,n_p=npp,n_q=nq)

NB=9   # base grid per axis per box (cheap; eps via local Lipschitz)
MAXD=int(sys.argv[2]) if len(sys.argv)>2 else 22
EPS_FLOOR=float(sys.argv[3]) if len(sys.argv)>3 else 3e-4

import heapq
box=json.loads(sys.argv[1])
hr,pr,qr=tuple(box[0]),tuple(box[1]),tuple(box[2])

# work stack of (hr,pr,qr,depth); track regional min over leaves proven done
regional_min=np.inf; worst_leaf=None; n_leaves=0; n_open=0
stack=[(hr,pr,qr,0)]
all_cleared=True
import time; t0=time.time()
while stack:
    bh,bp,bq,dep=stack.pop()
    nh=NB if (bh[1]-bh[0])>1e-12 else 1
    npp=NB if (bp[1]-bp[0])>1e-12 else 1
    nq=NB if (bq[1]-bq[0])>1e-12 else 1
    lb,pt,wit,gmin,eps,Lm=cm(bh,bp,bq,nh,npp,nq)
    if lb>=TARGET:
        # this sub-box is fully cleared; contributes lb to regional min
        n_leaves+=1
        if lb<regional_min: regional_min=lb; worst_leaf=(bh,bp,bq,pt,eps)
        continue
    # not cleared: must subdivide unless we hit limits
    if dep>=MAXD or eps<EPS_FLOOR:
        # cannot refine further: this is a genuine leaf BELOW target
        n_leaves+=1; all_cleared=False
        if lb<regional_min: regional_min=lb; worst_leaf=(bh,bp,bq,pt,eps)
        continue
    # split longest axis
    spans=[bh[1]-bh[0],bp[1]-bp[0],bq[1]-bq[0]]
    ax=int(np.argmax(spans))
    if ax==0:
        m=0.5*(bh[0]+bh[1]); stack+=[((bh[0],m),bp,bq,dep+1),((m,bh[1]),bp,bq,dep+1)]
    elif ax==1:
        m=0.5*(bp[0]+bp[1]); stack+=[(bh,(bp[0],m),bq,dep+1),(bh,(m,bp[1]),bq,dep+1)]
    else:
        m=0.5*(bq[0]+bq[1]); stack+=[(bh,bp,(bq[0],m),dep+1),(bh,bp,(m,bq[1]),dep+1)]
    n_open+=1
    if n_open%2000==0:
        print('  ...explored %d splits, leaves=%d, regional_min so far=%.7f (%.0fs)'%(
            n_open,n_leaves,regional_min,time.time()-t0),flush=True)

print('DONE splits=%d leaves=%d'%(n_open,n_leaves),flush=True)
print('regional_floor (min over partition lower bounds) = %.7f'%regional_min)
print('all sub-boxes cleared TARGET=%.7f : %s'%(TARGET,all_cleared))
wb=worst_leaf
print('worst leaf: h=%s p=%s q=%s worst_pt=(%.4f,%.4f,%.4f) eps=%.3e'%(
    wb[0],wb[1],wb[2],wb[3][0],wb[3][1],wb[3][2],wb[4]))
print('vs TARGET: %+.4e'%(regional_min-TARGET))
