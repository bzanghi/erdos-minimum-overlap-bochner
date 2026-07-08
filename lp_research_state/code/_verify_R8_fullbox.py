"""Final: strict rigorous adaptive subdivision over the ENTIRE R8 box
h[0,0.08] p[0,1] q[0.05,1.0] using ONLY cover_min_over_box (saved duals,
grid+Lipschitz). High-q is included whole (no infeasibility exclusion needed:
cover grid_min rises with q, so it clears once eps shrinks; any genuinely
infeasible sub-box only helps). Refine until lb>=TARGET or eps<EPS_FLOOR.
A refined leaf with lb<TARGET => REFUTE."""
import sys, json; sys.path.insert(0,'.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box
TARGET=0.3802838; EPS_FLOOR=8e-4
hv=M.harvest_centers(); kept=[]
for c in hv[0]:
    if c.get('primal') is None: c=dict(c); c['primal']=float(c['dual_lb'])+1e-5
    kept.append(c)
def cover(box,nh,npp,nq):
    (h0,h1),(p0,p1),(q0,q1)=box
    return cover_min_over_box(kept,'primal_m1e5',(h0,h1),(p0,p1),(q0,q1),n_h=nh,n_p=npp,n_q=nq)
def res(box):
    (h0,h1),(p0,p1),(q0,q1)=box
    return (max(7,min(21,int((h1-h0)/0.004)+2)),
            max(9,min(49,int((p1-p0)/0.007)+2)),
            max(7,min(21,int((q1-q0)/0.004)+2)))
def split(box):
    (h0,h1),(p0,p1),(q0,q1)=box
    ax=max([('h',h1-h0),('p',p1-p0),('q',q1-q0)],key=lambda t:t[1])[0]
    if ax=='h': m=.5*(h0+h1); return[((h0,m),(p0,p1),(q0,q1)),((m,h1),(p0,p1),(q0,q1))]
    if ax=='p': m=.5*(p0+p1); return[((h0,h1),(p0,m),(q0,q1)),((h0,h1),(m,p1),(q0,q1))]
    m=.5*(q0+q1); return[((h0,h1),(p0,p1),(q0,m)),((h0,h1),(p0,p1),(m,q1))]
leaves=[]; lowest=[None,np.inf]
import sys
sys.setrecursionlimit(100000)
def rec(box,d):
    nh,npp,nq=res(box)
    lb,pt,wit,gmin,eps,Lm=cover(box,nh,npp,nq)
    if lb>=TARGET:
        leaves.append((box,'cleared',lb,gmin,eps,pt))
        if lb<lowest[1]: lowest[0],lowest[1]=box,lb
        return
    if eps<EPS_FLOOR or d>=30:
        leaves.append((box,'refined',lb,gmin,eps,pt))
        if lb<lowest[1]: lowest[0],lowest[1]=box,lb
        return
    for sb in split(box): rec(sb,d+1)
rec(((0.0,0.08),(0.0,1.0),(0.05,1.0)),0)
floor=min(l[2] for l in leaves)
below=[l for l in leaves if l[1]=='refined' and l[2]<TARGET]
print(f'FULL R8 box: n_leaves={len(leaves)} floor={floor:.7f} eps_floor={EPS_FLOOR}')
print(f'  refined-and-below-target: {len(below)}')
print(f'  binding leaf {lowest[0]} lb={lowest[1]:.7f} margin={lowest[1]-TARGET:+.2e}')
b=lowest[0]
json.dump({'floor':floor,'n_leaves':len(leaves),'n_below':len(below),
           'binding_box':b,'binding_lb':lowest[1],'eps_floor':EPS_FLOOR,'target':TARGET},
          open('_verify_R8_fullbox.json','w'),indent=1)
for l in below[:5]: print('  BELOW:',l[0],l[2])
print('wrote _verify_R8_fullbox.json')
