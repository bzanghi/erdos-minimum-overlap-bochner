"""Rigorous branch-and-bound certification of R9 full box to >= TARGET.
For each box: run cover_min_over_box (full 121 centers). If certified floor (grid-eps)
>= TARGET, the box is CLEARED (rigorous). Else split along the dimension whose
eps-contribution (L_axis * cell) is largest and recurse. The min certified floor over
all CLEARED leaves is a rigorous lower bound iff every leaf clears (or we track the
true min floor over all leaves). Returns the global min certified floor and binding box.
"""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import _r9_verify as V
from _fullspace_eval import cover_min_over_box, anchor_value

TARGET = V.TARGET

def axis_L(kept, hr, pr, qr):
    """Per-axis Lipschitz upper bounds (gh,gp,gq) max over centers, for split choice."""
    h0,h1=hr; p0,p1=pr; q0,q1=qr
    def gmax(lin,quad,lo,hi): return max(abs(lin+quad*lo),abs(lin+quad*hi))
    GH=GP=GQ=0.0
    for c in kept:
        d=c['duals']
        GH=max(GH,gmax(d['con_53'],-d['con_54'],h0,h1))
        GP=max(GP,gmax(d['con_512_pL']-d['con_512_pU'],-d['con_513'],p0,p1))
        GQ=max(GQ,gmax(d['con_512_qL']-d['con_512_qU'],-d['con_513'],q0,q1))
    return GH,GP,GQ

def certify(kept, hr, pr, qr, max_leaves=200000):
    stack=[(hr,pr,qr,0)]
    min_floor=np.inf; binding=None; n_leaves=0; max_depth=0; max_eps_cleared=0.0
    while stack:
        hr_,pr_,qr_,depth=stack.pop()
        max_depth=max(max_depth,depth)
        lb,pt,wit,gmin,eps,Lm=cover_min_over_box(kept,'primal_m1e5',hr_,pr_,qr_,n_h=15,n_p=31,n_q=11)
        if lb>=TARGET:
            n_leaves+=1
            if lb<min_floor: min_floor=lb; binding=(hr_,pr_,qr_,pt,wit,gmin,eps)
            max_eps_cleared=max(max_eps_cleared,eps)
            continue
        # not cleared. If grid_min itself < TARGET -> genuine deficit (REFUTE locus).
        if gmin < TARGET:
            return {'verdict':'REFUTE','floor':lb,'grid_min':gmin,'box':(hr_,pr_,qr_),
                    'worst_point':pt,'witness':wit,'eps':eps,'depth':depth}
        # grid_min >= TARGET but eps too big -> split widest eps-axis
        GH,GP,GQ=axis_L(kept,hr_,pr_,qr_)
        ch=(hr_[1]-hr_[0]); cp=(pr_[1]-pr_[0]); cq=(qr_[1]-qr_[0])
        contrib=[('h',GH*ch,hr_),('p',GP*cp,pr_),('q',GQ*cq,qr_)]
        ax=max(contrib,key=lambda t:t[1])[0]
        if n_leaves+len(stack) > max_leaves:
            return {'verdict':'INCONCLUSIVE','reason':'max_leaves','floor':min_floor}
        if ax=='h':
            m=0.5*(hr_[0]+hr_[1]); stack.append(((hr_[0],m),pr_,qr_,depth+1)); stack.append(((m,hr_[1]),pr_,qr_,depth+1))
        elif ax=='p':
            m=0.5*(pr_[0]+pr_[1]); stack.append((hr_,(pr_[0],m),qr_,depth+1)); stack.append((hr_,(m,pr_[1]),qr_,depth+1))
        else:
            m=0.5*(qr_[0]+qr_[1]); stack.append((hr_,pr_,(qr_[0],m),depth+1)); stack.append((hr_,pr_,(m,qr_[1]),depth+1))
    return {'verdict':'CONFIRM','floor':float(min_floor),'binding':binding,
            'n_leaves':n_leaves,'max_depth':max_depth,'max_eps_cleared':max_eps_cleared}
