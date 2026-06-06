"""R8 adaptive subdivision on the FEASIBLE region using ONLY the rigorous
cover_min_over_box (saved per-center duals; grid + Lipschitz). No SDP, no
poly_moment tail-bound. Recursively split sub-boxes whose lb < target,
take min over leaves of cover_min_over_box(subbox). A sub-box is a 'leaf'
once its lb >= target (cleared) or its eps_grid < EPS_FLOOR (refined to
resolution; report its true grid-min-based lb)."""
import sys, json; sys.path.insert(0,'.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box

TARGET = 0.3802838
EPS_FLOOR = 5e-4   # refine until per-subbox eps_grid < this

hv = M.harvest_centers()
centers = hv[0]
kept = []
for c in centers:
    if c.get('primal') is None:
        c = dict(c); c['primal'] = float(c['dual_lb']) + 1e-5
    kept.append(c)
print(f'loaded {len(kept)} centers', flush=True)

def cover(box, nh, npp, nq):
    (h0,h1),(p0,p1),(q0,q1)=box
    return cover_min_over_box(kept,'primal_m1e5',(h0,h1),(p0,p1),(q0,q1),n_h=nh,n_p=npp,n_q=nq)

def res_for(box):
    """pick grid resolution ~ so cell sizes are comparable & modest cost."""
    (h0,h1),(p0,p1),(q0,q1)=box
    nh = max(9, min(33, int((h1-h0)/0.002)+1))
    npp= max(9, min(81, int((p1-p0)/0.004)+1))
    nq = max(9, min(33, int((q1-q0)/0.002)+1))
    return nh,npp,nq

# Adaptive recursion. We split the longest axis (in eps-contribution sense).
def split(box):
    (h0,h1),(p0,p1),(q0,q1)=box
    # split the axis with the largest cell*Lipschitz contribution; approximate by
    # raw extent weighted: q & h carry L~7.7 (full grad), p smaller. Use extent.
    exts = [('h',h1-h0),('p',(p1-p0)),('q',(q1-q0))]
    ax = max(exts, key=lambda t:t[1])[0]
    if ax=='h':
        m=0.5*(h0+h1); return [((h0,m),(p0,p1),(q0,q1)),((m,h1),(p0,p1),(q0,q1))]
    if ax=='p':
        m=0.5*(p0+p1); return [((h0,h1),(p0,m),(q0,q1)),((h0,h1),(m,p1),(q0,q1))]
    m=0.5*(q0+q1); return [((h0,h1),(p0,p1),(q0,m)),((h0,h1),(p0,p1),(m,q1))]

leaves=[]
worst=[None, np.inf]  # [box, lb] -- the binding (lowest) cleared/refined leaf
def recurse(box, depth):
    nh,npp,nq = res_for(box)
    lb,pt,wit,gmin,eps,Lm = cover(box,nh,npp,nq)
    if lb >= TARGET:
        leaves.append((box,'cleared',lb,gmin,eps,pt)); 
        if lb<worst[1]: worst[0],worst[1]=box,lb
        return
    if eps < EPS_FLOOR or depth>=22:
        # fully refined; this leaf's TRUE lb stands (may be < target => REFUTE evidence)
        leaves.append((box,'refined',lb,gmin,eps,pt))
        if lb<worst[1]: worst[0],worst[1]=box,lb
        return
    for sb in split(box):
        recurse(sb, depth+1)

# FEASIBLE region: per grid_min monotonicity, binding is the q=0.05 edge.
# We subdivide the full FEASIBLE box h[0,0.08] p[0,1] q[0.05, QHI].
# QHI = upper end of feasible region (set generously; high-q clears trivially by
# grid_min anyway, and infeasible high-q only helps). Use QHI=0.30 (grid_min there
# already ~0.42, far above target, so cover clears regardless of feasibility).
FEAS_BOX = ((0.0,0.08),(0.0,1.0),(0.05,0.30))
print(f'recursing on FEASIBLE box {FEAS_BOX} target={TARGET} eps_floor={EPS_FLOOR}', flush=True)
recurse(FEAS_BOX, 0)

floor = min(l[2] for l in leaves)
n_refined_below = sum(1 for l in leaves if l[1]=='refined' and l[2] < TARGET)
n_cleared = sum(1 for l in leaves if l[1]=='cleared')
n_refined = sum(1 for l in leaves if l[1]=='refined')
print(f'\nn_leaves={len(leaves)}  cleared={n_cleared}  refined={n_refined}', flush=True)
print(f'refined leaves still BELOW target: {n_refined_below}', flush=True)
print(f'==> FEASIBLE certified floor (min over leaves of cover lb) = {floor:.7f}', flush=True)
print(f'    binding leaf: {worst[0]} lb={worst[1]:.7f}', flush=True)
# show the lowest few leaves
leaves.sort(key=lambda l:l[2])
print('\nLowest leaves:', flush=True)
for box,kind,lb,gmin,eps,pt in leaves[:8]:
    (h0,h1),(p0,p1),(q0,q1)=box
    print(f'  [{kind}] h[{h0:.4f},{h1:.4f}] p[{p0:.3f},{p1:.3f}] q[{q0:.4f},{q1:.4f}] lb={lb:.6f} grid_min={gmin:.6f} eps={eps:.2e} worst@({pt[0]:.3f},{pt[1]:.3f},{pt[2]:.3f})', flush=True)

json.dump({'floor':floor,'n_leaves':len(leaves),'n_refined_below':n_refined_below,
           'binding_box':worst[0],'binding_lb':worst[1],'eps_floor':EPS_FLOOR,'target':TARGET,
           'lowest_leaves':[{'box':l[0],'kind':l[1],'lb':l[2],'grid_min':l[3],'eps':l[4],'worst':l[5]} for l in leaves[:12]]},
          open('_verify_R8_adaptive.json','w'), indent=1)
print('wrote _verify_R8_adaptive.json', flush=True)
