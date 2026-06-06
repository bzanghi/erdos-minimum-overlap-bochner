"""Rigorous box-min LB of the augmented Phi-cover over R7, via ADAPTIVE
subdivision so the Lipschitz eps_grid is driven down where the cover is tight.

Loads: existing 12 core anchors + prior 11 halo + fresh R7 centers (this task).
All composed under the conservative primal_m1e5 anchor convention (anchor = the
stored 'primal' field minus 1e-5; for halo/fresh we set primal := conservative
anchor so the convention matches). Cover = max_c Phi_c is a valid mu-LB at every
point; box-min via fine grid + L_max*half_diag is rigorous.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _fullspace_eval import load_centers, cover_min_over_box, CORE_HEADLINE, WHITE_OUTSIDE_FLOOR

def build_combined():
    existing,_ = load_centers()
    hd = json.load(open(CODE.parent/"parallel_results"/"fullspace_stage2_halo_centers.json"))
    halo = hd["centers"]
    for c in halo: c["primal"] = c["dual_lb"]
    fp = CODE.parent/"parallel_results"/"fullspace_promote_R7.json"
    fresh = []
    if fp.exists():
        fresh = json.load(open(fp))["centers"]
        for c in fresh: c["primal"] = c["anchor"]
    return existing, halo, fresh, existing+halo+fresh

def adaptive_boxmin(centers, h_range, p_range, q_range, target,
                    depth=0, max_depth=6, base=(21,21,21)):
    """Return rigorous box-min LB over the (h,p,q) box. Subdivide the worst axis
    until either lb>=target (cleared) or eps_grid is negligible vs (grid_min-target)
    or max_depth reached. Returns (lb, grid_min_min, worst_pt, worst_wit, n_leaves)."""
    nh,npp,nq = base
    lb,pt,wit,gm,eps,L = cover_min_over_box(centers,"primal_m1e5",h_range,p_range,q_range,
                                            n_h=nh,n_p=npp,n_q=nq)
    if lb >= target or depth >= max_depth:
        return lb, gm, pt, wit, 1
    # subdivide the longest axis (by span) to shrink eps
    spans = [(h_range[1]-h_range[0], 'h'),(p_range[1]-p_range[0],'p'),(q_range[1]-q_range[0],'q')]
    spans.sort(reverse=True)
    axis = spans[0][1]
    def split(rng):
        m=0.5*(rng[0]+rng[1]); return (rng[0],m),(m,rng[1])
    worst_lb=np.inf; worst=None; tot=0
    if axis=='h': parts=[(a,p_range,q_range) for a in split(h_range)]
    elif axis=='p': parts=[(h_range,a,q_range) for a in split(p_range)]
    else: parts=[(h_range,p_range,a) for a in split(q_range)]
    for hr,pr,qr in parts:
        l,g,pp,ww,n = adaptive_boxmin(centers,hr,pr,qr,target,depth+1,max_depth,base)
        tot+=n
        if l<worst_lb: worst_lb=l; worst=(g,pp,ww)
    return worst_lb, worst[0], worst[1], worst[2], tot

if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    existing,halo,fresh,combined = build_combined()
    print(f"centers: existing={len(existing)} halo={len(halo)} fresh={len(fresh)} combined={len(combined)}")
    for c in fresh:
        print(f"  fresh {c['label']}: anchor={c['anchor']:.6f} resid={c['dual_resid']:.1e} status={c['status']}")
    H=(0.0,0.08); P=(0.0,1.0); Q=(-0.05,-0.025)
    # Coarse pass over full R7 to find weak p-bands, then adaptive on each.
    print("\n=== Adaptive box-min over full R7 (h[0,0.08] p[0,1] q[-0.05,-0.025]) ===")
    pedges=np.array([0.0,0.25,0.33,0.35,0.37,0.38,0.39,0.40,0.42,0.45,0.6,1.0])
    worst=np.inf; worst_box=None; below_floor=[]
    for i in range(len(pedges)-1):
        p0,p1=float(pedges[i]),float(pedges[i+1])
        lb,gm,pt,wit,nl = adaptive_boxmin(combined,H,(p0,p1),Q,CORE_HEADLINE,max_depth=7)
        tag = "OK>=CORE" if lb>=CORE_HEADLINE else (">=FLOOR" if lb>=WHITE_OUTSIDE_FLOOR else "<FLOOR!")
        if lb<worst: worst=lb; worst_box=(p0,p1,gm,pt,wit)
        if lb<WHITE_OUTSIDE_FLOOR: below_floor.append((p0,p1,lb))
        print(f"  p=[{p0:.3f},{p1:.3f}] lb={lb:.6f} ceil(gm)={gm:.6f} leaves={nl} worst@(h={pt[0]:.4f},p={pt[1]:.4f},q={pt[2]:.4f}) {tag}")
    print(f"\nR7 independent floor (worst sub-box lb) = {worst:.6f}")
    print(f"   worst box p=[{worst_box[0]},{worst_box[1]}] ceiling grid_min={worst_box[2]:.6f}")
    print(f"   clears 0.380000 independently: {worst>=WHITE_OUTSIDE_FLOOR}")
    print(f"   clears 0.380284 independently: {worst>=CORE_HEADLINE}")
    if below_floor: print(f"   BELOW FLOOR bands: {below_floor}")
