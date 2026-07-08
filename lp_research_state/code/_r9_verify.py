"""Independent verification of salvaged R9 floor (claimed 0.380367, target >=0.3802838).
Reuses harvest_centers + cover_min_over_box (rigorous grid+Lipschitz). Adds:
 (A) full-box cover, (B) feasibility map via SDP, (C) adaptive subdivision over feasible region.
"""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box

TARGET = 0.3802838

def load_kept():
    hv = M.harvest_centers()
    centers = hv[0]
    kept = []
    for c in centers:
        if c.get('primal') is None:
            if c.get('dual_lb') is None:
                continue
            c = dict(c); c['primal'] = float(c['dual_lb']) + 1e-5
        kept.append(c)
    return kept

def feasibility(h, p, q, N=4000, T=1600, R=10, bochner_n=20):
    """Solve build_problem_with_dual_handles at single point (q1=q2=q). Return status."""
    import cvxpy as cp
    from path_b_analytical import build_problem_with_dual_handles
    Omega, cons, H = build_problem_with_dual_handles(N, T, R, h, h, p, p, q, q, bochner_n=bochner_n)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception as e:
        return 'solver_error:' + str(e)[:60]
    return prob.status

def adaptive_floor(kept, hr, pr, qr, eps_target=5e-4, max_depth=6, base=(41,81,41), depth=0):
    """Recursive: cover_min_over_box on box; if eps_grid > eps_target and depth<max,
    split the box along its largest dimension and recurse, return min over children.
    Returns (floor, worst_point, eps_achieved, n_boxes)."""
    nh, npp, nq = base
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(kept, 'primal_m1e5', hr, pr, qr,
                                                    n_h=nh, n_p=npp, n_q=nq)
    if eps <= eps_target or depth >= max_depth:
        return lb, pt, eps, 1, gmin, Lm
    # split along widest of (h,p,q) measured in box units scaled by L contribution-ish:
    wh = hr[1]-hr[0]; wp = pr[1]-pr[0]; wq = qr[1]-qr[0]
    dims = [('h', wh), ('p', wp), ('q', wq)]
    dim = max(dims, key=lambda t: t[1])[0]
    best = (np.inf, None, eps, 0, gmin, Lm)
    if dim == 'h':
        mid = 0.5*(hr[0]+hr[1]); subs = [((hr[0],mid),pr,qr), ((mid,hr[1]),pr,qr)]
    elif dim == 'p':
        mid = 0.5*(pr[0]+pr[1]); subs = [(hr,(pr[0],mid),qr), (hr,(mid,pr[1]),qr)]
    else:
        mid = 0.5*(qr[0]+qr[1]); subs = [(hr,pr,(qr[0],mid)), (hr,pr,(mid,qr[1]))]
    tot = 0
    for (shr, spr, sqr) in subs:
        f, p2, e2, n2, gm2, lm2 = adaptive_floor(kept, shr, spr, sqr, eps_target, max_depth, base, depth+1)
        tot += n2
        if f < best[0]:
            best = (f, p2, e2, n2, gm2, lm2)
    return best[0], best[1], best[2], tot, best[4], best[5]
