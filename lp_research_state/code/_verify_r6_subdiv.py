"""Independent adaptive subdivision of R6 over a (feasible) sub-box.
Recursively splits the box into sub-boxes, taking the MIN over sub-boxes of
cover_min_over_box(subbox). Each sub-box has smaller half_diag => smaller
eps_grid => tighter (rigorous) bound. Refine until per-subbox eps_grid < tol.

This is the SAME rigorous machinery (grid + Lipschitz) the candidate relies on;
we re-run it INDEPENDENTLY on the feasible region only to find the true floor."""
import sys; sys.path.insert(0, '.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box


DROP = set()  # labels to drop from cover (set via CLI)


def load_kept():
    hv = M.harvest_centers()
    kept = []
    for c in hv[0]:
        if c.get('primal') is None:
            if c.get('dual_lb') is not None:
                c = dict(c); c['primal'] = float(c['dual_lb']) + 1e-5
            else:
                continue
        if c['label'] in DROP:
            continue
        kept.append(c)
    return kept


def cover_box(kept, hr, pr, qr, nh, npp, nq):
    return cover_min_over_box(kept, 'primal_m1e5', hr, pr, qr, n_h=nh, n_p=npp, n_q=nq)


def adaptive_min(kept, hr, pr, qr, eps_tol=5e-4, max_depth=14, base=21,
                 _depth=0, _stats=None):
    """Return (floor_lb, worst_pt, worst_eps, max_depth_reached). Recursively
    subdivide along the longest axis until eps_grid < eps_tol or max_depth."""
    if _stats is None:
        _stats = {'leaves': 0, 'max_depth': 0, 'min_eps': np.inf}
    nh = base if (hr[1] - hr[0]) > 1e-9 else 1
    npp = base if (pr[1] - pr[0]) > 1e-9 else 1
    nq = base if (qr[1] - qr[0]) > 1e-9 else 1
    lb, pt, wit, gmin, eps, Lm = cover_box(kept, hr, pr, qr, nh, npp, nq)
    if eps < eps_tol or _depth >= max_depth:
        _stats['leaves'] += 1
        _stats['max_depth'] = max(_stats['max_depth'], _depth)
        _stats['min_eps'] = min(_stats['min_eps'], eps)
        return lb, pt, eps, _stats
    # split along the axis with the largest span*Lipschitz contribution; simplest:
    # split the physically longest axis (after scaling p,q which carry the L_max).
    spans = [(hr[1] - hr[0]), (pr[1] - pr[0]), (qr[1] - qr[0])]
    ax = int(np.argmax(spans))
    best = (np.inf, None, None)
    if ax == 0:
        mid = 0.5 * (hr[0] + hr[1])
        boxes = [((hr[0], mid), pr, qr), ((mid, hr[1]), pr, qr)]
    elif ax == 1:
        mid = 0.5 * (pr[0] + pr[1])
        boxes = [(hr, (pr[0], mid), qr), (hr, (mid, pr[1]), qr)]
    else:
        mid = 0.5 * (qr[0] + qr[1])
        boxes = [(hr, pr, (qr[0], mid)), (hr, pr, (mid, qr[1]))]
    for (bh, bp, bq) in boxes:
        sub_lb, sub_pt, sub_eps, _ = adaptive_min(
            kept, bh, bp, bq, eps_tol, max_depth, base, _depth + 1, _stats)
        if sub_lb < best[0]:
            best = (sub_lb, sub_pt, sub_eps)
    return best[0], best[1], best[2], _stats


if __name__ == '__main__':
    import json
    import os
    drop_env = os.environ.get('R6_DROP', '')
    if drop_env:
        DROP.update(drop_env.split(','))
    kept = load_kept()
    print('usable centers:', len(kept), flush=True)
    # box passed as JSON: [[h0,h1],[p0,p1],[q0,q1]]
    box = json.loads(sys.argv[1])
    hr, pr, qr = tuple(box[0]), tuple(box[1]), tuple(box[2])
    eps_tol = float(sys.argv[2]) if len(sys.argv) > 2 else 5e-4
    base = int(sys.argv[3]) if len(sys.argv) > 3 else 21
    md = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    floor, pt, eps, stats = adaptive_min(kept, hr, pr, qr, eps_tol=eps_tol,
                                         max_depth=md, base=base)
    print('BOX h=%s p=%s q=%s' % (hr, pr, qr))
    print('  certified floor (min over subboxes) = %.7f' % floor)
    print('  worst point = (h=%.4f, p=%.4f, q=%.4f)' % pt)
    print('  worst-subbox eps_grid = %.4e' % eps)
    print('  leaves=%d max_depth=%d min_eps=%.2e' % (
        stats['leaves'], stats['max_depth'], stats['min_eps']))
    print('  vs 0.3802838: %+.4e' % (floor - 0.3802838))
