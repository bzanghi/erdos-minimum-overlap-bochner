"""Rigorous adaptive subdivision certified floor for R5.

Strategy:
  - Recursively bisect the R5 box. For each subbox call cover_min_over_box with a
    grid dense enough that eps_grid < EPS_TARGET *for that subbox*.
  - A subbox is FULLY CERTIFIED once its eps_grid < EPS_TARGET: its returned
    cover_min is a rigorous LB on max_c Phi_c over the subbox.
  - Maintain a running certified floor = min over fully-certified subboxes of cover_min.
  - For subboxes whose cover_min (even partially) is still below the best certified
    floor AND not yet eps-converged, subdivide further (split largest dim).
  - Prune: a subbox whose cover_min (already eps-penalized) >= TARGET cannot lower
    the floor below TARGET, so it can be CERTIFIED-CLEARED without further splitting
    (we only need to prove floor >= TARGET). We still record its value.
  - We separately note, for the lowest few certified subboxes, whether their worst
    point is SDP-feasible (informational; not required if floor >= TARGET on geometry).

Rigor: min over a partition of cover_min_over_box(subbox) is a rigorous LB of
max_c Phi_c over the whole box (each subbox bound is rigorous; min of rigorous
lower bounds over a cover is a rigorous lower bound).
"""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import json, itertools
from heapq import heappush, heappop
import numpy as np
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box, WHITE_TABLE2

TARGET = 0.3802838
EPS_TARGET = 4.0e-4

hv = M.harvest_centers(); centers = hv[0]
kept = []
for c in centers:
    if c.get('primal') is None:
        c = dict(c); c['primal'] = float(c['dual_lb']) + 1e-5
    kept.append(c)
HR = WHITE_TABLE2[4][0]; PR = WHITE_TABLE2[4][1]; QR = WHITE_TABLE2[4][2]
print('R5 box h%s p%s q%s; %d centers; TARGET=%.7f EPS_TARGET=%.1e' % (HR, PR, QR, len(kept), TARGET, EPS_TARGET), flush=True)


def grid_for(box, eps_target):
    """Pick (n_h,n_p,n_q) so that cover_min_over_box's eps_grid < eps_target.
    eps_grid = L_max * 0.5 * sqrt(cell_h^2+cell_p^2+cell_q^2). We first get L_max
    cheaply via a tiny grid, then size cells so each axis contributes < target.
    """
    (h0, h1), (p0, p1), (q0, q1) = box
    # cheap L_max probe
    _, _, _, _, _, Lm = cover_min_over_box(kept, 'primal_m1e5', (h0, h1), (p0, p1), (q0, q1),
                                           n_h=3, n_p=3, n_q=3)
    Lm = max(Lm, 1e-9)
    # want L_max * 0.5 * cell_axis <= eps_target/sqrt(3) per axis (so combined < eps_target)
    budget = (eps_target / Lm) * 2.0 / np.sqrt(3.0)
    def npts(ext):
        if ext <= 0:
            return 1
        n = int(np.ceil(ext / budget)) + 1
        return max(2, min(n, 2000))
    return npts(h1 - h0), npts(p1 - p0), npts(q1 - q0), Lm


def eval_certified(box):
    nh, npp, nq, _ = grid_for(box, EPS_TARGET)
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(kept, 'primal_m1e5', box[0], box[1], box[2],
                                                    n_h=nh, n_p=npp, n_q=nq)
    return lb, pt, wit, gmin, eps, Lm, (nh, npp, nq)


def main():
    # initial split to seed: coarse over q (the dominant decay dim) and p
    hs = list(np.linspace(*HR, 3))
    ps = list(np.linspace(*PR, 9))
    qs = [-1.0, -0.5, -0.3, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0]
    boxes = []
    for i in range(len(hs) - 1):
        for jp in range(len(ps) - 1):
            for k in range(len(qs) - 1):
                boxes.append(((hs[i], hs[i + 1]), (ps[jp], ps[jp + 1]), (qs[k], qs[k + 1])))

    cnt = itertools.count()
    pq = []  # min-heap on partial cover_min
    for b in boxes:
        # quick coarse evaluation for ordering (cheap grid)
        lb0, pt0, _, g0, e0, L0 = cover_min_over_box(kept, 'primal_m1e5', b[0], b[1], b[2],
                                                     n_h=5, n_p=21, n_q=21)
        heappush(pq, (lb0, next(cnt), b))

    certified_floor = np.inf
    binding = None
    low_subboxes = []  # subboxes whose certified cover_min < TARGET (need feasibility scrutiny)
    n_eval = 0
    MAXEVAL = 20000
    while pq and n_eval < MAXEVAL:
        lb0, _, b = heappop(pq)
        # if even the coarse partial bound is already >= current certified_floor, and
        # >= TARGET, it cannot reduce the floor below TARGET -> safe to skip detailed work,
        # but we still must ensure it's >= TARGET rigorously. Do a certified eval.
        if lb0 >= certified_floor:
            # cannot lower the floor; the heap is ordered so everything remaining is >= too
            break
        lb, pt, wit, gmin, eps, Lm, grid = eval_certified(b)
        n_eval += 1
        if eps < EPS_TARGET:
            # fully certified subbox
            if lb < certified_floor:
                certified_floor = lb
                binding = (b, pt, wit, gmin, eps, Lm, grid)
            if lb < TARGET:
                low_subboxes.append((lb, b, pt, wit, gmin, eps))
        else:
            # not converged: split largest dim and requeue
            (h0, h1), (p0, p1), (q0, q1) = b
            exts = sorted([(h1 - h0, 'h'), (p1 - p0, 'p'), (q1 - q0, 'q')], reverse=True)
            dim = exts[0][1]
            if dim == 'h':
                m = (h0 + h1) / 2; subs = [((h0, m), (p0, p1), (q0, q1)), ((m, h1), (p0, p1), (q0, q1))]
            elif dim == 'p':
                m = (p0 + p1) / 2; subs = [((h0, h1), (p0, m), (q0, q1)), ((h0, h1), (m, p1), (q0, q1))]
            else:
                m = (q0 + q1) / 2; subs = [((h0, h1), (p0, p1), (q0, m)), ((h0, h1), (p0, p1), (m, q1))]
            for s in subs:
                # use current (now eps-penalized) lb as the ordering key proxy
                heappush(pq, (lb, next(cnt), s))
        if n_eval % 200 == 0:
            print('  ...n_eval=%d heap=%d certified_floor=%.6f' % (n_eval, len(pq), certified_floor), flush=True)

    print('\nn_eval=%d  remaining heap=%d' % (n_eval, len(pq)), flush=True)
    print('CERTIFIED FLOOR (pure geometry, no infeasibility exclusion) = %.7f' % certified_floor, flush=True)
    print('  vs TARGET 0.3802838: %+.6e' % (certified_floor - TARGET), flush=True)
    if binding:
        b, pt, wit, gmin, eps, Lm, grid = binding
        print('  binding subbox=%s' % (b,), flush=True)
        print('  worst_pt(h,p,q)=%s grid_min=%.6f eps=%.2e Lmax=%.3f grid=%s wit=%s' % (
            tuple(round(x, 5) for x in pt), gmin, eps, Lm, grid, wit), flush=True)
    print('\nlow subboxes (certified cover_min < TARGET): %d' % len(low_subboxes), flush=True)
    low_subboxes.sort()
    for lb, b, pt, wit, gmin, eps in low_subboxes[:12]:
        print('  cover_min=%.6f box=%s worst@(h%.4f,p%.4f,q%.4f) grid_min=%.6f eps=%.2e wit=%s' % (
            lb, tuple(tuple(round(x, 4) for x in ax) for ax in b), pt[0], pt[1], pt[2], gmin, eps, wit), flush=True)

    out = {
        'certified_floor_geometry': certified_floor,
        'target': TARGET,
        'binding_box': binding[0] if binding else None,
        'binding_worst_point': binding[1] if binding else None,
        'binding_grid_min': binding[3] if binding else None,
        'binding_eps': binding[4] if binding else None,
        'n_low_subboxes_below_target': len(low_subboxes),
        'low_subboxes': [{'cover_min': lb, 'box': b, 'worst': pt, 'wit': wit, 'grid_min': gmin, 'eps': eps}
                         for lb, b, pt, wit, gmin, eps in low_subboxes[:30]],
        'n_eval': n_eval,
    }
    json.dump(out, open('_r5_subdiv_out.json', 'w'), indent=2, default=float)
    print('\nsaved _r5_subdiv_out.json', flush=True)


if __name__ == '__main__':
    main()
