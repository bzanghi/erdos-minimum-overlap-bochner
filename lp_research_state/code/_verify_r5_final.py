"""FINAL rigorous certified floor for R5. Corrected: refine ALL dims (incl. h) so
each tile's eps_grid < EPS_TARGET. Two phases:
  A) coarse q-tiles outside the low band, cleared if cover_min >= TARGET (dense enough).
  B) dense (p,q,h) tiling of the low band q in [-0.42, 0.20]; certified floor = min.
Rigorous: cover_min_over_box is a rigorous LB of max_c Phi_c on a box; min over a
partition of the box is a rigorous LB on the whole box.
"""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import json
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
HR, PR, QR = WHITE_TABLE2[4][0], WHITE_TABLE2[4][1], WHITE_TABLE2[4][2]


def cb(box, nh, npp, nq):
    return cover_min_over_box(kept, 'primal_m1e5', box[0], box[1], box[2], n_h=nh, n_p=npp, n_q=nq)


def grid_for(box, eps_target):
    (h0, h1), (p0, p1), (q0, q1) = box
    Lm = max(cb(box, 3, 3, 3)[5], 1e-9)
    budget = (eps_target / Lm) * 2.0 / np.sqrt(3.0)
    def n(ext):
        return 1 if ext <= 1e-15 else max(2, min(int(np.ceil(ext / budget)) + 1, 6000))
    return n(h1 - h0), n(p1 - p0), n(q1 - q0), Lm


def certify(box, eps_target=EPS_TARGET):
    nh, npp, nq, Lm = grid_for(box, eps_target)
    lb, pt, wit, gmin, eps, L = cb(box, nh, npp, nq)
    bump = 0
    while eps >= eps_target and bump < 3:
        nh = min(nh * 2, 6000); npp = min(npp * 2, 6000); nq = min(nq * 2, 6000)
        lb, pt, wit, gmin, eps, L = cb(box, nh, npp, nq)
        bump += 1
    return lb, pt, wit, gmin, eps, L


def main():
    print('R5 FINAL certify box h%s p%s q%s TARGET=%.7f EPS=%.1e' % (HR, PR, QR, TARGET, EPS_TARGET), flush=True)
    LOW_Q = (-0.42, 0.20)   # low band; outside this cover>>target (verified)

    # Phase A: clear q outside low band with dense-enough grid
    qA = [(-1.0, -0.7), (-0.7, -0.55), (-0.55, -0.48), (-0.48, -0.42),
          (0.20, 0.25), (0.25, 0.32), (0.32, 0.42), (0.42, 0.55), (0.55, 0.75), (0.75, 1.0)]
    floor_A = np.inf; A_rows = []
    for (q0, q1) in qA:
        box = (HR, PR, (q0, q1))
        lb, pt, wit, gmin, eps, L = certify(box, eps_target=2e-3)  # looser ok (cover >> target)
        floor_A = min(floor_A, lb); A_rows.append((q0, q1, lb, gmin, eps))
        print('  [A] q[%+.2f,%+.2f] cover_min=%.6f gmin=%.6f eps=%.2e %s' % (
            q0, q1, lb, gmin, eps, 'OK' if lb >= TARGET else '*** BELOW ***'), flush=True)

    # Phase B: dense low band. p-tiles 0.05 wide; q-tiles 0.02 wide.
    p_edges = np.linspace(0.0, 1.0, 21)
    qB_edges = np.arange(LOW_Q[0], LOW_Q[1] + 1e-9, 0.02)
    floor_B = np.inf; binding = None; low_below = []; n_tiles = 0
    for ip in range(len(p_edges) - 1):
        pb = (p_edges[ip], p_edges[ip + 1])
        row_min = np.inf
        for iq in range(len(qB_edges) - 1):
            box = (HR, pb, (qB_edges[iq], qB_edges[iq + 1]))
            # cheap pre-check: raw grid_min margin
            lb0, pt0, w0, g0, e0, L0 = cb(box, 9, 41, 41)
            if g0 > TARGET + 0.03:
                # big margin; certify lightly to a looser eps (still < margin)
                lb, pt, wit, gmin, eps, L = certify(box, eps_target=3e-3)
            else:
                lb, pt, wit, gmin, eps, L = certify(box, eps_target=EPS_TARGET)
            n_tiles += 1
            row_min = min(row_min, lb)
            if lb < floor_B:
                floor_B = lb; binding = (box, pt, wit, gmin, eps, L)
            if lb < TARGET:
                low_below.append((lb, box, pt, gmin, eps))
        print('  [B] p[%.2f,%.2f] row_min=%.6f floor_B=%.6f tiles=%d' % (pb[0], pb[1], row_min, floor_B, n_tiles), flush=True)

    certified = min(floor_A, floor_B)
    print('\n=== R5 CERTIFIED FLOOR (geometry) = %.7f ===' % certified, flush=True)
    print('    Phase A floor = %.7f ; Phase B floor = %.7f' % (floor_A, floor_B), flush=True)
    print('    vs TARGET 0.3802838: %+.6e' % (certified - TARGET), flush=True)
    if binding:
        box, pt, wit, gmin, eps, L = binding
        print('    binding tile=%s' % (box,), flush=True)
        print('    worst@(h%.4f,p%.4f,q%.4f) grid_min=%.6f eps=%.2e Lmax=%.3f wit=%s' % (
            pt[0], pt[1], pt[2], gmin, eps, L, wit), flush=True)
    print('    tiles below TARGET: %d' % len(low_below), flush=True)
    low_below.sort()
    for lb, box, pt, gmin, eps in low_below[:10]:
        print('      cover_min=%.6f box=%s worst@(h%.4f,p%.4f,q%.4f) gmin=%.6f eps=%.2e' % (
            lb, tuple(tuple(round(x, 4) for x in ax) for ax in box), pt[0], pt[1], pt[2], gmin, eps), flush=True)

    out = {'region': 5, 'certified_floor_geometry': certified, 'phaseA': floor_A, 'phaseB': floor_B,
           'target': TARGET, 'binding_tile': binding[0] if binding else None,
           'binding_worst': binding[1] if binding else None, 'binding_grid_min': binding[3] if binding else None,
           'binding_eps': binding[4] if binding else None, 'n_below_target': len(low_below),
           'low_below_target': [{'cover_min': lb, 'box': box, 'worst': pt} for lb, box, pt, g, e in low_below[:20]]}
    json.dump(out, open('_r5_final_out.json', 'w'), indent=2, default=float)
    print('\nsaved _r5_final_out.json', flush=True)


if __name__ == '__main__':
    main()
