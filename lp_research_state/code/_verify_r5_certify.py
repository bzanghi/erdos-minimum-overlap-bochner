"""Efficient rigorous certified floor for R5 via uniform tiling + per-tile eps control.

Two-phase, fully rigorous:
  Phase A (cheap clearance): tile the box coarsely in q (the dominant decay dim).
    For any q-tile whose cover_min_over_box (with a modest grid) is already
    >= TARGET, that tile is CERTIFIED-CLEARED (cover >= cover_min >= TARGET on it).
  Phase B (dense certification of the low band): for the remaining q-tiles where
    cover dips near TARGET, tile finely in (p,q) (and h) and size each tile's grid
    so eps_grid < EPS_TARGET. cover_min on each fine tile is rigorous; the min over
    all tiles is the certified floor on the low band.
  Certified floor over R5 = min(min over cleared tiles' cover_min, min over fine tiles).

All bounds use cover_min_over_box, which is a rigorous LB of max_c Phi_c over its box
(grid value minus Lipschitz cell error). Min over a partition is rigorous.
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


def lmax_of(box):
    return cb(box, 3, 3, 3)[5]


def grid_for(box, eps_target):
    (h0, h1), (p0, p1), (q0, q1) = box
    Lm = max(lmax_of(box), 1e-9)
    budget = (eps_target / Lm) * 2.0 / np.sqrt(3.0)
    def npts(ext):
        if ext <= 1e-15:
            return 1
        return max(2, min(int(np.ceil(ext / budget)) + 1, 4000))
    return npts(h1 - h0), npts(p1 - p0), npts(q1 - q0)


def certify_tile(box, eps_target=EPS_TARGET):
    """Refine grid until eps < eps_target; return rigorous cover_min for box."""
    nh, npp, nq = grid_for(box, eps_target)
    lb, pt, wit, gmin, eps, Lm = cb(box, nh, npp, nq)
    # safety: bump density if eps still above target (L_max underestimated on tiny probe)
    bump = 0
    while eps >= eps_target and bump < 4:
        nh = min(nh * 2, 4000); npp = min(npp * 2, 4000); nq = min(nq * 2, 4000)
        lb, pt, wit, gmin, eps, Lm = cb(box, nh, npp, nq)
        bump += 1
    return lb, pt, wit, gmin, eps, Lm, (nh, npp, nq)


def main():
    print('R5 certify: box h%s p%s q%s TARGET=%.7f EPS_TARGET=%.1e' % (HR, PR, QR, TARGET, EPS_TARGET), flush=True)

    # --- Phase A: coarse q-tiles, cheap clearance ---
    q_edges = np.array([-1.0, -0.7, -0.5, -0.4, -0.35, -0.3, -0.25, -0.2,
                        0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 1.0])
    # NOTE: gap [-0.2, 0.15] is the LOW BAND (handled in Phase B densely).
    cleared = []
    floor_cleared = np.inf
    for i in range(len(q_edges) - 1):
        q0, q1 = q_edges[i], q_edges[i + 1]
        if q0 >= -0.2 and q1 <= 0.15:
            continue  # low band, skip here
        box = (HR, PR, (q0, q1))
        # modest grid; these tiles have cover >> TARGET, just need cover_min >= TARGET
        lb, pt, wit, gmin, eps, Lm = cb(box, 9, 121, max(11, int((q1 - q0) / 0.01) + 2))
        cleared.append((lb, (q0, q1), gmin, eps, pt))
        floor_cleared = min(floor_cleared, lb)
        tag = 'CLEARED' if lb >= TARGET else '*** BELOW TARGET ***'
        print('  [A] q[%+.2f,%+.2f] cover_min=%.6f grid_min=%.6f eps=%.2e %s @(h%.3f,p%.3f,q%.3f)' % (
            q0, q1, lb, gmin, eps, tag, pt[0], pt[1], pt[2]), flush=True)

    # --- Phase B: dense certification of the LOW BAND q in [-0.2, 0.15] ---
    # tile (p,q); for each tile size grid for eps<EPS_TARGET. h kept whole (small, [0.08,0.1]).
    p_edges = np.linspace(0.0, 1.0, 21)        # 0.05-wide p-tiles
    qB_edges = np.linspace(-0.2, 0.15, 36)     # 0.01-wide q-tiles
    floor_band = np.inf
    binding = None
    low_below = []
    n_tiles = 0
    for ip in range(len(p_edges) - 1):
        for iq in range(len(qB_edges) - 1):
            box = (HR, (p_edges[ip], p_edges[ip + 1]), (qB_edges[iq], qB_edges[iq + 1]))
            # quick pre-check with cheap grid; if already comfortably above target, certify light
            lb0, pt0, _, g0, e0, L0 = cb(box, 5, 41, 41)
            if g0 - 0.0 > TARGET + 0.02 and lb0 > TARGET:
                # plenty of margin; lb0 already certified >= TARGET
                val = lb0; pt = pt0; gmin = g0; eps = e0; grid = (5, 41, 41)
            else:
                val, pt, wit, gmin, eps, Lm, grid = certify_tile(box)
            n_tiles += 1
            if val < floor_band:
                floor_band = val; binding = (box, pt, gmin, eps, grid)
            if val < TARGET:
                low_below.append((val, box, pt, gmin, eps))
        print('  [B] p[%.2f,%.2f] done; running floor_band=%.6f (tiles=%d)' % (
            p_edges[ip], p_edges[ip + 1], floor_band, n_tiles), flush=True)

    certified = min(floor_cleared, floor_band)
    print('\n=== R5 CERTIFIED FLOOR (geometry; min over all tiles) = %.7f ===' % certified, flush=True)
    print('    Phase A cleared-tiles floor = %.7f' % floor_cleared, flush=True)
    print('    Phase B low-band floor       = %.7f' % floor_band, flush=True)
    print('    vs TARGET 0.3802838: %+.6e' % (certified - TARGET), flush=True)
    if binding:
        box, pt, gmin, eps, grid = binding
        print('    binding tile=%s' % (box,), flush=True)
        print('    worst@(h%.4f,p%.4f,q%.4f) grid_min=%.6f eps=%.2e grid=%s' % (pt[0], pt[1], pt[2], gmin, eps, grid), flush=True)
    print('    low tiles below TARGET: %d' % len(low_below), flush=True)
    low_below.sort()
    for val, box, pt, gmin, eps in low_below[:10]:
        print('      cover_min=%.6f box=%s worst@(h%.4f,p%.4f,q%.4f) grid_min=%.6f eps=%.2e' % (
            val, tuple(tuple(round(x, 4) for x in ax) for ax in box), pt[0], pt[1], pt[2], gmin, eps), flush=True)

    out = {'region': 5, 'certified_floor_geometry': certified,
           'phaseA_cleared_floor': floor_cleared, 'phaseB_lowband_floor': floor_band,
           'target': TARGET, 'binding_tile': binding[0] if binding else None,
           'binding_worst': binding[1] if binding else None,
           'binding_grid_min': binding[2] if binding else None,
           'binding_eps': binding[3] if binding else None,
           'n_low_below_target': len(low_below)}
    json.dump(out, open('_r5_certify_out.json', 'w'), indent=2, default=float)
    print('\nsaved _r5_certify_out.json', flush=True)


if __name__ == '__main__':
    main()
