"""Rigorous box-min LB of the augmented Phi-cover over GATE region R6, via ADAPTIVE
subdivision (same mechanism as _eval_r7_box.py) so the Lipschitz eps_grid is driven
down LOCALLY where the cover is tight.

R6 = White Table-2 region 6: h(E(M)) in [0,0.08], p(c1) in [0,1], q(d1) in [-1,-0.05].
width_class WIDE. The binding locus is the q=-0.05 boundary strip adjacent to the core
(p~0.33-0.45), where the program value is intrinsically ~0.381 so the cover ceiling
(grid_min) is ~0.3808 -- a thin ~5e-4 headroom over the 0.380284 target. A single global
grid in cover_min_over_box cannot drive eps_grid below that headroom (the wide q-axis
dominates half_diag and the global L_max is inflated by the deep-q centers). Adaptive
subdivision recomputes L_max per sub-box and shrinks eps where it matters.

NOTE on the infeasible deep-q region (|d1|>~0.45-0.55, verified infeasible at N=2500 and
N=5000): no admissible function lives there, so it contributes nothing to mu. Our fresh
centers (built with steep con_513) happen to extrapolate the cover UP to 0.5-0.86 there,
so those sub-boxes clear the target with huge margin anyway -- the adaptive method handles
them trivially. The honest binding region is the feasible shallow-q strip.

Composes: 12 core anchors + fresh R6 centers (this task), conservative primal_m1e5 anchor.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _fullspace_eval import load_centers, cover_min_over_box, CORE_HEADLINE, WHITE_OUTSIDE_FLOOR

PR = CODE.parent / "parallel_results"
PROMOTE = PR / "fullspace_promote_R6.json"
# Centers dropped to keep L_max controlled (spiky con_513, cover <2% of grid):
DROP_DEFAULT = ("R6_c2_p05_q05", "R6_c3_p05_q20")


def build_combined(drop=DROP_DEFAULT):
    existing, _ = load_centers()
    fresh = []
    if PROMOTE.exists():
        fresh = [c for c in json.load(open(PROMOTE))["centers"]
                 if c["label"] not in drop and c.get("dual_lb") is not None]
        # conservative anchor convention: anchor_value(.,primal_m1e5) = primal - 1e-5;
        # the stored 'primal' is the fresh solve primal -> matches directly.
    return existing, fresh, existing + fresh


def adaptive_boxmin(centers, h_range, p_range, q_range, target,
                    depth=0, max_depth=8, base=(21, 21, 21)):
    """Rigorous box-min LB over the (h,p,q) box. Subdivide the longest axis until
    lb>=target (cleared) or max_depth. Returns (lb, grid_min, worst_pt, worst_wit, n_leaves)."""
    nh, npp, nq = base
    lb, pt, wit, gm, eps, L = cover_min_over_box(
        centers, "primal_m1e5", h_range, p_range, q_range, n_h=nh, n_p=npp, n_q=nq)
    if lb >= target or depth >= max_depth:
        return lb, gm, pt, wit, 1
    spans = [(h_range[1]-h_range[0], 'h'), (p_range[1]-p_range[0], 'p'),
             (q_range[1]-q_range[0], 'q')]
    spans.sort(reverse=True)
    axis = spans[0][1]
    def split(rng):
        m = 0.5 * (rng[0] + rng[1]); return (rng[0], m), (m, rng[1])
    worst_lb = np.inf; worst = None; tot = 0
    if axis == 'h':
        parts = [(a, p_range, q_range) for a in split(h_range)]
    elif axis == 'p':
        parts = [(h_range, a, q_range) for a in split(p_range)]
    else:
        parts = [(h_range, p_range, a) for a in split(q_range)]
    for hr, pr, qr in parts:
        l, g, pp, ww, n = adaptive_boxmin(centers, hr, pr, qr, target,
                                          depth+1, max_depth, base)
        tot += n
        if l < worst_lb:
            worst_lb = l; worst = (g, pp, ww)
    return worst_lb, worst[0], worst[1], worst[2], tot


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    existing, fresh, combined = build_combined()
    print(f"centers: core={len(existing)} fresh(R6, spiky dropped)={len(fresh)} combined={len(combined)}")
    for c in fresh:
        print(f"  fresh {c['label']}: anchor={c['primal']-1e-5:.6f} resid={c['dual_resid']:.1e}")
    H = (0.0, 0.08); P = (0.0, 1.0); Q = (-1.0, -0.05)
    # Seed: split q into bands (the hard axis) and p into bands (binding near core-p).
    # Coarsely partition; adaptive_boxmin refines each.
    qedges = [-1.0, -0.55, -0.45, -0.30, -0.20, -0.12, -0.08, -0.05]
    pedges = [0.0, 0.25, 0.33, 0.45, 0.6, 1.0]
    print("\n=== Adaptive box-min over R6 (seeded q x p bands) ===")
    worst = np.inf; worst_box = None; below_floor = []; below_target = []
    for qi in range(len(qedges)-1):
        q0, q1 = qedges[qi], qedges[qi+1]
        for pi in range(len(pedges)-1):
            p0, p1 = pedges[pi], pedges[pi+1]
            lb, gm, pt, wit, nl = adaptive_boxmin(combined, H, (p0, p1), (q0, q1),
                                                  CORE_HEADLINE, max_depth=8)
            if lb < worst:
                worst = lb; worst_box = (p0, p1, q0, q1, gm, pt, wit)
            if lb < WHITE_OUTSIDE_FLOOR:
                below_floor.append((p0, p1, q0, q1, lb))
            if lb < CORE_HEADLINE:
                below_target.append((p0, p1, q0, q1, lb, gm))
            tag = "OK>=CORE" if lb >= CORE_HEADLINE else (
                ">=FLOOR" if lb >= WHITE_OUTSIDE_FLOOR else "<FLOOR!")
            print(f"  q[{q0:.3f},{q1:.3f}] p[{p0:.3f},{p1:.3f}] lb={lb:.6f} "
                  f"ceil(gm)={gm:.6f} leaves={nl} {tag}")
    print(f"\nR6 independent floor (worst sub-box lb) = {worst:.6f}")
    print(f"   worst box q[{worst_box[2]},{worst_box[3]}] p[{worst_box[0]},{worst_box[1]}] "
          f"ceiling grid_min={worst_box[4]:.6f} worst@(h={worst_box[5][0]:.4f},"
          f"p={worst_box[5][1]:.4f},q={worst_box[5][2]:.4f})")
    print(f"   clears 0.380000 independently: {worst >= WHITE_OUTSIDE_FLOOR}")
    print(f"   clears 0.380284 independently: {worst >= CORE_HEADLINE}")
    if below_floor:
        print(f"   BELOW 0.380000 bands: {[(round(b[0],3),round(b[1],3),round(b[2],3),round(b[3],3),round(b[4],6)) for b in below_floor]}")
    if below_target:
        print(f"   below 0.380284 bands (count={len(below_target)}):")
        for b in below_target:
            print(f"      q[{b[2]:.3f},{b[3]:.3f}] p[{b[0]:.3f},{b[1]:.3f}] lb={b[4]:.6f} ceil={b[5]:.6f}")
