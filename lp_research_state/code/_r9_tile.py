"""Rigorous full-box R9 certification by uniform tiling.
Tile h(0,0.08) x p(0,1) x q(0.025,0.05) into sub-boxes; run cover_min_over_box on each
with FULL 121-center set; min over tiles is a rigorous box-min lower bound. Sub-box sizes
chosen so eps_grid < EPS_TARGET given global L_max~7.72. Reports the binding tile."""
import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import time
import _r9_verify as V
from _fullspace_eval import cover_min_over_box

TARGET = V.TARGET

def tile_certify(kept, hr, pr, qr, nh_tiles, np_tiles, nq_tiles,
                 ng_h=21, ng_p=21, ng_q=11, verbose_below=None):
    """Partition box into nh*np*nq tiles; cover_min on each; return global min tile.
    ng_* = grid points per tile per axis."""
    h_edges = np.linspace(hr[0], hr[1], nh_tiles+1)
    p_edges = np.linspace(pr[0], pr[1], np_tiles+1)
    q_edges = np.linspace(qr[0], qr[1], nq_tiles+1)
    gmin_overall = np.inf  # min envelope grid value (true floor proxy)
    floor_overall = np.inf # min certified (grid - eps) over tiles
    worst_floor = None
    worst_gmin = None
    max_eps = 0.0
    n = 0
    for i in range(nh_tiles):
        for j in range(np_tiles):
            for k in range(nq_tiles):
                shr = (h_edges[i], h_edges[i+1])
                spr = (p_edges[j], p_edges[j+1])
                sqr = (q_edges[k], q_edges[k+1])
                lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
                    kept, 'primal_m1e5', shr, spr, sqr, n_h=ng_h, n_p=ng_p, n_q=ng_q)
                n += 1
                max_eps = max(max_eps, eps)
                if gmin < gmin_overall:
                    gmin_overall = gmin; worst_gmin = (pt, wit, gmin, eps, Lm)
                if lb < floor_overall:
                    floor_overall = lb; worst_floor = (pt, wit, gmin, eps, Lm)
                if verbose_below is not None and lb < verbose_below:
                    print(f'   tile h{shr} p{spr} q{sqr}: floor={lb:.7f} gmin={gmin:.7f} eps={eps:.3e} wit={wit} worst@{pt}', flush=True)
    return floor_overall, worst_floor, gmin_overall, worst_gmin, max_eps, n
