import warnings; warnings.filterwarnings('ignore')
import sys, json; sys.path.insert(0, '.')
import numpy as np
import _fs_recon as R
from _fullspace_eval import cover_min_over_box, phi_center, anchor_value, CORE_HEADLINE, WHITE_OUTSIDE_FLOOR

NEED = ('con_53', 'con_54', 'con_512_pL', 'con_512_pU', 'con_512_qL', 'con_512_qU', 'con_513')
core, src = R.load_core(); halo = R.load_halo()
core_halo = core + halo
fr = json.load(open('../parallel_results/fullspace_promote_R17.json'))
fresh = [{'label': c['label'], 'h_c': c['h_c'], 'p_c': c['p_c'], 'q1': c['q1'], 'q2': c['q2'],
          'primal': c['primal'], 'duals': {k: float(c['duals'][k]) for k in NEED}} for c in fr['centers']]
combined = core_halo + fresh
print('core_halo=%d fresh=%d combined=%d' % (len(core_halo), len(fresh), len(combined)))

hr = (0.0, 0.06); pr = (0.33, 0.45); qr = (0.02, 0.025)


def show(name, cs, n_h=41, n_p=81, n_q=41):
    lb, pt, wit, gm, eps, Lm = cover_min_over_box(cs, 'primal_m1e5', hr, pr, qr, n_h=n_h, n_p=n_p, n_q=n_q)
    print('%-44s phi_min=%.7f grid_min=%.7f eps=%.2e Lmax=%.4f worst=(h=%.5f,p=%.5f,q=%.5f) wit=%s c380284=%s'
          % (name, lb, gm, eps, Lm, pt[0], pt[1], pt[2], wit, lb >= CORE_HEADLINE))
    return lb, gm, eps, Lm, pt, wit


print('--- single-box cover ---')
show('R17 corehalo (baseline)', core_halo)
show('R17 combined (41x81x41)', combined)
show('R17 combined (121x241x121)', combined, n_h=121, n_p=241, n_q=121)


# ---------------- Rigorous SUBDIVIDED box-min cover ----------------
# Split R17 into a grid of sub-boxes. For each sub-box, call cover_min_over_box
# (which uses the LOCAL L_max over that sub-box's corners and the sub-box half-diag
# for eps_grid). Take the min over sub-boxes. This is rigorous: each sub-box LB is a
# valid box-min LB; min over a partition is the global box-min LB. Subdivision shrinks
# the per-box half-diagonal => eps_grid shrinks, while grid_min within each sub-box
# is unchanged. Concentrate subdivision in p (the binding dimension).
def subdivided_min(cs, nh_sub, np_sub, nq_sub, gh=21, gp=21, gq=11, verbose_worst=True):
    h_edges = np.linspace(hr[0], hr[1], nh_sub + 1)
    p_edges = np.linspace(pr[0], pr[1], np_sub + 1)
    q_edges = np.linspace(qr[0], qr[1], nq_sub + 1)
    gmin_all = np.inf
    box_lb_min = np.inf
    worst = None
    for i in range(nh_sub):
        for jx in range(np_sub):
            for k in range(nq_sub):
                sb_h = (h_edges[i], h_edges[i + 1])
                sb_p = (p_edges[jx], p_edges[jx + 1])
                sb_q = (q_edges[k], q_edges[k + 1])
                lb, pt, wit, gm, eps, Lm = cover_min_over_box(
                    cs, 'primal_m1e5', sb_h, sb_p, sb_q, n_h=gh, n_p=gp, n_q=gq)
                gmin_all = min(gmin_all, gm)
                if lb < box_lb_min:
                    box_lb_min = lb
                    worst = (lb, gm, eps, Lm, pt, wit, sb_h, sb_p, sb_q)
    return box_lb_min, gmin_all, worst


print('\n--- SUBDIVIDED box-min cover (divide-and-conquer eps_grid) ---')
for (nh, npb, nq) in [(2, 12, 1), (3, 24, 1), (3, 48, 2), (6, 60, 5)]:
    blb, gall, w = subdivided_min(combined, nh, npb, nq)
    lb, gm, eps, Lm, pt, wit, sh, sp, sq = w
    print('subdiv h%d p%d q%d: box_min_LB=%.7f (min grid=%.7f) c380284=%s | worst sub-box h%s p%s q%s eps=%.2e Lmax=%.3f wit=%s @ (h=%.5f,p=%.5f,q=%.5f)'
          % (nh, npb, nq, blb, gall, blb >= CORE_HEADLINE,
             tuple(round(x, 4) for x in sh), tuple(round(x, 4) for x in sp), tuple(round(x, 5) for x in sq),
             eps, Lm, wit, pt[0], pt[1], pt[2]))
