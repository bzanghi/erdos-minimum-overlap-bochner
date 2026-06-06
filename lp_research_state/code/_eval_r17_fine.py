import warnings; warnings.filterwarnings('ignore')
import sys, json; sys.path.insert(0, '.')
import numpy as np
import _fs_recon as R
from _fullspace_eval import cover_min_over_box, CORE_HEADLINE

NEED = ('con_53', 'con_54', 'con_512_pL', 'con_512_pU', 'con_512_qL', 'con_512_qU', 'con_513')
core, src = R.load_core(); halo = R.load_halo()
core_halo = core + halo
fr = json.load(open('../parallel_results/fullspace_promote_R17.json'))
fresh = [{'label': c['label'], 'h_c': c['h_c'], 'p_c': c['p_c'], 'q1': c['q1'], 'q2': c['q2'],
          'primal': c['primal'], 'duals': {k: float(c['duals'][k]) for k in NEED}} for c in fr['centers']]
combined = core_halo + fresh
print('combined=%d centers' % len(combined))

hr = (0.0, 0.06); pr = (0.33, 0.45); qr = (0.02, 0.025)


def subdivided_min(cs, nh_sub, np_sub, nq_sub, gh=15, gp=15, gq=9):
    h_edges = np.linspace(hr[0], hr[1], nh_sub + 1)
    p_edges = np.linspace(pr[0], pr[1], np_sub + 1)
    q_edges = np.linspace(qr[0], qr[1], nq_sub + 1)
    box_lb_min = np.inf; gmin_all = np.inf; eps_at_worst = None; worst = None
    for i in range(nh_sub):
        for jx in range(np_sub):
            for k in range(nq_sub):
                sb_h = (h_edges[i], h_edges[i + 1]); sb_p = (p_edges[jx], p_edges[jx + 1]); sb_q = (q_edges[k], q_edges[k + 1])
                lb, pt, wit, gm, eps, Lm = cover_min_over_box(cs, 'primal_m1e5', sb_h, sb_p, sb_q, n_h=gh, n_p=gp, n_q=gq)
                gmin_all = min(gmin_all, gm)
                if lb < box_lb_min:
                    box_lb_min = lb; worst = (lb, gm, eps, Lm, pt, wit, sb_h, sb_p, sb_q)
    return box_lb_min, gmin_all, worst


print('\n--- finer subdivision until clear ---')
results = []
for (nh, npb, nq, gp) in [(20, 360, 10, 15), (24, 480, 12, 17), (30, 600, 15, 21), (40, 800, 20, 21)]:
    blb, gall, w = subdivided_min(combined, nh, npb, nq, gp=gp)
    lb, gm, eps, Lm, pt, wit, sh, sp, sq = w
    cleared = blb >= CORE_HEADLINE
    results.append((nh, npb, nq, blb, gall, eps, cleared))
    print('subdiv h%d p%d q%d (gp=%d): box_min_LB=%.8f (min grid=%.8f) MARGIN=%+.2e c380284=%s | eps=%.2e Lmax=%.3f wit=%s @ (h=%.5f,p=%.5f,q=%.5f)'
          % (nh, npb, nq, gp, blb, gall, blb - CORE_HEADLINE, cleared, eps, Lm, wit, pt[0], pt[1], pt[2]))
    if cleared:
        print('   >>> CLEARED 0.380284 at this subdivision <<<')
        break
