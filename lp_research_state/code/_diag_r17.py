import warnings; warnings.filterwarnings('ignore')
import sys, json; sys.path.insert(0, '.')
import numpy as np
import _fs_recon as R
from _fullspace_eval import phi_center, anchor_value
from path_b_analytical import dual_objective_shift

NEED = ('con_53', 'con_54', 'con_512_pL', 'con_512_pU', 'con_512_qL', 'con_512_qU', 'con_513')
core, src = R.load_core(); halo = R.load_halo()
fr = json.load(open('../parallel_results/fullspace_promote_R17.json'))
fresh = [{'label': c['label'], 'h_c': c['h_c'], 'p_c': c['p_c'], 'q1': c['q1'], 'q2': c['q2'],
          'primal': c['primal'], 'duals': {k: float(c['duals'][k]) for k in NEED}} for c in fr['centers']]
cde3 = [c for c in core if c['label'] == 'cde_n30_iter3'][0]

for (h, p, q) in [(0.0, 0.3915, 0.02438), (0.0, 0.3915, 0.025), (0.003, 0.39, 0.025)]:
    print('=== point (h=%g,p=%g,q=%g) ===' % (h, p, q))
    a = anchor_value(cde3, 'primal_m1e5')
    print('  cde_n30_iter3            anchor=%.7f Phi=%.7f' % (a, phi_center(cde3, a, h, p, q)))
    for c in fresh:
        a = anchor_value(c, 'primal_m1e5')
        print('  %-28s anchor=%.7f Phi=%.7f' % (c['label'], a, phi_center(c, a, h, p, q)))

c = fresh[1]; center = {k: c[k] for k in ('h_c', 'p_c', 'q1', 'q2')}; d = c['duals']
h, p, q = 0.0, 0.3915, 0.02438
print('--- shift decomposition fresh c0.390 at worst pt (h=0,p=0.3915,q=0.02438) ---')
print('  53: Drhs=%g lam=%.4f term=%g' % (h - c['h_c'], d['con_53'], d['con_53'] * (h - c['h_c'])))
print('  54: Drhs=%g lam=%.4f term=%g' % ((h**2 - c['h_c']**2) / 2, d['con_54'], -d['con_54'] * ((h**2 - c['h_c']**2) / 2)))
print('  pL: Drhs=%g lam=%.4f  pU: Drhs=%g lam=%.4f' % (p - c['p_c'], d['con_512_pL'], p - c['p_c'], d['con_512_pU']))
print('  qL: Drhs=%g lam=%.4f  qU: Drhs=%g lam=%.4f' % (q - c['q1'], d['con_512_qL'], q - c['q2'], d['con_512_qU']))
print('  513: q-part Drhs=%g lam=%.4f' % (-0.5 * (q * q - max(c['q1']**2, c['q2']**2)), d['con_513']))
print('  total shift=%g' % dual_objective_shift(h, p, q, q, center, d))
