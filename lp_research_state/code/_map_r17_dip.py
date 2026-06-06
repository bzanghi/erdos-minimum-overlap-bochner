import warnings; warnings.filterwarnings('ignore')
import sys, json; sys.path.insert(0, '.')
import numpy as np
import _fs_recon as R
from _fullspace_eval import phi_center, anchor_value, CORE_HEADLINE

NEED = ('con_53', 'con_54', 'con_512_pL', 'con_512_pU', 'con_512_qL', 'con_512_qU', 'con_513')
core, src = R.load_core(); halo = R.load_halo()
core_halo = core + halo


def cover_at(cs, h, p, q):
    return max(phi_center(c, anchor_value(c, 'primal_m1e5'), h, p, q) for c in cs)


# Map the cover along p at the worst h=0, q=0.025 (worst |q|) and find where it dips below target
print('Baseline corehalo cover along p (h=0, q=0.025): dip extent below %.6f' % CORE_HEADLINE)
ps = np.linspace(0.33, 0.45, 241)
below = []
for p in ps:
    v = cover_at(core_halo, 0.0, p, 0.025)
    if v < CORE_HEADLINE:
        below.append((p, v))
if below:
    lo = min(b[0] for b in below); hi = max(b[0] for b in below)
    worst = min(below, key=lambda b: b[1])
    print('  dip p-range = [%.4f, %.4f]  (width %.4f)  worst Phi=%.7f at p=%.4f'
          % (lo, hi, hi - lo, worst[1], worst[0]))
else:
    print('  NO dip at q=0.025 along h=0')

# Also check at h=0, scanning q in [0.02,0.025] for the worst over the dip band
print('Worst over (h=0, p in dip band, q in [0.02,0.025]):')
worst = 1e9; wp = None
for p in np.linspace(0.385, 0.397, 61):
    for q in np.linspace(0.02, 0.025, 11):
        v = cover_at(core_halo, 0.0, p, q)
        if v < worst:
            worst = v; wp = (p, q)
print('  worst=%.7f at p=%.4f q=%.4f' % (worst, wp[0], wp[1]))
