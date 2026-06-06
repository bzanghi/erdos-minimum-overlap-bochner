import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import _fs_recompute as M

hv = M.harvest_centers()
centers = hv[0]
kept = []
for c in centers:
    if c.get('primal') is None:
        c = dict(c); c['primal'] = float(c['dual_lb']) + 1e-5
    kept.append(c)

# On a small subbox, compute per-center gradient contributions to L_max
h0, h1 = 0.08, 0.09; p0, p1 = 0.3, 0.4; q0, q1 = -0.15, -0.1
def gmax(lin_const, quad_coeff, lo, hi):
    return max(abs(lin_const + quad_coeff * lo), abs(lin_const + quad_coeff * hi))
rows = []
for c in kept:
    d = c['duals']
    gh = gmax(d['con_53'], -d['con_54'], h0, h1)
    gp = gmax(d['con_512_pL'] - d['con_512_pU'], -d['con_513'], p0, p1)
    gq = gmax(d['con_512_qL'] - d['con_512_qU'], -d['con_513'], q0, q1)
    L = float(np.sqrt(gh * gh + gp * gp + gq * gq))
    rows.append((L, c['label'], gh, gp, gq, d))
rows.sort(reverse=True)
print('Top centers driving L_max on small subbox h[0.08,0.09] p[0.3,0.4] q[-0.15,-0.1]:')
for L, lab, gh, gp, gq, d in rows[:8]:
    print('  L=%8.3f %-30s gh=%7.3f gp=%9.3f gq=%7.3f' % (L, lab, gh, gp, gq))
    print('       con_53=%.3f con_54=%.3f con_512_pL=%.3f con_512_pU=%.3f con_513=%.3f con_512_qL=%.3f con_512_qU=%.3f' % (
        d['con_53'], d['con_54'], d['con_512_pL'], d['con_512_pU'], d['con_513'], d['con_512_qL'], d['con_512_qU']))
print()
print('L_max over all = %.4f driven by %s' % (rows[0][0], rows[0][1]))
hd = 0.5 * np.sqrt((h1 - h0)**2 + (p1 - p0)**2 + (q1 - q0)**2)
print('half_diag(box)=%.4f -> eps=%.4f' % (hd, rows[0][0] * hd))

# Now: which center is the WITNESS (argmax of cover) at the worst point, and what is
# ITS gradient? The eps penalty should really only need the Lipschitz const of the
# ACTIVE envelope near the min, but cover_min_over_box uses max over ALL centers.
print()
print('=== gradient of the WINNING center near worst point matters less; but check ===')
# Evaluate which center wins at worst point (0.08, 0.324, -0.108)
from _fullspace_eval import phi_center_grid, anchor_value
hw, pw, qw = 0.08, 0.324, -0.108
best = (-1e9, None)
for c in kept:
    a = anchor_value(c, 'primal_m1e5')
    F = float(phi_center_grid(c, a, np.array([[hw]]), np.array([[pw]]), qw)[0, 0])
    if F > best[0]:
        best = (F, c['label'])
print('winner at worst pt (%.3f,%.3f,%.3f): %s with Phi=%.6f' % (hw, pw, qw, best[1], best[0]))
