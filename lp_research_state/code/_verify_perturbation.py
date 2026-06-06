"""
Empirical rigor test of the Path-B ellipse LB direction.

White's argument: a dual-feasible point y extracted at center theta_c gives
   Phi(theta) = V_c + (b(theta) - b(theta_c)) . y
which is a VALID lower bound on the true SDP optimum V(theta) for every theta,
because y stays dual-feasible (parameters live only in the rhs b).

This is the one claim the two static codebases can't self-check (they share the
sign convention). Here we test it empirically: build the ellipse at a center,
then actually re-solve the SDP at perturbed (h,p) and confirm

        true_V(h', p')  >=  Phi(h', p')        (LB holds)
        true_V(h', p')  ~=  Phi(h', p')        (linearization is exact, tight)

A sign error would make Phi OVER-shoot true_V -> LB invalid.

Cheap scale (N=3000, bochner_n=20); the math is scale-independent.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import solve_and_extract_duals, find_ellipse_h_p

N, T, R, BN = 3000, 1200, 10, 20
h_c, p_c, q1, q2 = 0.004, 0.3875, -0.02, 0.02   # row4 center

print(f"solving center row4 (h={h_c}, p={p_c}) at N={N} bochner_n={BN} ...")
cen = solve_and_extract_duals(N, T, R, h_c, p_c, q1, q2, BN)
print(f"  V_c = {cen['value']:.8f}  ({cen['status']}, {cen['time']:.1f}s)")
ell = find_ellipse_h_p(cen, cen["duals"], q1, q2, target=0.379005)

def phi(h, p):
    return (ell["V_c"] + ell["const_q"]
            + ell["A_h2"]*h*h + ell["A_h1"]*h + ell["A_h0"]
            + ell["A_p2"]*p*p + ell["A_p1"]*p + ell["A_p0"])

# perturbed points spread across White's (h,p) box
pts = [(0.000, 0.390), (0.020, 0.400), (0.004, 0.420),
       (0.030, 0.370), (0.000, 0.450), (0.012, 0.39015)]
print(f"\n{'(h, p)':>20} {'true_V':>12} {'Phi(LB)':>12} {'true-Phi':>12}  verdict")
worst = 1e9
for (h, p) in pts:
    r = solve_and_extract_duals(N, T, R, h, p, q1, q2, BN)
    tv = r["value"]; pv = phi(h, p); gap = tv - pv
    worst = min(worst, gap)
    print(f"{('('+format(h,'.3f')+', '+format(p,'.3f')+')'):>20} "
          f"{tv:>12.7f} {pv:>12.7f} {gap:>+12.2e}  "
          f"{'LB OK' if gap > -2e-6 else '*** LB VIOLATED ***'}")
print(f"\nworst (true_V - Phi) = {worst:+.2e}")
print("LB direction VALID" if worst > -2e-6 else "LB DIRECTION FAILS")
