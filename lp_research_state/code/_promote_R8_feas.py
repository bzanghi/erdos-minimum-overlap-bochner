"""
PRO-38 R8 — feasibility-boundary sweep in q (=d1).

The augmented program enforces |d| <= 2/pi (path_b_analytical.py:115), so d1=d[0]
cannot exceed 2/pi ~ 0.6366. Above that the SDP is infeasible: there is NO function
f:[-1,1]->[0,1] with first sine coefficient that large, so those (E(M),c1,d1) points
are unrealizable and vacuously covered by ANY bound. This sweep maps the exact feasible
q-boundary and confirms CLARABEL's infeasible certificate is stable (not solver noise),
so we know which q-sub-strip actually needs fresh centers.

Cheap config (feasibility only needs a coarse solve): N=2000, T=800, bn=20, pm_k=20.
We sweep q at a few p values (the |d|<=2/pi box is p-independent, but check anyway).
"""
from __future__ import annotations
import sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints


def feas_status(h_c, p_c, q, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q, q, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver="CLARABEL", verbose=False)
    return prob.status, (float(prob.value) if prob.value is not None and np.isfinite(prob.value) else None), time.time()-t0


def main():
    N, T, R, bn, pm_k = 2000, 800, 10, 20, 20
    print(f"=== R8 feasibility sweep (N={N} T={T} bn={bn}) ; 2/pi={2/np.pi:.6f} ===\n", flush=True)
    for p_c in [0.0, 0.5, 1.0]:
        print(f"-- p_c={p_c} --", flush=True)
        for q in [0.30, 0.50, 0.60, 0.63, 0.635, 0.637, 0.64, 0.65, 0.70, 0.80, 1.0]:
            st, val, dt = feas_status(0.0, p_c, q, N, T, R, bn, pm_k)
            print(f"  q={q:.3f}: {st:22s} val={val} ({dt:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
