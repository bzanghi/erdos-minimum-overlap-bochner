import sys, warnings; warnings.filterwarnings('ignore'); sys.path.insert(0, '.')
import numpy as np, cvxpy as cp
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints


def feas(h, p, q, N=2000, T=800, R=10, bn=20, pmk=20):
    Om, cons, H = build_problem_with_dual_handles(N, T, R, h, h, p, p, q, q, bochner_n=bn)
    pc, _ = build_even_moment_nonneg_constraints(H['c'], H['d'], T, k_max=pmk)
    cons.extend(pc)
    prob = cp.Problem(cp.Minimize(Om), cons)
    try:
        prob.solve(solver="CLARABEL", verbose=False)
        st = prob.status
        v = float(prob.value) if (prob.value is not None and np.isfinite(prob.value)) else None
    except Exception as e:
        st = "SOLVER_ERROR:" + type(e).__name__
        v = None
    return st, v


if __name__ == "__main__":
    for p in [0.0, 0.5, 1.0]:
        for q in [0.3, 0.35, 0.40, 0.42, 0.44, 0.45, 0.5, 0.55, 0.6, 0.65]:
            st, v = feas(0.0, p, q)
            print(f"RESULT p={p} q={q:.3f} {st:30s} val={v}", flush=True)
