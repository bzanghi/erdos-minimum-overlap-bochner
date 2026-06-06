"""R8 low-q value map: primal at the binding low-q edge across p, at h=0 and h=0.08.
This is the zone where the cover is tightest (high-q is easy: primal rises steeply;
high-q corners are infeasible). Cheap config; just need primal magnitude + boundary."""
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
        st = "SOLVER_ERROR"; v = None
    return st, v


if __name__ == "__main__":
    # low-q binding edge: q in {0.05,0.08,0.12}, p across [0,1], h at both extremes
    for h in [0.0, 0.08]:
        for p in [0.0, 0.25, 0.5, 0.75, 1.0]:
            row = []
            for q in [0.05, 0.08, 0.12, 0.18]:
                st, v = feas(h, p, q)
                row.append(f"q{q}={'INF' if v is None else format(v,'.4f')}({st[:4]})")
            print(f"RESULT h={h} p={p}: " + "  ".join(row), flush=True)
