"""Independent feasibility probe for R6 box: solve build_problem_with_dual_handles
at a (h,p,q) point and report status. Sequential, moderate N."""
import sys; sys.path.insert(0, '.')
import time
import cvxpy as cp
from path_b_analytical import build_problem_with_dual_handles


def feas(N, T, R, h, p, q, bochner_n=20):
    t0 = time.time()
    Omega, cons, handles = build_problem_with_dual_handles(
        N, T, R, h, h, p, p, q, q, bochner_n=bochner_n)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception as e:
        return ('SOLVE_ERROR:%s' % type(e).__name__, None, time.time() - t0)
    val = prob.value if prob.value is not None else float('nan')
    return (prob.status, val, time.time() - t0)


if __name__ == '__main__':
    import json
    pts = json.loads(sys.argv[1])  # list of [h,p,q]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    T = int(sys.argv[3]) if len(sys.argv) > 3 else 1600
    bn = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    for (h, p, q) in pts:
        st, val, dt = feas(N, T, 10, h, p, q, bochner_n=bn)
        vs = ('%.6f' % val) if (val is not None and val == val) else 'none/inf'
        print('h=%.4f p=%.4f q=%.4f -> %-20s val=%s (%.1fs)' % (h, p, q, st, vs, dt), flush=True)
