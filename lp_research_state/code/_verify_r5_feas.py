import sys; sys.path.insert(0, '.')
import warnings; warnings.filterwarnings('ignore')
import time
import numpy as np
import cvxpy as cp
from path_b_analytical import build_problem_with_dual_handles

def feas(h, p, q, N=4000, T=1600, R=10, bn=20):
    Omega, cons, H = build_problem_with_dual_handles(N, T, R, h, h, p, p, q, q, bochner_n=bn)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    try:
        prob.solve(solver='CLARABEL', verbose=False)
    except Exception as e:
        return ('ERROR:' + str(e)[:40], None, time.time() - t0)
    return (prob.status, None if prob.value is None else float(prob.value), time.time() - t0)

if __name__ == '__main__':
    import json
    pts = eval(sys.argv[1]) if len(sys.argv) > 1 else [
        (0.08, 0.324, -0.11),   # binding worst point of cover
        (0.1, 0.0, -1.0),       # candidate's claimed worst corner
        (0.08, 0.0, -1.0),
        (0.1, 1.0, 1.0),
        (0.08, 0.5, -1.0),
        (0.1, 0.5, 1.0),
        (0.09, 0.3, 0.0),
    ]
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    out = []
    for (h, p, q) in pts:
        st, val, el = feas(h, p, q, N=N)
        out.append({'h': h, 'p': p, 'q': q, 'status': st, 'val': val, 't': round(el, 1)})
        print('h=%.4f p=%.4f q=%+.4f -> %-20s val=%s  (%.1fs)' % (
            h, p, q, st, ('%.5f' % val) if val is not None else 'None', el), flush=True)
    print('JSON:' + json.dumps(out))
