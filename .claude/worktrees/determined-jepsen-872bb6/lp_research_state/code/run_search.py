"""
Run multistart refinement at larger N using the fast evaluator.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from fast_eval import smooth_max_grad_fast, J_all
from evaluator import bound_float
from refiner import project_to_constraints
import time, json, sys


def penalty_obj(h, N, β, λ):
    sm, g = smooth_max_grad_fast(h, β)
    target = N / 2.0
    s = float(h.sum())
    pen = 0.5 * λ * (s - target) ** 2
    pen_g = λ * (s - target) * np.ones(N)
    return sm + pen, g + pen_g


def refine_fast(h0, β_schedule, λ=200.0, iters=600, symmetric=False):
    h = h0.copy()
    N = len(h)
    if symmetric:
        h = 0.5 * (h + h[::-1])
    h = project_to_constraints(h)
    half = (N + 1) // 2
    for β in β_schedule:
        if symmetric:
            def to_full(v):
                f = np.empty(N)
                f[:half] = v
                f[N - half:] = v[::-1]
                if N % 2 == 1:
                    f[half - 1] = v[half - 1]
                return f
            def obj(v):
                f = to_full(v)
                val, gf = penalty_obj(f, N, β, λ)
                gv = gf[:half] + gf[N - half:][::-1]
                if N % 2 == 1:
                    gv[half - 1] /= 2
                return val, gv
            v0 = h[:half].copy()
            res = minimize(obj, v0, jac=True, method="L-BFGS-B",
                           bounds=[(0.0, 1.0)] * half,
                           options={"maxiter": iters, "ftol": 1e-14, "gtol": 1e-12})
            h = to_full(res.x)
        else:
            res = minimize(penalty_obj, h, args=(N, β, λ),
                           jac=True, method="L-BFGS-B",
                           bounds=[(0.0, 1.0)] * N,
                           options={"maxiter": iters, "ftol": 1e-14, "gtol": 1e-12})
            h = res.x
        h = project_to_constraints(h)
    return h


def random_starts(N, n_starts, seed=0):
    rng = np.random.default_rng(seed)
    starts = []
    xs = np.linspace(1.0/N, 2 - 1.0/N, N)
    starts.append(np.minimum(xs, 2-xs)*(N/2)/np.minimum(xs, 2-xs).sum())
    starts.append(0.5*np.ones(N))
    h = np.zeros(N); h[N//4:3*N//4] = 1; starts.append(h)
    starts.append(0.5*(1 - np.cos(np.pi*xs)))
    # 0/1 with two plateaus at width 0.6 each (heuristic)
    h = np.where((xs > 0.3) & (xs < 1.7), 1.0, 0.0); starts.append(h*(N/2)/h.sum())
    # symmetric 'M' shape
    h = np.where((xs > 0.5) & (xs < 1.5), 0.7, 0.4); starts.append(project_to_constraints(h))
    for _ in range(n_starts):
        starts.append(project_to_constraints(0.5 + 0.4*rng.standard_normal(N)))
    return starts


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    n_starts = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    schedule = (10.0, 40.0, 200.0, 1000.0, 5000.0, 20000.0)
    starts = random_starts(N, n_starts, seed=42)
    best = (1e9, None, None)
    t0 = time.time()
    for k, h0 in enumerate(starts):
        t1 = time.time()
        h = refine_fast(h0, schedule, symmetric=True)
        b, j, s = bound_float(h)
        dt = time.time() - t1
        marker = ""
        if b < best[0]:
            best = (b, h.copy(), (j, s))
            marker = "  *"
        print(f"  start {k:2d}  bound = {b:.7f}  argmax j={j} ({s})  t={dt:5.1f}s{marker}")
    print(f"\nBEST (N={N}, symmetric): {best[0]:.7f}  argmax j={best[2][0]} ({best[2][1]})")
    print(f"total time: {time.time()-t0:.1f}s")
    out = {
        "N": N,
        "bound": best[0],
        "argmax_j": best[2][0],
        "argmax_sign": best[2][1],
        "h": best[1].tolist(),
    }
    with open(f"best_N{N}.json", "w") as f:
        json.dump(out, f)
    print(f"saved best_N{N}.json")
