"""
Basin-hopping refinement: alternate L-BFGS with random kicks to escape
local minima.  Uses the fast evaluator.

Also: a sequential LP refinement step using scipy.optimize.linprog,
linearizing the bilinear form around the current iterate within a trust
region.
"""
from __future__ import annotations
import numpy as np, json, time, sys
from scipy.optimize import minimize, linprog
from fast_eval import smooth_max_grad_fast, J_all
from evaluator import bound_float
from refiner import project_to_constraints


def penalty_obj(h, N, β, λ):
    sm, g = smooth_max_grad_fast(h, β)
    target = N / 2.0
    s = float(h.sum())
    pen = 0.5 * λ * (s - target) ** 2
    pen_g = λ * (s - target) * np.ones(N)
    return sm + pen, g + pen_g


def lbfgs_anneal(h0, schedule, λ, iters, symmetric):
    h = h0.copy()
    N = len(h)
    if symmetric:
        h = 0.5 * (h + h[::-1])
    h = project_to_constraints(h)
    half = (N + 1) // 2
    for β in schedule:
        if symmetric:
            def to_full(v):
                f = np.empty(N); f[:half] = v; f[N-half:] = v[::-1]
                if N % 2 == 1: f[half-1] = v[half-1]
                return f
            def obj(v):
                f = to_full(v)
                val, gf = penalty_obj(f, N, β, λ)
                gv = gf[:half] + gf[N-half:][::-1]
                if N % 2 == 1: gv[half-1] /= 2
                return val, gv
            v0 = h[:half].copy()
            res = minimize(obj, v0, jac=True, method="L-BFGS-B",
                           bounds=[(0.,1.)]*half,
                           options={"maxiter": iters, "ftol": 1e-15, "gtol": 1e-13})
            h = to_full(res.x)
        else:
            res = minimize(penalty_obj, h, args=(N, β, λ), jac=True,
                           method="L-BFGS-B", bounds=[(0.,1.)]*N,
                           options={"maxiter": iters, "ftol": 1e-15, "gtol": 1e-13})
            h = res.x
        h = project_to_constraints(h)
    return h


def slp_step(h, trust=0.05):
    """One step of sequential LP: linearize bilinear constraints at h.
       Solve  min c  s.t.  for each (j, σ):
           J^σ_j(h) + ∇J^σ_j(h) · (h' - h) ≤ c
           0 ≤ h' ≤ 1, sum h' = N/2, |h' - h|_∞ ≤ trust
    """
    N = len(h)
    Δ = 2.0 / N
    # Build constraint matrix A_ub h' + (-1)*c ≤ b_ub  (one row per (j, σ)).
    rows = []
    bs = []
    Jp_now, Jm_now = J_all(h)
    for j in range(1, N):
        # ∇J^+_j: from earlier analysis,
        #   (∂J^+_j/∂h_k) = Δ [1_{k>j} - h_{k-j}·1_{k>j} - h_{k+j}·1_{k+j ≤ N}]
        gp = np.zeros(N)
        for k in range(1, N + 1):
            v = 0.0
            if k > j:
                v += 1.0 - h[k - j - 1]   # h_{k-j} 1-indexed → 0-idx k-j-1
            if k + j <= N:
                v += -h[k + j - 1]
            gp[k - 1] = Δ * v
        # constraint: gp · h' - c ≤ -(J_now - gp · h)
        rows.append(np.concatenate([gp, [-1.0]]))
        bs.append(-(Jp_now[j - 1] - float(gp @ h)))
        # ∇J^-_j: (∂J^-_j/∂h_k) = Δ [1_{k+j ≤ N} - h_{k-j}·1_{k>j} - h_{k+j}·1_{k+j ≤ N}]
        gm = np.zeros(N)
        for k in range(1, N + 1):
            v = 0.0
            if k + j <= N:
                v += 1.0 - h[k + j - 1]
            if k > j:
                v += -h[k - j - 1]
            gm[k - 1] = Δ * v
        rows.append(np.concatenate([gm, [-1.0]]))
        bs.append(-(Jm_now[j - 1] - float(gm @ h)))
    A_ub = np.array(rows)
    b_ub = np.array(bs)
    # Equality: sum h' = N/2  → row of ones with c-coeff 0.
    A_eq = np.zeros((1, N + 1)); A_eq[0, :N] = 1.0
    b_eq = np.array([N / 2.0])
    # Bounds: 0 ≤ h'_i ≤ 1, c free.  Trust region: |h'_i - h_i| ≤ trust.
    bounds = [(max(0.0, h[i] - trust), min(1.0, h[i] + trust)) for i in range(N)] + [(None, None)]
    c_obj = np.zeros(N + 1); c_obj[-1] = 1.0
    res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method="highs")
    if not res.success:
        return h, None
    h_new = res.x[:N]
    h_new = project_to_constraints(h_new)
    return h_new, res.x[-1]


def basin_hop(N, n_starts=8, n_kicks=4, schedule=(50., 200., 1000., 5000., 20000., 100000.),
              symmetric=True, slp_iters=10, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    starts = []
    xs = np.linspace(1./N, 2-1./N, N)
    starts.append(np.minimum(xs, 2-xs)*(N/2)/np.minimum(xs, 2-xs).sum())
    starts.append(np.where((xs > 0.3) & (xs < 1.7), 1., 0.) * (N/2) / max(1, ((xs>0.3)&(xs<1.7)).sum()))
    starts.append(0.5*(1 - np.cos(np.pi*xs)))
    for _ in range(n_starts - 3):
        starts.append(project_to_constraints(0.5 + 0.4*rng.standard_normal(N)))

    best = (1e9, None, None)
    t0 = time.time()
    for k, h0 in enumerate(starts):
        h = lbfgs_anneal(h0, schedule, λ=400., iters=600, symmetric=symmetric)
        b, j, s = bound_float(h)
        if verbose:
            print(f"  init {k}: bound = {b:.7f}")
        # SLP polishing
        for it in range(slp_iters):
            h_try, c_try = slp_step(h, trust=0.03)
            if c_try is None: break
            b2, _, _ = bound_float(h_try)
            if b2 < b - 1e-7:
                h = h_try
                b = b2
            else:
                break
        # Random kicks
        for kick in range(n_kicks):
            h_kick = h + 0.05 * rng.standard_normal(N)
            if symmetric: h_kick = 0.5 * (h_kick + h_kick[::-1])
            h_kick = project_to_constraints(np.clip(h_kick, 0, 1))
            h2 = lbfgs_anneal(h_kick, schedule, λ=400., iters=400, symmetric=symmetric)
            b2, _, _ = bound_float(h2)
            if b2 < b - 1e-7:
                h = h2; b = b2
        b_final, j_final, s_final = bound_float(h)
        if b_final < best[0]:
            best = (b_final, h.copy(), (j_final, s_final))
        if verbose:
            print(f"  init {k} final: bound = {b_final:.7f}  argmax j={j_final}({s_final})")
    return best


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    out = basin_hop(N, n_starts=6, n_kicks=2, slp_iters=15, seed=7)
    print(f"\nBEST (N={N}, symmetric): {out[0]:.7f}  argmax j={out[2][0]} ({out[2][1]})")
    json.dump({"N": N, "bound": out[0], "argmax_j": int(out[2][0]), "argmax_sign": out[2][1],
               "h": out[1].tolist()}, open(f"basinhop_N{N}.json","w"))
