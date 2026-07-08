"""Pretest A+B: active-set contiguity of Together's h*, and a rigorous
SLP step with exact quadratic-error accounting.

Convention (together_loader.py): h on [0,2], n=600 cells, width 2/n,
sum(h) = n/2 = 300, g = 1-h on-support, zero-extended off-support.
M_j = (2/n) * C_j,  C_j = sum_i h_i g_{i+j},  j = -(n-1)..(n-1).
M is EXACTLY quadratic in h:
  M_j(h+d) = M_j(h) + L_j(d) - (2/n) sum_i d_i d_{i+j}
so |true - linear| <= (2/n)*||d||_2^2 for every j.
"""
import json, numpy as np
from scipy.optimize import linprog

d = json.load(open('lp_research_state/data/together_f_star.json'))
h = np.array(d['together']['values'], dtype=float)
n = len(h)
w = 2.0 / n
print(f"n={n}, sum={h.sum():.12f}, min={h.min():.3e}, max={h.max():.12f}")

def corr_all(h):
    """C_j for j=-(n-1)..(n-1); returns (lags, C)."""
    g = 1.0 - h
    lags = np.arange(-(n - 1), n)
    C = np.empty(2 * n - 1)
    for idx, j in enumerate(lags):
        if j >= 0:
            C[idx] = h[:n - j] @ g[j:]
        else:
            C[idx] = h[-j:] @ g[:n + j]
    return lags, C

lags, C = corr_all(h)
M = w * C
Mmax = M.max()
print(f"\nM_max = {Mmax:.15f}  at lag(s) {lags[M > Mmax - 1e-13].tolist()}")

# ---- Pretest A: active-set contiguity ----
print("\n=== Pretest A: active-set structure ===")
for tol in [1e-12, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6]:
    act = lags[M > Mmax - tol]
    # contiguous runs
    runs = []
    if len(act):
        s = p = act[0]
        for a in act[1:]:
            if a == p + 1:
                p = a
            else:
                runs.append((s, p)); s = p = a
        runs.append((s, p))
    runsummary = runs if len(runs) <= 8 else runs[:4] + ['...'] + runs[-3:]
    print(f"tol={tol:.0e}: |S|={len(act):4d}  runs={len(runs):3d}  {runsummary}")

# ---- Pretest B: one SLP step with exact quadratic accounting ----
print("\n=== Pretest B: SLP step, exact quadratic bound ===")
g = 1.0 - h
# gradient rows: dM_j/dh_k = w * ( g[k+j] (in range) - h[k-j] (in range) )
# Build A_ub (2n-1) x (n+1): w*Lrow_j . delta - u <= -M_j  with vars (delta, u)
def build_lp(h, r):
    g = 1.0 - h
    m = 2 * n - 1
    A = np.zeros((m, n + 1))
    b = np.empty(m)
    lags_, C_ = corr_all(h)
    for idx, j in enumerate(lags_):
        row = np.zeros(n)
        # g[k+j] term
        if j >= 0:
            row[:n - j] += g[j:]
        else:
            row[-j:] += g[:n + j]
        # -h[k-j] term
        if j >= 0:
            row[j:] -= h[:n - j]
        else:
            row[:n + j] -= h[-j:]
        A[idx, :n] = w * row
        A[idx, n] = -1.0
        b[idx] = -w * C_[idx]
    lo = np.maximum(-r, -h)
    hi = np.minimum(r, 1.0 - h)
    bounds = [(lo[k], hi[k]) for k in range(n)] + [(None, None)]
    A_eq = np.zeros((1, n + 1)); A_eq[0, :n] = 1.0
    c = np.zeros(n + 1); c[n] = 1.0
    return linprog(c, A_ub=A, b_ub=b, A_eq=A_eq, b_eq=[0.0],
                   bounds=bounds, method='highs')

for r in [1e-3, 1e-4, 1e-5]:
    res = build_lp(h, r)
    if not res.success:
        print(f"r={r:.0e}: LP FAILED: {res.message}"); continue
    delta = res.x[:n]; u_pred = res.x[n]
    quad_bound = w * (delta @ delta)
    h2 = h + delta
    _, C2 = corr_all(h2)
    Mnew = w * C2.max()
    print(f"r={r:.0e}: pred_max={u_pred:.12f} (gain {Mmax-u_pred:+.3e})  "
          f"quad_bnd={quad_bound:.3e}  TRUE new max={Mnew:.12f} "
          f"(true gain {Mmax-Mnew:+.3e})  feas: sum={h2.sum():.9f} "
          f"min={h2.min():.2e} max={h2.max():.12f}")
