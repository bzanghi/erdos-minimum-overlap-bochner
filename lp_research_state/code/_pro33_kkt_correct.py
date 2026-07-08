"""Pretest D: extract dual multipliers gamma from the correct minimax LP
at n=600 and test the KKT functional equation with THEM (not with a
threshold-classified LP as PRO-23 did).

KKT stationarity for min_h max_j M_j(h) s.t. 0<=h<=1, sum h = n/2:
exists gamma_j >= 0 (sum=1 on active set), lambda in R, mu_lo, mu_hi >= 0:
   sum_j gamma_j * dM_j/dh_k = lambda + mu_hi_k - mu_lo_k
with complementary slackness. On strictly interior cells:
   sum_j gamma_j * w*(g[k+j] - h[k-j]) = lambda.
Note dM_j/dh_k = w*(g[k+j] - h[k-j]) with in-range indicators, and
g = 1-h, so this is a slightly different equation than PRO-23's
"sum gamma [h(x+t)+h(x-t)] = kappa" — PRO-23 used the symmetrized form.
Check both.
"""
import json, numpy as np
from scipy.optimize import linprog

d = json.load(open('lp_research_state/data/together_f_star.json'))
h = np.array(d['together']['values'], dtype=float)
n = len(h); w = 2.0 / n
g = 1.0 - h

def corr_all(h):
    n = len(h); g = 1.0 - h
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

r = 1e-4
m = 2 * n - 1
A = np.zeros((m, n + 1)); b = np.empty(m)
for idx, j in enumerate(lags):
    row = np.zeros(n)
    if j >= 0:
        row[:n - j] += g[j:]
        row[j:] -= h[:n - j]
    else:
        row[-j:] += g[:n + j]
        row[:n + j] -= h[-j:]
    A[idx, :n] = w * row
    A[idx, n] = -1.0
    b[idx] = -w * C[idx]
lo = np.maximum(-r, -h); hi = np.minimum(r, 1.0 - h)
bounds = [(lo[k], hi[k]) for k in range(n)] + [(None, None)]
A_eq = np.zeros((1, n + 1)); A_eq[0, :n] = 1.0
c = np.zeros(n + 1); c[n] = 1.0
res = linprog(c, A_ub=A, b_ub=b, A_eq=A_eq, b_eq=[0.0],
              bounds=bounds, method='highs')
print("LP status:", res.status, res.message[:60])
print("pred max:", res.x[n], " gain:", M.max() - res.x[n])

gam = -np.asarray(res.ineqlin.marginals)  # duals of <= constraints, >=0
lam = -np.asarray(res.eqlin.marginals)[0]
print(f"\ngamma: sum={gam.sum():.6f}, nnz(>1e-8)={int((gam>1e-8).sum())}")
print("top-10 gamma lags:", lags[np.argsort(-gam)[:10]].tolist())
print("gamma symmetric? ||gam(j)-gam(-j)||_inf =",
      np.abs(gam - gam[::-1]).max())

# stationarity residual on interior cells using LP duals directly:
# grad_k = sum_j gam_j * w * (g[k+j] - h[k-j])   (in-range)
gradsum = np.zeros(n)
for idx, j in enumerate(lags):
    if gam[idx] < 1e-14:
        continue
    row = np.zeros(n)
    if j >= 0:
        row[:n - j] += g[j:]
        row[j:] -= h[:n - j]
    else:
        row[-j:] += g[:n + j]
        row[:n + j] -= h[-j:]
    gradsum += gam[idx] * w * row

interior = (h > 1e-6) & (h < 1 - 1e-6)
lower = h <= 1e-6
upper = h >= 1 - 1e-6
print(f"\ncells: interior={interior.sum()}, lower={lower.sum()}, upper={upper.sum()}")
res_int = np.abs(gradsum[interior] - lam)
print(f"lambda (eq dual) = {lam:.6e}")
print(f"interior stationarity residual: max={res_int.max():.3e}, "
      f"mean={res_int.mean():.3e}")
# sign conditions: lower-active cells need gradsum >= lam (pushing up costs),
# upper-active need gradsum <= lam
viol_lo = np.maximum(0, lam - gradsum[lower])
viol_hi = np.maximum(0, gradsum[upper] - lam)
print(f"lower-cell violation: max={viol_lo.max() if lower.any() else 0:.3e}")
print(f"upper-cell violation: max={viol_hi.max() if upper.any() else 0:.3e}")

# Also test PRO-23's symmetrized form with these gammas:
# sum gamma_j [h(x+j) + h(x-j)] ?= const on interior
conv = np.zeros(n)
hz = np.concatenate([np.zeros(n - 1), h, np.zeros(n - 1)])
for idx, j in enumerate(lags):
    if gam[idx] < 1e-14:
        continue
    conv += gam[idx] * (hz[n - 1 + j: 2 * n - 1 + j])
kappa = conv[interior].mean()
print(f"\nPRO-23-form: sum gam h(x+t): interior spread "
      f"max-dev={np.abs(conv[interior]-kappa).max():.3e} around {kappa:.4f}")
