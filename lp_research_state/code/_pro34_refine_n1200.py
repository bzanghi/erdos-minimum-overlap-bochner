"""Pretest C: refine Together's h* from n=600 to n=1200 by cell-doubling
(value-preserving, M-preserving on the doubled lag grid), then run the
correct SLP (all lags, box + mass constraints, exact quadratic accounting)
iteratively to see whether the new degrees of freedom admit descent.
"""
import json, numpy as np
from scipy.optimize import linprog

d = json.load(open('lp_research_state/data/together_f_star.json'))
h600 = np.array(d['together']['values'], dtype=float)

def corr_all(h):
    n = len(h)
    g = 1.0 - h
    lags = np.arange(-(n - 1), n)
    C = np.empty(2 * n - 1)
    for idx, j in enumerate(lags):
        if j >= 0:
            C[idx] = h[:n - j] @ g[j:]
        else:
            C[idx] = h[-j:] @ g[:n + j]
    return lags, C

def Mmax_of(h):
    n = len(h)
    _, C = corr_all(h)
    return (2.0 / n) * C.max()

M600 = Mmax_of(h600)
h = np.repeat(h600, 2)
n = len(h)
w = 2.0 / n
M0 = Mmax_of(h)
print(f"n=600 M={M600:.15f}   n=1200 doubled M={M0:.15f}  (diff {M0-M600:+.2e})")

def slp_step(h, r):
    n = len(h); w = 2.0 / n
    g = 1.0 - h
    lags_, C_ = corr_all(h)
    m = 2 * n - 1
    A = np.zeros((m, n + 1)); b = np.empty(m)
    for idx, j in enumerate(lags_):
        row = np.zeros(n)
        if j >= 0:
            row[:n - j] += g[j:]
            row[j:] -= h[:n - j]
        else:
            row[-j:] += g[:n + j]
            row[:n + j] -= h[-j:]
        A[idx, :n] = w * row
        A[idx, n] = -1.0
        b[idx] = -w * C_[idx]
    lo = np.maximum(-r, -h); hi = np.minimum(r, 1.0 - h)
    bounds = [(lo[k], hi[k]) for k in range(n)] + [(None, None)]
    A_eq = np.zeros((1, n + 1)); A_eq[0, :n] = 1.0
    c = np.zeros(n + 1); c[n] = 1.0
    res = linprog(c, A_ub=A, b_ub=b, A_eq=A_eq, b_eq=[0.0],
                  bounds=bounds, method='highs')
    if not res.success:
        return None, None
    return res.x[:n], res.x[n]

# one diagnostic step at a few radii
Mcur = M0
for r in [1e-3, 1e-4]:
    delta, u = slp_step(h, r)
    if delta is None:
        print(f"r={r:.0e}: LP failed"); continue
    quad = w * (delta @ delta)
    Mnew = Mmax_of(np.clip(h + delta, 0.0, 1.0))
    print(f"r={r:.0e}: pred gain {M0-u:+.3e}  quad_bnd {quad:.2e}  "
          f"true gain {M0-Mnew:+.3e}")

# iterate with adaptive trust region
print("\n--- iterative SLP at n=1200 ---")
hbest = h.copy(); Mbest = M0
hcur = h.copy(); r = 5e-4
for it in range(40):
    delta, u = slp_step(hcur, r)
    if delta is None:
        print(f"it{it}: LP failed at r={r:.1e}"); break
    trial = np.clip(hcur + delta, 0.0, 1.0)
    # restore exact mass by tiny uniform shift on interior cells
    excess = trial.sum() - n / 2.0
    interior = (trial > 1e-9) & (trial < 1 - 1e-9)
    if interior.any():
        trial[interior] -= excess / interior.sum()
        trial = np.clip(trial, 0.0, 1.0)
    Mt = Mmax_of(trial)
    pred_gain = Mmax_of(hcur) - u
    act_gain = Mmax_of(hcur) - Mt
    if Mt < Mmax_of(hcur) - 1e-14:
        hcur = trial; note = "ACCEPT"
        if Mt < Mbest:
            Mbest = Mt; hbest = trial.copy()
        r = min(r * 1.5, 2e-3)
    else:
        note = "reject"; r *= 0.4
    if it % 5 == 0 or note == "ACCEPT":
        print(f"it{it:3d} r={r:.1e} pred_gain={pred_gain:+.2e} "
              f"act_gain={act_gain:+.2e} M={Mt:.15f} {note}")
    if r < 1e-8:
        print("trust region collapsed"); break

print(f"\nFINAL: M_best(n=1200) = {Mbest:.15f}")
print(f"       vs Together UB    0.380870310586220  (improvement {M600-Mbest:+.3e})")
np.save('/private/tmp/claude-501/-Users-benzanghi-Documents-Claude-Projects-Erdos--claude-worktrees-creative-problem-solving-5279bd/9d445af9-366e-4d48-b6bd-b702477204eb/scratchpad/h1200_best.npy', hbest)
