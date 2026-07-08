import json, numpy as np, os, time
np.random.seed(3)
os.chdir('/Users/benzanghi/Documents/Claude/Projects/Erdos/.claude/worktrees/creative-problem-solving-5279bd')
d = json.load(open('lp_research_state/data/together_f_star.json'))
h = np.repeat(np.array(d['together']['values'], dtype=float), 2)
n = len(h); w = 2.0/n; g = 1.0-h
lags = np.arange(-(n-1), n)
C = np.empty(2*n-1)
for idx, j in enumerate(lags):
    C[idx] = h[:n-j] @ g[j:] if j >= 0 else h[-j:] @ g[:n+j]
M = w*C; Mmax = M.max()
act = np.where(M > Mmax - 1e-9)[0]
free = np.where((h > 1e-9) & (h < 1 - 1e-9))[0]
G = np.zeros((len(act) + 1, len(free)))
for ci, idx in enumerate(act):
    j = int(lags[idx])
    row = np.zeros(n)
    if j >= 0:
        row[:n-j] += g[j:]; row[j:] -= h[:n-j]
    else:
        row[-j:] += g[:n+j]; row[:n+j] -= h[-j:]
    G[ci] = row[free]
G[-1] = 1.0
u_, s_, vt = np.linalg.svd(G, full_matrices=True)
rank = (s_ > s_[0]*1e-10).sum()
Nsp = vt[rank:].T
P = np.zeros((n, Nsp.shape[1])); P[free] = Nsp
K = P.shape[1]; m = len(act)
Bs = np.empty((m, K, K))
for ci, idx in enumerate(act):
    j = abs(int(lags[idx]))
    SP = np.zeros_like(P)
    if j == 0:
        SP = -w * P
    else:
        SP[:n-j] += P[j:]; SP[j:] += P[:n-j]; SP *= -w/2
    Bs[ci] = P.T @ SP
scale = np.abs(Bs).max()
Bn = Bs / scale
Bflat = Bn.reshape(m, K*K)
print(f"K={K} m={m}", flush=True)

def q_of(x):
    return Bflat @ np.outer(x, x).ravel()

best = (np.inf, None)
t0 = time.time()
for restart in range(30):
    x = np.random.randn(K); x /= np.linalg.norm(x)
    for beta in [50, 200, 1000, 5000]:
        for it in range(400):
            q = q_of(x)
            mq = q.max()
            wts = np.exp(beta*(q-mq)); wts /= wts.sum()
            Aw = (wts[:,None]*Bflat).sum(0).reshape(K,K)
            grad = 2*(Aw @ x)
            grad -= (grad @ x)*x
            gn = np.linalg.norm(grad)
            if gn < 1e-14: break
            x = x - 0.05*grad/gn
            x /= np.linalg.norm(x)
    q = q_of(x)
    if q.max() < best[0]:
        best = (q.max(), x.copy())
        print(f"restart {restart}: new best max_j Q_j = {q.max():.6e} (scaled)  t={time.time()-t0:.0f}s", flush=True)
print(f"\nbest max_j Q_j over sphere (scaled): {best[0]:.6e}  (orig {best[0]*scale:.3e})", flush=True)
x = best[1]
np.save('/private/tmp/claude-501/-Users-benzanghi-Documents-Claude-Projects-Erdos--claude-worktrees-creative-problem-solving-5279bd/9d445af9-366e-4d48-b6bd-b702477204eb/scratchpad/descent_x.npy', x)
if best[0] < 0:
    delta = P @ x
    print(f"DESCENT CANDIDATE: ||delta||_inf={np.abs(delta).max():.3f}", flush=True)
    def corr_max(hh):
        gg = 1.0-hh; b = -1.0
        for j in range(-(n-1), n):
            v = hh[:n-j] @ gg[j:] if j >= 0 else hh[-j:] @ gg[:n+j]
            if v > b: b = v
        return w*b
    for t in [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]:
        for sgn in (1,-1):
            h2 = h + sgn*t*delta
            if h2.min() < -1e-15 or h2.max() > 1+1e-15:
                print(f"t={sgn*t:+.0e}: infeasible", flush=True); continue
            Mn = corr_max(np.clip(h2,0,1))
            print(f"t={sgn*t:+.0e}: M={Mn:.15f} gain={Mmax-Mn:+.3e}", flush=True)
else:
    print("no negative-curvature common direction found by smoothed minimax", flush=True)
