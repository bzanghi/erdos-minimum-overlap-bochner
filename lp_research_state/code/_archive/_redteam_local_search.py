"""Local-search refinement of Together's h via successive LP linearization.
Repeats the SLP step in Together's stated method but starting from their h."""
import sys
sys.path.insert(0, '/tmp/together_repo/erdos-minimum-overlap/solutions')
import numpy as np
import cvxpy as cp
from together_ai_2026 import h_values

h = np.array(h_values, dtype=np.float64)
n = len(h)
print(f"n = {n}")

def M_true(h):
    return np.max(np.correlate(h, 1-h, mode='full'))/len(h)*2

def A_k_and_grad(h, k, n):
    """A_k(h) = sum_t h[t]*h[t+k]; grad_t = h[t+k] + h[t-k] (zeros outside range)."""
    if k >= 0:
        A = float(np.sum(h[:n-k]*h[k:]))
        grad = np.zeros(n)
        grad[:n-k] += h[k:]
        grad[k:] += h[:n-k]
    else:
        kk = -k
        A = float(np.sum(h[:n-kk]*h[kk:]))
        grad = np.zeros(n)
        grad[:n-kk] += h[kk:]
        grad[kk:] += h[:n-kk]
    return A, grad

def slp_step(h, n, delta=0.005, tol_active=1e-4):
    """One step of sequential LP: linearize A_k around h for shifts active within tol."""
    conv = np.correlate(h, 1-h, mode='full')
    mx = np.max(conv)
    shifts = np.where(conv > mx - tol_active*n/2)[0] - (n-1)
    # filter k!=0 (k=0 is always largest by trivial reasoning? no actually)
    shifts = [int(k) for k in shifts if k != 0]
    print(f"  {len(shifts)} active shifts (tol={tol_active})")
    x = cp.Variable(n)
    t = cp.Variable()
    cons = [cp.sum(x) == n/2.0, x >= 0, x <= 1, cp.norm(x - h, 'inf') <= delta]
    for k in shifts:
        A_h, grad = A_k_and_grad(h, k, n)
        # we want min over k of A_k(x) >= t  --->  A_h + grad·(x-h) >= t
        cons.append(t <= A_h + grad @ (x - h))
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver='CLARABEL')
    return x.value, t.value, prob.status

# Run SLP iterations
M0 = M_true(h)
print(f"M(h_0) = {M0:.15f}")

current = h.copy()
for it in range(8):
    print(f"--- Iter {it} ---")
    new, t_val, status = slp_step(current, n, delta=0.002, tol_active=2e-4)
    if new is None:
        print("LP failed:", status)
        break
    M_new = M_true(new)
    print(f"  status={status}  predicted t={t_val:.10f}  true M={M_new:.15f}  delta={M_new - M_true(current):.3e}")
    if M_new < M_true(current) - 1e-12:
        current = new
    else:
        print("  no improvement; stopping")
        break

print(f"\nFinal M(h_*) = {M_true(current):.15f}")
print(f"Improvement vs baseline: {M0 - M_true(current):.3e}")
print(f"||h_* - h_0||_inf = {np.max(np.abs(current - h)):.3e}")
print(f"||h_* - h_0||_2   = {np.linalg.norm(current - h):.3e}")
