"""PRO-4: upsample Together's 600-cell h* to finer discretization (1200, 2400)
and run LP-minimax descent. The hypothesis: with more degrees of freedom,
the discrete optimum may go below 0.3808703.

Upsampling preserves M at the *even* shifts (which include the active set
of the 600-cell problem). The new ODD shifts may have larger or smaller M,
so the initial max may be the same, larger, or smaller than 0.3808703.
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def compute_M(h: np.ndarray, L: float) -> np.ndarray:
    n = len(h)
    corr = np.correlate(h, 1.0 - h, mode="full")
    return L * corr[n - 1 : 2 * n - 1]


def upsample(h: np.ndarray, factor: int = 2) -> np.ndarray:
    """Each cell h[i] → factor copies."""
    return np.repeat(h, factor)


def gradient_M_at(j: int, h: np.ndarray, L: float) -> np.ndarray:
    n = len(h)
    if j == 0:
        return L * (1.0 - 2.0 * h)
    g = np.zeros(n)
    g[: n - j] += 1.0 - h[j:]
    g[j:] -= h[: n - j]
    return L * g


def lp_descent(h, L, active_tol=1e-7, step_inf=1.0, box_eps=1e-9):
    n = len(h)
    Ms = compute_M(h, L)
    M_max = Ms.max()
    S = np.where(M_max - Ms <= active_tol)[0]
    if len(S) == 0:
        S = np.array([int(np.argmax(Ms))])

    G = np.zeros((len(S), n))
    for k, j in enumerate(S):
        G[k] = gradient_M_at(int(j), h, L)

    n_var = n + 1
    c = np.zeros(n_var); c[-1] = 1.0
    A_ub = np.zeros((len(S), n_var))
    A_ub[:, :n] = G; A_ub[:, n] = -1.0
    b_ub = np.zeros(len(S))
    A_eq = np.zeros((1, n_var)); A_eq[0, :n] = 1.0
    b_eq = np.array([0.0])
    bounds = []
    for i in range(n):
        lo, up = -step_inf, +step_inf
        if h[i] <= box_eps:
            lo = 0.0
        if h[i] >= 1.0 - box_eps:
            up = 0.0
        bounds.append((lo, up))
    bounds.append((None, None))
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=bounds, method="highs", options={"presolve": True})
    if not res.success:
        return None, None, len(S), M_max
    return res.x[:n], res.x[n], len(S), M_max


def line_search(h, d, L, alpha_max=1.0):
    """Find largest improving α."""
    M_init = compute_M(h, L).max()
    target_sum = h.sum()
    # Box-feasibility cap
    eps = 1e-14
    a_box = alpha_max
    for i in range(len(h)):
        if d[i] > eps:
            a_box = min(a_box, (1.0 - h[i]) / d[i])
        elif d[i] < -eps:
            a_box = min(a_box, -h[i] / d[i])
    a_box = max(0.0, min(a_box, alpha_max))
    alpha = a_box
    for _ in range(60):
        if alpha < 1e-16:
            break
        h_new = np.clip(h + alpha * d, 0.0, 1.0)
        s = h_new.sum()
        if abs(s - target_sum) > 1e-9:
            h_new = h_new * (target_sum / s)
            h_new = np.clip(h_new, 0.0, 1.0)
        M_new = compute_M(h_new, L).max()
        if M_new < M_init - 1e-15:
            return alpha, M_new, h_new
        alpha *= 0.5
    return 0.0, M_init, h.copy()


def refine(h0, L, *, max_iters=300, active_tol_init=1e-5, verbose=True):
    h = h0.copy()
    target_sum = float(h0.sum())
    M_init = compute_M(h, L).max()
    history = [M_init]
    active_tol = active_tol_init
    no_improve = 0
    t0 = time.time()
    for it in range(max_iters):
        d, t_star, S_size, M_max = lp_descent(h, L, active_tol=active_tol)
        if d is None or t_star is None or t_star >= -1e-14:
            if active_tol > 1e-13:
                active_tol = max(active_tol * 0.3, 1e-13)
                continue
            if verbose:
                print(f"  it={it:3d} t*={t_star} converged")
            break
        alpha, M_new, h_new = line_search(h, d, L, alpha_max=1.0)
        if alpha == 0:
            if active_tol < active_tol_init * 10:
                # Widen active set
                active_tol = min(active_tol * 5, active_tol_init * 10)
                no_improve += 1
                if no_improve > 10:
                    if verbose:
                        print(f"  it={it:3d} no improve after widening")
                    break
                continue
            if active_tol > 1e-13:
                active_tol = max(active_tol * 0.3, 1e-13)
                continue
            break
        h = h_new
        no_improve = 0
        history.append(M_new)
        if verbose and (it % 10 == 0 or it < 5):
            print(
                f"  it={it:3d} |S|={S_size:4d} M={M_new:.10f} "
                f"Δ={history[-1]-history[-2]:+.3e} α={alpha:.2e} "
                f"t*={t_star:+.2e} tol={active_tol:.0e} t={time.time()-t0:.1f}s"
            )
    return h, history


def main(factor: int = 2, max_iters: int = 100):
    with open(DATA) as f:
        data = json.load(f)
    h_orig = np.array(data["together"]["values"], dtype=np.float64)
    n_orig = len(h_orig)
    L_orig = 2.0 / n_orig
    M_orig = compute_M(h_orig, L_orig).max()
    print(f"Original: n={n_orig}, M = {M_orig:.10f}")

    h_up = upsample(h_orig, factor=factor) / float(factor)
    # Wait: we need to PRESERVE sum_h = (n_new)/2. Let me think.
    # n_orig=600, sum=300, L=1/300, ∫h = L*sum = 1.
    # If we upsample by factor=2: each cell becomes 2 cells of width L/2.
    # We want integral preserved: ∫h_new = L_new * sum(h_new) = 1.
    # If we keep h_new[2i] = h_new[2i+1] = h_orig[i], then sum_new = 2*sum_orig = 600.
    # L_new = L_orig/2 = 1/600.
    # ∫h_new = (1/600) * 600 = 1.  GOOD.
    # And constraint Σh_new = n_new/2 = 1200/2 = 600. GOOD.
    # So we should NOT divide by factor; just repeat.
    h_up = upsample(h_orig, factor=factor)
    n_new = len(h_up)
    L_new = 2.0 / n_new
    sum_new = h_up.sum()
    expected_sum = n_new / 2
    print(f"Upsampled: n={n_new}, sum={sum_new}, expected={expected_sum}, L={L_new}")
    M_up = compute_M(h_up, L_new).max()
    print(f"Upsampled M (before refine) = {M_up:.10f}  Δ vs orig = {M_up - M_orig:+.3e}")

    print(f"\nRefining at n={n_new}...")
    h_refined, history = refine(h_up, L_new, max_iters=max_iters,
                                 active_tol_init=1e-5, verbose=True)
    M_final = compute_M(h_refined, L_new).max()
    print(f"\nUpsampled+refined M = {M_final:.10f}")
    print(f"Δ vs Together (0.380871): {M_final - M_orig:+.3e}")
    if M_final < M_orig - 1e-9:
        print(f"✓ IMPROVEMENT! Δ = {M_orig - M_final:.3e}")
    else:
        print(f"✗ No improvement (Δ = {M_final - M_orig:+.3e})")

    return M_orig, M_final


if __name__ == "__main__":
    import sys
    factor = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    main(factor=factor, max_iters=iters)
