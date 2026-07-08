"""PRO-4: refine Together's h* via subgradient/LP minimax descent.

The KKT residual in PRO-23 (ε ≈ 7.6e-3) showed Together's h* is not a tight
optimum. Attempt to push the UB μ ≤ M(h*) below 0.380871 via local descent.

Method
------
At each iteration:
  1. Compute M(jL) for all j ∈ [0, n-1].
  2. Active set S = {j : M(jL) ≥ max_M - tol}.
  3. Gradient g[j, :] = ∂M(jL)/∂h for j ∈ S.
     ∂M(jL)/∂h[a] = L*[1_{a+j<n}*(1 - h[a+j]) - 1_{a-j≥0}*h[a-j]]
  4. Solve LP: min t s.t.  g[j, :] · d ≤ t for j ∈ S
                          Σ d = 0
                          d_i ≥ 0 if h_i ≤ ε_lo (lower-active)
                          d_i ≤ 0 if h_i ≥ 1-ε_up (upper-active)
                          ||d||_∞ ≤ 1
  5. Backtracking line search on α; h ← h + α*d.
  6. Repeat until t* ≥ -eps_descent.

We report new max_j M(jL) and compare to 0.380871.
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"
OUT = ROOT / "lp_research_state" / "data" / "pro4_refined_h.json"


def compute_M(h: np.ndarray, L: float) -> np.ndarray:
    """M(jL) = L * Σ_{i=0..n-1-j} h[i]*(1 - h[i+j]) for j = 0..n-1."""
    n = len(h)
    # Use np.correlate which matches Together's convention:
    # corr[k] for k = -(n-1)..(n-1); positive lags are indices n-1..2n-2
    corr = np.correlate(h, 1.0 - h, mode="full")
    # Lag j (>= 0): corr[n-1+j]
    return L * corr[n - 1 : 2 * n - 1]


def gradient_M_at(j: int, h: np.ndarray, L: float) -> np.ndarray:
    """∂M(jL)/∂h[a] for a = 0..n-1.

    M(jL) = L*Σ_{i=0..n-1-j} h[i]*(1 - h[i+j])
          = L*Σ h[i] - L*Σ h[i]*h[i+j]
    So:
        ∂/∂h[a] = L * [1_{a≤n-1-j} - 1_{a≤n-1-j}*h[a+j] - 1_{a≥j}*h[a-j]]
    """
    n = len(h)
    g = np.zeros(n)
    # +1 from the linear part: a ∈ [0, n-1-j]
    if j == 0:
        # Special: M(0) = L*Σ h[i](1-h[i]); ∂/∂h[a] = L*(1 - 2*h[a])
        return L * (1.0 - 2.0 * h)
    g[: n - j] += 1.0
    # -h[a+j] from a in [0, n-1-j]
    g[: n - j] -= h[j:]
    # -h[a-j] from a in [j, n-1]
    g[j:] -= h[: n - j]
    return L * g


def descent_step(
    h: np.ndarray, L: float,
    active_tol: float = 1e-7,
    box_eps: float = 1e-8,
    step_inf_bound: float = 1.0,
):
    """Compute one descent direction via LP. Returns (d, t_star, info)."""
    n = len(h)
    Ms = compute_M(h, L)
    M_max = Ms.max()
    S = np.where(M_max - Ms <= active_tol)[0]
    if len(S) == 0:
        S = np.array([int(np.argmax(Ms))])

    # Build gradient matrix G of shape (|S|, n).
    G = np.zeros((len(S), n))
    for idx, j in enumerate(S):
        G[idx] = gradient_M_at(int(j), h, L)

    # LP variables: x = (d_0, ..., d_{n-1}, t).
    # min t  s.t.  G·d - t ≤ 0,  Σd = 0,  d ∈ box,  ||d||_∞ ≤ step_inf_bound, t free.
    n_var = n + 1
    c = np.zeros(n_var); c[-1] = 1.0

    # Inequality: G·d - t ≤ 0
    A_ub = np.zeros((len(S), n_var))
    A_ub[:, :n] = G
    A_ub[:, n] = -1.0
    b_ub = np.zeros(len(S))

    # Equality: Σ d = 0  (and we leave t free)
    A_eq = np.zeros((1, n_var)); A_eq[0, :n] = 1.0
    b_eq = np.array([0.0])

    # Bounds on d: box-respecting
    bounds = []
    for i in range(n):
        lo = -step_inf_bound
        up = +step_inf_bound
        if h[i] <= box_eps:
            lo = 0.0  # can't go below 0
        if h[i] >= 1.0 - box_eps:
            up = 0.0  # can't go above 1
        bounds.append((lo, up))
    bounds.append((None, None))  # t free

    res = linprog(
        c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
        method="highs", options={"presolve": True},
    )
    if not res.success:
        return None, None, {"status": res.message, "S_size": len(S), "M_max": M_max}

    d = res.x[:n]
    t_star = res.x[n]
    return d, t_star, {
        "S_size": len(S),
        "M_max": M_max,
        "min_d": float(d.min()),
        "max_d": float(d.max()),
        "norm_d": float(np.linalg.norm(d)),
    }


def line_search(h: np.ndarray, d: np.ndarray, L: float, *, alpha_max: float = 0.1):
    """Backtracking line search: find largest α ∈ (0, alpha_max] such that
    h + α*d ∈ [0,1]^n and max_j M(h + α*d) < max_j M(h)."""
    M_init = compute_M(h, L).max()

    # Cap α by box constraints: h + αd ∈ [0,1]
    eps = 1e-12
    a_box = alpha_max
    for i in range(len(h)):
        if d[i] > eps:
            a_box = min(a_box, (1.0 - h[i]) / d[i])
        elif d[i] < -eps:
            a_box = min(a_box, -h[i] / d[i])
    a_box = max(0.0, min(a_box, alpha_max))

    alpha = a_box
    best_alpha = 0.0
    best_M = M_init
    for _ in range(40):  # 40 halvings
        if alpha < 1e-14:
            break
        h_new = h + alpha * d
        # Clip to handle numerical drift
        h_new = np.clip(h_new, 0.0, 1.0)
        # Project sum back (rare; only if clip changed anything)
        s = h_new.sum()
        if abs(s - h.sum()) > 1e-9:
            h_new = h_new * (h.sum() / s)
            h_new = np.clip(h_new, 0.0, 1.0)
        M_new = compute_M(h_new, L).max()
        if M_new < best_M - 1e-15:
            best_M = M_new
            best_alpha = alpha
            return best_alpha, best_M  # accept first improving step
        alpha *= 0.5
    return best_alpha, best_M


def refine(
    h0: np.ndarray, L: float,
    *,
    max_iters: int = 200,
    active_tol_init: float = 1e-6,
    descent_tol: float = 1e-12,
    verbose: bool = True,
):
    h = h0.copy()
    M_history = [compute_M(h, L).max()]
    active_tol = active_tol_init
    log = []
    t0 = time.time()
    for it in range(max_iters):
        d, t_star, info = descent_step(h, L, active_tol=active_tol)
        if d is None or t_star is None or t_star >= -descent_tol:
            if active_tol > 1e-12:
                active_tol = max(active_tol * 0.1, 1e-12)
                if verbose:
                    print(f"  it={it:3d} no descent; tightening active_tol → {active_tol:.0e}")
                continue
            if verbose:
                print(f"  it={it:3d} converged: t*={t_star}, info={info}")
            break
        # Line search
        alpha, M_new = line_search(h, d, L, alpha_max=1.0)
        if alpha == 0:
            if active_tol > 1e-12:
                active_tol = max(active_tol * 0.5, 1e-12)
                if verbose:
                    print(f"  it={it:3d} line-search fail; tightening active_tol → {active_tol:.0e}")
                continue
            if verbose:
                print(f"  it={it:3d} line-search fail at minimum tol — stopping.")
            break
        h_new = np.clip(h + alpha * d, 0.0, 1.0)
        s = h_new.sum()
        if abs(s - h0.sum()) > 1e-9:
            h_new = h_new * (h0.sum() / s)
            h_new = np.clip(h_new, 0.0, 1.0)
        h = h_new
        M_history.append(M_new)
        if verbose:
            print(
                f"  it={it:3d} |S|={info['S_size']:4d} M={M_new:.10f} "
                f"Δ={M_history[-1] - M_history[-2]:+.3e} α={alpha:.3e} "
                f"t*={t_star:+.3e} sum={h.sum():.6f} t={time.time()-t0:.1f}s"
            )
    return h, M_history, log


def main():
    with open(DATA) as f:
        data = json.load(f)
    h0 = np.array(data["together"]["values"], dtype=np.float64)
    n = len(h0)
    L = 2.0 / n
    M_init = compute_M(h0, L).max()
    print(f"Initial M = {M_init:.10f} (Together reports 0.380871)")

    h_refined, M_hist, log = refine(
        h0, L,
        max_iters=200, active_tol_init=1e-6,
        descent_tol=1e-14, verbose=True,
    )
    M_final = compute_M(h_refined, L).max()
    print(f"\nFinal M  = {M_final:.10f}")
    print(f"ΔM        = {M_final - M_init:+.3e}")
    print(f"|Δh|_inf  = {np.abs(h_refined - h0).max():.3e}")
    print(f"sum(h)    = {h_refined.sum():.10f}")

    if M_final < M_init - 1e-9:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        out = {
            "n": n,
            "L": L,
            "M_init": M_init,
            "M_final": M_final,
            "delta_M": M_final - M_init,
            "h": h_refined.tolist(),
            "M_history": list(map(float, M_hist)),
        }
        OUT.write_text(json.dumps(out, indent=2))
        print(f"\nSaved refined h to: {OUT}")
    else:
        print("\nNo improvement; not saving.")


if __name__ == "__main__":
    main()
