"""PRO-4 v2: refine Together's h* via smoothed minimax (log-sum-exp) descent.

The LP-based steepest-descent (v1) found descent directions but kept being
defeated by shifts that were just below the active threshold — they overtook
the active max after any nontrivial step.

This version uses a smooth surrogate:
    f_τ(h) = (1/τ) log Σ_j exp(τ * M(jL))
which → max_j M(jL) as τ → ∞ and has smooth gradient. We do projected
gradient descent on f_τ with τ ramped up gradually, projecting onto:
  - Σ h_i = constant (the sum)
  - h_i ∈ [0, 1] (box)
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"
OUT = ROOT / "lp_research_state" / "data" / "pro4_refined_h.json"


def compute_M(h: np.ndarray, L: float) -> np.ndarray:
    n = len(h)
    corr = np.correlate(h, 1.0 - h, mode="full")
    return L * corr[n - 1 : 2 * n - 1]


def gradient_M_matrix(h: np.ndarray, L: float, J: np.ndarray) -> np.ndarray:
    """Gradient of M(jL) w.r.t. h, stacked over j ∈ J.

    Returns matrix G of shape (|J|, n) with G[k] = ∂M(J[k]*L)/∂h.
    ∂M(jL)/∂h[a] = L * [1_{a+j<n}*(1 - h[a+j]) - 1_{a-j≥0}*h[a-j]]   for j ≥ 1
                 = L * (1 - 2*h[a])                                     for j = 0
    """
    n = len(h)
    G = np.zeros((len(J), n))
    for k, j_ in enumerate(J):
        j = int(j_)
        if j == 0:
            G[k] = L * (1.0 - 2.0 * h)
        else:
            row = np.zeros(n)
            row[: n - j] += 1.0 - h[j:]
            row[j:] -= h[: n - j]
            G[k] = L * row
    return G


def smoothed_gradient(h: np.ndarray, L: float, tau: float, top_k: int = None):
    """Return ∇f_τ(h) where f_τ = (1/τ) log Σ_j exp(τ*M(jL)).
    Also returns max_j M, the active mask of high-weight shifts, and weights.
    """
    Ms = compute_M(h, L)
    M_max = Ms.max()
    # Softmax weights (stable)
    z = tau * (Ms - M_max)
    w = np.exp(z)
    w_sum = w.sum()
    w /= w_sum
    # Restrict to top-k shifts for efficiency (after softmax)
    if top_k is not None and top_k < len(w):
        idx = np.argsort(w)[::-1][:top_k]
        # Renormalize
        w_top = w[idx]
        w_top = w_top / w_top.sum()
        G = gradient_M_matrix(h, L, idx)
        grad = w_top @ G
        return grad, M_max, idx, w_top
    else:
        # Full gradient sum (n^2 work; OK for n=600 if rare)
        J = np.arange(len(Ms))
        G = gradient_M_matrix(h, L, J)
        grad = w @ G
        return grad, M_max, J, w


def project_step(h_proposed: np.ndarray, target_sum: float) -> np.ndarray:
    """Project onto {sum = target_sum, 0 ≤ h ≤ 1} by simple clip-rescale.
    Iterative: clip to [0,1], rescale by missing mass spread over free cells.
    """
    h = np.clip(h_proposed, 0.0, 1.0)
    for _ in range(30):
        s = h.sum()
        diff = target_sum - s
        if abs(diff) < 1e-13:
            break
        # Cells with room to move in the needed direction
        if diff > 0:
            mask = h < 1.0 - 1e-15
        else:
            mask = h > 1e-15
        if not mask.any():
            break
        room = (1.0 - h) if diff > 0 else (-h)
        room_total = room[mask].sum()
        if abs(room_total) < 1e-15:
            break
        h[mask] += room[mask] * (diff / room_total)
        h = np.clip(h, 0.0, 1.0)
    return h


def refine(h0, L, *,
           taus=(1e4, 5e4, 2e5, 1e6, 5e6, 2e7),
           iters_per_tau=200,
           lr_init=1e-4,
           verbose=True,
           top_k=None):
    h = h0.copy()
    target_sum = float(h0.sum())
    M_init = compute_M(h, L).max()
    history = [M_init]
    t0 = time.time()
    for tau in taus:
        lr = lr_init
        best_M = compute_M(h, L).max()
        no_improve = 0
        for it in range(iters_per_tau):
            grad, M_max_curr, _, _ = smoothed_gradient(h, L, tau, top_k=top_k)
            # Project gradient onto Σd = 0 (sum-preserving)
            grad_proj = grad - grad.mean()
            # Zero gradient at boundary cells (so direction stays feasible)
            grad_proj = np.where((h <= 1e-12) & (grad_proj > 0), 0.0, grad_proj)
            grad_proj = np.where((h >= 1.0 - 1e-12) & (grad_proj < 0), 0.0, grad_proj)
            # Re-project sum
            free_mask = (h > 1e-12) & (h < 1.0 - 1e-12)
            if free_mask.any():
                grad_proj = grad_proj.copy()
                grad_proj[free_mask] -= grad_proj[free_mask].mean()
            g_norm = np.linalg.norm(grad_proj)
            if g_norm < 1e-14:
                if verbose:
                    print(f"  τ={tau:.0e} it={it:3d} ||grad||≈0; stopping inner loop")
                break
            d = -grad_proj / g_norm
            # Try a step
            accepted = False
            for _ in range(30):
                h_new = h + lr * d
                h_new = project_step(h_new, target_sum)
                M_new = compute_M(h_new, L).max()
                if M_new < best_M - 1e-15:
                    accepted = True
                    h = h_new
                    best_M = M_new
                    history.append(M_new)
                    no_improve = 0
                    lr *= 1.2  # grow step on success
                    break
                lr *= 0.5
            if not accepted:
                no_improve += 1
                if no_improve > 30:
                    break
                lr *= 0.5
            if verbose and (it % 20 == 0 or it < 5):
                print(
                    f"  τ={tau:.0e} it={it:3d} M={best_M:.10f} "
                    f"|grad|={g_norm:.2e} lr={lr:.2e} t={time.time()-t0:.1f}s"
                )
        if verbose:
            print(f"  τ={tau:.0e} END M={best_M:.10f} ΔM_total={best_M-M_init:+.3e}")
    return h, history


def main():
    with open(DATA) as f:
        data = json.load(f)
    h0 = np.array(data["together"]["values"], dtype=np.float64)
    n = len(h0)
    L = 2.0 / n
    M_init = compute_M(h0, L).max()
    print(f"Initial M = {M_init:.10f} (Together reports 0.380871)")

    h_refined, history = refine(
        h0, L,
        taus=(1e4, 1e5, 1e6, 1e7),
        iters_per_tau=100, lr_init=1e-3,
        verbose=True, top_k=400,
    )
    M_final = compute_M(h_refined, L).max()
    print(f"\nFinal M  = {M_final:.10f}")
    print(f"ΔM       = {M_final - M_init:+.3e}")
    print(f"|Δh|_inf = {np.abs(h_refined - h0).max():.3e}")
    print(f"sum(h)   = {h_refined.sum():.10f}")

    if M_final < M_init - 1e-9:
        out = {
            "n": n, "L": L,
            "M_init": M_init, "M_final": M_final,
            "delta_M": M_final - M_init,
            "h": h_refined.tolist(),
            "M_history": list(map(float, history)),
        }
        OUT.write_text(json.dumps(out, indent=2))
        print(f"\nSaved refined h to: {OUT}")
    else:
        print("\nNo improvement; not saving.")


if __name__ == "__main__":
    main()
