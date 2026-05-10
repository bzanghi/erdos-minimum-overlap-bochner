"""
Phase 2 — lower bound on C_5.

What is rigorously proved here:

  Trivial pigeonhole bound:  M(n)/n ≥ n/(4n − 2)  for all n ≥ 2.
                            Hence C_5 ≥ 1/4.

What is sketched but not solved:

  Sequential McCormick relaxation of the discrete LP.  Using the bilinear
  decomposition J^σ_j(h) = (linear) − (autocorrelation term), one can replace
  z_{i,j} = h_i h_{i+j} by its McCormick envelope and solve a finite LP for a
  rigorous lower bound on inf_h sup_t J(t).  We show below that the *upper*
  McCormick z ≤ min(h_i, h_{i+j}) yields only the trivial 0 bound on uniform
  λ_j, and that nontrivial bounds require the full Fourier-LP framework of
  White (Acta Arith. 2023).
"""
from __future__ import annotations
import numpy as np


def trivial_lower_bound(n: int) -> float:
    """M(n)/n ≥ n/(4n−2)."""
    return n / (4 * n - 2)


def mccormick_lp_uniform_lambda(N: int) -> tuple[float, np.ndarray]:
    """
    Compute inf_h Σ_j (2/N)(S^+_j(h) − Σ_i min(h_i, h_{i+j}))  with λ uniform.

    Since min is concave, this is a convex problem in h ∈ [0,1]^N with
    Σh = N/2.  Solved with scipy convex tools (or a simple gradient descent).

    Returns (objective_value, optimal_h).  This gives a (typically weak)
    lower bound on inf_h sup_t J(t) ≥ Σ_j λ_j J^+_j(h).
    """
    from scipy.optimize import minimize

    def obj(h):
        # objective = sum_{j=1}^{N-1} (2/N) (S^+_j − sum_i min(h_i, h_{i+j}))
        Δ = 2.0 / N
        total = 0.0
        for j in range(1, N):
            s_plus = h[j:].sum()
            mins = np.minimum(h[:N-j], h[j:]).sum()
            total += Δ * (s_plus - mins)
        return total

    def grad(h):
        # Use numeric gradient; min(a,b) has subgradient ≈ 1_{a < b} on a, etc.
        Δ = 2.0 / N
        g = np.zeros(N)
        for j in range(1, N):
            # ∂S^+_j/∂h_k = 1_{k > j-1}  (k from j to N-1, 0-indexed)
            g[j:] += Δ
            mask_a_smaller = h[:N-j] < h[j:]
            # ∂min(h_i, h_{i+j})/∂h_i = 1_{h_i < h_{i+j}} (else 0; use 1/2 for ties for stability)
            g[:N-j] += -Δ * np.where(mask_a_smaller, 1.0, np.where(h[:N-j] > h[j:], 0.0, 0.5))
            g[j:]   += -Δ * np.where(mask_a_smaller, 0.0, np.where(h[:N-j] > h[j:], 1.0, 0.5))
        return g

    # Project onto {sum = N/2, 0 ≤ h ≤ 1}
    def project(h):
        target = N / 2.0
        lo, hi = -2.0, 2.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            s = np.clip(h - mid, 0.0, 1.0).sum()
            if s > target: lo = mid
            else: hi = mid
        return np.clip(h - 0.5*(lo+hi), 0.0, 1.0)

    h = 0.5 * np.ones(N)
    h = project(h)
    res = minimize(obj, h, jac=grad, method="L-BFGS-B",
                   bounds=[(0., 1.)] * N,
                   constraints=None,
                   options={"maxiter": 500, "ftol": 1e-12})
    h = project(res.x)
    return obj(h), h


if __name__ == "__main__":
    print("Trivial lower bound  M(n)/n ≥ n/(4n−2):")
    for n in [10, 100, 1000, 10000]:
        print(f"  n={n:5d}: bound = {trivial_lower_bound(n):.6f}")
    print(f"  n→∞:   bound → 0.250000  ⇒  C_5 ≥ 1/4")

    # Verify against brute force:
    import json
    table = json.load(open("Mn_brute.json"))
    print("\nVerified consistency with brute force (M(n)/n ≥ trivial bound):")
    for row in table:
        n, M = row["n"], row["M"]
        tb = trivial_lower_bound(n)
        ok = "✓" if M / n >= tb else "✗"
        print(f"  n={n:2d}: M(n)/n = {M/n:.5f}  ≥  {tb:.5f}  {ok}")

    print("\nMcCormick-relaxation LP (uniform λ_j) — illustrative:")
    for N in [40, 80, 120]:
        val, _ = mccormick_lp_uniform_lambda(N)
        # The LP gives a lower bound on Σ λ_j J^+_j(h), which is ≤ sup_t J(t).
        # Comparing against avg-J bound (3/2 - c_h)/2 for c_h=1: avg = 1/4.
        print(f"  N={N:3d}: LP value = {val:.6f}  (this is Σ λ_j J^+_j, NOT directly C_5)")

    print("\nNOTE: This relaxation, with uniform λ, recovers a value ≈ 0 for")
    print("symmetric h, so it does NOT improve the trivial 1/4 bound on C_5.")
    print("White's 0.379005 lower bound uses a specific Fourier-analytic dual")
    print("that goes well beyond uniform-λ McCormick.")
