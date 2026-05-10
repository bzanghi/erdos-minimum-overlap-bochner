"""
Reproduce White (Acta Arith. 2023, arXiv:2201.05704) Section 4 simplified LP.

Variables: Ω, w_1, ..., w_N
Minimize:  Ω
Subject to:
   (4.1)  0 ≤ w_j ≤ Ω,                     1 ≤ j ≤ N
   (4.2)  Σ_j w_j = N/4
   (4.3)  Σ_j α^-_{j,2m} w_j ≤ 0,          1 ≤ m ≤ R
   (4.4)  L^3 Σ_j (j-1)^2 w_j ≤ 1/3
where L = 2/N and α^-_{j,m} = cos(πmL(j-1/2)/2) − πmL/4.

Proves: ∥M∥_∞ ≥ Ω*  for all even M  (and so µ ≥ Ω*).

White reports for N = 80000, R = 20:  Ω* ≈ 0.375169.
"""
import numpy as np
import cvxpy as cp
import time, json, sys


def alpha_minus(j_arr, m, L):
    """α^-_{j,m} = cos(πmL(j-1/2)/2) − πmL/4."""
    x_mid = (j_arr - 0.5) * L           # midpoint of [(j-1)L, jL]
    return np.cos(np.pi * m * x_mid / 2) - np.pi * m * L / 4


def solve_simplified_LP(N: int, R: int, solver: str = "SCS", verbose: bool = False):
    L = 2.0 / N
    j = np.arange(1, N + 1)

    w = cp.Variable(N, nonneg=True)
    Omega = cp.Variable()

    constraints = [
        w <= Omega,
        cp.sum(w) == N / 4,
        L**3 * cp.sum((j - 1) ** 2 * w) <= 1.0 / 3,
    ]
    # (4.3): R constraints on even cosine coefficients (m = 2, 4, ..., 2R)
    for m_half in range(1, R + 1):
        m = 2 * m_half
        a_minus = alpha_minus(j, m, L)
        constraints.append(a_minus @ w <= 0)

    prob = cp.Problem(cp.Minimize(Omega), constraints)
    t0 = time.time()
    prob.solve(solver=solver, verbose=verbose)
    dt = time.time() - t0
    return float(prob.value), w.value, dt


if __name__ == "__main__":
    print("White simplified LP (Section 4) — symmetric case, target ≈ 0.375")
    print("=" * 70)
    for (N, R) in [(500, 5), (2000, 10), (10000, 15), (20000, 20)]:
        try:
            val, _, dt = solve_simplified_LP(N, R, solver="CLARABEL", verbose=False)
            print(f"  N = {N:6d}  R = {R:3d}   Ω* = {val:.7f}   t = {dt:5.1f}s")
        except Exception as e:
            print(f"  N = {N:6d}  R = {R:3d}   ERROR: {e}")
