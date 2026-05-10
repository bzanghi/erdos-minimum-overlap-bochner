"""
Hankel-PSD constraint probe.

For f >= 0 on [-1, 1] with m_k = ∫ x^k f(x) dx, the Hausdorff moment theorem says
H_n[i,j] = m_{i+j} is PSD for all n, and the interval-positivity matrix
A_n[i,j] = m_{i+j} - m_{i+j+2} is also PSD on [-1, 1].

We compute m_k from the LP-optimal (c*, d*) via Fourier expansion:
    m_k = (1/2) ⟨x^k, 1⟩ + Σ_{j=1..T} (c_j α_j^(k) + d_j β_j^(k))
where α_j^(k) = ∫_{-1}^{1} x^k cos(πjx) dx, β_j^(k) = ∫_{-1}^{1} x^k sin(πjx) dx.

These integrals have CLOSED FORMS (integration by parts k times). For low k (≤6),
the coefficients α, β decay as O(1/(πj)^{k+1}) or so, making truncation tails small.

Goal of this probe: solve the current Bochner-augmented LP, reconstruct truncated
m_k from (c*, d*), compute eigenvalues of H_n. If H_n has a strongly negative
min eigenvalue, the Hankel constraint cuts. If marginal, it might still be useful
combined with a tail bound.
"""
from __future__ import annotations
import argparse
import sys
import time
from pathlib import Path

import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem


def fourier_coeffs_of_xk(k_max: int, j_max: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute α_j^(k) = ∫_{-1}^{1} x^k cos(πjx) dx and β_j^(k) = ∫_{-1}^{1} x^k sin(πjx) dx,
    plus α_0^(k) = ∫_{-1}^{1} x^k dx.

    Returns:
      alpha0:  shape (k_max+1,)      α_0^(k)
      alpha:   shape (k_max+1, j_max)  α_j^(k) for j = 1..j_max
      beta:    shape (k_max+1, j_max)  β_j^(k) for j = 1..j_max

    Uses recurrence from integration by parts:
      ∫ x^k cos(πjx) = [x^k sin(πjx)/(πj)] - k/(πj) ∫ x^(k-1) sin(πjx) dx
      ∫ x^k sin(πjx) = [-x^k cos(πjx)/(πj)] + k/(πj) ∫ x^(k-1) cos(πjx) dx
    Evaluated over [-1, 1]: x^k at x=±1 gives (-1)^k or (-1)^{k+1} times 1.
    sin(πj·(±1)) = 0; cos(πj·(±1)) = (-1)^j.

    Boundary term in ∫ x^k cos = (x^k sin(πjx)/(πj)) from -1 to 1 = 0.
    Boundary term in ∫ x^k sin = (-x^k cos(πjx)/(πj)) from -1 to 1
        = -(1^k · (-1)^j - (-1)^k · (-1)^j)/(πj)
        = -(-1)^j (1 - (-1)^k)/(πj)
        = 0 if k even, -2(-1)^j/(πj) if k odd.
    """
    j = np.arange(1, j_max + 1)
    alpha = np.zeros((k_max + 1, j_max))
    beta = np.zeros((k_max + 1, j_max))
    alpha0 = np.zeros(k_max + 1)

    # k = 0
    alpha0[0] = 2.0  # ∫_{-1}^{1} 1 dx
    alpha[0, :] = 0.0  # ∫ cos(πjx) dx = 0
    beta[0, :] = 0.0   # ∫ sin(πjx) dx = 0 (odd function)

    inv_pij = 1.0 / (np.pi * j)
    sign_j = (-1.0) ** j

    for kk in range(1, k_max + 1):
        # α_j^(k) = boundary - (k/(πj)) β_j^(k-1).  boundary = 0 for cos.
        alpha[kk, :] = -kk * inv_pij * beta[kk - 1, :]
        # β_j^(k) = boundary + (k/(πj)) α_j^(k-1)
        if kk % 2 == 1:
            boundary = -2 * sign_j * inv_pij
        else:
            boundary = np.zeros_like(j, dtype=float)
        beta[kk, :] = boundary + kk * inv_pij * alpha[kk - 1, :]
        # alpha0^(k) = ∫ x^k dx = (1 - (-1)^(k+1))/(k+1) = 2/(k+1) if k even, else 0
        alpha0[kk] = 2.0 / (kk + 1) if kk % 2 == 0 else 0.0

    return alpha0, alpha, beta


def moments_from_cd(c_star: np.ndarray, d_star: np.ndarray, k_max: int) -> np.ndarray:
    """m_k = (1/2) α_0^(k) + Σ_j (c_j α_j^(k) + d_j β_j^(k)) for k = 0..k_max."""
    T = len(c_star)
    alpha0, alpha, beta = fourier_coeffs_of_xk(k_max, T)
    m = 0.5 * alpha0 + (alpha @ c_star) + (beta @ d_star)
    return m


def build_hankel(m: np.ndarray, n: int) -> np.ndarray:
    """H_n[i,j] = m_{i+j} for i, j in [0, n]. Requires m of length ≥ 2n+1."""
    H = np.empty((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(n + 1):
            H[i, j] = m[i + j]
    return H


def build_interval_pos(m: np.ndarray, n: int) -> np.ndarray:
    """A_n[i,j] = m_{i+j} - m_{i+j+2} for i, j in [0, n-1].  Interval positivity on [-1,1].
    Requires m of length ≥ 2n+1."""
    A = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = m[i + j] - m[i + j + 2]
    return A


def estimate_tail_bound(k: int, T_lp: int, j_far: int = 100000) -> float:
    """Estimate the maximum possible value of |Σ_{j>T_lp} (c_j α_j^(k) + d_j β_j^(k))|
    given the LP bounds |c_j|, |d_j| ≤ 2/π.

    Returns an upper bound on the truncation error in m_k.
    """
    alpha0, alpha, beta = fourier_coeffs_of_xk(k, j_far)
    # Tail = j_max - T_lp:
    tail_alpha = np.abs(alpha[k, T_lp:]).sum()
    tail_beta = np.abs(beta[k, T_lp:]).sum()
    bound_c = (2.0 / np.pi) * tail_alpha
    bound_d = (2.0 / np.pi) * tail_beta
    return float(bound_c + bound_d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, default=4)
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--T", type=int, default=800)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=15)
    ap.add_argument("--n_hankel", type=int, default=6)
    args = ap.parse_args()

    ROWS = {
        1: (0.015, 0.381,  -0.02, 0.02),
        2: (0.015, 0.385,  -0.02, 0.02),
        3: (0.020, 0.375,  -0.02, 0.02),
        4: (0.004, 0.3875, -0.02, 0.02),
        5: (0.000, 0.4,    -0.02, 0.02),
        6: (0.000, 0.381,  -0.02, 0.02),
        7: (0.030, 0.375,  -0.02, 0.02),
    }
    h, p, qm, qp = ROWS[args.row]
    print(f"# Hankel probe at row {args.row}  N={args.N} T={args.T} bochner_n={args.bochner_n}")
    print(f"# Hankel level n = {args.n_hankel}  (matrix size {args.n_hankel+1}×{args.n_hankel+1})")
    print()

    t0 = time.time()
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        args.N, args.T, args.R, h, h, p, p, qm, qp,
        bochner_n=args.bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  LP solved in {time.time() - t0:.1f}s, status = {prob.status}")
    print(f"  Ω* = {prob.value:.10f}")

    c_star = np.array(c.value); d_star = np.array(d.value)

    # 1) Compute truncated moments
    k_max = 2 * args.n_hankel + 2
    m_trunc = moments_from_cd(c_star, d_star, k_max)
    print(f"\n  Truncated moments m_0..m_{k_max}:")
    for k in range(min(k_max + 1, 12)):
        tail = estimate_tail_bound(k, args.T)
        print(f"    m_{k} = {m_trunc[k]:+.6f}    tail bound: ±{tail:.2e}")

    # 2) Build Hankel matrix and check PSDness
    H = build_hankel(m_trunc, args.n_hankel)
    eigs_H = np.linalg.eigvalsh(H)
    print(f"\n  Hankel H_{args.n_hankel} eigenvalues:")
    print(f"    min = {eigs_H[0]:+.6e}")
    print(f"    max = {eigs_H[-1]:+.6e}")
    print(f"    all = {eigs_H}")

    A = build_interval_pos(m_trunc, args.n_hankel)
    eigs_A = np.linalg.eigvalsh(A)
    print(f"\n  Interval-pos A_{args.n_hankel} eigenvalues:")
    print(f"    min = {eigs_A[0]:+.6e}")
    print(f"    max = {eigs_A[-1]:+.6e}")
    print(f"    all = {eigs_A}")

    # 3) Verdict
    print(f"\n  === verdict ===")
    if eigs_H[0] < -1e-4:
        print(f"  ✓ HANKEL CUTS: min eig {eigs_H[0]:+.2e} dominates plausible tail ~1e-6")
    elif eigs_H[0] < -1e-6:
        print(f"  ~ marginal: min eig {eigs_H[0]:+.2e}, may be tail-bound noise")
    else:
        print(f"  ✗ HANKEL DOES NOT CUT at this scale")

    if eigs_A[0] < -1e-4:
        print(f"  ✓ INTERVAL POS CUTS: min eig {eigs_A[0]:+.2e}")
    elif eigs_A[0] < -1e-6:
        print(f"  ~ interval pos marginal")
    else:
        print(f"  ✗ INTERVAL POS DOES NOT CUT")


if __name__ == "__main__":
    main()
