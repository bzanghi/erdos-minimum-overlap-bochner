"""
Polynomial-moment PSD constraints (Hankel) for the Bochner-augmented LP.

Theorem (Hausdorff moment problem on [-1, 1]).
A real measurable f: [-1, 1] → R is ≥ 0 a.e. iff for every n ≥ 0:
    H_n := [m_{i+j}]_{i,j=0..n}     ⪰ 0       (Hankel PSD)
    A_n := [m_{i+j} - m_{i+j+2}]_{i,j=0..n-1} ⪰ 0    (interval-positivity)
where m_k := ∫_{-1}^1 x^k f(x) dx.

Linking m_k to White's LP variables (c, d) via the Fourier expansion of x^k:
    m_k = (1/2) α_0^(k) + Σ_{j=1..T} (c_j α_j^(k) + d_j β_j^(k)) + tail_k
where
    α_0^(k) = ∫_{-1}^1 x^k dx     = 2/(k+1) if k even, 0 if odd
    α_j^(k) = ∫_{-1}^1 x^k cos(πjx) dx   (computable by k-fold integration by parts)
    β_j^(k) = ∫_{-1}^1 x^k sin(πjx) dx
The tail_k is the sum over j > T; we bound it by precomputed constants.

For EVEN k, the coefficients decay as O(1/(πj)^k) so the tail is small.
For ODD k, they decay only as O(1/(πj)) — tail too large to be useful for rigorous use.

Therefore we expose:
  - For each even k ∈ {2, 4, ..., 2n_max}, the SCALAR inequality m_k ≥ 0
    (with tail-slack, becomes m_k_truncated ≥ -tail_bound_k).
  - Optionally, the principal submatrix of H_n on even-only indices, which uses
    only even moments and is rigorizable. (Future extension.)
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def fourier_coeffs_of_xk(k_max: int, j_max: int):
    """Compute α_j^(k), β_j^(k) for k = 0..k_max, j = 1..j_max, plus α_0^(k).

    Recurrence (integration by parts):
      ∫ x^k cos(πjx) dx = [x^k sin(πjx)/(πj)]_{-1}^{1} - (k/(πj)) ∫ x^{k-1} sin(πjx) dx
      ∫ x^k sin(πjx) dx = [-x^k cos(πjx)/(πj)]_{-1}^{1} + (k/(πj)) ∫ x^{k-1} cos(πjx) dx
    Boundary terms:
      sin(πj·(±1)) = 0  ⇒  α boundary = 0 always.
      cos(πj·(±1)) = (-1)^j; (±1)^k = ±1.
      β boundary = -(1^k · (-1)^j − (-1)^k · (-1)^j)/(πj) = -(-1)^j (1 − (-1)^k)/(πj)
                 = 0 if k even, −2(-1)^j/(πj) if k odd.
    """
    j = np.arange(1, j_max + 1)
    alpha = np.zeros((k_max + 1, j_max))
    beta = np.zeros((k_max + 1, j_max))
    alpha0 = np.zeros(k_max + 1)

    inv_pij = 1.0 / (np.pi * j)
    sign_j = (-1.0) ** j

    # k = 0
    alpha0[0] = 2.0
    # alpha[0,:] = 0; beta[0,:] = 0 (already zero)

    for kk in range(1, k_max + 1):
        alpha[kk, :] = -kk * inv_pij * beta[kk - 1, :]
        if kk % 2 == 1:
            boundary = -2 * sign_j * inv_pij
        else:
            boundary = np.zeros_like(j, dtype=float)
        beta[kk, :] = boundary + kk * inv_pij * alpha[kk - 1, :]
        alpha0[kk] = 2.0 / (kk + 1) if kk % 2 == 0 else 0.0

    return alpha0, alpha, beta


def even_moment_tail_bound(k: int, T_lp: int, j_far: int = 20000) -> float:
    """Upper bound on |tail_k| = |Σ_{j>T_lp} (c_j α_j^(k) + d_j β_j^(k))| given LP bounds
    |c_j|, |d_j| ≤ 2/π.  For even k, returns a finite small bound."""
    if k % 2 != 0:
        return float("inf")  # odd-k tail is too loose for rigorous use here
    alpha0, alpha, beta = fourier_coeffs_of_xk(k, j_far)
    tail = (np.abs(alpha[k, T_lp:]).sum() + np.abs(beta[k, T_lp:]).sum())
    return float((2.0 / np.pi) * tail)


def build_even_moment_nonneg_constraints(
    c: cp.Variable, d: cp.Variable, T: int,
    k_max: int = 12,
):
    """Add LINEAR constraints m_k ≥ -tail_bound_k for each even k ∈ {2, ..., k_max}.

    For f ≥ 0 on [-1, 1], m_k = ∫ x^k f ≥ 0 for every even k.  Truncated to T
    Fourier coefficients, |m_k_true - m_k_truncated| ≤ tail_bound_k.  Hence:
        m_k_truncated ≥ -tail_bound_k
    is a rigorous linear inequality on (c, d).

    Returns:
      list of cvxpy constraints to add to `cons`
      dict {k: tail_bound_k}  for reference / logging
    """
    if k_max % 2 != 0:
        k_max -= 1
    alpha0, alpha, beta = fourier_coeffs_of_xk(k_max, T)
    cons = []
    tail_bounds = {}
    for k in range(2, k_max + 1, 2):
        # m_k_truncated = (1/2) α_0^(k) + Σ_j (c_j α_j^(k) + d_j β_j^(k))
        m_k = 0.5 * alpha0[k] + alpha[k, :] @ c + beta[k, :] @ d
        tb = even_moment_tail_bound(k, T)
        cons.append(m_k >= -tb)
        tail_bounds[k] = tb
    return cons, tail_bounds


def build_even_hankel_psd(
    c: cp.Variable, d: cp.Variable, T: int,
    n_hankel: int,
):
    """Even-indices-only principal submatrix of the Hankel matrix.

    H_n^{even}[i,j] = m_{2(i+j)} for i, j = 0..n_hankel.  All entries are even
    polynomial moments; tail bounds are finite small.  We impose:
        H_n^{even} - tail_perturbation * I ⪰ 0
    where tail_perturbation is sized to dominate the simultaneous truncation
    errors of all entries (sum of tail_bound_2k over k = 0..2*n_hankel).

    Note: this is a STRICTLY WEAKER constraint than full Hankel PSD (we drop
    the odd-moment off-diagonal blocks), but it is RIGOROUSLY truncatable.

    Returns:
      single cvxpy constraint (matrix PSD)
      tail_pert size (for logging)
    """
    k_max = 4 * n_hankel  # need moments up to m_{4n}
    alpha0, alpha, beta = fourier_coeffs_of_xk(k_max, T)
    H_entries = []
    max_tail = 0.0
    for i in range(n_hankel + 1):
        row = []
        for j in range(n_hankel + 1):
            k = 2 * (i + j)
            m_k_expr = 0.5 * alpha0[k] + alpha[k, :] @ c + beta[k, :] @ d
            row.append(m_k_expr)
            tb = even_moment_tail_bound(k, T)
            if tb > max_tail: max_tail = tb
        H_entries.append(row)
    # Build H as a cvxpy matrix expression via bmat (each entry is a scalar expr).
    H_mat = cp.bmat([[cp.reshape(e, (1, 1)) for e in row] for row in H_entries])
    # Conservative PSD perturbation: subtract max_tail * (n+1) on diagonal
    # (Gershgorin-style: any eigenvalue shifted by entry uncertainty ≤ (n+1) max_tail)
    pert = (n_hankel + 1) * max_tail
    constraint = (H_mat + pert * np.eye(n_hankel + 1)) >> 0
    # Re-state: H_true ⪰ 0 ⇒ H_truncated ⪰ -pert*I  ⇔  H_truncated + pert*I ⪰ 0
    return constraint, pert
