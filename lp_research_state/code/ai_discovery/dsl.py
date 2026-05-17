"""Constraint DSL: parametric families of candidate convex constraints.

Each family is `family(theta) -> constraint_fn`, where `constraint_fn` is
called with the cvxpy variables of build_problem and returns a list of
extra cvxpy constraints to add.

VALIDITY PRINCIPLE: every family must produce constraints satisfied by
EVERY feasible f. If invalid, the SDP becomes infeasible or produces a
false lower bound. We document validity per family and verify by
solving with the family added (checking status != infeasible).

== Implemented Families ==

F2. cell_envelope_high_freq(m_high)
    Add the existing (W.1) cell-envelope cosine constraint at *higher*
    frequencies m ∈ (2R, 2R+m_high]. Validity: same as the existing
    family. This isn't truly novel but tests whether more lags help.

F3. test_function_T5pk(k)
    Generalize T5p (which uses test function φ(x) = 1 - cos(πx)) to
    φ_k(x) = 1 - cos(πkx) for integer k ≥ 1. φ_k ≥ 0 always.
    Constraint: c^T Q_k c + d^T Q_k d ≤ 1/2, where
    Q_k = I_T - (1/2)(E_{+k} + E_{-k}) and E_{±k} is the matrix with
    1's on the ±k-th off-diagonal.

    Validity: Q_k is PSD with eigenvalues 1 - cos(kπj/(T+1)) ∈ [0, 2].
    Proof: φ_k ≥ 0 + ∫(f - f²) φ_k = (1/2) - (c.T Q_k c + d.T Q_k d) ≥ 0
    (following the derivation in white_full_convex.py:221-226 for k=1).

    Genuinely NEW family. T5p is k=1; we add k=2, 3, ..., k_max.

F4. test_function_sumcos(theta)
    More general test function: φ(x) = sum_k (theta_k · (1 - cos(πkx)))
    for nonneg coefficients theta. φ ≥ 0 for any nonneg theta (sum of
    nonneg).
    Constraint: c^T (sum_k theta_k Q_k) c + d^T (sum_k theta_k Q_k) d ≤
                                                              (sum_k theta_k)/2.

    Validity: linear combination of F3's at varying k, scaled by sum.

F5. test_function_fejer(n_fejer)
    φ(x) = K_n(x) = Fejér kernel = (1/(n+1)) · |sum_{k=0}^n e^{iπkx}|².
    Always nonneg (it's a squared modulus / (n+1)). Has compact Fourier
    support: K_n = sum_{|k|≤n} (1 - |k|/(n+1)) cos(πkx).
    Constraint built from F4's coefficients.

F6. f_squared_times_one_minus(theta)
    Test f²(1 - f) ≥ 0 (since f ≤ 1) against constant 1: ∫f² - ∫f³ ≥ 0.
    The third moment ∫f³ is CUBIC, not directly convex. Drop unless
    rewritable as SOC.

F7. cell_envelope_complement(theta)
    Apply (W.1)-style constraint to 1-f instead of f. Since 1-f ∈ [0, 1]
    when f ∈ [0, 1], we get a parallel SDP with new variables 1-w, 1-v.
    Many constraints transfer. NOT yet implemented.
"""
from __future__ import annotations
import cvxpy as cp
import numpy as np
from typing import Callable, Dict


# ============================================================================
# F3: test_function_T5pk -- generalize T5p to higher k
# ============================================================================

def build_Qk(T: int, k: int) -> np.ndarray:
    """Q_k = I_T - (1/2)(E_{+k} + E_{-k})."""
    Q = np.eye(T)
    if k > 0 and k < T:
        Q -= 0.5 * np.eye(T, k=k)
        Q -= 0.5 * np.eye(T, k=-k)
    return Q


def family_T5pk(k: int) -> Callable:
    """Test function φ_k(x) = 1 - cos(πkx) for f² ≤ f.

    Constraint: c^T Q_k c + d^T Q_k d ≤ 1/2.
    """
    if k < 1:
        raise ValueError("k must be >= 1")

    def constraint_fn(Omega, w, v, c, d, eps, dlt, cfg):
        T = cfg["T"]
        if k >= T:
            return []  # out of Fourier range
        Q = build_Qk(T, k)
        return [cp.quad_form(c, cp.psd_wrap(Q)) + cp.quad_form(d, cp.psd_wrap(Q)) <= 0.5]

    return constraint_fn


# ============================================================================
# F4: test_function_sumcos -- linear combination of T5pk's
# ============================================================================

def family_T5p_sumcos(theta: np.ndarray) -> Callable:
    """φ(x) = sum_k theta_k · (1 - cos(πkx)), theta_k ≥ 0."""
    theta = np.asarray(theta, dtype=float)
    theta = np.abs(theta)  # enforce nonneg
    K = len(theta)

    def constraint_fn(Omega, w, v, c, d, eps, dlt, cfg):
        T = cfg["T"]
        # Q = sum_k theta_k · Q_k = (sum_k theta_k) I - sum_k theta_k · (E_{+k} + E_{-k})/2
        total = float(theta.sum())
        if total <= 0:
            return []
        Q = np.zeros((T, T))
        for k_idx in range(K):
            if theta[k_idx] == 0:
                continue
            kk = k_idx + 1  # k starts from 1
            if kk >= T:
                continue
            Q += theta[k_idx] * build_Qk(T, kk)
        # RHS = (sum_k theta_k) / 2
        rhs = total / 2.0
        return [cp.quad_form(c, cp.psd_wrap(Q)) + cp.quad_form(d, cp.psd_wrap(Q)) <= rhs]

    return constraint_fn


# ============================================================================
# F2: cell_envelope_high_freq
# ============================================================================

def family_cell_env_high_freq(m_high: int) -> Callable:
    """Add (W.1) cosine cell-envelope at m = 2R+1, ..., 2R+m_high (even m only).

    Same shape as the existing family at line 182 of white_full_convex.py,
    just with larger m.
    """
    def constraint_fn(Omega, w, v, c, d, eps, dlt, cfg):
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from white_full_convex import cos_cell_bounds_exact
        N = cfg["N"]
        T = cfg["T"]
        R = cfg["R"]
        L = 2.0 / N
        j = np.arange(1, N + 1)
        cons = []
        for m in range(2 * R + 1, 2 * R + m_high + 1):
            sin_pi_half_m = np.sin(np.pi * m / 2)
            if m % 2 != 0:
                continue  # odd m needs eps which we don't extend
            half = m // 2
            if half - 1 >= T:
                continue
            am = 0.5 * c[half - 1]
            bm = 0.5 * d[half - 1]
            a_minus, _ = cos_cell_bounds_exact(j, m, L)
            lhs = (L / 2) * (a_minus @ (w + v))
            rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
            cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)
        return cons

    return constraint_fn


# ============================================================================
# F5: fejer kernel
# ============================================================================

def family_fejer(n_fejer: int) -> Callable:
    """φ(x) = K_n(x) (Fejér kernel), positive trig polynomial of degree n.

    K_n(x) = sum_{k=-n}^n (1 - |k|/(n+1)) cos(πkx).
    Tests f² ≤ f against K_n: ∫(f - f²) K_n ≥ 0.

    Coefficients: theta_k = (1 - k/(n+1)) for k=1..n. theta_0 = 1 but
    we only use k≥1 coefficients (the k=0 part contributes to RHS).
    """
    theta = np.array([(1.0 - k / (n_fejer + 1)) for k in range(1, n_fejer + 1)])
    return family_T5p_sumcos(theta)


# ============================================================================
# Registry
# ============================================================================

FAMILIES: Dict[str, Callable] = {
    "F3_T5p_k2":  lambda: family_T5pk(2),
    "F3_T5p_k3":  lambda: family_T5pk(3),
    "F3_T5p_k4":  lambda: family_T5pk(4),
    "F3_T5p_k5":  lambda: family_T5pk(5),
    "F3_T5p_k6":  lambda: family_T5pk(6),
    "F3_T5p_k8":  lambda: family_T5pk(8),
    "F4_sumcos_uniform_5":   lambda: family_T5p_sumcos(np.ones(5)),
    "F4_sumcos_uniform_10":  lambda: family_T5p_sumcos(np.ones(10)),
    "F4_sumcos_decay_1overK_5":  lambda: family_T5p_sumcos(np.array([1/k for k in range(1, 6)])),
    "F4_sumcos_decay_1overK_10": lambda: family_T5p_sumcos(np.array([1/k for k in range(1, 11)])),
    "F5_fejer_5":  lambda: family_fejer(5),
    "F5_fejer_10": lambda: family_fejer(10),
    "F5_fejer_20": lambda: family_fejer(20),
    "F2_cell_env_high_m4":  lambda: family_cell_env_high_freq(4),
    "F2_cell_env_high_m8":  lambda: family_cell_env_high_freq(8),
}


def list_families():
    return list(FAMILIES.keys())
