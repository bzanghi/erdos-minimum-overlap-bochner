"""PRO-26 Phase 2a v1: Chebyshev ansatz + BFGS optimization on Together's h*.

Strategy:
1. Project Together's h* onto a Chebyshev polynomial basis of degree K on [0, 2].
2. Re-evaluate h from the truncated basis (≤ K coefficients) — this is the smooth ansatz.
3. Optimize via scipy.optimize.minimize (BFGS) on the K+1 coefficients.
4. Objective: max_j M(jL) using a soft-max (LogSumExp) surrogate.
5. Constraints: ∫h = 1 enforced by adjusting constant term; 0 ≤ h ≤ 1 via clipping (soft penalty).

The hope: smooth perturbations can escape Together's discrete local min that PRO-4's
600-DoF descent could not.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def load_together_h() -> np.ndarray:
    with open(DATA) as f:
        return np.array(json.load(f)["together"]["values"], dtype=np.float64)


def cheb_basis(n: int, K: int):
    """Return (B, x_grid) where B[k, i] is the value of T_k at x_i (i = cell midpoint
    in [0, 2], mapped to [-1, 1])."""
    # Cell midpoints in [0, 2]
    L = 2.0 / n
    x_phys = (np.arange(n) + 0.5) * L          # in [L/2, 2 - L/2]
    x_norm = (x_phys - 1.0)                     # in [-1+L/2, 1-L/2] → essentially [-1, 1]
    # Chebyshev T_0..T_K evaluated at x_norm
    B = np.zeros((K + 1, n))
    for k in range(K + 1):
        c = np.zeros(K + 1)
        c[k] = 1.0
        B[k] = cheb.chebval(x_norm, c)
    return B, x_norm


def project_onto_cheb(h: np.ndarray, K: int) -> np.ndarray:
    """Least-squares projection of h onto T_0..T_K."""
    B, _ = cheb_basis(len(h), K)
    # Solve B^T a = h  (where a is the (K+1)-vector of coefficients)
    a, *_ = np.linalg.lstsq(B.T, h, rcond=None)
    return a


def h_from_coefs(a: np.ndarray, n: int, K: int) -> np.ndarray:
    """Reconstruct h from Chebyshev coefficients (K+1 of them)."""
    B, _ = cheb_basis(n, K)
    return a @ B  # shape (n,)


def compute_M(h: np.ndarray, L: float) -> np.ndarray:
    """BUG FIXED 2026-05-18: include all 2n-1 lags (positive + negative).
    For asymmetric h, M(t) ≠ M(-t); BFGS exploited the slicing bug to find
    bad h with low positive-lag M but huge negative-lag M."""
    corr = np.correlate(h, 1.0 - h, mode="full")
    return L * corr


def soft_max(values: np.ndarray, tau: float) -> float:
    """Stable LogSumExp soft-max."""
    M_max = values.max()
    return M_max + np.log(np.sum(np.exp(tau * (values - M_max)))) / tau


def objective_with_constraints(a: np.ndarray, n: int, K: int, L: float,
                                target_sum: float, tau: float,
                                box_penalty: float) -> float:
    """Compute soft-max M plus penalties for sum constraint and box violation."""
    h = h_from_coefs(a, n, K)

    # Adjust constant coefficient to enforce sum exactly
    # ∫h = L * Σh; if Σh ≠ target_sum, shift by (target_sum - Σh) / n on every cell
    # That's equivalent to adjusting the constant term: shift a[0] by ... hmm
    # The basis function T_0 = 1, so B[0, i] = 1 for all i. So adding δ to a[0]
    # adds δ to every h_i. To make Σh = target_sum, set δ = (target_sum - Σh) / n.
    shift = (target_sum - h.sum()) / n
    h = h + shift

    # Box penalty: 0 ≤ h ≤ 1
    box_violation = np.maximum(0, -h).sum() + np.maximum(0, h - 1.0).sum()

    # M values
    Ms = compute_M(h, L)
    M_smax = soft_max(Ms, tau)

    return M_smax + box_penalty * box_violation


def evaluate_true_max(a: np.ndarray, n: int, K: int, L: float,
                       target_sum: float) -> tuple[float, np.ndarray]:
    """Reconstruct h and return (max_M, h_clipped)."""
    h = h_from_coefs(a, n, K)
    shift = (target_sum - h.sum()) / n
    h = h + shift
    h_clipped = np.clip(h, 0.0, 1.0)
    # Adjust again after clipping
    if abs(h_clipped.sum() - target_sum) > 1e-9:
        scale = target_sum / h_clipped.sum()
        h_clipped = np.clip(h_clipped * scale, 0.0, 1.0)
    Ms = compute_M(h_clipped, L)
    return float(Ms.max()), h_clipped


def main():
    print("=" * 78)
    print("PRO-26 Phase 2a v1: Chebyshev ansatz + BFGS on Together's h*")
    print("=" * 78)

    h_init = load_together_h()
    n = len(h_init)
    L = 2.0 / n
    target_sum = float(h_init.sum())
    M_init = compute_M(h_init, L).max()
    print(f"n = {n}, L = {L}, sum_h = {target_sum}, M(h_init) = {M_init:.10f}\n")

    # Try various K (Chebyshev degree)
    for K in [10, 20, 30, 50, 100, 200]:
        print(f"--- K = {K} ({K+1} Chebyshev coefficients) ---")
        # Projection error
        a0 = project_onto_cheb(h_init, K)
        h_proj = h_from_coefs(a0, n, K)
        proj_err = np.abs(h_proj - h_init).max()
        # After clipping + sum-correction
        M_proj, _ = evaluate_true_max(a0, n, K, L, target_sum)
        print(f"  Proj error ‖h_proj - h_init‖_inf = {proj_err:.3e}")
        print(f"  M(proj) = {M_proj:.10f}  Δ vs h_init = {M_proj - M_init:+.3e}")

        # Optimize
        tau = 1e6
        box_penalty = 100.0
        res = minimize(
            objective_with_constraints, a0,
            args=(n, K, L, target_sum, tau, box_penalty),
            method="BFGS",
            options={"maxiter": 200, "gtol": 1e-10, "disp": False},
        )
        M_opt, h_opt = evaluate_true_max(res.x, n, K, L, target_sum)
        print(f"  After BFGS ({res.nit} iters): M = {M_opt:.10f}  Δ vs init = {M_opt - M_init:+.3e}")
        print(f"  ‖h_opt - h_init‖_inf = {np.abs(h_opt - h_init).max():.3e}\n")


if __name__ == "__main__":
    main()
