"""PRO-29: Spectral reformulation of μ = inf_h sup_t M(h, t) as a min-max
eigenvalue-like problem.

Reformulation
-------------
Let `(T_t f)(x) := f(x + t)` denote translation. Define the operator
    A_t := I - T_t  (on L²[0, 2] with zero extension)

Then M(h, t) = ⟨h, A_t h⟩ — quadratic form in h.

For h orthogonal to the constant function (away from the integral constraint),
sup over unit-norm h of ⟨h, A_t h⟩ = λ_max(symmetric part of A_t) = λ_max((A_t + A_t^T)/2).

But our constraints are ‖h‖_∞ ≤ 1 and ⟨h, 1⟩ = 1, not ‖h‖_2 = 1. So the
quadratic-form-eigenvalue analogy is approximate.

This prototype:
1. Builds A_t as an N×N matrix (cell-shift discretization, n=600 cells).
2. Computes the spectrum at a sweep of shifts t.
3. Examines whether the spectrum gives any useful bound on μ.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
from scipy import linalg as spla

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "lp_research_state" / "data" / "together_f_star.json"


def shift_matrix(n: int, j: int) -> np.ndarray:
    """T_j: cells [(i+j) mod n] → cell i. With zero extension (not periodic),
    use a non-periodic shift: T_j[i, i+j] = 1 if i+j < n, else 0.
    """
    T = np.zeros((n, n))
    for i in range(n - j):
        T[i, i + j] = 1.0
    return T


def A_t(n: int, j: int) -> np.ndarray:
    """A_j := I - T_j."""
    return np.eye(n) - shift_matrix(n, j)


def M_via_matrix(h: np.ndarray, A: np.ndarray, L: float) -> float:
    """M = L * h^T A h."""
    return L * float(h @ A @ h)


def M_via_correlate(h: np.ndarray, L: float, j: int) -> float:
    """Reference: M(jL) using np.correlate (Together convention)."""
    n = len(h)
    if j == 0:
        # M(0) = L * Σ h_i * (1 - h_i)
        return L * float(np.sum(h * (1 - h)))
    return L * float((h[: n - j] * (1 - h[j:])).sum())


def verify_M_matches():
    """Sanity check: matrix form == correlate form."""
    with open(DATA) as f:
        h = np.array(json.load(f)["together"]["values"])
    n = len(h)
    L = 2.0 / n
    print(f"--- Sanity check M_via_matrix == M_via_correlate ---")
    for j in [0, 33, 100, 200, 500]:
        A = A_t(n, j)
        # M via matrix: L * h @ (I - T_j) @ h
        # = L * (h@h - h@T_j h) = L * (Σh_i² - Σ h_i h_{i+j})
        # But Together convention is M = L * Σ h_i (1 - h_{i+j}) = L*Σh_i - L*Σh_i h_{i+j}
        # The difference: h@h vs Σh_i is L²[norm] vs L¹[mass]
        # So these are DIFFERENT quantities. We need to use the L¹-version.
        M_corr = M_via_correlate(h, L, j)
        M_mat = M_via_matrix(h, A, L)
        print(f"  j={j:4d}: M_correlate = {M_corr:.10f}, M_matrix(L²) = {M_mat:.10f}, diff = {M_mat - M_corr:+.3e}")


def spectral_sweep():
    """Compute spectrum of A_j across a sweep of shifts."""
    n = 200  # smaller for fast spectral analysis
    print(f"\n--- Spectral sweep at n = {n} ---")
    print(f"{'j':>4} {'λ_max(A_j)':>12} {'λ_min(A_j)':>12} {'σ_max(A_j)':>12}")
    results = []
    for j in [0, 10, 20, 33, 50, 100, 150, 199]:
        A = A_t(n, j)
        eigvals = np.linalg.eigvals(A)
        re_eigvals = eigvals.real
        # Symmetric part
        A_sym = (A + A.T) / 2
        sym_eigvals = np.linalg.eigvalsh(A_sym)
        results.append({
            "j": j,
            "lambda_max_A": float(re_eigvals.max()),
            "lambda_min_A": float(re_eigvals.min()),
            "lambda_max_sym": float(sym_eigvals.max()),
            "lambda_min_sym": float(sym_eigvals.min()),
        })
        print(f"  {j:>4d} {sym_eigvals.max():>12.6f} {sym_eigvals.min():>12.6f} {sym_eigvals.max() - sym_eigvals.min():>12.6f}")

    return results


def compare_spectral_bound_to_M():
    """For Together's h*, compare the spectral bound to actual M.

    The Rayleigh-quotient bound: for h with ‖h‖_2 = h_norm, ⟨h, A_t h⟩ ≤ λ_max(A_t)·h_norm².
    Need to express this in our normalization (L¹ unit mass).
    """
    print(f"\n--- Rayleigh-quotient bound vs Together's actual M ---")
    with open(DATA) as f:
        h = np.array(json.load(f)["together"]["values"])
    n = len(h)
    L = 2.0 / n
    h_l2 = np.sqrt((h ** 2).sum())  # discrete L² norm (with cell weight)
    print(f"  h_l2 (no L weighting) = {h_l2:.6f}")
    print(f"  ‖h‖_2 (with L weighting) = √(L · Σh²) = {np.sqrt(L * (h**2).sum()):.6f}")
    print(f"  ‖h‖_1 = L · Σh = {L * h.sum()}")

    print(f"\n  Actual sup_j M(jL) = {max(M_via_correlate(h, L, j) for j in range(n)):.10f}")

    for j in [33]:
        A = A_t(n, j)
        A_sym = (A + A.T) / 2
        eigs = np.linalg.eigvalsh(A_sym)
        # ⟨h, A_j h⟩ ≤ λ_max · ‖h‖²
        lambda_max = eigs.max()
        ub_via_spectral = L * lambda_max * (h ** 2).sum()
        actual = M_via_correlate(h, L, j)
        print(f"\n  j={j}: λ_max(sym A_j) = {lambda_max:.6f}")
        print(f"    spectral UB on M (Rayleigh):  L · λ_max · ‖h‖² = {ub_via_spectral:.10f}")
        print(f"    actual M(jL):                                    {actual:.10f}")
        print(f"    bound is {ub_via_spectral / actual:.2f}× actual (loose)")


def what_might_work():
    """Outline what an actual useful spectral attack would look like."""
    print(f"\n--- Why naive spectral doesn't give μ directly ---")
    print("""
    The issue: ⟨h, A_j h⟩ = ‖h‖² · λ_max only when h is UNCONSTRAINED unit vector.
    Our h has constraints ⟨h, 1⟩ = 1 and 0 ≤ h ≤ 1, which restrict h to a polytope
    smaller than the unit L² ball. The Rayleigh quotient bound is therefore very loose.

    The min-over-h of sup_t ⟨h, A_t h⟩ requires solving a constrained min-max QP at
    each h-step, which is exactly what cvxpy/CLARABEL does in the SDP. Spectral analysis
    of A_t alone doesn't shortcut this.

    BUT: spectral analysis can be useful for:
    1. Identifying the EIGENVECTOR for the binding shift t (information about optimal h structure)
    2. Building a spectral PRECONDITIONER for the constrained optimization
    3. Lower bound via Weyl-type inequalities (spectral interlacing)

    For now: prototype shows that the Rayleigh bound at the binding shift t=33 is
    ~ {bound}×actual — clearly the constrained problem is much tighter than the
    unconstrained spectrum.
    """)


def main():
    print("=" * 78)
    print("PRO-29: Spectral reformulation prototype")
    print("=" * 78)
    verify_M_matches()
    spectral_sweep()
    compare_spectral_bound_to_M()
    what_might_work()


if __name__ == "__main__":
    main()
