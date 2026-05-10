"""
Phase 2 — Verified evaluator for step-function upper bounds on C_5.

For h piecewise-constant on N intervals of width Δ = 2/N, h_i ∈ [0,1],
sum h_i = N/2 (so ∫h = 1), define for j = 1,...,N-1:

   J^+_j(h) = Δ · sum_{i=1}^{N-j} h_{i+j}(1 - h_i)         [shift t = +jΔ]
   J^-_j(h) = Δ · sum_{i=1}^{N-j} h_i (1 - h_{i+j})        [shift t = −jΔ]

The continuous bound is
   B(h) = max_{j ∈ [1,N-1]} max(J^+_j(h), J^-_j(h))

(both signs needed because the discrete max_k takes k ∈ Z\{0}, both signs).
For h symmetric about x=1 (h_i = h_{N+1-i}) we have J^+_j = J^-_j.

The continuous sup_{t∈(0,2)} J(t) is attained at some t = jΔ because the
restricted autocorrelation t ↦ ∫_0^{2−t} h(x)h(x+t)dx is piecewise-linear in
t between grid points (h piecewise-constant on uniform grid). So checking
j = 1,...,N−1 suffices.
"""
from __future__ import annotations
import numpy as np
from fractions import Fraction
from typing import Sequence


def J_plus(h: np.ndarray, j: int) -> float:
    """Δ · sum_{i=1}^{N-j} h_{i+j}(1 - h_i).  Vectorized."""
    N = len(h)
    return float(2.0 / N) * float(np.dot(h[j:], 1.0 - h[: N - j]))


def J_minus(h: np.ndarray, j: int) -> float:
    """Δ · sum_{i=1}^{N-j} h_i(1 - h_{i+j}).  Vectorized."""
    N = len(h)
    return float(2.0 / N) * float(np.dot(h[: N - j], 1.0 - h[j:]))


def J_two_sided(h: np.ndarray, j: int) -> float:
    return max(J_plus(h, j), J_minus(h, j))


def bound_float(h: np.ndarray) -> tuple[float, int, str]:
    """Returns (max-bound, j*, sign)."""
    N = len(h)
    best = (-1.0, -1, "")
    for j in range(1, N):
        jp, jm = J_plus(h, j), J_minus(h, j)
        if jp > best[0]:
            best = (jp, j, "+")
        if jm > best[0]:
            best = (jm, j, "-")
    return best


def J_plus_exact(h: Sequence[Fraction], j: int) -> Fraction:
    N = len(h)
    Δ = Fraction(2, N)
    return Δ * sum(h[i + j] * (1 - h[i]) for i in range(N - j))


def J_minus_exact(h: Sequence[Fraction], j: int) -> Fraction:
    N = len(h)
    Δ = Fraction(2, N)
    return Δ * sum(h[i] * (1 - h[i + j]) for i in range(N - j))


def bound_exact(h: Sequence[Fraction]) -> tuple[Fraction, int, str]:
    N = len(h)
    best = (Fraction(-1), -1, "")
    for j in range(1, N):
        jp, jm = J_plus_exact(h, j), J_minus_exact(h, j)
        if jp > best[0]:
            best = (jp, j, "+")
        if jm > best[0]:
            best = (jm, j, "-")
    return best


def check_constraints_exact(h: Sequence[Fraction]) -> tuple[bool, str]:
    N = len(h)
    Δ = Fraction(2, N)
    if not all(Fraction(0) <= hi <= Fraction(1) for hi in h):
        return False, "values out of [0,1]"
    if Δ * sum(h) != Fraction(1):
        return False, f"normalization {Δ * sum(h)} ≠ 1"
    return True, "ok"


# --- Sanity checks ---------------------------------------------------------
def _heaviside_h(N: int) -> np.ndarray:
    h = np.zeros(N)
    h[N // 2 :] = 1.0
    return h


def _rectangle_h(N: int) -> np.ndarray:
    assert N % 4 == 0
    h = np.zeros(N)
    h[N // 4 : 3 * N // 4] = 1.0
    return h


if __name__ == "__main__":
    print("Sanity checks:")
    for N in [40, 200, 1000]:
        h = _heaviside_h(N)
        b, j, s = bound_float(h)
        print(
            f"  heaviside N={N}: bound = {b:.6f} at j={j} ({s}); "
            f"discrete expectation: 1.0"
        )
    for N in [40, 200, 1000]:
        h = _rectangle_h(N)
        b, j, s = bound_float(h)
        print(
            f"  rectangle N={N}: bound = {b:.6f} at j={j} ({s}); "
            f"discrete expectation: 0.5"
        )
    # Constant h = 1/2 (∫h = 1):
    for N in [40, 200, 1000]:
        h = 0.5 * np.ones(N)
        b, j, s = bound_float(h)
        print(f"  constant 1/2 N={N}: bound = {b:.6f} at j={j} ({s})")
