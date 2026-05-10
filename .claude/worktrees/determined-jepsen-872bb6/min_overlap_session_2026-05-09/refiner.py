"""
Phase 2 — Local refinement of step-function upper bounds on C_5.

Solve  min_h  max_{j ∈ [1, N-1], σ ∈ {+,-}}  J^σ_j(h)
       s.t.   sum h_i = N/2,    0 ≤ h_i ≤ 1.

Smooth max with logsumexp, project, L-BFGS-B with box bounds + Lagrangian
for the sum constraint, β-annealing, multistart.

Optionally restrict to symmetric h (h_i = h_{N+1-i}) — under symmetry J^+=J^-.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from evaluator import bound_float, J_plus, J_minus


def J_array(h: np.ndarray) -> np.ndarray:
    """Returns the 2(N-1) values [J^+_1, ..., J^+_{N-1}, J^-_1, ..., J^-_{N-1}]."""
    N = len(h)
    Δ = 2.0 / N
    out = np.empty(2 * (N - 1))
    for j in range(1, N):
        out[j - 1] = Δ * float(np.dot(h[j:], 1.0 - h[: N - j]))  # J^+
        out[N - 1 + j - 1] = Δ * float(np.dot(h[: N - j], 1.0 - h[j:]))  # J^-
    return out


def smooth_max_grad(h: np.ndarray, β: float):
    """Returns (smoothed_max(J), gradient w.r.t. h)."""
    N = len(h)
    Δ = 2.0 / N
    Js = np.empty(2 * (N - 1))
    grads = np.zeros((2 * (N - 1), N))
    for j in range(1, N):
        a = h[: N - j]      # h_i, i = 1..N-j
        b = h[j:]           # h_{i+j}
        # J^+: sum b_i (1 - a_i)
        Js[j - 1] = Δ * float(np.dot(b, 1.0 - a))
        grads[j - 1, j:] += Δ * (1.0 - a)
        grads[j - 1, : N - j] += -Δ * b
        # J^-: sum a_i (1 - b_i)
        Js[N - 1 + j - 1] = Δ * float(np.dot(a, 1.0 - b))
        grads[N - 1 + j - 1, : N - j] += Δ * (1.0 - b)
        grads[N - 1 + j - 1, j:] += -Δ * a
    m = float(Js.max())
    w = np.exp(β * (Js - m))
    Z = float(w.sum())
    sm = m + np.log(Z) / β
    g = (w @ grads) / Z
    return sm, g


def project_to_constraints(h: np.ndarray) -> np.ndarray:
    """Project onto {h: 0 ≤ h_i ≤ 1, sum h_i = N/2}."""
    N = len(h)
    target = N / 2.0
    lo, hi = -2.0, 2.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        s = np.clip(h - mid, 0.0, 1.0).sum()
        if s > target:
            lo = mid
        else:
            hi = mid
    return np.clip(h - 0.5 * (lo + hi), 0.0, 1.0)


def penalty_objective(h_flat: np.ndarray, N: int, β: float, λ: float):
    sm, g = smooth_max_grad(h_flat, β)
    sum_h = float(h_flat.sum())
    target = N / 2.0
    pen = 0.5 * λ * (sum_h - target) ** 2
    pen_g = λ * (sum_h - target) * np.ones(N)
    return sm + pen, g + pen_g


def refine(
    h0: np.ndarray,
    β_schedule=(20.0, 80.0, 300.0, 1500.0, 8000.0),
    λ: float = 200.0,
    iters: int = 400,
    symmetric: bool = False,
    verbose: bool = False,
) -> np.ndarray:
    h = h0.copy()
    N = len(h)
    if symmetric:
        h = 0.5 * (h + h[::-1])
    h = project_to_constraints(h)
    for β in β_schedule:
        if symmetric:
            # parameterize by first half v ∈ R^{N/2} (or R^{(N+1)/2})
            half = (N + 1) // 2

            def to_full(v):
                f = np.empty(N)
                f[: len(v)] = v
                # mirror
                f[N - len(v) :] = v[::-1]
                if N % 2 == 0:
                    pass  # already filled
                return f

            def obj(v):
                f = to_full(v)
                val, gf = penalty_objective(f, N, β, λ)
                gv = gf[: half] + gf[N - half :][::-1]
                if N % 2 == 1:
                    gv[half - 1] /= 2  # middle counted twice
                return val, gv

            v0 = h[:half].copy()
            res = minimize(
                obj,
                v0,
                jac=True,
                method="L-BFGS-B",
                bounds=[(0.0, 1.0)] * half,
                options={"maxiter": iters, "ftol": 1e-13, "gtol": 1e-11},
            )
            h = to_full(res.x)
        else:
            res = minimize(
                penalty_objective,
                h,
                args=(N, β, λ),
                jac=True,
                method="L-BFGS-B",
                bounds=[(0.0, 1.0)] * N,
                options={"maxiter": iters, "ftol": 1e-13, "gtol": 1e-11},
            )
            h = res.x
        h = project_to_constraints(h)
        if verbose:
            b, j, s = bound_float(h)
            print(f"   β={β:>7.1f}  bound={b:.6f}  argmax j={j} ({s})")
    return h


def random_starts(N: int, n_starts: int, seed: int = 0) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    starts: list[np.ndarray] = []
    starts.append(0.5 * np.ones(N))
    xs_mid = np.linspace(1.0 / N, 2 - 1.0 / N, N)
    triangle = np.minimum(xs_mid, 2.0 - xs_mid)
    starts.append(project_to_constraints(triangle * (N / 2.0) / triangle.sum()))
    rect = np.zeros(N)
    rect[N // 4 : 3 * N // 4] = 1.0
    starts.append(rect.copy())
    # Trapezoid: 0 on outer ramp, c on plateau, smooth.
    h = np.where((xs_mid > 0.5) & (xs_mid < 1.5), 0.7, 0.3)
    starts.append(project_to_constraints(h))
    # Cosine bump
    h = 0.5 * (1 - np.cos(np.pi * xs_mid))
    starts.append(project_to_constraints(h))
    # Random
    for _ in range(n_starts):
        h = 0.5 + 0.4 * rng.standard_normal(N)
        starts.append(project_to_constraints(h))
    return starts


if __name__ == "__main__":
    N = 80
    starts = random_starts(N, n_starts=10, seed=2)
    best = (1e9, None)
    for k, h0 in enumerate(starts):
        h = refine(h0, symmetric=True)
        b, j, s = bound_float(h)
        if b < best[0]:
            best = (b, h.copy())
        print(f"start {k:2d}: bound = {b:.6f}  argmax j={j} ({s})")
    print(f"\nbest (symmetric, N={N}): {best[0]:.6f}")
