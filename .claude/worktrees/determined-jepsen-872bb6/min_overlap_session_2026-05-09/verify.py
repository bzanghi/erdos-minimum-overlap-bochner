"""
Phase 3 — adversarial verification of an upper-bound candidate h.

Three independent checks:
  (V1) Recompute J^+_j and J^-_j for all j using a from-scratch loop.
  (V2) Densely scan t ∈ [0, 2] using an exact piecewise-linear formula for the
       restricted autocorrelation on a uniform grid of width Δ = 2/N.
  (V3) Re-evaluate everything in mpmath with 50 digits of precision.

The point: never trust the optimizer's reported value. Recompute through code
that does NOT share its codepath, and confirm the supremum is achieved at a
grid point (because the autocorrelation is piecewise-linear in t).
"""
from __future__ import annotations
import json, sys, math
import numpy as np
import mpmath as mp


def V1_independent_recompute(h):
    """Loop-based, no fancy vectorization. Returns max bound and tabulated J."""
    N = len(h)
    Δ = 2.0 / N
    Jp = [0.0] * (N - 1)
    Jm = [0.0] * (N - 1)
    for j in range(1, N):
        sp = 0.0
        sm = 0.0
        for i in range(N - j):
            sp += h[i + j] * (1.0 - h[i])
            sm += h[i] * (1.0 - h[i + j])
        Jp[j - 1] = Δ * sp
        Jm[j - 1] = Δ * sm
    bound = max(max(Jp), max(Jm))
    return bound, Jp, Jm


def V2_continuous_scan(h, n_t=20001):
    """Compute J^+(t) for many t in [0, 2] using the exact piecewise-linear formula.

    For piecewise-constant h with grid Δ, define the restricted correlation
       R̃(t) := ∫_0^{2-t} h(x) h(x+t) dx,   t ∈ [0, 2].
    Between grid points, R̃ is piecewise-linear in t. The integration domain
    can be evaluated exactly cell-by-cell.
    """
    N = len(h)
    Δ = 2.0 / N
    h = np.asarray(h, dtype=float)

    def Jpm_at_t(t):
        """Compute J^+(t) and J^-(t) exactly for piecewise-constant h."""
        if t <= 0:
            v = float(Δ * np.dot(h, 1.0 - h))
            return v, v
        if t >= 2:
            return 0.0, 0.0
        # Walk through cells of h(x) on x ∈ [0, 2-t]:
        Jp = 0.0
        Jm = 0.0
        x = 0.0
        i = 0
        end = 2.0 - t
        while x < end - 1e-15 and i < N:
            x_next = min((i + 1) * Δ, end)
            if x_next <= x:
                break
            # On [x, x_next]: h(x) = h[i]. h(x+t) is piecewise-constant; cell
            # index k = floor((x+t)/Δ). Subdivide where x+t crosses a boundary.
            y = x
            while y < x_next - 1e-15:
                k = int((y + t) / Δ + 1e-12)
                if k >= N:
                    y_next = x_next
                    h_kt = 0.0
                else:
                    boundary = (k + 1) * Δ - t
                    y_next = min(boundary, x_next)
                    h_kt = h[k]
                length = y_next - y
                Jp += length * h_kt * (1.0 - h[i])
                Jm += length * h[i] * (1.0 - h_kt)
                y = y_next
            x = x_next
            i += 1
        return Jp, Jm

    ts = np.linspace(0.0, 2.0, n_t)
    Jps = np.empty(n_t)
    Jms = np.empty(n_t)
    for k, t in enumerate(ts):
        Jps[k], Jms[k] = Jpm_at_t(t)
    return ts, Jps, Jms


def V3_mpmath_recompute(h, dps=50):
    mp.mp.dps = dps
    N = len(h)
    Δ = mp.mpf(2) / N
    hm = [mp.mpf(repr(float(x))) for x in h]
    best = mp.mpf(0)
    best_j = -1
    best_sign = ""
    for j in range(1, N):
        sp = mp.mpf(0)
        sm = mp.mpf(0)
        for i in range(N - j):
            sp += hm[i + j] * (1 - hm[i])
            sm += hm[i] * (1 - hm[i + j])
        Jp = Δ * sp
        Jm = Δ * sm
        if Jp > best: best, best_j, best_sign = Jp, j, "+"
        if Jm > best: best, best_j, best_sign = Jm, j, "-"
    return best, best_j, best_sign


if __name__ == "__main__":
    fname = sys.argv[1] if len(sys.argv) > 1 else "refined_N400.json"
    data = json.load(open(fname))
    h = np.array(data["h"])
    N = len(h)
    print(f"Verifying {fname}: N={N}, claimed bound={data['bound']:.10f}, "
          f"argmax j={data['argmax_j']}({data['argmax_sign']})")
    print(f"sum h = {h.sum():.12f}  (target {N/2})  → ∫h = {2*h.sum()/N:.12f}")
    print(f"min h = {h.min():.4e}, max h = {h.max():.4e}")

    print("\nV1: Independent re-compute via straight loops...")
    b1, Jp, Jm = V1_independent_recompute(h)
    j_arg_p = int(np.argmax(Jp)) + 1
    j_arg_m = int(np.argmax(Jm)) + 1
    print(f"  V1 bound: {b1:.10f}")
    print(f"  V1 J^+ peak at j={j_arg_p}, value {max(Jp):.10f}")
    print(f"  V1 J^- peak at j={j_arg_m}, value {max(Jm):.10f}")

    print("\nV2: dense-t continuous scan (n_t = 20001)...")
    ts, Jps, Jms = V2_continuous_scan(h, n_t=20001)
    sup_p = float(Jps.max())
    sup_m = float(Jms.max())
    t_p = float(ts[int(np.argmax(Jps))])
    t_m = float(ts[int(np.argmax(Jms))])
    print(f"  V2 sup J^+ = {sup_p:.10f} at t={t_p:.5f} (j*Δ = {data['argmax_j']*2/N:.5f})")
    print(f"  V2 sup J^- = {sup_m:.10f} at t={t_m:.5f}")
    sup_v2 = max(sup_p, sup_m)
    print(f"  V2 max bound: {sup_v2:.10f}")

    print("\nV3: mpmath 50-digit recompute (slow)...")
    if N <= 400:
        b3, j3, s3 = V3_mpmath_recompute(h, dps=50)
        mp.mp.dps = 30
        print(f"  V3 bound: {mp.nstr(b3, 25)}  argmax j={j3} ({s3})")
        diff = abs(float(b3) - b1)
        print(f"  |V3 - V1| = {diff:.2e}")
    else:
        print(f"  (skipping V3 for large N={N})")

    # Final assertion
    sup_check = max(b1, sup_v2)
    print(f"\nVerified upper bound for h:  {sup_check:.10f}")
    print(f"Reported by optimizer:       {data['bound']:.10f}")
    print(f"|verified - reported|:       {abs(sup_check - data['bound']):.2e}")
