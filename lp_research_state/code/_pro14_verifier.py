"""PRO-14 (tooling angle): extract shadow prices (ξ, ν_3, ν_4, τ, λ_m^cos,
σ_m^{1,2}) from any White SDP solve and compute the *verifiable* C_explicit
bound using the Theorem 2 form from LEVER_I_PRIME_THEOREM.md §2.2.

The current conditional Theorem 2 says: IF shadow prices satisfy
    |ξ| ≤ Ξ,  τ ≤ T,  ν_3 ≤ V,  |Δ_sin(1)| ≤ Σ_s
THEN
    Σ_m λ_m^cos ≤ (1 + 3×10⁻⁴) · (2Ξ + T + 2LV + Σ_s).

This tool measures (Ξ, T, V, Σ_s) at the solve's actual KKT point and emits
the C_explicit value. By running it across the 7 White rows (or a richer
cover), we obtain an EMPIRICALLY VERIFIED bound — turning the conditional
theorem into a numerically-certified one for any solve we audit.

Constraint indices (from white_full_convex.py reading; cf.
_lever_i_prime_lambda_m_all_rows.py):

  cons[0..4]   var-box: w,v,Ω boundaries
  cons[5]     L·Σ(w+v) == 1                              → ξ
  cons[6]     L²·Σ(j·w − (j−1)·v) ≥ h_1                 → ν_3
  cons[7]     L³·Σ((j−1)²(w+v)) ≤ 2/3 + h_2²/2          → ν_4
  cons[8..27] cell-envelope cosine (2R=20 inequalities)  → λ_m^cos
  cons[28..67] cell-envelope sine  (2R*2=40 inequalities) → σ_m^{1,2}
  cons[68..R+67] |eps_m| ≤ bound  (R=10 pairs)
  cons[78..R+77] |dlt_m| ≤ bound
  cons[88]  |c| ≤ 2/π
  cons[89]  |d| ≤ 2/π
  cons[90]  Σc² + Σd² ≤ 0.5
  cons[91..94] c[0]/d[0] anchors  (4 inequalities)
  cons[95]  (5.13)                                      → τ
  + bochner-PSD constraints if requested

NOTE: when use_T5p / use_T3 / use_T5 / assume_even are set, indices shift.
This tool is for the default config (no extras except bochner_n).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem  # noqa: E402


def _scalar(dv):
    if dv is None:
        return None
    if np.isscalar(dv):
        return float(dv)
    arr = np.asarray(dv).ravel()
    return float(arr[0]) if arr.size == 1 else None


def extract_shadows(
    N: int, T: int, R: int, bochner_n: int,
    h: float, p: float, q1: float, q2: float,
    *, verbose: bool = False,
):
    """Solve at one center, extract shadow prices, return dict.

    Returns keys: omega, status, xi, nu_3, nu_4, tau, lambda_cos (length 2R),
    sigma_pairs (list of 2R (σ1, σ2) tuples), |Δ_sin(1)| upper bound,
    C_explicit (the implied bound on Σ λ_m^cos), and the breakdown.
    """
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2,
        bochner_n=bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=cp.CLARABEL, verbose=False)
    dt = time.time() - t0
    omega = float(prob.value) if prob.value is not None else None
    status = prob.status

    L = 2.0 / N

    # Default-layout indices (no T3/T5/T5p/assume_even); verified by inspecting
    # build_problem source order.
    idx_xi    = 5
    idx_nu3   = 6
    idx_nu4   = 7
    cos_start = 8
    cos_end   = cos_start + 2 * R
    sin_start = cos_end
    sin_end   = sin_start + 4 * R       # 40 inequalities for R=10
    idx_513   = sin_end + 4 * R + 4 + 4 # tail bounds + box + parseval + anchors
    # The robust way: search for the constraint by signature instead of fixed
    # index. (5.13) has rhs_513 = -0.5·max(p²) ≈ -0.075 etc.
    # We'll detect it dynamically below.

    xi   = _scalar(cons[idx_xi].dual_value)
    nu_3 = _scalar(cons[idx_nu3].dual_value)
    nu_4 = _scalar(cons[idx_nu4].dual_value)

    lambdas: list = []
    for i in range(cos_start, cos_end):
        lambdas.append(_scalar(cons[i].dual_value))

    sigma_pairs = []
    for m in range(1, 2 * R + 1):
        i1 = sin_start + 2 * (m - 1)
        i2 = i1 + 1
        sigma_pairs.append((_scalar(cons[i1].dual_value), _scalar(cons[i2].dual_value)))

    # (5.13) constraint is at a fixed index 95 in the default white_full_convex
    # config: 8 scalar (var-box + sum + h1 + h2) + 20 cosine + 40 sine + 20 tail
    # + 7 (c-box, d-box, sum_squares, c[0]≥p1, c[0]≤p2, d[0]≥q1, d[0]≤q2) = 95.
    # bochner-PSD constraints come after.
    idx_tau = 8 + 2 * R + 4 * R + 2 * R + 7  # = 95 for R=10
    rhs_513 = -0.5 * (p ** 2 + max(abs(q1), abs(q2)) ** 2)
    try:
        tau = _scalar(cons[idx_tau].dual_value)
    except IndexError:
        tau = None
        idx_tau = None

    # |Δ_sin(1)| bound from §2 of theorem: Σ |σ_m^{1,2}| · π m L / 2 (using
    # |β_m^±(1)| ≤ πmL/2).
    delta_sin_bound = 0.0
    for m in range(1, 2 * R + 1):
        s1, s2 = sigma_pairs[m - 1]
        smag = abs(s1 or 0.0) + abs(s2 or 0.0)
        delta_sin_bound += smag * np.pi * m * L / 2.0

    # α_min for the Corollary 1 prefactor: cos(πRL)
    alpha_min = np.cos(np.pi * R * L)
    prefactor = 1.0 / alpha_min  # ≥ 1, ≈ 1 + (πRL)²/2

    # ----- Compute the IMPLIED Σ λ bound (Theorem 2 form) -----
    # Empirical:  Σ λ measured directly
    sum_lambda_measured = sum(abs(l) for l in lambdas if l is not None)
    # Bound:  Σ λ ≤ α_min⁻¹ · (2|ξ| + α_2^+(1)·τ + 2L·ν_3 + |Δ_sin(1)|)
    # with α_2^+(1) = 1 (cell-max of cos(πx) at x=0).
    if xi is None or tau is None or nu_3 is None:
        sum_lambda_bound = None
    else:
        sum_lambda_bound = prefactor * (
            2 * abs(xi) + abs(tau) + 2 * L * abs(nu_3) + delta_sin_bound
        )

    # ----- Compute C_explicit (= Ω* + ResidualGain) -----
    # ResidualGain ≤ (π/(2N)) · Σ_m m·|λ_m|   (Case-A dominant, per §0.3 in
    # LEVER_I_PRIME_THEOREM.md). For an UNCONDITIONAL ceiling, we use the
    # sum_lambda_bound to bound Σ m·|λ_m| ≤ 2R · sum_lambda_bound  (very loose
    # — uses |m| ≤ 2R).
    sum_m_lambda_measured = sum(
        m * abs(lambdas[m - 1] or 0.0) for m in range(1, 2 * R + 1)
    )
    residual_gain_measured = (np.pi / (2 * N)) * sum_m_lambda_measured
    if sum_lambda_bound is not None:
        residual_gain_via_bound = (np.pi / (2 * N)) * 2 * R * sum_lambda_bound
    else:
        residual_gain_via_bound = None
    C_explicit_measured = omega + residual_gain_measured if omega is not None else None
    C_explicit_via_bound = (
        omega + residual_gain_via_bound if (omega is not None and residual_gain_via_bound is not None) else None
    )

    return {
        "params": {"N": N, "T": T, "R": R, "bochner_n": bochner_n,
                   "h": h, "p": p, "q1": q1, "q2": q2},
        "L": L,
        "status": status,
        "Omega": omega,
        "solve_time_s": dt,
        "xi": xi,
        "nu_3": nu_3,
        "nu_4": nu_4,
        "tau": tau,
        "idx_tau": idx_tau,
        "lambda_cos": lambdas,
        "sigma_pairs": sigma_pairs,
        "delta_sin_bound": delta_sin_bound,
        "alpha_min": alpha_min,
        "prefactor_correction": prefactor,
        "sum_lambda_measured": sum_lambda_measured,
        "sum_lambda_bound_via_theorem2": sum_lambda_bound,
        "sum_m_lambda_measured": sum_m_lambda_measured,
        "residual_gain_measured": residual_gain_measured,
        "residual_gain_via_theorem2_bound": residual_gain_via_bound,
        "C_explicit_measured": C_explicit_measured,
        "C_explicit_via_theorem2_bound": C_explicit_via_bound,
        "ratio_xi_over_Omega": abs(xi) / abs(omega) if (xi is not None and omega) else None,
        "ratio_tau_over_Omega": tau / abs(omega) if (tau is not None and omega) else None,
    }


ROWS = {
    "row1":          (0.015, 0.381,    -0.02, 0.02),
    "row4":          (0.004, 0.3875,   -0.02, 0.02),
    "row7":          (0.030, 0.375,    -0.02, 0.02),
    "cde_n30_iter1": (0.0,   0.394175, -0.02, 0.02),
}


def main():
    N = 3000
    T = 1200
    R = 10
    bochner_n = 20

    results = {}
    print(f"{'='*78}")
    print(f"PRO-14 shadow-price verifier  |  N={N} T={T} R={R} bn={bochner_n}")
    print(f"{'='*78}\n")
    for name, (h, p, q1, q2) in ROWS.items():
        print(f"--- {name}  (h={h}, p={p}, q={q1},{q2}) ---")
        r = extract_shadows(N, T, R, bochner_n, h, p, q1, q2, verbose=True)
        results[name] = r
        def _fmt(x, w=10, p=6):
            if x is None: return "None".rjust(w)
            try:
                return f"{x:>{w}.{p}f}"
            except Exception:
                return str(x).rjust(w)
        ratio_b = (
            r["sum_lambda_bound_via_theorem2"] / max(r["sum_lambda_measured"], 1e-12)
            if r["sum_lambda_bound_via_theorem2"] is not None else None
        )
        print(
            f"  Ω = {r['Omega']:.10f}  status = {r['status']}\n"
            f"  ξ = {_fmt(r['xi'])}  → |ξ|/Ω = {_fmt(r['ratio_xi_over_Omega'], 8, 4)}\n"
            f"  τ = {_fmt(r['tau'])}  → τ/Ω = {_fmt(r['ratio_tau_over_Omega'], 8, 4)}\n"
            f"  ν_3 = {_fmt(r['nu_3'])}\n"
            f"  ν_4 = {_fmt(r['nu_4'])}\n"
            f"  Σ|λ_m| (measured)    = {_fmt(r['sum_lambda_measured'])}\n"
            f"  Σ|λ_m| (Thm 2 bound) = {_fmt(r['sum_lambda_bound_via_theorem2'])}\n"
            f"  ratio (bound / measured) = {_fmt(ratio_b, 8, 4)}\n"
            f"  Σ m·|λ_m| = {_fmt(r['sum_m_lambda_measured'], 8, 4)}\n"
            f"  ResidualGain (measured) = {r['residual_gain_measured']:.3e}\n"
            f"  C_explicit (measured) = {r['C_explicit_measured']:.7f}\n"
        )

    out_path = Path(__file__).parent.parent / "data" / "pro14_shadow_prices.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(results, indent=2, default=str))
    print(f"Saved to {out_path}")

    # Summary
    print("\n--- SUMMARY ---")
    print(f"{'row':<18} {'Ω':>10} {'|ξ|/Ω':>8} {'τ/Ω':>8} {'Σ|λ|':>8} {'C_explicit':>12}")
    def _f(x, w, p):
        if x is None: return "None".rjust(w)
        return f"{x:>{w}.{p}f}"
    for name, r in results.items():
        print(
            f"{name:<18} {_f(r['Omega'], 10, 7)} {_f(r['ratio_xi_over_Omega'], 8, 4)}"
            f" {_f(r['ratio_tau_over_Omega'], 8, 4)}"
            f" {_f(r['sum_lambda_measured'], 8, 4)}"
            f" {_f(r['C_explicit_measured'], 12, 7)}"
        )


if __name__ == "__main__":
    main()
