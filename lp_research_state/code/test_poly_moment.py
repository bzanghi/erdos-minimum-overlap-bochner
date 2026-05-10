"""
Test: does adding the even-polynomial-moment-nonneg constraint to the
Bochner-augmented LP raise Ω* at row 4?
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem
from poly_moment import build_even_moment_nonneg_constraints

ROWS = {
    4: (0.004, 0.3875, -0.02, 0.02),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, default=4)
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--T", type=int, default=800)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=15)
    ap.add_argument("--k_max", type=int, default=10)
    args = ap.parse_args()
    h, p, qm, qp = ROWS[args.row]

    print(f"=== Test polynomial-moment-nonneg at row {args.row} ===")
    print(f"  N={args.N} T={args.T} R={args.R} bochner_n={args.bochner_n} k_max={args.k_max}\n")

    # Baseline
    t0 = time.time()
    Omega1, w, v, c, d, eps, dlt, cons1 = build_problem(
        args.N, args.T, args.R, h, h, p, p, qm, qp, bochner_n=args.bochner_n)
    prob1 = cp.Problem(cp.Minimize(Omega1), cons1)
    prob1.solve(solver="CLARABEL", verbose=False)
    t_base = time.time() - t0
    print(f"[Baseline]  Ω* = {prob1.value:.10f}   ({t_base:.1f}s, {prob1.status})")

    # With poly-moment constraints
    t0 = time.time()
    Omega2, w2, v2, c2, d2, eps2, dlt2, cons2 = build_problem(
        args.N, args.T, args.R, h, h, p, p, qm, qp, bochner_n=args.bochner_n)
    pm_cons, tb = build_even_moment_nonneg_constraints(c2, d2, args.T, k_max=args.k_max)
    cons2 += pm_cons
    print(f"  Added {len(pm_cons)} poly-moment constraints (k=2,4,...,{args.k_max}).")
    print(f"  Tail bounds: {tb}")
    prob2 = cp.Problem(cp.Minimize(Omega2), cons2)
    prob2.solve(solver="CLARABEL", verbose=False)
    t_aug = time.time() - t0
    print(f"[+PolyMoment] Ω* = {prob2.value:.10f}   ({t_aug:.1f}s, {prob2.status})")

    delta = prob2.value - prob1.value
    print(f"\nΔΩ* = {delta:+.7f}")
    if delta > 1e-5:
        print("✓ Polynomial-moment constraint CUTS — new valid lever discovered")
    elif delta > 1e-7:
        print("~ Marginal cut (numeric noise possible)")
    else:
        print("✗ No cut from polynomial-moment constraint at this LP optimum")

    # Diagnose: which of the added constraints are tight at the augmented optimum?
    c2_s = np.array(c2.value); d2_s = np.array(d2.value)
    from poly_moment import fourier_coeffs_of_xk
    alpha0, alpha, beta = fourier_coeffs_of_xk(args.k_max, args.T)
    print("\n  m_k at augmented optimum vs tail bound:")
    for k in range(2, args.k_max + 1, 2):
        m_k = 0.5 * alpha0[k] + alpha[k, :] @ c2_s + beta[k, :] @ d2_s
        slack = m_k - (-tb[k])  # how far above the lower bound
        print(f"    m_{k} = {m_k:+.6e}   tail = ±{tb[k]:.2e}   slack from constraint = {slack:+.2e}")


if __name__ == "__main__":
    main()
