"""
Test: does adding the even-only Hankel-PSD constraint to the Bochner-augmented
LP cut more than the scalar m_k ≥ -tail constraints?
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import cvxpy as cp
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem
from poly_moment import (build_even_moment_nonneg_constraints,
                         build_even_hankel_psd, fourier_coeffs_of_xk,
                         even_moment_tail_bound)

ROWS = {4: (0.004, 0.3875, -0.02, 0.02)}


def solve(N, T, R, h, p, q1, q2, bochner_n, kind, *, pm_k_max=14, hankel_n=4):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2, bochner_n=bochner_n)
    if kind in ("pm", "both"):
        pmc, _ = build_even_moment_nonneg_constraints(c, d, T, k_max=pm_k_max)
        cons += pmc
    if kind in ("hankel", "both"):
        hcons, m_var, tails = build_even_hankel_psd(c, d, T, n_hankel=hankel_n)
        cons += hcons
        print(f"  Hankel n={hankel_n}: tight slack encoding ({len(hcons)} constraints, "
              f"tails: {[f'{t:.1e}' for t in tails[:5]]}...)")
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time(); prob.solve(solver="CLARABEL", verbose=False); dt = time.time() - t0
    return {"value": float(prob.value), "status": prob.status, "time": dt,
            "c": np.array(c.value), "d": np.array(d.value)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, default=4)
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=30)
    ap.add_argument("--pm_k_max", type=int, default=14)
    ap.add_argument("--hankel_n", type=int, default=4)
    args = ap.parse_args()
    h, p, q1, q2 = ROWS[args.row]

    print(f"=== Hankel-PSD test at row {args.row} ===")
    print(f"  N={args.N} T={args.T} bochner_n={args.bochner_n} "
          f"pm_k_max={args.pm_k_max} hankel_n={args.hankel_n}\n")

    print("[1] poly_moment scalars only (baseline)")
    r1 = solve(args.N, args.T, args.R, h, p, q1, q2, args.bochner_n, "pm",
               pm_k_max=args.pm_k_max)
    print(f"    Ω* = {r1['value']:.10f}  ({r1['status']}, {r1['time']:.1f}s)\n")

    print(f"[2] + even-Hankel PSD (level n={args.hankel_n}, matrix {args.hankel_n+1}×{args.hankel_n+1}, uses m_0..m_{4*args.hankel_n})")
    r2 = solve(args.N, args.T, args.R, h, p, q1, q2, args.bochner_n, "both",
               pm_k_max=args.pm_k_max, hankel_n=args.hankel_n)
    print(f"    Ω* = {r2['value']:.10f}  ({r2['status']}, {r2['time']:.1f}s)\n")

    delta = r2['value'] - r1['value']
    print(f"=== ΔΩ* from adding Hankel-PSD = {delta:+.7f} ===")
    if delta > 1e-5: print("✓ Hankel-PSD CUTS even when scalar constraints already present")
    elif delta > 1e-7: print("~ Marginal cut")
    else: print("✗ Hankel-PSD does not add cutting power beyond scalar m_k ≥ -tail")

    # Show eigenvalue of Hankel at baseline (no Hankel constraint) optimum
    c_star = r1['c']; d_star = r1['d']
    alpha0, alpha, beta = fourier_coeffs_of_xk(4 * args.hankel_n, args.T)
    m = np.array([0.5 * alpha0[k] + alpha[k, :] @ c_star + beta[k, :] @ d_star
                  for k in range(4 * args.hankel_n + 1)])
    H = np.array([[m[2 * (i + j)] for j in range(args.hankel_n + 1)]
                   for i in range(args.hankel_n + 1)])
    eigs = np.linalg.eigvalsh(H)
    print(f"\n  Even-Hankel eigs at baseline optimum: min = {eigs[0]:.4e}  max = {eigs[-1]:.4e}")
    print(f"  Diagnostic: if min eig < -pert, Hankel could cut. (pert from constraint construction.)")


if __name__ == "__main__":
    main()
