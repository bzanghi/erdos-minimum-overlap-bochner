"""
Constraint Discovery Engine — Phase 1 probe.

Solve the Bochner-augmented White LP at a chosen row, reconstruct f̃ from the
LP optimum, and run a violation panel: which physical f-properties does the
reconstructed f̃ violate? Each violation is a hint that a corresponding convex
constraint would cut.

Usage:
    python3 probe.py [--row 4] [--N 2000] [--T 800] [--R 10] [--n 15]
"""
from __future__ import annotations
import argparse
import json
import sys
import time
from pathlib import Path

import cvxpy as cp
import numpy as np

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

from white_full_convex import build_problem  # noqa: E402
from dual_extractor import solve_with_dual_extraction  # noqa: E402

ROWS = {
    1: (0.015, 0.381,  -0.02, 0.02),
    2: (0.015, 0.385,  -0.02, 0.02),
    3: (0.020, 0.375,  -0.02, 0.02),
    4: (0.004, 0.3875, -0.02, 0.02),
    5: (0.000, 0.4,    -0.02, 0.02),
    6: (0.000, 0.381,  -0.02, 0.02),
    7: (0.030, 0.375,  -0.02, 0.02),
}


def reconstruct_f(c_star: np.ndarray, d_star: np.ndarray, x: np.ndarray) -> np.ndarray:
    """f̃(x) = 1/2 + Σ_{k=1..T} c_k cos(πkx) + d_k sin(πkx)."""
    T = len(c_star)
    k = np.arange(1, T + 1)
    # x: shape (M,)  →  πkx: shape (T, M)
    arg = np.pi * np.outer(k, x)
    return 0.5 + c_star @ np.cos(arg) + d_star @ np.sin(arg)


def reconstruct_M(c_star: np.ndarray, d_star: np.ndarray, x: np.ndarray) -> np.ndarray:
    """M̃(x) = (f̃ * g̃)(x) where g̃ = 1 − f̃, both viewed as 2-periodic.

    By Plancherel, the autocorrelation-style Fourier coeffs of M relate to
    |f̂|² minus a Fourier-pair structure. For diagnostic purposes we just
    convolve numerically on a dense grid.
    """
    # 2-periodic extension: f̃ is supported on [-1, 1].
    # Use a long dense grid and FFT-convolve for fast eval.
    M = len(x)
    f = reconstruct_f(c_star, d_star, x)
    g = 1.0 - f
    # Cyclic convolution scaled to interval length 2.
    dx = x[1] - x[0]
    Mtilde = np.real(np.fft.ifft(np.fft.fft(f) * np.fft.fft(g))) * dx
    return Mtilde


def bochner_min_eig(c_star: np.ndarray, d_star: np.ndarray, n: int, sign: int = +1) -> float:
    """Min eigenvalue of M_n(f̃) (sign=+1) or M_n(1−f̃) (sign=−1).

    f̂(0) = 1/2; f̂(k) = (c_k − i d_k)/2 for k ≥ 1. (1−f)̂ flips sign on k ≠ 0.
    Build (n+1)×(n+1) Hermitian Toeplitz; return min real eigenvalue.
    """
    T = len(c_star)
    n = min(n, T)
    # fhat[k] for k in 0..n
    fhat = np.empty(n + 1, dtype=complex)
    fhat[0] = 0.5 if sign == +1 else 0.5  # for (1-f), f̂(0) flips to 1/2 too (1 − 1/2)
    for k in range(1, n + 1):
        fhat[k] = sign * (c_star[k - 1] - 1j * d_star[k - 1]) / 2.0
    # Toeplitz: T_{jk} = fhat[j-k] for j-k≥0, conj for j-k<0
    T_mat = np.empty((n + 1, n + 1), dtype=complex)
    for j in range(n + 1):
        for kk in range(n + 1):
            ell = j - kk
            T_mat[j, kk] = fhat[ell] if ell >= 0 else np.conj(fhat[-ell])
    eigs = np.linalg.eigvalsh(T_mat)
    return float(eigs[0])


def f_integral_squares(c_star: np.ndarray, d_star: np.ndarray) -> tuple[float, float]:
    """Exact ∫_{-1}^{1} f̃² dx using Parseval, vs. ∫ f̃ = 1.

    f̃(x) = 1/2 + Σ c_k cos(πkx) + d_k sin(πkx). On [-1, 1] of length 2:
    ∫ 1·1 = 2,   ∫ cos² = 1,   ∫ sin² = 1,   cross-terms zero.
    So ∫ f̃² = (1/2)² · 2 + Σ (c_k² + d_k²) = 1/2 + Σ (c_k² + d_k²).
    ∫ f̃ = (1/2) · 2 = 1.
    Note: this is the formula even if f̃ goes outside [0,1] — it's a moment
    statement, not a pointwise one.
    """
    norm_sq = 0.5 + float(np.sum(c_star ** 2) + np.sum(d_star ** 2))
    norm_1 = 1.0
    return norm_sq, norm_1


def probe_row(row: int, N: int, T: int, R: int, bochner_n: int,
              extra_args: dict | None = None) -> dict:
    h, p, qm, qp = ROWS[row]
    extra = extra_args or {}
    print(f"# row {row}  (h, p, q) = ({h}, {p}, [{qm}, {qp}])  "
          f"N={N} T={T} R={R} bochner_n={bochner_n} extra={extra}")
    t0 = time.time()
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, qm, qp,
        bochner_n=bochner_n, **extra,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        res = solve_with_dual_extraction(prob)
    except Exception as e:
        print(f"  dual extraction failed ({type(e).__name__}: {e}); "
              "falling back to non-verbose solve")
        val = prob.solve(solver="CLARABEL")
        res = {"reported_value": val, "status": prob.status, "rigorous_dual_LB": val}
    elapsed = time.time() - t0
    print(f"  solved in {elapsed:.1f}s  status={res['status']}")
    c_star = np.array(c.value)
    d_star = np.array(d.value)
    print(f"  Ω* = {res['reported_value']:.10g}  (rigorous_dual_LB = "
          f"{res.get('rigorous_dual_LB')})")

    # Diagnostic 1: pointwise f̃ extremes on a dense grid
    xs = np.linspace(-1, 1, 8001)
    f_vals = reconstruct_f(c_star, d_star, xs)
    f_min = float(f_vals.min()); f_max = float(f_vals.max())
    arg_f_min = float(xs[np.argmin(f_vals)]); arg_f_max = float(xs[np.argmax(f_vals)])

    # Diagnostic 2: ∫ f̃² (Parseval) vs ∫ f̃
    Lsq, L1 = f_integral_squares(c_star, d_star)

    # Diagnostic 3: M_n(f̃) eigenvalues for n past bochner_n
    bochner_eigs = {}
    for n_probe in [bochner_n, bochner_n + 5, bochner_n + 10, 2 * bochner_n, T - 1]:
        if n_probe <= 0 or n_probe > T - 1: continue
        bochner_eigs[n_probe] = {
            "f":   bochner_min_eig(c_star, d_star, n_probe, sign=+1),
            "1-f": bochner_min_eig(c_star, d_star, n_probe, sign=-1),
        }

    # Diagnostic 4: pointwise M̃ extremes
    Mtilde = reconstruct_M(c_star, d_star, xs)
    M_min = float(Mtilde.min()); M_max = float(Mtilde.max())

    # Diagnostic 5: ∫ f̃(1-f̃) dx  (this is ∫ M̃(0) ·   ...wait, this is just
    # a single number. ∫ f(1-f) = ∫ f - ∫ f² = L1 - Lsq.)
    f_one_minus_f = L1 - Lsq

    print()
    print("  === violation panel ===")
    print(f"  pointwise f̃ range:        [{f_min:.6f}, {f_max:.6f}]  "
          f"(at x={arg_f_min:.4f}, x={arg_f_max:.4f})")
    print(f"    violation f̃ < 0:        {max(0, -f_min):.6f}")
    print(f"    violation f̃ > 1:        {max(0, f_max - 1):.6f}")
    print(f"  ∫f̃² = {Lsq:.6f}    ∫f̃ = {L1:.6f}    ∫f̃(1-f̃) = {f_one_minus_f:.6f}")
    print(f"    violation ∫f̃² > ∫f̃:    {max(0, Lsq - L1):.6f}")
    print(f"  pointwise M̃ range:        [{M_min:.6f}, {M_max:.6f}]")
    print(f"    violation M̃ < 0:        {max(0, -M_min):.6f}")
    print("  M_n(f̃) min eigenvalues past bochner_n level:")
    for n_p, ev in bochner_eigs.items():
        flag_f = "⚠️" if ev["f"] < -1e-8 else "  "
        flag_g = "⚠️" if ev["1-f"] < -1e-8 else "  "
        print(f"    n={n_p:4d}:  M_n(f̃) min eig = {ev['f']:+.6e} {flag_f}  "
              f"M_n(1-f̃) min eig = {ev['1-f']:+.6e} {flag_g}")

    return {
        "row": row,
        "params": {"N": N, "T": T, "R": R, "bochner_n": bochner_n, **extra},
        "Omega_star": res["reported_value"],
        "rigorous_dual_LB": res.get("rigorous_dual_LB"),
        "status": res["status"],
        "f_min": f_min, "f_max": f_max,
        "arg_f_min": arg_f_min, "arg_f_max": arg_f_max,
        "integral_f": L1, "integral_f2": Lsq,
        "integral_f_one_minus_f": f_one_minus_f,
        "M_min": M_min, "M_max": M_max,
        "bochner_eigs": bochner_eigs,
        "elapsed_s": elapsed,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--row", type=int, default=4)
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--T", type=int, default=800)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--n", type=int, default=15, help="bochner_n level")
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    result = probe_row(args.row, args.N, args.T, args.R, args.n)
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2, default=float))
        print(f"\n→ saved to {args.out}")


if __name__ == "__main__":
    main()
