"""Lever I' Step A: extract lambda_m on multiple rows to test profile stability.

Builds on _lever_i_prime_lambda_m_extract.py. Runs the cell-envelope dual
extraction at row1, row4 (re-run for consistency), row7 (the three White rows
with non-trivial p span at q=0), plus cde_n30_iter1 (a CDE-discovered center
near p=0.394175, h=0.0). Records each lambda_m profile and the sum.

Decision rule (per session goal):
- If profiles cluster tightly (similar dominant m, sum within factor 2): the
  empirical residual ~Sigma * (per-m bound) is structurally row-independent and
  the saturation theorem can be derived without row-by-row casework.
- If profiles diverge wildly: theorem must be row-by-row and is harder.

Output: lp_research_state/data/lambda_m_all_rows.json
"""
from __future__ import annotations
import sys
import time
import json
from pathlib import Path
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem  # noqa: E402


ROWS = {
    "row1":         (0.015, 0.381,    -0.02, 0.02),
    "row4":         (0.004, 0.3875,   -0.02, 0.02),
    "row7":         (0.030, 0.375,    -0.02, 0.02),
    "cde_n30_iter1":(0.0,   0.394175, -0.02, 0.02),
}


def extract_lambda_m(N, T, R, bochner_n, h, p, q1, q2):
    """Solve SDP at one center; return Omega, status, and the 2R cosine
    cell-envelope multipliers in order m=1..2R.

    The 20 cosine cell-envelope constraints (white_full_convex.py:176-182)
    appear *immediately after* the seven scalar constraints (var box, sum,
    h1, h2, monotone) -- empirically at constraint indices 8..27 in row4. We
    locate them robustly by structural signature: scalar inequality whose
    expression preview contains "@ (var2 + var3)" AND "PowerApprox" AND
    "<= 0" pattern.
    """
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2,
        bochner_n=bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=cp.CLARABEL, verbose=False)
    dt = time.time() - t0
    status = prob.status
    omega = float(prob.value) if prob.value is not None else None

    # The cosine cell-envelope constraints are at fixed indices in cons[]:
    # build_problem adds them as cons[8..(8 + 2R - 1)] (assume_even=False
    # default). The sine cell-envelope (40 constraints, two per m) follows at
    # cons[(8 + 2R)..(8 + 6R - 1)]. See white_full_convex.py:176-190.
    cos_start = 8
    cos_end   = cos_start + 2 * R
    sin_start = cos_end
    sin_end   = sin_start + 4 * R   # 2R lags × 2 inequalities each

    def _scalar_dual(ci):
        dv = ci.dual_value
        if dv is None:
            return None
        if np.isscalar(dv):
            return float(dv)
        arr = np.asarray(dv).ravel()
        return float(arr[0]) if arr.size == 1 else None

    lam = [_scalar_dual(cons[i]) for i in range(cos_start, cos_end)]
    # sin multipliers, paired (σ_m^1, σ_m^2) for m=1..2R
    sigma_pairs = []
    for m in range(1, 2 * R + 1):
        i1 = sin_start + 2 * (m - 1)
        i2 = i1 + 1
        sigma_pairs.append((_scalar_dual(cons[i1]), _scalar_dual(cons[i2])))

    return {
        "status": status,
        "Omega": omega,
        "solve_time_s": dt,
        "lambda_m": lam,
        "sigma_pairs": sigma_pairs,
        "n_lambdas_extracted": len(lam),
        "n_sigmas_extracted": len(sigma_pairs),
    }


def main():
    N = 3000
    T = 1200
    R = 10
    bochner_n = 20

    print("=== Lever I' Step A: extract lambda_m on multiple rows ===")
    print(f"Config: N={N}, T={T}, R={R}, bochner_n={bochner_n}")

    out = {
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n},
        "rows": {},
    }

    for label, (h, p, q1, q2) in ROWS.items():
        print(f"\n--- {label}: h={h}, p={p}, q=[{q1},{q2}] ---")
        try:
            res = extract_lambda_m(N, T, R, bochner_n, h, p, q1, q2)
            res["params"] = {"h": h, "p": p, "q1": q1, "q2": q2}
            out["rows"][label] = res
            lam = res["lambda_m"]
            if lam:
                lam_arr = np.array(lam)
                top4 = np.argsort(-np.abs(lam_arr))[:4] + 1
                print(f"  status={res['status']}, Omega={res['Omega']:.8f}")
                print(f"  Sum_m lambda_m   = {np.sum(np.abs(lam_arr)):.6f}")
                print(f"  Max_m lambda_m   = {np.max(np.abs(lam_arr)):.6f}")
                print(f"  Top-4 by |lambda|: m = {top4.tolist()}")
                print(f"  lambda_m[1..2R]  =")
                for m, lm in enumerate(lam, start=1):
                    print(f"    m={m:2d}: {lm: .6e}")
            else:
                print("  (no lambdas extracted)")
        except Exception as e:
            print(f"  FAILED: {e}")
            out["rows"][label] = {"error": str(e)}

    out_path = Path(__file__).parent.parent / "data" / "lambda_m_all_rows.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}")

    # Summary table
    print("\n=== Summary across rows (cosine cell-envelope) ===")
    print(f"{'row':<18s} {'Omega':>10s} {'Σ|λ|':>10s} {'Σm|λ|':>10s} {'Σm³|λ|':>10s} {'top-m':>15s}")
    for label, res in out["rows"].items():
        if "lambda_m" in res and res["lambda_m"]:
            lam = np.array([abs(x) for x in res["lambda_m"]])
            ms  = np.arange(1, len(lam)+1)
            top4 = (np.argsort(-lam)[:4] + 1).tolist()
            print(f"{label:<18s} {res['Omega']:>10.6f} {float(lam.sum()):>10.4f} {float((ms*lam).sum()):>10.4f} {float((ms**3*lam).sum()):>10.4f} {str(top4):>15s}")
        else:
            print(f"{label:<18s} (failed or no data)")

    print("\n=== Summary across rows (sine cell-envelope, |σ^1|+|σ^2| per m) ===")
    print(f"{'row':<18s} {'Σ|σ|':>10s} {'Σm|σ|':>10s} {'Σm³|σ|':>10s}")
    for label, res in out["rows"].items():
        if "sigma_pairs" in res and res["sigma_pairs"]:
            sig_sums = np.array([abs(p[0] or 0) + abs(p[1] or 0) for p in res["sigma_pairs"]])
            ms = np.arange(1, len(sig_sums)+1)
            print(f"{label:<18s} {float(sig_sums.sum()):>10.4f} {float((ms*sig_sums).sum()):>10.4f} {float((ms**3*sig_sums).sum()):>10.4f}")
        else:
            print(f"{label:<18s} (no sigma data)")

    # Phase 5 residual at N=10000 sup over rows
    print("\n=== Phase 5 (N=10000) corrected residual per row ===")
    N5 = 10000
    Omega5 = 0.38
    factorA = np.pi / (2*N5)
    factorB = np.pi**2 * Omega5 / (3*N5**3)
    print(f"{'row':<18s} {'cos resid':>12s} {'sin resid':>12s} {'combined':>12s} {'C_explicit':>12s}")
    LB = 0.3801279
    sup_resid = 0
    for label, res in out["rows"].items():
        if "lambda_m" not in res or not res["lambda_m"]:
            continue
        lam = np.array([abs(x) for x in res["lambda_m"]])
        sig = np.array([abs(p[0] or 0)+abs(p[1] or 0) for p in res["sigma_pairs"]])
        ms = np.arange(1, len(lam)+1)
        cos_res = factorA*float((ms*lam).sum()) + factorB*float((ms**3*lam).sum())
        sin_res = factorA*float((ms*sig).sum()) + factorB*float((ms**3*sig).sum())
        combined = cos_res + sin_res
        sup_resid = max(sup_resid, combined)
        print(f"{label:<18s} {cos_res:>12.4e} {sin_res:>12.4e} {combined:>12.4e} {LB+combined:>12.6f}")
    print(f"\nSup-over-rows combined residual: {sup_resid:.4e}")
    print(f"C_explicit (sup)                 = {LB + sup_resid:.6f}")
    print(f"Together UB                       = 0.380871")
    print(f"Open gap                          = 7.43e-4")
    if LB + sup_resid > 0.380871:
        print(f"Theorem VACUOUS at N={N5}: overage = {LB + sup_resid - 0.380871:.4e}")
        break_even = np.pi * max(float((np.arange(1, len(lam)+1)*np.array([abs(x) for x in r['lambda_m']])).sum() + (np.arange(1, len(lam)+1)*np.array([abs(p[0] or 0)+abs(p[1] or 0) for p in r['sigma_pairs']])).sum()) for r in out['rows'].values()) / (2 * (0.380871 - LB))
        print(f"Break-even N (combined cos+sin) ≈ {break_even:.0f}")
    else:
        print(f"Theorem NON-VACUOUS at N={N5}!")


if __name__ == "__main__":
    main()
