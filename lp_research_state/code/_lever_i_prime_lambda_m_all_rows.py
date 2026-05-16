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
    # default; with assume_even=True, they shift by 3). See
    # white_full_convex.py:140-182.
    start = 8
    end = start + 2 * R
    lam = []
    indices = []
    for i in range(start, end):
        ci = cons[i]
        dv = ci.dual_value
        if dv is None:
            lam.append(None)
            indices.append(i)
            continue
        if np.isscalar(dv):
            lam.append(float(dv))
        else:
            arr = np.asarray(dv).ravel()
            lam.append(float(arr[0]) if arr.size == 1 else None)
        indices.append(i)
    # Sanity check: the constraint at each index should have "@ (var2 + var3)"
    # in its stringification.
    sanity_ok = all(("@ (var2 + var3)" in str(cons[i])) for i in indices)
    if not sanity_ok:
        print(f"  WARNING: cell-envelope index detection failed sanity check.")
        for i in indices:
            print(f"    cons[{i}]: {str(cons[i])[:120]}")

    return {
        "status": status,
        "Omega": omega,
        "solve_time_s": dt,
        "lambda_m": lam,
        "constraint_indices": indices,
        "n_lambdas_extracted": len(lam),
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
    print("\n=== Summary across rows ===")
    print(f"{'row':<18s} {'Omega':>10s} {'Sum|lam|':>10s} {'Max|lam|':>10s} {'top-m':>15s}")
    for label, res in out["rows"].items():
        if "lambda_m" in res and res["lambda_m"]:
            lam = np.array(res["lambda_m"])
            sumlam = float(np.sum(np.abs(lam)))
            maxlam = float(np.max(np.abs(lam)))
            top4 = (np.argsort(-np.abs(lam))[:4] + 1).tolist()
            print(f"{label:<18s} {res['Omega']:>10.6f} {sumlam:>10.4f} {maxlam:>10.4f} {str(top4):>15s}")
        else:
            print(f"{label:<18s} (failed or no data)")


if __name__ == "__main__":
    main()
