"""Lever I' Step D: scaled lambda_m, sigma_m extraction.

Runs the cosine and sine cell-envelope dual extraction at user-specified
(N, row) or (N, list-of-rows), incrementally appending to
lp_research_state/data/lambda_m_scaled.json.

Usage:
    python3 _lever_i_prime_lambda_m_scaled.py <N> <row_label> [more_row_labels...]
    python3 _lever_i_prime_lambda_m_scaled.py 15000 row7
    python3 _lever_i_prime_lambda_m_scaled.py 15000 row7 row4 row1 cde_n30_iter1

Output file: lp_research_state/data/lambda_m_scaled.json
Each (N, row_label) entry is keyed as "N=<N>_row=<label>" and includes:
  Omega, status, solve_time_s, lambda_m, sigma_pairs, params
Reruns are skipped if the entry already exists (idempotent).
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

DATA_PATH = Path(__file__).parent.parent / "data" / "lambda_m_scaled.json"


def load_existing():
    if DATA_PATH.exists():
        return json.loads(DATA_PATH.read_text())
    return {"runs": {}}


def save(blob):
    DATA_PATH.write_text(json.dumps(blob, indent=2, default=str))


def extract(N, T, R, bochner_n, h, p, q1, q2):
    """Solve and pull both the 2R cosine and 4R sine cell-envelope duals."""
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h, h, p, p, q1, q2,
        bochner_n=bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=cp.CLARABEL, verbose=False)
    dt = time.time() - t0

    cos_start = 8
    cos_end = cos_start + 2 * R
    sin_start = cos_end
    sin_end = sin_start + 4 * R

    def _scalar(ci):
        dv = ci.dual_value
        if dv is None:
            return None
        if np.isscalar(dv):
            return float(dv)
        arr = np.asarray(dv).ravel()
        return float(arr[0]) if arr.size == 1 else None

    lam = [_scalar(cons[i]) for i in range(cos_start, cos_end)]
    sigma_pairs = []
    for m in range(1, 2 * R + 1):
        i1 = sin_start + 2 * (m - 1)
        i2 = i1 + 1
        sigma_pairs.append((_scalar(cons[i1]), _scalar(cons[i2])))

    return {
        "status": prob.status,
        "Omega": float(prob.value) if prob.value is not None else None,
        "solve_time_s": dt,
        "lambda_m": lam,
        "sigma_pairs": sigma_pairs,
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n},
        "params": {"h": h, "p": p, "q1": q1, "q2": q2},
    }


def _summary(res):
    """Return (Σ m|λ|, Σ m³|λ|, Σ m|σ|, Σ m³|σ|, ResidualGain at config N, C_explicit)."""
    LB = 0.3801279
    Omega = 0.38
    if res.get("lambda_m") is None:
        return None
    lam = np.array([abs(x or 0) for x in res["lambda_m"]])
    sig = np.array([abs(p[0] or 0) + abs(p[1] or 0) for p in res["sigma_pairs"]])
    ms = np.arange(1, len(lam) + 1)
    Sml = float((ms * lam).sum())
    Sm3l = float((ms**3 * lam).sum())
    Sms = float((ms * sig).sum())
    Sm3s = float((ms**3 * sig).sum())
    N = res["config"]["N"]
    cosR = np.pi / (2 * N) * Sml + np.pi**2 * Omega / (3 * N**3) * Sm3l
    sinR = np.pi / (2 * N) * Sms + np.pi**2 * Omega / (3 * N**3) * Sm3s
    combined = cosR + sinR
    return Sml, Sm3l, Sms, Sm3s, cosR, sinR, combined, LB + combined


def main():
    if len(sys.argv) < 3:
        print("Usage: _lever_i_prime_lambda_m_scaled.py <N> <row_label> [more_rows...]")
        sys.exit(2)
    N = int(sys.argv[1])
    row_labels = sys.argv[2:]
    T = 1200
    R = 10
    bochner_n = 20

    blob = load_existing()
    runs = blob.setdefault("runs", {})

    for label in row_labels:
        if label not in ROWS:
            print(f"  Unknown row label {label!r}; skipping")
            continue
        key = f"N={N}_row={label}"
        if key in runs and runs[key].get("Omega") is not None:
            print(f"  [{key}] already done (Omega={runs[key]['Omega']:.6f}); skipping")
            continue
        h, p, q1, q2 = ROWS[label]
        print(f"  [{key}] solving (h={h}, p={p}, q=[{q1},{q2}]), N={N}, T={T}, R={R}, bochner_n={bochner_n}")
        try:
            res = extract(N, T, R, bochner_n, h, p, q1, q2)
        except Exception as e:
            print(f"  [{key}] FAILED: {type(e).__name__}: {e}")
            runs[key] = {"error": str(e), "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n},
                         "params": {"h": h, "p": p, "q1": q1, "q2": q2}}
            save(blob)
            continue

        runs[key] = res
        save(blob)

        s = _summary(res)
        if s is not None:
            Sml, Sm3l, Sms, Sm3s, cosR, sinR, combined, Cexpl = s
            print(f"  [{key}] DONE status={res['status']}, Omega={res['Omega']:.8f}, time={res['solve_time_s']:.1f}s")
            print(f"         Σ m·λ = {Sml:.4f}, Σ m³·λ = {Sm3l:.2f}")
            print(f"         Σ m·σ = {Sms:.4f}, Σ m³·σ = {Sm3s:.2f}")
            print(f"         cos resid = {cosR:.4e}, sin resid = {sinR:.4e}, combined = {combined:.4e}")
            print(f"         C_explicit at N={N}: {Cexpl:.6f}  (vs Together UB 0.380871)")
            if Cexpl < 0.380871:
                print(f"         *** SATURATION THEOREM NON-VACUOUS *** (margin {0.380871 - Cexpl:.4e})")
            else:
                print(f"         theorem still vacuous (overage {Cexpl - 0.380871:.4e})")


if __name__ == "__main__":
    main()
