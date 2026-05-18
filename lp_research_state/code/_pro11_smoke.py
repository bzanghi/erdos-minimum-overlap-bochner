"""PRO-11 smoke test: row 4 at small scale, CLARABEL vs SDPA-GMP via our serializer."""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
from sdpa_runner import solve_with_sdpa_gmp
import cvxpy as cp


def run(N, T, R, bochner_n, timeout=1800, label=""):
    print(f"\n=== {label} :  N={N} T={T} R={R} bochner_n={bochner_n} ===")
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bochner_n
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)

    print("--- CLARABEL ---")
    t0 = time.time()
    cl = solve_with_dual_extraction(prob)
    t_cl = time.time() - t0
    print(f"  reported value          = {cl['reported_value']!r}")
    print(f"  rigorous dual LB        = {cl['rigorous_dual_LB']!r}")
    print(f"  status                  = {cl['status']}")
    print(f"  time                    = {cl['time']:.2f}s  (wall {t_cl:.2f}s)")

    # Build a FRESH problem instance for SDPA (cvxpy caches state otherwise)
    Omega2, w2, v2, c2, d2, eps2, dlt2, cons2 = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bochner_n
    )
    prob2 = cp.Problem(cp.Minimize(Omega2), cons2)

    print("--- SDPA-GMP via custom serializer ---")
    sd = solve_with_sdpa_gmp(prob2, timeout_sec=timeout, verbose=False)
    print(f"  phase                   = {sd['phase']}")
    print(f"  primal_obj              = {sd['primal_obj']!r}")
    print(f"  dual_obj                = {sd['dual_obj']!r}")
    print(f"  duality_gap             = {sd['duality_gap']!r}")
    print(f"  rigorous_dual_LB        = {sd['rigorous_dual_LB']!r}")
    print(f"  precision_digits        = {sd['precision_digits']!r}")
    print(f"  iterations              = {sd['iterations']!r}")
    print(f"  runtime                 = {sd['runtime_sec']:.2f}s "
          f"(serialize={sd['serialize_sec']:.2f}s, solve={sd['solve_sec']:.2f}s)")
    print(f"  m                       = {sd['m']}, blocks = {len(sd['block_structure'])}")
    print(f"  dat-s path              = {sd['sdpa_dat_s_path']}")

    if cl["rigorous_dual_LB"] is not None and sd["dual_obj"] is not None:
        diff_dual = sd["dual_obj"] - cl["rigorous_dual_LB"]
        print(f"\n  >> SDPA dual_obj - CLARABEL rigorous LB = {diff_dual:.3e}")
        if cl["reported_value"] is not None and sd["primal_obj"] is not None:
            diff_pri = sd["primal_obj"] - cl["reported_value"]
            print(f"  >> SDPA primal_obj - CLARABEL reported   = {diff_pri:.3e}")
    return cl, sd


if __name__ == "__main__":
    # Smoke: very small
    run(N=200, T=80, R=10, bochner_n=4, label="SMOKE")
