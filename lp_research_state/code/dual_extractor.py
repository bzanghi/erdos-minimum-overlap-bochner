"""
Dual-objective extractor for CLARABEL on the Bochner-augmented SDP.

*** READ THIS BEFORE QUOTING ANY NUMBER FROM THIS MODULE AS A BOUND. ***

`rigorous_dual_LB` is a MISNOMER kept for backward compatibility.  It is the
solver's dual objective at a near-converged iteration, with NO correction for
dual infeasibility.  A dual objective lower-bounds the primal optimum only when
the dual point is exactly FEASIBLE; CLARABEL's iterate is feasible only up to
its dual residual, and nothing here subtracts a compensating margin.  So this
number is a very good *estimate* of a lower bound, and is fine for ranking
configurations or steering a search — but it is NOT a certificate.

For a certificate use `_jansson_verify.jansson_lower_bound`, which computes a
true a-posteriori bound (Jansson-Chaykin-Keil) in directed-rounding interval
arithmetic, including the dual-defect and cone-distance penalties this module
omits.  At the production centers the two differ by ~1e-5 to ~1e-3, in both
directions, so the distinction is not academic.

TWO DEFECTS WERE FOUND AND FIXED HERE (2026-07-25):
 1. Column mis-parse.  CLARABEL's iteration table is
        iter  pcost  dcost  gap  pres  dres  k/t  mu  step
    The old regex captured five columns and named the fifth `dual_residual` —
    but the fifth column is `pres`, the PRIMAL residual.  The eligibility gate
    was therefore applied to the wrong quantity.  The regex now captures both
    `pres` and `dres`, and the gate uses `dres`.
 2. No feasibility margin.  The old docstring conceded that "for a strict bound
    we would need ZERO dual residual ... we can absorb that into a margin", but
    no margin was ever absorbed.  Rather than invent one, the return value is
    now explicitly flagged `is_certificate: False`.

The cvxpy interface returns `prob.value` at termination with status
`optimal_inaccurate` when the gap is in [tol, reduced_tol]; that status is a
labelling artifact at these gaps and is not itself a reason to distrust the run.
"""
from __future__ import annotations
import re, io, contextlib, time
from typing import Optional


# iter  pcost  dcost  gap  pres  dres  k/t  mu  step
CLARABEL_ITER_RE = re.compile(
    r"^\s*(\d+)\s+([+\-]?\d+\.\d+e[+\-]?\d+)\s+([+\-]?\d+\.\d+e[+\-]?\d+)\s+"
    r"([+\-]?\d+\.\d+e[+\-]?\d+)\s+([+\-]?\d+\.\d+e[+\-]?\d+)\s+"
    r"([+\-]?\d+\.\d+e[+\-]?\d+)\s+"
)


def parse_clarabel_iterations(out: str):
    """Rows of iter, primal_obj, dual_obj, gap, primal_residual, dual_residual.

    `dual_residual` is CLARABEL's `dres` column (the 6th), not `pres`.
    """
    rows = []
    for line in out.splitlines():
        m = CLARABEL_ITER_RE.match(line)
        if not m:
            continue
        rows.append({
            "iter": int(m.group(1)),
            "primal_obj": float(m.group(2)),
            "dual_obj": float(m.group(3)),
            "gap": float(m.group(4)),
            "primal_residual": float(m.group(5)),
            "dual_residual": float(m.group(6)),
        })
    return rows


def best_dual_lower_bound(rows, max_dual_residual: float = 1e-4):
    """Highest dual objective among iterations whose DUAL residual is small.

    NOT a certificate — see the module docstring.  The returned value is the
    dual objective of a point that is dual-feasible only to within
    `dual_residual`, so it can exceed the true optimum by an amount this
    function does not bound.  Use `_jansson_verify.jansson_lower_bound` when a
    number has to survive scrutiny.

    Returns (best_dual_obj, dual_residual_at_that_iter, iter_index, n_eligible).
    """
    eligible = [r for r in rows if r["dual_residual"] <= max_dual_residual]
    if not eligible:
        return None, None, None, 0
    # dual_obj strengthens monotonically along CLARABEL's IPM, so the last
    # eligible iteration carries the tightest value.
    last = eligible[-1]
    return last["dual_obj"], last["dual_residual"], last["iter"], len(eligible)


def solve_with_dual_extraction(prob, solver: str = "CLARABEL"):
    """Solve the cvxpy Problem and extract a rigorous dual lower bound."""
    import cvxpy as cp
    buf = io.StringIO()
    t0 = time.time()
    with contextlib.redirect_stdout(buf):
        prob.solve(solver=solver, verbose=True)
    elapsed = time.time() - t0
    output = buf.getvalue()
    rows = parse_clarabel_iterations(output)
    best_dual, residual, it, n_elig = best_dual_lower_bound(rows)
    return {
        "status": prob.status,
        "reported_value": (float(prob.value) if prob.value is not None else None),
        # name kept for backward compatibility; NOT a certificate (see module docstring)
        "rigorous_dual_LB": best_dual,
        "is_certificate": False,
        "certificate_hint": "use _jansson_verify.jansson_lower_bound for a real bound",
        "dual_residual_at_LB": residual,
        "primal_residual_at_LB": (rows[it]["primal_residual"]
                                  if it is not None and it < len(rows) else None),
        "best_iter": it,
        "n_eligible_iters": n_elig,
        "n_iters_total": rows[-1]["iter"] if rows else 0,
        "time": elapsed,
        "raw_iterations": rows,
    }


# ----- self-test -----------------------------------------------------------
if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from white_full_convex import build_problem
    import cvxpy as cp

    # Run a known case: row1, N=1500, T=600, R=10, bochner_n=10
    N, T, R = 1500, 600, 10
    h, p, qm, qp = 0.015, 0.381, -0.02, 0.02
    Omega, w, v, c, d, eps_v, dlt, cons = build_problem(
        N, T, R, h, h, p, p, qm, qp, bochner_n=10
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)

    res = solve_with_dual_extraction(prob)
    print("Dual extraction self-test:")
    print(f"  status                  = {res['status']}")
    print(f"  reported value          = {res['reported_value']:.10f}")
    print(f"  rigorous dual LB        = {res['rigorous_dual_LB']:.10f}"
          if res['rigorous_dual_LB'] is not None else "  rigorous dual LB        = (none)")
    print(f"  dual residual at LB iter = {res['dual_residual_at_LB']:.2e}"
          if res['dual_residual_at_LB'] is not None else "")
    print(f"  iter at LB              = {res['best_iter']}/{res['n_iters_total']}")
    print(f"  eligible iters          = {res['n_eligible_iters']}")
    if res["rigorous_dual_LB"] is not None and res["reported_value"] is not None:
        print(f"  reported - rigorous     = {res['reported_value'] - res['rigorous_dual_LB']:.2e}")
