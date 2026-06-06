"""Lever I' Step 1: extract λ_m (cell-envelope dual multipliers) from a single SDP solve
to test the O(1/m²) decay hypothesis from LEVER_I_PRIME_POC.md.

If the decay holds, the residual-enumeration saturation theorem becomes tractable.
If it fails, we need a different multiplier bound (or the theorem is not viable).

This solves at row 4 (binding row) at moderate scale (N=3000, T=1200) to save time,
since dual decay is a structural property and should not depend sensitively on N, T.

Output: lp_research_state/data/lambda_m_extracted.json
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from white_full_convex import (
    build_problem,  # use the standard builder to get the full constraint stack
)


def solve_and_extract_lambdas(N=3000, T=1200, R=10, bochner_n=20,
                               h1=0.004, h2=0.004, p1=0.3875, p2=0.3875,
                               q1=-0.02, q2=0.02):
    """Solve row 4 SDP at modest scale, return per-constraint dual values keyed by index.

    The cell-envelope constraints in white_full_convex.py:176-200 add 2*R = 20 inequalities
    indexed by m = 1..2R. We solve, extract every constraint's dual_value, and
    cross-reference by inspecting the cvxpy expression source.
    """
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, q1, q2,
        bochner_n=bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=cp.CLARABEL, verbose=False)
    dt = time.time() - t0
    print(f"Solve done: status={prob.status}, Omega={prob.value:.10f}, time={dt:.1f}s")

    # Collect every constraint's dual_value
    duals = []
    for i, c_i in enumerate(cons):
        dv = c_i.dual_value
        if dv is None:
            duals.append((i, None, str(c_i)[:200]))
        else:
            # Try to convert to scalar
            try:
                if np.isscalar(dv):
                    duals.append((i, float(dv), str(c_i)[:200]))
                else:
                    arr = np.asarray(dv).ravel()
                    duals.append((i, [float(x) for x in arr[:50]], str(c_i)[:200]))
            except Exception:
                duals.append((i, "(unparseable)", str(c_i)[:200]))
    return prob.status, float(prob.value), dt, duals


def main():
    print("=== Lever I' Step 1: extract λ_m from cell-envelope constraints ===")
    status, value, dt, duals = solve_and_extract_lambdas()
    out = {
        "config": {"N": 3000, "T": 1200, "R": 10, "bochner_n": 20,
                   "row": "row4", "h_c": 0.004, "p_c": 0.3875, "q1": -0.02, "q2": 0.02},
        "status": status,
        "Omega": value,
        "solve_time_s": dt,
        "n_constraints": len(duals),
        "duals": [{"idx": i, "dual": dv, "expr_preview": expr}
                  for (i, dv, expr) in duals],
    }
    out_path = Path(__file__).parent.parent / "data" / "lambda_m_extracted.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"Wrote {out_path}")

    # Now: identify which constraints are the cell-envelope ones.
    # Pattern: lines 176-200 of white_full_convex.py have constraints involving
    # `cos_cell_bounds_*` or `sin_cell_bounds_*` results, applied as inequalities
    # against am, bm expressions. They appear consecutively in the cons list
    # AFTER the (5.3), (5.4) constraints and BEFORE the eps/dlt tail bounds.
    print("\n=== Searching for cell-envelope constraint indices ===")
    cell_envelope_candidates = []
    for i, (idx, dv, expr) in enumerate(duals):
        # Heuristic: cell-envelope constraints contain w + v in their LHS as a
        # dot product, plus quadratic-in-(am, bm) terms.
        if "+" in expr and "@" in expr:
            cell_envelope_candidates.append((idx, dv, expr[:120]))
    print(f"Found {len(cell_envelope_candidates)} candidate constraints")
    for idx, dv, expr in cell_envelope_candidates[:25]:
        if isinstance(dv, float):
            print(f"  idx={idx}: dual={dv:.6e}  expr={expr}")
        else:
            print(f"  idx={idx}: dual=<array>  expr={expr}")


if __name__ == "__main__":
    main()
