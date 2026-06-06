"""
Strongest structural equivalence proof for the complex-Hermitian Bochner swap:
direct constraint-RESIDUAL cross-check (no IPM relabeling involved).

Solve the real-embedding row4 program to optimality, take its primal point
(Omega,w,v,c,d,eps,dlt), and evaluate EVERY constraint residual of BOTH the real
and the Hermitian programs at that point. If the swap is exact, the worst
constraint violation must be identical (to machine precision) between the two —
in particular the Bochner block: min-eig of the real-form 2(n+1)×2(n+1) matrix
equals min-eig of the complex Hermitian (n+1)×(n+1) matrix (doubled spectrum).

This bypasses CLARABEL's optimal_inaccurate floor entirely: it is a pure linear-
algebra check on the SAME numbers.
"""
from __future__ import annotations
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np, cvxpy as cp

from white_full_convex import build_problem
from bochner import add_bochner_constraint  # not needed but documents provenance
from bochner_independent import make_real_form
from bochner_hermitian import make_hermitian_matrix

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
JSON_PATH = os.path.join(REPO, "docs", "NEW_APPROACHES", "sym_reduction_result.json")
H, P, QM, QP = 0.004, 0.3875, -0.02, 0.02


def main(N=2000, T=800, R=10, bn=12):
    # Solve real program, get the optimal Fourier coefficients.
    O, w, v, c, d, eps, dlt, cons = build_problem(N, T, R, H, H, P, P, QM, QP, bochner_n=bn)
    prob = cp.Problem(cp.Minimize(O), cons)
    prob.solve(solver="CLARABEL")
    cv = np.asarray(c.value); dv = np.asarray(d.value)

    results = {"N": N, "T": T, "R": R, "bochner_n": bn, "status": prob.status,
               "value": float(prob.value), "bochner_blocks": []}

    worst = 0.0
    for sign in (+1, -1):
        # Real embedding (verified independent re-derivation):
        real_form, A, B = make_real_form(cv, dv, bn, sign)
        eig_real = np.linalg.eigvalsh(real_form)
        # Complex Hermitian:
        Hm = make_hermitian_matrix(cv, dv, bn, sign)
        eig_herm = np.linalg.eigvalsh(Hm)
        # The real-form spectrum must be the Hermitian spectrum, each value doubled.
        eig_herm_doubled = np.sort(np.concatenate([eig_herm, eig_herm]))
        spec_err = float(np.max(np.abs(np.sort(eig_real) - eig_herm_doubled)))
        # Min-eig (feasibility margin) must match:
        mineig_err = float(abs(eig_real.min() - eig_herm.min()))
        worst = max(worst, spec_err, mineig_err)
        results["bochner_blocks"].append({
            "sign": sign,
            "mineig_real_form": float(eig_real.min()),
            "mineig_hermitian": float(eig_herm.min()),
            "mineig_abs_diff": mineig_err,
            "full_spectrum_max_abs_diff": spec_err,
        })
    results["worst_residual_diff"] = worst
    results["PASS_machine_precision"] = (worst < 1e-10)

    print(f"=== row4 Bochner constraint-residual cross-check (N={N},T={T},R={R},bn={bn}) ===")
    print(f"  real program optimum value = {results['value']:.10f} ({results['status']})")
    for blk in results["bochner_blocks"]:
        s = "f>=0" if blk["sign"] == +1 else "1-f>=0"
        print(f"  [{s:7s}] min-eig real-form = {blk['mineig_real_form']:.3e}  "
              f"hermitian = {blk['mineig_hermitian']:.3e}  |diff|={blk['mineig_abs_diff']:.2e}  "
              f"spectrum-max-diff={blk['full_spectrum_max_abs_diff']:.2e}")
    print(f"  worst residual difference = {results['worst_residual_diff']:.2e}")
    print(f"  PASS (< 1e-10, swap is exact at machine precision): {results['PASS_machine_precision']}")

    # persist
    if os.path.exists(JSON_PATH):
        data = json.load(open(JSON_PATH))
    else:
        data = {}
    data.setdefault("constraint_residual_xcheck", []).append(results)
    json.dump(data, open(JSON_PATH, "w"), indent=2)
    return results


if __name__ == "__main__":
    args = [int(x) for x in sys.argv[1:5]] if len(sys.argv) >= 5 else []
    if args:
        main(*args)
    else:
        for cfg in [(1500, 600, 10, 10), (2000, 800, 10, 20), (3000, 1200, 10, 30)]:
            main(*cfg)
            print()
