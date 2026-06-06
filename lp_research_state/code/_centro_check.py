"""
PRO-49 centrosymmetric split — verification scratchpad.

(A) CONE INVENTORY: prove CLARABEL actually receives the half-size PSD cone
    (the test the complex-Hermitian route FAILED — cvxpy lowered it back to 2(n+1)).
(B) EXACT-EQUIVALENCE at the binding row4 optimum, via Bochner-block eigenvalue
    spectra + constraint residuals (NOT objective digits — CLARABEL nondeterminism
    caps objective agreement at ~7-9 digits per the prior finding).

Moderate scale only (N<=5000) to coexist with the co-running production L2 solve.
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

from white_full_convex import build_problem
from white_full_convex_centro import build_problem_centro
from bochner_independent import make_real_form
from bochner_centro import make_centro_block

# row4 = binding row: (h, p, q) = (0.004, 0.3875, ±0.02)
ROW4 = dict(h1=0.004, h2=0.004, p1=0.3875, p2=0.3875, q1=-0.02, q2=0.02)


def psd_sides(prob):
    """Return sorted list of PSD cone side-lengths CLARABEL actually receives."""
    data, chain, inv = prob.get_problem_data(solver=cp.CLARABEL)
    dims = data["dims"]
    # cvxpy ConeDims: .psd is a list of side lengths
    return sorted(list(dims.psd)), data["A"].shape, int(data["A"].nnz)


def cone_inventory(N, T, R, bn):
    print(f"\n=== (A) CONE INVENTORY  N={N} T={T} R={R} bochner_n={bn} ===")
    Omega_r, *_, cons_r = build_problem(N, T, R, **ROW4, bochner_n=bn)
    prob_r = cp.Problem(cp.Minimize(Omega_r), cons_r)
    sides_r, shape_r, nnz_r = psd_sides(prob_r)

    Omega_c, *_, cons_c = build_problem_centro(N, T, R, **ROW4, bochner_n=bn)
    prob_c = cp.Problem(cp.Minimize(Omega_c), cons_c)
    sides_c, shape_c, nnz_c = psd_sides(prob_c)

    print(f"  REAL form  PSD sides = {sides_r}   A={shape_r} nnz={nnz_r}")
    print(f"  CENTRO     PSD sides = {sides_c}   A={shape_c} nnz={nnz_c}")
    print(f"  expected: real has two {2*(bn+1)}-blocks; centro has two {bn+1}-blocks")
    return sides_r, sides_c


def exact_equivalence(N, T, R, bn):
    """Solve both at row4; compare Bochner-block min-eigs + full spectra at the
    optimum (machine-precision, IPM-independent), plus objective for context."""
    print(f"\n=== (B) EXACT-EQUIVALENCE at row4  N={N} T={T} R={R} bochner_n={bn} ===")
    res = {}
    for label, builder in (("real", build_problem), ("centro", build_problem_centro)):
        Omega, w, v, c, d, eps, dlt, cons = builder(N, T, R, **ROW4, bochner_n=bn)
        prob = cp.Problem(cp.Minimize(Omega), cons)
        prob.solve(solver="CLARABEL")
        res[label] = dict(
            status=prob.status,
            value=(float(prob.value) if prob.value is not None else None),
            c=np.asarray(c.value), d=np.asarray(d.value),
        )
        print(f"  [{label}] status={prob.status} value={prob.value:.10f}")

    # Reconstruct Bochner blocks at EACH solver's own optimum and at the OTHER's,
    # to show (i) feasibility carries over, (ii) spectra match for shared (c,d).
    print("  -- Bochner block min-eigenvalues at each optimum (both signs) --")
    worst_share = 0.0
    for src in ("real", "centro"):
        cv, dv = res[src]["c"], res[src]["d"]
        for sign in (+1, -1):
            RF, _, _ = make_real_form(cv, dv, bn, sign)
            Bk = make_centro_block(cv, dv, bn, sign)
            eRF = np.linalg.eigvalsh(RF)
            eBk = np.linalg.eigvalsh(Bk)
            # RF spectrum should equal Bk spectrum doubled
            spec_err = np.max(np.abs(np.sort(eRF) - np.sort(np.concatenate([eBk, eBk]))))
            worst_share = max(worst_share, spec_err)
            print(f"    ({src} opt, sign={sign:+d}) min_eig RF={eRF.min():.3e}  "
                  f"min_eig Bk={eBk.min():.3e}  |spec_RF - 2x spec_Bk|={spec_err:.2e}")
    print(f"  max |RF spectrum - Bk spectrum doubled| over all = {worst_share:.2e}")
    obj_diff = abs(res["real"]["value"] - res["centro"]["value"])
    cd_diff = max(np.max(np.abs(res["real"]["c"] - res["centro"]["c"])),
                  np.max(np.abs(res["real"]["d"] - res["centro"]["d"])))
    print(f"  objective real-vs-centro diff = {obj_diff:.2e}  "
          f"(expect ~1e-7..1e-9 = CLARABEL solve-to-solve noise, NOT an error)")
    print(f"  optimal (c,d) real-vs-centro max diff = {cd_diff:.2e}")
    return res


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=2000)
    ap.add_argument("--T", type=int, default=800)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=20)
    ap.add_argument("--mode", choices=["cone", "equiv", "both"], default="both")
    a = ap.parse_args()
    if a.mode in ("cone", "both"):
        cone_inventory(a.N, a.T, a.R, a.bn)
    if a.mode in ("equiv", "both"):
        exact_equivalence(a.N, a.T, a.R, a.bn)
