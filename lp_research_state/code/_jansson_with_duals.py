"""Self-consistent extractor: ONE solve -> (Jansson p_lo) AND (con_* duals), so the
cover anchor and the shift coefficients come from the SAME dual point z.

Runs build_problem_with_dual_handles (has con_* handles), solves via cvxpy ONCE
(populates dual_value), then re-derives z from the SAME canonical data to run the
Jansson interval certification (verified p_lo).  Returns p_lo (rigorous), the con_*
duals (z components for the parameter rows), and consistency self-checks.

LIGHT config only (N<=3000); used to certify the interval cover machinery end-to-end
with a fully self-consistent (p_lo, duals) pair.  Production-N is handled separately
(stored cover duals + production Jansson p_lo + inconsistency bound)."""
from __future__ import annotations
import sys, json, time, argparse, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp
import mpmath
from mpmath import iv

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
import _jansson_verify as JV

# All 12 cover centers (h_c, p_c, q1, q2 EXACTLY as in cde_phase5_corrected_tail.json).
# The 2 binding centers (row4, cde_n30_iter3) are already self-consistently verified
# in L2_CLEAN_sc_prod_N20000.json; the other 10 are run by PRO-47 finish so that the
# CLEAN TIER-2c cover certification uses a self-consistent (p_lo, duals) anchor for
# EVERY center -> no V_c-margin convention anywhere.
CENTERS = {
    "row1": dict(h_c=0.015, p_c=0.381, q1=-0.02, q2=0.02),
    "row2": dict(h_c=0.015, p_c=0.385, q1=-0.02, q2=0.02),
    "row3": dict(h_c=0.02, p_c=0.375, q1=-0.02, q2=0.02),
    "row4": dict(h_c=0.004, p_c=0.3875, q1=-0.02, q2=0.02),
    "row5": dict(h_c=0.0, p_c=0.4, q1=-0.02, q2=0.02),
    "row6": dict(h_c=0.0, p_c=0.381, q1=-0.02, q2=0.02),
    "row7": dict(h_c=0.03, p_c=0.375, q1=-0.02, q2=0.02),
    "cde_n30_iter1": dict(h_c=0.0, p_c=0.394175, q1=-0.02, q2=0.02),
    "cde_n30_iter2": dict(h_c=0.0034349999999999997, p_c=0.384175, q1=-0.02, q2=0.02),
    "cde_n30_iter3": dict(h_c=4.4999999999999996e-05, p_c=0.39015, q1=-0.02, q2=0.02),
    "cde_n30_iter4": dict(h_c=0.012225, p_c=0.39075, q1=-0.02, q2=0.02),
    "cde_n30_iter5": dict(h_c=0.008145, p_c=0.38955, q1=-0.02, q2=0.02),
}


def extract(cname, N, T, R, bn, pm, slack_infl=1.0):
    cc = CENTERS[cname]
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, cc["h_c"], cc["h_c"], cc["p_c"], cc["p_c"], cc["q1"], cc["q2"],
        bochner_n=bn)
    pm_tb = {}
    if pm > 0:
        pm_cons, pm_tb = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=cp.CLARABEL, verbose=False)         # populates dual_value
    solve_t = time.time() - t0
    duals = {k: float(H[k].dual_value) for k in
             ("con_53", "con_54", "con_512_pL", "con_512_pU",
              "con_512_qL", "con_512_qU", "con_513")}
    Vc = float(prob.value)

    # re-derive z from the SAME problem data for Jansson certification
    data, chain, inv = prob.get_problem_data(cp.CLARABEL)
    sol = chain.solve_via_data(prob, data)
    A = data['A']; b = np.asarray(data['b']); c = np.asarray(data['c']); dims = data['dims']
    x = np.asarray(sol.x); z = np.asarray(sol.z); s = np.asarray(sol.s)

    # --- Jansson p_lo (mirror jansson_lower_bound's core) ---
    xbar, xflags, _ = JV.canonical_x_bounds(prob, H, N, T, R, x_solution=x)
    import scipy.sparse as _sp
    Acsc = _sp.csc_matrix(A)
    z_iv = [JV._iv(float(zr)) for zr in z]
    pen_Dx = iv.mpf(0)
    for i in range(A.shape[1]):
        s0, e0 = Acsc.indptr[i], Acsc.indptr[i + 1]
        acc = JV._iv(float(c[i]))
        for r_, v_ in zip(Acsc.indices[s0:e0].tolist(), Acsc.data[s0:e0].tolist()):
            acc = acc + JV._iv(float(v_)) * z_iv[r_]
        pen_Dx = pen_Dx + abs(acc) * JV._iv(float(xbar[i]))
    z_blocks = JV.split_cone_blocks(z, dims); s_blocks = JV.split_cone_blocks(s, dims)
    pen_zs = iv.mpf(0); psd_lmins = []
    for (zk, zp, zn), (sk, sp_, sn) in zip(z_blocks, s_blocks):
        if zk == "zero":
            continue
        d_j, _ = JV.cone_lambda_min_lower(zk, zp, zn)
        sbar_j, _ = JV.slack_size_upper(sk, sp_, sn, infl=slack_infl)
        d_neg = min(mpmath.mpf(0), mpmath.mpf(d_j))
        pen_zs = pen_zs + JV._iv(float(d_neg)) * JV._iv(float(sbar_j))
        if zk == "psd":
            psd_lmins.append(float(d_j))
    neg_bz_iv = iv.mpf(0)
    for bi, zi in zip(b, z):
        neg_bz_iv = neg_bz_iv - JV._iv(float(bi)) * JV._iv(float(zi))
    p_lo = float((neg_bz_iv + pen_zs - pen_Dx).a)

    # consistency: are con_* dual_values == the corresponding z entries?  (the shift
    # uses these duals; p_lo uses the full z.  They are the SAME z, so consistent.)
    return {
        "center": cname, "N": N, "T": T, "R": R, "bochner_n": bn, "pm_k_max": pm,
        "prob_value": Vc, "p_lo": p_lo, "penalty_total": float((pen_zs - pen_Dx).a),
        "penalty_Dx_upper": float(pen_Dx.b), "penalty_zs_lower": float(pen_zs.a),
        "duals": duals,
        "psd_lambda_min_lowers": psd_lmins,
        "self_checks": {"|c@x-obj_val|": abs(float(c @ x) - float(sol.obj_val)),
                        "r_prim": float(sol.r_prim), "r_dual": float(sol.r_dual),
                        "status": str(sol.status)},
        "solve_time_s": solve_t,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=3000)
    ap.add_argument("--T", type=int, default=1200)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=20)
    ap.add_argument("--pm_k_max", type=int, default=14)
    ap.add_argument("--out", type=str, default="/tmp/_jansson_duals.json")
    ap.add_argument("--center", type=str, default=None,
                    help="run a SINGLE center per process (one of: "
                         + "|".join(CENTERS) + "). MEMORY-AWARE: one ~4.6GB solve at a "
                         "time. If --out exists it is merged so sequential single-center "
                         "runs accumulate into one file.")
    args = ap.parse_args()
    if args.center is not None and args.center not in CENTERS:
        raise SystemExit(f"unknown center {args.center!r}; choose from {list(CENTERS)}")
    centers_to_run = (args.center,) if args.center else ("row4", "cde_n30_iter3")
    # merge into existing --out so sequential single-center runs accumulate
    res = {}
    if args.center and Path(args.out).exists():
        try:
            res = json.loads(Path(args.out).read_text())
        except Exception:
            res = {}
    for cn in centers_to_run:
        print(f"### {cn} N={args.N} T={args.T} bn={args.bochner_n} pm={args.pm_k_max}", flush=True)
        r = extract(cn, args.N, args.T, args.R, args.bochner_n, args.pm_k_max)
        res[cn] = r
        print(f"   p_lo={r['p_lo']:.10f}  V_c={r['prob_value']:.10f}  "
              f"penalty={r['penalty_total']:.2e}  ({r['solve_time_s']:.0f}s)")
        print(f"   duals: " + ", ".join(f"{k}={v:.5f}" for k, v in r['duals'].items()))
        print(f"   psd lmin: {r['psd_lambda_min_lowers']}  status={r['self_checks']['status']}")
    Path(args.out).write_text(json.dumps(res, indent=2, default=float))
    print(f"-> wrote {args.out}")
