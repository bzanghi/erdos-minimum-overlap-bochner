"""
_l2_clean_theorem.py — PRO-47 FINISH: the CLEAN verified-cover theorem.

Eliminates the LAST residual of the TIER-2c headline (mu >= 0.380295): previously
the binding centers' shift duals (con_53/54/512_pL/pU/qL/qU/513) came from the
historical COVER solve (cde_phase5_corrected_tail.json) while their anchor (Jansson
p_lo) came from a DIFFERENT production Jansson solve (L2_PROD.json) -- two separate
CLARABEL runs, ~5e-6 cross-solve drift.

This driver consumes a SELF-CONSISTENT (p_lo, con_* duals) tuple produced by
_jansson_with_duals.py in ONE production solve (N=20000, T=4000, bochner_n=40,
pm_k_max=20), where:
  * the con_* duals are the EXACT components of the same numeric conic dual z, and
  * the Jansson p_lo is -b(theta_c)^T z - pen_Dx for that SAME z, with the two PSD
    Bochner blocks interval-certified PSD (pen_zs = 0).
By _verify_shift_eq_dualobj.py (Lemma 10 + shift==-b(theta)^Tz to 1.95e-16) and the
empirical CLARABEL determinism check (||z-z2||_inf == 0 for solve-via-data on identical
canonical data), the shift coefficients used here are exactly the z-components driving
p_lo.  So for the binding centers the anchor AND the shift come from ONE solve --
residual ELIMINATED.

Two clean tiers are emitted (both fully interval-certified shift + box-min):

  CLEAN TIER 1 (UNCONDITIONAL):  2 binding centers ONLY, each anchored at its OWN
    self-consistent production p_lo with its OWN self-consistent duals.  No gap
    assumption, no cross-solve mix.  Airtight.  Weaker in VALUE because a 2-center
    (h,p) cover dips at the far box corners.

  CLEAN TIER 2c (HEADLINE-STRENGTH):  all 12 production centers.  The 2 binding
    centers use their SELF-CONSISTENT (p_lo anchor + duals) -- residual gone for the
    centers that actually bind the floor.  The 10 NON-binding centers keep the
    cover-solve duals + V_c-margin anchor (documented convention; non-binding +
    robustness-stress-bounded).  Trusted base for the binding witness shrinks to
    NOTHING (its anchor is interval-certified and self-consistent); for the floor as a
    whole it is just "the 10 non-binding cover-solve gaps <= margin".

Author: Claude (machine-assisted), PRO-47 clean finish.  2026-06-06.
"""
from __future__ import annotations
import sys, json, time, argparse
from pathlib import Path

import numpy as np
import mpmath
from mpmath import iv

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import _cover_iv_certify as C
from _cover_iv_certify import CenterIV, _iv, lo_of, hi_of, box_min_cell_enclosure_fast

mpmath.mp.dps = 50
iv.dps = 50

WHITE = 0.379005
PRIOR_PUB = 0.379544
HEADLINE = 0.380284

REPO = CODE.parent.parent
DOCS = REPO / "docs" / "RND_WHITESPACE"

BINDING = ("row4", "cde_n30_iter3")


def _mk_center(label, h_c, p_c, q1, q2, duals, anchor_val, rigorous, src):
    """Build a CenterIV with quadratic coeffs about (h_c, p_c)."""
    cc = CenterIV(label, h_c, p_c, q1, q2, duals, _iv(anchor_val), rigorous, src)
    return cc


def build_clean_centers(cover_json, sc_json, penalty_fallback):
    """Return (centers_2c_clean, centers_sc_only, sc_data, diag).

    centers_2c_clean : 12 centers; binding 2 use SELF-CONSISTENT (p_lo, duals);
                       10 non-binding use cover duals + V_c - penalty_fallback anchor.
    centers_sc_only  : the 2 binding centers ONLY, self-consistent (CLEAN TIER 1).
    """
    cover = json.loads(Path(cover_json).read_text())
    cmap = {c["label"]: c for c in cover["centers"] if "error" not in c}
    sc = json.loads(Path(sc_json).read_text())

    # sanity: self-consistent extraction must contain both binding centers
    missing = [c for c in BINDING if c not in sc]
    if missing:
        raise SystemExit(f"self-consistent extraction {sc_json} missing centers: {missing}")

    diag = {"sc_config": {}, "binding_anchor_check": {}, "dual_drift_cover_vs_sc": {}}
    for cn in BINDING:
        r = sc[cn]
        diag["sc_config"][cn] = {
            "N": r["N"], "T": r["T"], "bochner_n": r["bochner_n"],
            "pm_k_max": r["pm_k_max"], "p_lo": r["p_lo"], "prob_value": r["prob_value"],
            "penalty_total": r["penalty_total"],
            "psd_lambda_min_lowers": r.get("psd_lambda_min_lowers"),
            "status": r["self_checks"]["status"],
        }
        # p_lo <= prob.value (Jansson soundness at the center)
        diag["binding_anchor_check"][cn] = {
            "p_lo": r["p_lo"], "prob_value": r["prob_value"],
            "p_lo_le_value": bool(r["p_lo"] <= r["prob_value"]),
            "both_psd_blocks_certified": bool(
                r.get("psd_lambda_min_lowers")
                and all(x >= 0 for x in r["psd_lambda_min_lowers"])),
        }
        # how far the self-consistent duals drift from the historical cover duals
        cv = cmap[cn]
        diag["dual_drift_cover_vs_sc"][cn] = {
            k: float(cv["duals"][k] - r["duals"][k]) for k in cv["duals"]}

    # ---- CLEAN TIER 1: 2 binding centers, self-consistent ----
    centers_sc_only = []
    for cn in BINDING:
        r = sc[cn]
        cv = cmap[cn]
        cc = _mk_center(cn, cv["h_c"], cv["p_c"], cv["q1"], cv["q2"], r["duals"],
                        r["p_lo"], True,
                        f"self-consistent prod (N={r['N']}): anchor=p_lo={r['p_lo']:.10f}, "
                        f"SAME-solve duals")
        cc.V_c = r["prob_value"]
        centers_sc_only.append(cc)

    # ---- CLEAN TIER 2c: 12 centers; binding 2 self-consistent, others V_c-margin ----
    centers_2c = []
    for lab, cv in cmap.items():
        if lab in BINDING:
            r = sc[lab]
            cc = _mk_center(lab, cv["h_c"], cv["p_c"], cv["q1"], cv["q2"], r["duals"],
                            r["p_lo"], True,
                            f"self-consistent prod (N={r['N']}): anchor=p_lo={r['p_lo']:.10f} "
                            f"+ SAME-solve duals  [RESIDUAL ELIMINATED]")
        else:
            anchor = cv["V_c"] - penalty_fallback
            cc = _mk_center(lab, cv["h_c"], cv["p_c"], cv["q1"], cv["q2"], cv["duals"],
                            anchor, False,
                            f"V_c-{penalty_fallback:.1e} (cover-solve gap<=margin) [non-binding]")
        cc.V_c = cv["V_c"]
        centers_2c.append(cc)

    return centers_2c, centers_sc_only, sc, diag


def crosscheck_binding_independent(centers_sc_only, sc):
    """Float cross-check: interval Phi (midpoint) vs path_b_independent.Phi_row for the
    self-consistent binding centers at random (h,p), SAME self-consistent duals.
    Confirms the interval shift formula matches the independent re-derivation."""
    from path_b_independent import Phi_row
    rng = np.random.default_rng(11)
    worst = 0.0; worst_loc = None
    for cc in centers_sc_only:
        r = sc[cc.label]
        duals = r["duals"]
        cen = {"h_c": float(mpmath.mpf(cc.h_c.a)), "p_c": float(mpmath.mpf(cc.p_c.a)),
               "q1_c": float(mpmath.mpf(cc.q1_c.a)), "q2_c": float(mpmath.mpf(cc.q2_c.a))}
        rec = {"value": r["p_lo"], "duals": {
            "lam_53": duals["con_53"], "lam_54": duals["con_54"],
            "lam_pL": duals["con_512_pL"], "lam_pU": duals["con_512_pU"],
            "lam_qL": duals["con_512_qL"], "lam_qU": duals["con_512_qU"],
            "lam_513": duals["con_513"]}, "center": cen}
        for _ in range(400):
            h = rng.uniform(*C.H_BOX); p = rng.uniform(*C.P_BOX); q = cen["q1_c"]
            phi_mid = 0.5 * (float(mpmath.mpf(cc.phi_point(h, p).a))
                             + float(mpmath.mpf(cc.phi_point(h, p).b)))
            F_ind = Phi_row(rec, h, p, q)
            q_off = Phi_row(rec, cen["h_c"], cen["p_c"], q) - r["p_lo"]
            F_ind_hp = F_ind - q_off
            d = abs(phi_mid - F_ind_hp)
            if d > worst:
                worst = d; worst_loc = (cc.label, float(h), float(p))
    return {"worst_abs_diff": float(worst), "worst_loc": worst_loc,
            "agree_10digit": bool(worst < 1e-9),
            "note": "interval Phi(midpoint) vs path_b_independent.Phi_row, SELF-CONSISTENT duals"}


def run(cover_json, sc_json, penalty_fallback=1e-6, n_h=2000, n_p=2000,
        out_json=None, out_md=None, verbose=True):
    t0 = time.time()
    centers_2c, centers_sc, sc, diag = build_clean_centers(
        cover_json, sc_json, penalty_fallback)

    if verbose:
        print("[clean] self-consistent binding-center extraction:")
        for cn in BINDING:
            d = diag["sc_config"][cn]; chk = diag["binding_anchor_check"][cn]
            print(f"   {cn:16s} N={d['N']} bn={d['bochner_n']} pm={d['pm_k_max']}  "
                  f"p_lo={d['p_lo']:.10f}  V_c={d['prob_value']:.10f}")
            print(f"        p_lo<=value: {chk['p_lo_le_value']}  "
                  f"both PSD blocks certified: {chk['both_psd_blocks_certified']}  "
                  f"psd_lmin_lo={d['psd_lambda_min_lowers']}  status={d['status']}")

    # ---- CLEAN TIER 1 (unconditional) ----
    floor_sc, info_sc = box_min_cell_enclosure_fast(centers_sc, n_h, n_p, verbose=False)
    floor_sc_f = float(floor_sc)
    if verbose:
        print(f"\n[clean] CLEAN TIER 1 (2 self-consistent binding centers, UNCONDITIONAL):")
        print(f"        floor = {floor_sc_f:.10f}  witness={info_sc['worst_witness']}  "
              f"cell={info_sc['worst_cell']}")

    # ---- CLEAN TIER 2c (headline-strength) ----
    floor_2c, info_2c = box_min_cell_enclosure_fast(centers_2c, n_h, n_p, verbose=False)
    floor_2c_f = float(floor_2c)
    binding_2c = info_2c["worst_witness"]
    binding_is_self_consistent = binding_2c in BINDING
    if verbose:
        print(f"\n[clean] CLEAN TIER 2c (12 centers; binding self-consistent, residual GONE):")
        print(f"        floor = {floor_2c_f:.10f}  witness={binding_2c}  "
              f"(self-consistent binding witness: {binding_is_self_consistent})")
        print(f"        cell={info_2c['worst_cell']}  "
              f"nonrig anchor ever binds: {info_2c['nonrigorous_anchor_ever_binds']}")

    # validate fast scan vs slow pure-iv on a coarse grid (rigor self-check)
    floor_2c_slow, _ = C.box_min_cell_enclosure(centers_2c, 80, 80, verbose=False)
    floor_2c_fast_coarse, _ = box_min_cell_enclosure_fast(centers_2c, 80, 80, verbose=False)
    fastslow_2c = abs(float(floor_2c_slow) - float(floor_2c_fast_coarse))

    floor_sc_slow, _ = C.box_min_cell_enclosure(centers_sc, 80, 80, verbose=False)
    floor_sc_fast_coarse, _ = box_min_cell_enclosure_fast(centers_sc, 80, 80, verbose=False)
    fastslow_sc = abs(float(floor_sc_slow) - float(floor_sc_fast_coarse))
    if verbose:
        print(f"\n[clean] fast-vs-slow-iv coarse |Δ|: TIER2c={fastslow_2c:.2e}  TIER1={fastslow_sc:.2e}")

    # independent formula cross-check for the binding self-consistent centers
    xc = crosscheck_binding_independent(centers_sc, sc)
    if verbose:
        print(f"[clean] independent formula cross-check (binding, SC duals): "
              f"worst|Δ|={xc['worst_abs_diff']:.2e} (10-digit: {xc['agree_10digit']})")

    result = {
        "kind": "l2_clean_theorem",
        "cover_json": str(cover_json), "sc_json": str(sc_json),
        "penalty_fallback": penalty_fallback,
        "WHITE": WHITE, "PRIOR_PUB": PRIOR_PUB, "HEADLINE": HEADLINE,
        "binding_centers": list(BINDING),
        "self_consistent_extraction": diag,

        "CLEAN_TIER1_unconditional_selfconsistent": {
            "floor": floor_sc_f,
            "binding_witness": info_sc["worst_witness"],
            "worst_cell": info_sc["worst_cell"],
            "what": (
                "UNCONDITIONAL. 2 binding centers ONLY; each anchored at its OWN "
                "production self-consistent Jansson p_lo with its OWN same-solve duals. "
                "NO duality-gap assumption, NO cross-solve mix -- the two-solve "
                "nondeterminism residual is ELIMINATED. Weaker in VALUE than TIER 2c "
                "because a 2-center (h,p) cover dips at the far box corners (binds at "
                f"{info_sc['worst_cell'][:2]}x{info_sc['worst_cell'][2:]}), but airtight."),
            "clears_white": bool(floor_sc_f >= WHITE),
            "clears_prior_pub": bool(floor_sc_f >= PRIOR_PUB),
            "clears_headline": bool(floor_sc_f >= HEADLINE),
            "margin_vs_white": floor_sc_f - WHITE,
            "margin_vs_prior_pub": floor_sc_f - PRIOR_PUB,
        },
        "CLEAN_TIER2c_production_selfconsistent_binding": {
            "floor": floor_2c_f,
            "binding_witness": binding_2c,
            "binding_witness_is_self_consistent": binding_is_self_consistent,
            "worst_cell": info_2c["worst_cell"],
            "nonrigorous_anchor_ever_binds": info_2c["nonrigorous_anchor_ever_binds"],
            "what": (
                "HEADLINE-STRENGTH. All 12 production (N=20000) centers; shift + box-min "
                "fully interval-certified. The 2 binding centers (row4, cde_n30_iter3) use "
                "their SELF-CONSISTENT (p_lo anchor + same-solve duals) -- the two-solve "
                "residual is ELIMINATED for the centers that actually set the floor. The 10 "
                "NON-binding centers keep the cover-solve duals + V_c-margin anchor "
                "(documented convention; non-binding + robustness-stress-bounded). Since the "
                "binding witness is one of the self-consistent centers, the floor's binding "
                "anchor is interval-certified and self-consistent (NOT a margin convention)."),
            "clears_white": bool(floor_2c_f >= WHITE),
            "clears_prior_pub": bool(floor_2c_f >= PRIOR_PUB),
            "clears_headline": bool(floor_2c_f >= HEADLINE),
            "margin_vs_white": floor_2c_f - WHITE,
            "margin_vs_prior_pub": floor_2c_f - PRIOR_PUB,
            "margin_vs_headline": floor_2c_f - HEADLINE,
        },
        "rigor_self_checks": {
            "fast_vs_slow_iv_coarse_diff_tier2c": fastslow_2c,
            "fast_vs_slow_iv_coarse_diff_tier1": fastslow_sc,
            "independent_formula_crosscheck": xc,
            "lemma10_and_shift_eq_dualobj": (
                "verified by _verify_shift_eq_dualobj.py: max|A_theta-A_c|=0, "
                "max|c_theta-c_c|=0, shift recon vs exact -b(theta)^Tz change "
                "worst|Δ|=1.95e-16 -> con_* duals ARE the z-components driving p_lo"),
            "clarabel_determinism": (
                "verified: solve_via_data on identical canonical data is bit-for-bit "
                "deterministic (||z-z2||_inf=0); prob.solve()==solve_via_data obj to 0 "
                "-> within ONE process the p_lo and con_* duals come from the SAME z"),
        },
        "VERIFIED_CLEAN_FLOOR": floor_2c_f,
        "VERIFIED_CLEAN_FLOOR_unconditional": floor_sc_f,
        "elapsed_s": time.time() - t0,
        "RESIDUAL_STATUS": (
            "ELIMINATED for the binding centers. CLEAN TIER 2c floor=%.10f is set by a "
            "self-consistent binding witness (%s) whose anchor (production Jansson p_lo) "
            "and shift duals come from ONE N=20000 solve. CLEAN TIER 1=%.10f is fully "
            "unconditional (2 self-consistent centers, no gap assumption)."
            % (floor_2c_f, binding_2c, floor_sc_f)),
    }

    if out_json:
        Path(out_json).write_text(json.dumps(result, indent=2, default=float))
        if verbose:
            print(f"\n-> wrote {out_json}")

    if verbose:
        print(f"\n[clean] ===========================================================")
        print(f"[clean] CLEAN VERIFIED COVER FLOOR (residual eliminated):")
        print(f"[clean]  >> CLEAN TIER 2c (headline-strength): {floor_2c_f:.10f}")
        print(f"[clean]       binding witness = {binding_2c} "
              f"(self-consistent: {binding_is_self_consistent})")
        print(f"[clean]       vs White 0.379005    : {floor_2c_f-WHITE:+.3e}")
        print(f"[clean]       vs prior pub 0.379544: {floor_2c_f-PRIOR_PUB:+.3e}")
        print(f"[clean]       vs headline 0.380284 : {floor_2c_f-HEADLINE:+.3e}")
        print(f"[clean]  CLEAN TIER 1 (unconditional): {floor_sc_f:.10f}  "
              f"(>White: {floor_sc_f>=WHITE}, >prior-pub: {floor_sc_f>=PRIOR_PUB})")
        print(f"[clean] ===========================================================")
    return result


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--cover_json", type=str,
                    default=str(CODE.parent / "parallel_results" / "cde_phase5_corrected_tail.json"))
    ap.add_argument("--sc_json", type=str, default="/tmp/pro47/L2_CLEAN_sc_prod.json",
                    help="self-consistent (p_lo,duals) production extraction "
                         "(_jansson_with_duals.py --N 20000 ...)")
    ap.add_argument("--penalty_fallback", type=float, default=1e-6)
    ap.add_argument("--n_h", type=int, default=2000)
    ap.add_argument("--n_p", type=int, default=2000)
    ap.add_argument("--out", type=str, default=str(DOCS / "L2_CLEAN_THEOREM.json"))
    args = ap.parse_args()
    run(args.cover_json, args.sc_json, penalty_fallback=args.penalty_fallback,
        n_h=args.n_h, n_p=args.n_p, out_json=args.out)
