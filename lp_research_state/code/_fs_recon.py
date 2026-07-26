"""
PRO-38 — RECONSTRUCT the current best FULL-SPACE certified floor (pure evaluation).

NO heavy SDP solves. Loads the UNION of all saved dual-FEASIBLE Phi-centers and
re-runs the Stage-1 Phi-cover evaluator (_fullspace_eval.cover_min_over_box) over
White's 18 Table-2 outside regions, with the corrected mside_sin_coeff=4.0 baked
into path_b_analytical (verified by reproduce_core_headline ~= 0.380284).

CENTER UNION (dedup by label):
  (A) 12 CORE anchors      -- parallel_results/phase5_N20K_bn40_dualext.json
                              (single-point centers h_c,h_c,p_c,p_c,q1,q2; anchor
                               = primal-1e-5, conservative; the core headline set)
  (B) 11 HALO centers      -- parallel_results/fullspace_stage2_halo_centers.json
                              (genuine single-point Phi-centers solved at h_c,h_c,
                               p_c,p_c,qlo,qhi; producer set primal:=dual_lb so
                               anchor primal-1e-5 = dual_lb-1e-5, conservative)
  (C) stage2 box-LP leaves -- parallel_results/fullspace_stage2_centers.json
                              CLEARED leaves ONLY (verdict=='cleared'; the other
                              114 leaves are infeasible -> empty boxes, all-zero
                              duals, primal=inf, NOT centers). See VALIDITY note.

VALIDITY of treating a box-LP leaf as a single-point Phi-center
--------------------------------------------------------------
build_problem_with_dual_handles attaches the duals to a program whose RHS are:
  con_53 : h1 (LOWER h)      con_54 : 2/3 + h2^2/2 (UPPER h)
  con_512_pL : p1   con_512_pU : p2   con_512_qL : q1   con_512_qU : q2
  con_513 : -0.5*(max(p1^2,p2^2) + max(q1^2,q2^2))
For a SINGLE-POINT center h1=h2=h_c, p1=p2=p_c, so dual_objective_shift (which uses
ONE h_c, ONE p_c) is exactly consistent. For a BOX leaf (h1<h2, p1<p2, q1<q2) the
shift formula's single-(h_c,p_c) assumption does NOT match the program the duals
came from, so Phi reconstructed at the box midpoint is NOT a validated global LB.
We therefore EMPIRICALLY TEST each box leaf: rebuild Phi with center = the box's
(h1->h_c via h1, etc.) and check whether anchor + shift recovers the leaf's
dual_LB AT the box-defining params. Only leaves passing the self-consistency test
(|Phi(params) - dual_LB| < tol) are admitted as Phi-centers; the rest are reported
but EXCLUDED from the cover (we never feed unvalidated duals into a global LB).

Independent of the Phi-cover, the stage2 box-LP track ALSO certifies regions by its
own (stronger) mechanism: a region is box-LP-certified >= TARGET iff every leaf is
infeasible OR dual_LB-margin >= TARGET. We report that verdict separately so the
final floor can use the best available mechanism per region WITHOUT conflating the
two. EVERY claim distinguishes 'independently certified by ours' from 'reliant on
White's published 0.38'.
"""
from __future__ import annotations
import json, os, sys, warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

from path_b_analytical import dual_objective_shift
import _fullspace_eval as FE
from _fullspace_eval import (
    cover_min_over_box, reproduce_core_headline, anchor_value,
    WHITE_TABLE2, CORE_516, WHITE_OUTSIDE_FLOOR, CORE_HEADLINE,
    WHITE_STRIP_BOUND, phi_center,
)

PR = CODE.parent / "parallel_results"
# $LP_DUALEXT overrides, as in `_fullspace_eval.load_centers` -- so pointing the
# pipeline at a re-anchored core needs one env var, not an edit per module.
DUALEXT = Path(os.environ.get("LP_DUALEXT")
               or PR / "phase5_N20K_bn40_dualext.json")
DUALEXT_FALLBACK = PR / "cde_phase5_corrected_tail.json"
HALO = PR / "fullspace_stage2_halo_centers.json"
STAGE2 = PR / "fullspace_stage2_centers.json"

NEED_KEYS = ("con_53", "con_54", "con_512_pL", "con_512_pU",
             "con_512_qL", "con_512_qU", "con_513")


# ---------------------------------------------------------------------------
# loaders -> uniform center dict {label,h_c,p_c,q1,q2,primal,duals}
# anchor_value(c,'primal_m1e5') = primal - 1e-5  (so set primal := conservative anchor + 1e-5)
# ---------------------------------------------------------------------------
def _valid_duals(d):
    return (isinstance(d, dict) and all(k in d for k in NEED_KEYS)
            and all(d[k] is not None for k in NEED_KEYS))


def load_core():
    src = DUALEXT if DUALEXT.exists() else DUALEXT_FALLBACK
    data = json.load(open(src))
    out = []
    for c in data["centers"]:
        # corrected_tail format uses V_c / V_c_rigorous instead of primal/dual_lb
        primal = c.get("primal")
        if primal is None:
            primal = c.get("V_c")
        out.append({"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                    "q1": c["q1"], "q2": c["q2"], "primal": primal,
                    "duals": {k: float(c["duals"][k]) for k in NEED_KEYS},
                    "src": src.name, "kind": "core_singlepoint"})
    return out, src.name


def load_halo():
    if not HALO.exists():
        return []
    data = json.load(open(HALO))
    out = []
    for c in data.get("centers", []):
        if not _valid_duals(c.get("duals")):
            continue
        # halo producer convention: conservative anchor = dual_lb - 1e-5.
        # anchor_value(.,'primal_m1e5') subtracts 1e-5 from 'primal', so we set
        # primal := dual_lb (then anchor = dual_lb - 1e-5). Exactly matches the
        # producer's own `c['primal'] = c['dual_lb']` patch.
        anchor_primal = c.get("dual_lb", c.get("primal"))
        out.append({"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                    "q1": c["q1"], "q2": c["q2"], "primal": anchor_primal,
                    "duals": {k: float(c["duals"][k]) for k in NEED_KEYS},
                    "src": HALO.name, "kind": "halo_singlepoint"})
    return out


def load_stage2_leaves():
    """Return CLEARED box-LP leaves as candidate centers.

    box = [h1,h2,p1,p2,q1,q2]. We expose BOTH the box-defining params (for the
    self-consistency test) and a candidate single-point representation:
      h_c := box midpoint in h, p_c := box midpoint in p, q-range := [q1,q2].
    anchor (conservative) = dual_LB - 1e-5 (the leaf's certify margin). primal is
    set so anchor_value(.,'primal_m1e5') = dual_LB - 1e-5.
    """
    if not STAGE2.exists():
        return []
    data = json.load(open(STAGE2))
    out = []
    for reg in data.get("regions", []):
        rid = reg.get("region")
        for lf in reg.get("leaves", []):
            if lf.get("verdict") != "cleared":
                continue
            if lf.get("dual_LB") is None or not _valid_duals(lf.get("duals")):
                continue
            b = lf["box"]  # [h1,h2,p1,p2,q1,q2]
            h1, h2, p1, p2, q1, q2 = b
            out.append({
                "label": f"s2_R{rid}_h{h1:.3f}_{h2:.3f}_p{p1:.3f}_{p2:.3f}_q{q1:.3f}_{q2:.3f}",
                "region": rid,
                "box": b,
                # candidate single-point representation (midpoints)
                "h_c": 0.5 * (h1 + h2), "p_c": 0.5 * (p1 + p2),
                "q1": q1, "q2": q2,
                "primal": lf["dual_LB"],          # anchor = dual_LB - 1e-5
                "dual_LB": lf["dual_LB"], "dual_resid": lf.get("dual_resid"),
                "duals": {k: float(lf["duals"][k]) for k in NEED_KEYS},
                "src": STAGE2.name, "kind": "stage2_box_leaf",
            })
    return out


def _box_dep_terms(duals, h53, h54, p_pL, p_pU, q_qL, q_qU, maxp2, maxq2):
    """The (h,p,q)-dependent part of the dual objective D for fixed duals at the
    given RHS values. The RHS-independent dual mass K cancels in all differences."""
    return ( duals["con_53"] * h53
             - duals["con_54"] * (2.0 / 3 + h54 ** 2 / 2)
             + duals["con_512_pL"] * p_pL
             - duals["con_512_pU"] * p_pU
             + duals["con_512_qL"] * q_qL
             - duals["con_512_qU"] * q_qU
             + duals["con_513"] * (-0.5 * (maxp2 + maxq2)) )


def box_leaf_validity_offset(leaf):
    """Rigorous validity test for using a box-LP leaf's duals as a single-point
    Phi-center anchored at its dual_LB.

    A box leaf solved the augmented LP with RHS: con_53=h1, con_54 uses h2,
    con_512_pL=p1, pU=p2, qL=q1, qU=q2, con_513 uses max(p1^2,p2^2)+max(q1^2,q2^2).
    Its dual_LB is the dual objective at those (box) RHS constants, i.e. dual_LB ~=
    D_box + K. The duals are a feasible dual point (White Appendix II: params enter
    ONLY the objective, never dual feasibility), so D(query)+K is a valid LB on mu
    at EVERY query (h,p,q).

    We register the center at the box midpoint and reconstruct Phi via
    dual_objective_shift, which computes D_sp(query) - D_sp(center) where D_sp treats
    (h_c,p_c) as a single point. Hence for any query:
        Phi(query) - TrueD(query)
          = [anchor + D_sp(query) - D_sp(center)] - [D_sp(query) + K]
          = anchor - K - D_sp(center)
          = (dual_LB - margin) - K - D_mid_true
          = (D_box - D_mid_true) - margin          [since dual_LB ~= D_box + K]
    This is a CONSTANT (independent of query: the D_sp(query) terms cancel exactly).
    Phi is a valid global LB iff this constant <= 0, i.e. offset := D_box - D_mid_true
    <= 0. We return offset; admit the leaf iff offset <= +tol (we keep tol tiny; the
    -1e-5 anchor margin gives additional slack we do NOT rely on for the sign).
    """
    du = leaf["duals"]
    h1, h2, p1, p2, q1, q2 = leaf["box"]
    maxp2 = max(p1 ** 2, p2 ** 2); maxq2 = max(q1 ** 2, q2 ** 2)
    D_box = _box_dep_terms(du, h1, h2, p1, p2, q1, q2, maxp2, maxq2)
    hc = 0.5 * (h1 + h2); pc = 0.5 * (p1 + p2)
    # single-point reconstruction's reference (q range kept as the box's [q1,q2]):
    D_mid = _box_dep_terms(du, hc, hc, pc, pc, q1, q2, maxp2, maxq2)
    return D_box - D_mid


def test_box_leaf_consistency(leaf, tol=1e-9):
    """Admit a box leaf as a valid conservative Phi-center iff its validity offset
    (D_box - D_mid_true) <= tol (i.e. <= 0). Returns (offset, admit)."""
    offset = box_leaf_validity_offset(leaf)
    return offset, (offset <= tol)


def width_class(hr, pr, qr):
    hw = hr[1] - hr[0]; pw = pr[1] - pr[0]; qw = qr[1] - qr[0]
    if hw > 0.3 or qw > 0.3:
        return "wide"
    # narrow: every spanned dim small & near core
    near_core = (hr[0] >= 0.0 and hr[1] <= 0.08
                 and pr[0] >= 0.30 and pr[1] <= 0.50
                 and abs(qr[0]) <= 0.06 and abs(qr[1]) <= 0.06)
    small = (hw <= 0.08 and pw <= 0.15 and qw <= 0.06)
    if small and near_core:
        return "narrow"
    return "moderate"


def stage2_region_verdict():
    """Per-region box-LP certification verdict from the stage2 file (independent of
    Phi). Returns {rid: {certified, min_cleared_cert, n_cleared, n_infeasible,
    n_residual}}."""
    if not STAGE2.exists():
        return {}
    data = json.load(open(STAGE2))
    out = {}
    for reg in data.get("regions", []):
        out[reg["region"]] = {
            "box_lp_certified_ge_target": bool(reg.get("certified_ge_target")),
            "target": reg.get("target"),
            "min_cleared_cert": reg.get("min_cleared_cert"),
            "n_cleared": reg.get("n_cleared"), "n_infeasible": reg.get("n_infeasible"),
            "n_residual": reg.get("n_residual"),
        }
    return out


def main():
    print("=" * 78)
    print("PRO-38 full-space floor RECONSTRUCTION (pure eval of saved duals)")
    print("=" * 78)

    core, core_src = load_core()
    halo = load_halo()
    leaves = load_stage2_leaves()
    print(f"\n[load] core single-point centers : {len(core):3d}  (from {core_src})")
    print(f"[load] halo single-point centers : {len(halo):3d}  (from {HALO.name})")
    print(f"[load] stage2 CLEARED box-leaves : {len(leaves):3d}  (from {STAGE2.name})")

    # ---- test box leaves for admissibility as conservative Phi-centers
    admitted_leaves, rejected_leaves = [], []
    for lf in leaves:
        resid, ok = test_box_leaf_consistency(lf)
        lf["_phi_at_center_resid"] = resid
        (admitted_leaves if ok else rejected_leaves).append(lf)
    if leaves:
        offs = [lf["_phi_at_center_resid"] for lf in leaves]
        print(f"[test] box-leaf validity offset (D_box - D_mid_true) must be <=0 "
              f"to be a conservative Phi-center:")
        print(f"       admitted {len(admitted_leaves)}/{len(leaves)} "
              f"(worst offset = {max(offs):.3e}, all<=0 => {max(offs) <= 1e-9})")

    # ---- build the UNION, dedupe by label ----
    union = {}
    sources_present = []
    for grp, tag in [(core, "core"), (halo, "halo"), (admitted_leaves, "stage2_leaves")]:
        if grp:
            sources_present.append(tag)
        for c in grp:
            union.setdefault(c["label"], c)
    centers = list(union.values())
    center_sources = sorted({c["src"] for c in centers})
    print(f"\n[union] total deduped centers : {len(centers)}  "
          f"(core={sum(c['kind']=='core_singlepoint' for c in centers)}, "
          f"halo={sum(c['kind']=='halo_singlepoint' for c in centers)}, "
          f"stage2_leaves={sum(c['kind']=='stage2_box_leaf' for c in centers)})")
    print(f"[union] sources: {center_sources}")

    # ============================================================
    # (2) sanity: reproduce core headline (CANONICAL = 12 core centers)
    # ============================================================
    # The canonical core reproduction is the 12-center cover over the core box
    # (the project headline 0.380284). Box leaves have large gradients and only
    # INFLATE the Lipschitz cell-error there (they don't belong to the single-point
    # core cover); we therefore reproduce the headline with the core centers and
    # report the union value separately for transparency.
    core_only = [c for c in centers if c["kind"] == "core_singlepoint"]
    core_rep = reproduce_core_headline(core_only, "primal_m1e5")
    core_ok = abs(core_rep["rigorous_LB"] - CORE_HEADLINE) < 5e-5
    print("\n" + "-" * 78)
    print("[CORE 5.16] reproduce_core_headline (12 core centers, 4001x4001, coeff 4.0)")
    print(f"   rigorous_LB = {core_rep['rigorous_LB']:.10f}  "
          f"(grid_min={core_rep['grid_min']:.8f}, eps_grid={core_rep['eps_grid']:.2e})")
    print(f"   vs CORE_HEADLINE {CORE_HEADLINE}: {core_rep['rigorous_LB']-CORE_HEADLINE:+.2e}"
          f"  -> {'REPRODUCED' if core_ok else 'CHECK'}")
    print(f"   binding @ {core_rep['binding_point']} witness={core_rep['witness']}")
    core_rep_union = reproduce_core_headline(centers, "primal_m1e5")
    print(f"   (union-of-all LB = {core_rep_union['rigorous_LB']:.10f}; "
          f"union vs core-only delta = {core_rep_union['rigorous_LB']-core_rep['rigorous_LB']:+.2e} "
          f"-- box-leaf Lipschitz inflation, not used for headline)")

    # ============================================================
    # (3) per-region best INDEPENDENT floor over each of the 18 Table-2 regions
    # ============================================================
    # Three mutually-independent, each-rigorous lower-bound mechanisms exist per
    # region; the best certified floor is the MAX of those that apply (max of valid
    # LBs is a valid LB). Reported per region as ours_phi_min:
    #   (a) Phi-cover over CORE+HALO single-point centers (the canonical Phi track;
    #       NOT degraded by box-leaf Lipschitz inflation),
    #   (b) Phi-cover over the FULL UNION (helps the wide regions where box leaves
    #       dominate; can be looser than (a) on near-core regions due to the global
    #       L_max inflation -- so we keep both and take the max),
    #   (c) the STAGE-2 BOX-LP certificate (independent, stronger mechanism: a region
    #       is box-LP-certified >= 0.380284 iff every leaf is infeasible or its
    #       dual_LB-margin >= 0.380284; we use min_cleared_cert as the rigorous floor
    #       when certified).
    # ours_phi_min := max(a, b, c_if_certified). This is the region's INDEPENDENTLY
    # certified floor (NO White number). white_reliant := max(ours_phi_min, White).
    core_halo = [c for c in centers if c["kind"] in ("core_singlepoint", "halo_singlepoint")]
    s2v = stage2_region_verdict()

    def cover(cset, hr, pr, qr):
        n_h = 81 if (hr[1]-hr[0]) > 0.05 else 41
        n_p = 161 if (pr[1]-pr[0]) > 0.2 else 81
        n_q = 81 if (qr[1]-qr[0]) > 0.1 else 41
        return cover_min_over_box(cset, "primal_m1e5", hr, pr, qr,
                                  n_h=n_h, n_p=n_p, n_q=n_q)

    per_region = []
    gate_regions = []
    indep_floor = np.inf
    white_reliant_floor = np.inf

    print("\n" + "-" * 78)
    print("PER-REGION best INDEPENDENT floor = max(corehalo-Phi, union-Phi, box-LP cert)")
    print("-" * 78)
    for (hr, pr, qr, wbound) in WHITE_TABLE2:
        idx = WHITE_TABLE2.index((hr, pr, qr, wbound)) + 1
        wc = width_class(hr, pr, qr)

        lb_ch, pt_ch, wit_ch, gm_ch, eps_ch, Lm_ch = cover(core_halo, hr, pr, qr)
        lb_un, pt_un, wit_un, gm_un, eps_un, Lm_un = cover(centers, hr, pr, qr)

        s2 = s2v.get(idx)
        box_lp_cert = bool(s2 and s2["box_lp_certified_ge_target"])
        box_lp_floor = (float(s2["min_cleared_cert"]) if (box_lp_cert and
                        s2.get("min_cleared_cert") is not None) else None)

        # best independent Phi (max of the two cover variants) and its provenance
        if lb_ch >= lb_un:
            phi_best, phi_pt, phi_wit, phi_src = float(lb_ch), pt_ch, wit_ch, "corehalo"
        else:
            phi_best, phi_pt, phi_wit, phi_src = float(lb_un), pt_un, wit_un, "union"

        # best independent floor overall (Phi-best vs box-LP certificate)
        ours = phi_best
        mech = f"phi_{phi_src}"
        worst_corner = {"h": phi_pt[0], "p": phi_pt[1], "q": phi_pt[2]}
        worst_wit = phi_wit
        if box_lp_floor is not None and box_lp_floor > ours:
            ours = box_lp_floor
            mech = "box_lp"
            worst_corner = {"h": None, "p": None, "q": None}  # box-LP min over leaves
            worst_wit = "stage2_box_lp_min_cleared_cert"

        white_floor = WHITE_STRIP_BOUND if abs(wbound - WHITE_STRIP_BOUND) < 1e-9 \
            else WHITE_OUTSIDE_FLOOR
        white_reliant = max(ours, white_floor)
        indep_floor = min(indep_floor, ours)
        white_reliant_floor = min(white_reliant_floor, white_reliant)

        rec = {
            "region": idx, "width_class": wc,
            "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
            "ours_phi_min": ours,                       # best independent floor
            "ours_mechanism": mech,
            "phi_corehalo_min": float(lb_ch), "phi_union_min": float(lb_un),
            "box_lp_certified": box_lp_cert, "box_lp_floor": box_lp_floor,
            "worst_corner": worst_corner, "worst_witness": worst_wit,
            "phi_corehalo_worst": {"h": pt_ch[0], "p": pt_ch[1], "q": pt_ch[2],
                                   "witness": wit_ch, "L_max": float(Lm_ch),
                                   "eps_grid": float(eps_ch)},
            "white_floor": white_floor,
            "white_reliant_floor": white_reliant,
            "clears_380000_indep": bool(ours >= WHITE_OUTSIDE_FLOOR),
            "clears_380284_indep": bool(ours >= CORE_HEADLINE),
            "stage2_box_lp": s2,
        }
        per_region.append(rec)

        tag = f"  [{mech}]"
        if box_lp_cert:
            tag += " box-LP-CERT>=0.380284"
        if ours < CORE_HEADLINE:
            gate_regions.append({
                "region": idx,
                "worst_corner": worst_corner,
                "shortfall_380284": float(CORE_HEADLINE - ours),
                "shortfall_380000": float(max(0.0, WHITE_OUTSIDE_FLOOR - ours)),
                "width_class": wc,
                "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
                "stage2_box_lp_certified": box_lp_cert,
                "ours_mechanism": mech,
            })
        print(f"[R{idx:2d}] {wc:8s} h{hr} p{pr} q{qr}")
        print(f"      corehalo-Phi={lb_ch:.6f}  union-Phi={lb_un:.6f}  "
              f"box-LP={'%.6f' % box_lp_floor if box_lp_floor is not None else 'n/a':>8}"
              f"  -> ours={ours:.6f}{tag}")
        print(f"      white_reliant=max(ours,{white_floor:.5f})={white_reliant:.6f}"
              f"   indep_clears_0.380284={rec['clears_380284_indep']}"
              f"  indep_clears_0.380000={rec['clears_380000_indep']}\n")

    # ============================================================
    # (4) summary
    # ============================================================
    print("=" * 78)
    print("FULL-SPACE FLOOR SUMMARY")
    print("=" * 78)
    core_floor = core_rep["rigorous_LB"]
    indep_floor_with_core = min(indep_floor, core_floor)
    white_reliant_with_core = min(white_reliant_floor, core_floor)
    print(f"core (5.16) Phi floor (12 centers) : {core_floor:.7f}")
    print(f"INDEPENDENTLY-certified floor (ours only, NO White number anywhere):")
    print(f"   min over 18 regions of ours_phi_min = {indep_floor:.7f}")
    print(f"   combined with core                  = {indep_floor_with_core:.7f}")
    print(f"WHITE-RELIANT floor (max(ours, White 0.380000; 0.37925 for R18)):")
    print(f"   min over 18 regions                 = {white_reliant_floor:.7f}")
    print(f"   combined with core                  = {white_reliant_with_core:.7f}")
    gate_ids = [g["region"] for g in gate_regions]
    print(f"\nGATE regions (ours independent floor < {CORE_HEADLINE}): {gate_ids}")
    box_cert = [r["region"] for r in per_region if r["box_lp_certified"]]
    print(f"stage2 box-LP independently certifies >=0.380284: regions {box_cert}")
    indep_clear_380 = [r["region"] for r in per_region if not r["clears_380000_indep"]]
    print(f"regions whose INDEPENDENT floor is still < 0.380000 (rely on White): "
          f"{indep_clear_380}")

    out = {
        "centers_loaded": len(centers),
        "center_sources": center_sources,
        "n_core": sum(c["kind"] == "core_singlepoint" for c in centers),
        "n_halo": sum(c["kind"] == "halo_singlepoint" for c in centers),
        "n_stage2_leaves": sum(c["kind"] == "stage2_box_leaf" for c in centers),
        "core_reproduced_LB": core_floor,
        "core_ok": core_ok,
        "independently_certified_floor": float(indep_floor_with_core),
        "independently_certified_floor_regions_only": float(indep_floor),
        "white_reliant_floor": float(white_reliant_with_core),
        "white_reliant_floor_regions_only": float(white_reliant_floor),
        "per_region": per_region,
        "gate_regions": gate_regions,
        "stage2_box_lp_certified_regions": box_cert,
        "regions_indep_below_380000": indep_clear_380,
        "rejected_box_leaves": len(rejected_leaves),
    }
    OUT = PR / "fullspace_recon.json"
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nsaved -> {OUT}")
    return out


if __name__ == "__main__":
    main()
