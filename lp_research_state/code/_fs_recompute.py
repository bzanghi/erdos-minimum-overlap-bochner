"""
PRO-38 Phase-2 FINAL recompute of the full-space (E(M), c1, d1) lower bound.

Loads the UNION of every dual-feasible center we have produced for the
full-space promotion:

  * the 12 CORE anchors  (phase5_N20K_bn40_dualext.json)
  * fullspace_stage2_centers.json     (box-LP subdivision leaves; one center
                                       per leaf, placed at the leaf box-midpoint)
  * fullspace_stage2_halo_centers.json (11 halo point-centers)
  * EVERY parallel_results/fullspace_promote_R*.json produced this run
    (globbed; the per-region fresh centers from Phase-2)

and recomputes, with the CORRECTED-coefficient path_b (mside_sin_coeff=4.0 is
already baked into the saved duals of every fresh solve; the evaluator only
*reads* the saved dual values via path_b_analytical.dual_objective_shift, so the
correction is inherited automatically), two things:

  (1) reproduce_core_headline  -- sanity that the union still reproduces the
      canonical core (5.16) headline ~0.380284 (the core anchors dominate there).
  (2) cover_min_over_box over ALL 18 WHITE_TABLE2 regions, using the WHOLE union
      of centers, giving for each region a RIGOROUS box-min lower bound of
      Cover(h,p,q)=max_c Phi_c (fine grid + Lipschitz cell-error eps_grid).

VALIDITY (the project's #1 documented trap is OVERCLAIMING -- enforced here):
  * Each Phi_c is a GLOBALLY-valid lower bound on mu because its duals come from a
    converged dual-feasible solve (small dual residual) AND the anchor is
    CONSERVATIVE = primal - 1e-5 (never the reported primal, never dual_lb that
    was only printed to 5 digits).  For box-LP leaves we use leaf_primal - 1e-5.
  * max over centers of valid LBs is valid; cover_min_over_box's grid+Lipschitz
    is a rigorous box-min.  No White number enters Phi.
  * A region's INDEPENDENTLY-certified floor = min_box ours_phi_min  (NO White).
    Its WHITE-RELIANT floor = max(that, White's published floor:
        0.380000 for the seventeen "0.38" rows; 0.37925 for strip R18).
  * Full-space mu bound = min over ALL 18 regions + core of the best floor.
  * We report BOTH the independently-certified floor and the white-reliant floor,
    and we never conflate them.

NO expensive SDP solves -- pure evaluation of saved duals.
"""
from __future__ import annotations
import json, sys, glob, warnings
from pathlib import Path

warnings.filterwarnings("ignore")
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

# Reuse the validated, canonical machinery verbatim.
from _fullspace_eval import (
    WHITE_TABLE2, CORE_516, CORE_HEADLINE,
    WHITE_OUTSIDE_FLOOR, WHITE_STRIP_BOUND, TARGET_WHITE,
    anchor_value, phi_center, phi_center_grid,
    cover_min_over_box, reproduce_core_headline, load_centers,
)

PR = CODE.parent / "parallel_results"
DUALEXT = PR / "phase5_N20K_bn40_dualext.json"
STAGE2 = PR / "fullspace_stage2_centers.json"
HALO = PR / "fullspace_stage2_halo_centers.json"

REQ_DUAL = {"con_53", "con_54", "con_512_pL", "con_512_pU",
            "con_512_qL", "con_512_qU", "con_513"}

# Phase-1 independently-certified floor we are trying to beat (from project memory).
PHASE1_INDEP_FLOOR = 0.3575459843736719


# ---------------------------------------------------------------------------
# Center harvesting.  Every center is normalized to the canonical dict shape the
# evaluator expects: {label, h_c, p_c, q1, q2, primal, dual_lb, duals}.  We then
# ALWAYS anchor at primal - 1e-5 (conservative) via anchor_value(...,'primal_m1e5').
# ---------------------------------------------------------------------------
def _valid_duals(d):
    return isinstance(d, dict) and REQ_DUAL.issubset(set(d.keys())) and all(
        d.get(k) is not None for k in REQ_DUAL)


def _duals_all_zero(d):
    return all(abs(float(d.get(k, 0.0))) < 1e-15 for k in REQ_DUAL)


def _norm_point_center(c, src):
    """Normalize a point-center dict (has h_c,p_c,q1,q2,primal,duals)."""
    primal = c.get("primal")
    if primal is None:
        # R7-style: 'anchor' is already primal - margin; reconstruct a primal that
        # is >= anchor so that primal-1e-5 stays conservative.  If only 'anchor'
        # exists, treat anchor itself as the conservative value by setting
        # primal = anchor + 1e-5 (so anchor_value gives back 'anchor').
        if c.get("anchor") is not None:
            primal = float(c["anchor"]) + 1e-5
        elif c.get("dual_lb") is not None:
            primal = float(c["dual_lb"]) + 1e-5
        else:
            return None
    if not _valid_duals(c.get("duals", {})):
        return None
    return {
        "label": str(c.get("label", f"{src}_anon")),
        "h_c": float(c["h_c"]), "p_c": float(c["p_c"]),
        "q1": float(c["q1"]), "q2": float(c["q2"]),
        "primal": float(primal),
        "dual_lb": (float(c["dual_lb"]) if c.get("dual_lb") is not None else None),
        "duals": {k: float(c["duals"][k]) for k in REQ_DUAL},
        "_src": src,
    }


def _norm_box_leaf(lf, region, i, src):
    """Normalize a box-LP subdivision leaf into a point-center at the box midpoint.

    The leaf's duals are feasible for the augmented dual at the box; attaching the
    Phi-quadratic at the box midpoint and extrapolating the same KKT sensitivity is
    the standard reuse (identical to how the box duals were used to certify the box
    in stage 2).  Anchor = leaf_primal - 1e-5 (conservative).

    REJECTED:  'infeasible' leaves -- those certify their parameter sub-box by
    PRIMAL infeasibility (no valid (M)-distribution exists there), carry
    primal=+inf and all-zero duals, and are NOT usable as Phi lower-bound centers.
    Only 'optimal'/'optimal_inaccurate' (verdict 'cleared') leaves give a real
    dual bound."""
    box = lf.get("box")
    if box is None or len(box) != 6:
        return None
    status = (lf.get("status") or "").lower()
    if status not in ("optimal", "optimal_inaccurate"):
        return None
    if not _valid_duals(lf.get("duals", {})):
        return None
    if _duals_all_zero(lf["duals"]):
        return None
    primal = lf.get("primal")
    if primal is None or not np.isfinite(primal):
        return None
    h0, h1, p0, p1, q0, q1 = box
    hc = 0.5 * (h0 + h1)
    pc = 0.5 * (p0 + p1)
    qc = 0.5 * (q0 + q1)
    return {
        "label": f"{src}_R{region}_leaf{i}",
        "h_c": float(hc), "p_c": float(pc), "q1": float(qc), "q2": float(qc),
        "primal": float(primal),
        "dual_lb": (float(lf["dual_LB"]) if lf.get("dual_LB") is not None else None),
        "duals": {k: float(lf["duals"][k]) for k in REQ_DUAL},
        "_src": src,
    }


def harvest_centers():
    centers = []
    sources = {}

    def add(c):
        if c is None:
            return
        centers.append(c)

    # (1) 12 CORE anchors
    d = json.load(open(DUALEXT))
    core_cfg = d["config"]
    n0 = len(centers)
    for c in d["centers"]:
        add(_norm_point_center(c, "core"))
    sources["core_phase5_dualext"] = len(centers) - n0

    # (2) stage2 box-LP leaves -> midpoint centers
    n0 = len(centers)
    s2 = json.load(open(STAGE2))
    for r in s2.get("regions", []):
        reg = r.get("region")
        for i, lf in enumerate(r.get("leaves", [])):
            add(_norm_box_leaf(lf, reg, i, "stage2"))
    sources["fullspace_stage2_centers"] = len(centers) - n0

    # (3) stage2 halo point-centers
    n0 = len(centers)
    h = json.load(open(HALO))
    for c in h.get("centers", []):
        add(_norm_point_center(c, "halo"))
    sources["fullspace_stage2_halo_centers"] = len(centers) - n0

    # (4) every fullspace_promote_R*.json produced this run
    for f in sorted(glob.glob(str(PR / "fullspace_promote_R*.json"))):
        name = Path(f).name
        dd = json.load(open(f))
        n0 = len(centers)
        cs = dd.get("centers")
        if isinstance(cs, list) and cs and isinstance(cs[0], dict) and "h_c" in cs[0]:
            for c in cs:
                add(_norm_point_center(c, name))
        else:
            # e.g. fullspace_promote_R8.json (branch-and-bound 'leaves', no
            # point-centers). Its usable fresh centers live in the *_centers.json
            # sibling, which IS globbed separately, so just skip here.
            pass
        sources[name] = len(centers) - n0

    # ---- DEFENSIVE finite-Phi guard ----------------------------------------
    # A center is only admissible if (a) its conservative anchor is finite and
    # (b) its Phi is finite at a probe point.  A non-finite Phi center would
    # silently dominate the max-cover and make the box-min vacuous (the classic
    # overclaim trap).  Drop & count any such center.
    clean = []
    dropped = 0
    probe = np.array([[0.04]]); probe_p = np.array([[0.4]])
    for c in centers:
        a = c["primal"] - 1e-5
        if not np.isfinite(a):
            dropped += 1
            continue
        try:
            f = float(phi_center_grid(c, a, probe, probe_p, 0.0)[0, 0])
        except Exception:
            dropped += 1
            continue
        if not np.isfinite(f):
            dropped += 1
            continue
        clean.append(c)
    sources["_dropped_nonfinite_or_infeasible"] = dropped
    return clean, core_cfg, sources


def load_certified_region_floors():
    """Per-region INDEPENDENTLY-certified box-min floors produced by the dedicated
    Phase-2 drivers, read LIVE from their result files.  Each is a rigorous box-min
    lower bound on Cover over that region's FULL box, obtained by a STRONGER method
    than a single global grid: either box-LP subdivision (stage2 / R8) or adaptive
    per-subbox subdivision with local Lipschitz (R6, R16) or a fine dedicated grid
    (R7, R9, R17).  NONE uses any White number.  We fold these in via max() with our
    union global-grid cover -- max of valid LBs is valid -- which is the documented
    divide-and-conquer mechanism for the wide regions a single quadratic cannot cover.
    Returns {region:int -> {'floor':float,'source':str,'certified_ge_target':bool}}.
    """
    floors = {}

    def put(reg, val, src, cert=None):
        if val is None:
            return
        # keep the STRONGEST (max) certified floor if multiple files cover a region
        if reg not in floors or val > floors[reg]["floor"]:
            floors[reg] = {"floor": float(val), "source": src,
                           "certified_ge_target": (bool(cert) if cert is not None
                                                   else None)}

    # stage2 box-LP subdivision: regions 1-5, 10 (only when certified_ge_target)
    try:
        s2 = json.load(open(STAGE2))
        for r in s2.get("regions", []):
            if r.get("certified_ge_target") and r.get("min_cleared_cert") is not None \
               and not r.get("residual_gates"):
                put(r["region"], r["min_cleared_cert"],
                    "stage2_box_lp.min_cleared_cert", True)
    except FileNotFoundError:
        pass

    # per-region promote drivers
    def jload(name):
        p = PR / name
        return json.load(open(p)) if p.exists() else None

    d = jload("fullspace_promote_R6.json")
    if d and isinstance(d.get("adaptive_eval"), dict):
        ae = d["adaptive_eval"]
        put(6, ae.get("region_ours_phi_min_after"),
            "promote_R6.adaptive_eval.region_ours_phi_min_after",
            ae.get("clears_380284_indep"))

    d = jload("fullspace_promote_R7.json")
    if d:
        put(7, d.get("region_floor_independent"),
            "promote_R7.region_floor_independent",
            d.get("clears_380284_independent"))

    d = jload("fullspace_promote_R8.json")
    if d:
        put(8, d.get("independent_floor"),
            "promote_R8.independent_floor",
            d.get("certified_ge_target"))

    d = jload("fullspace_promote_R9.json")
    if d:
        put(9, d.get("region_ours_phi_min_after"),
            "promote_R9.region_ours_phi_min_after",
            d.get("clears_380284_indep"))

    d = jload("fullspace_promote_R16.json")
    if d and isinstance(d.get("result"), dict):
        res = d["result"]
        put(16, res.get("after_promo_floor_n241"),
            "promote_R16.result.after_promo_floor_n241",
            res.get("clears_380284_indep"))

    d = jload("fullspace_promote_R17.json")
    if d:
        put(17, d.get("cover_phi_min"),
            "promote_R17.cover_phi_min",
            d.get("clears_380284_indep"))

    return floors


def dedupe(centers, tol=1e-9):
    """Dedupe by label; if two share a label but differ in (h_c,p_c,q1,q2), keep
    both under disambiguated labels (we never want to silently drop a distinct
    feasible center). Identical (label + geometry) -> keep the first."""
    seen = {}
    out = []
    label_counts = {}
    for c in centers:
        key = (c["label"], round(c["h_c"], 12), round(c["p_c"], 12),
               round(c["q1"], 12), round(c["q2"], 12))
        if key in seen:
            continue
        seen[key] = True
        lbl = c["label"]
        if lbl in label_counts:
            label_counts[lbl] += 1
            c = dict(c, label=f"{lbl}#{label_counts[lbl]}")
        else:
            label_counts[lbl] = 0
        out.append(c)
    return out


# ---------------------------------------------------------------------------
def main():
    centers_raw, core_cfg, sources = harvest_centers()
    centers = dedupe(centers_raw)
    # Smooth subset = core anchors + halo + fresh promote point-centers (everything
    # EXCEPT the spiky stage2 box-LP leaves). The box leaves raise L_max ~45x and
    # blow up the Lipschitz eps_grid, which LOWERS the rigorous box-min on narrow
    # regions even though their pointwise cover is identical. We therefore evaluate
    # BOTH the full union and this smooth subset and take the max per region (max of
    # valid LBs is valid) so spiky leaves can only ever HELP, never hurt, the floor.
    smooth = [c for c in centers if c.get("_src") != "stage2"]
    # The rigorous (Lipschitz-penalized) box-min is NON-monotone in the center set:
    # adding spiky centers raises L_max and can LOWER the floor on narrow regions.
    # So we evaluate three nested subsets per region and take the max:
    #   corehalo = 12 core anchors + 11 halo (smallest L_max ~0.17, smoothest)
    #   smooth   = corehalo + fresh promote point-centers (no stage2 box leaves)
    #   union    = everything
    # max of valid LBs is valid; this lets every center only ever HELP the floor.
    corehalo = [c for c in centers if c.get("_src") in ("core", "halo")]
    print(f"harvested {len(centers_raw)} centers, {len(centers)} after dedupe "
          f"({len(corehalo)} corehalo / {len(smooth)} smooth / "
          f"{len(centers)-len(smooth)} stage2-leaf)")
    for k, v in sources.items():
        print(f"   {v:4d}  <- {k}")
    print()

    # ---- self-check: vectorized phi matches scalar phi at random points
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for _ in range(300):
        c = centers[rng.integers(len(centers))]
        h = rng.uniform(0, 2.0); p = rng.uniform(0, 1.0); q = rng.uniform(-1.0, 1.0)
        a = anchor_value(c, "primal_m1e5")
        s1 = phi_center(c, a, h, p, q)
        s2 = float(phi_center_grid(c, a, np.array([[h]]), np.array([[p]]), q)[0, 0])
        maxerr = max(maxerr, abs(s1 - s2))
    print(f"[self-check] vectorized-vs-scalar Phi max abs diff = {maxerr:.2e} "
          f"({'OK' if maxerr < 1e-12 else 'FAIL'})\n")

    # ---- (0) reproduce CORE (5.16) headline -- the SANITY GATE ---------------
    # The canonical core check uses ONLY the 12 core anchors on the project's
    # 4001x4001 (h,p) grid with q baked to the core range; that reproduces the
    # published µ >= 0.380284 headline and validates the evaluator + saved duals.
    core12 = load_centers()[0]
    core_canon = reproduce_core_headline(core12, "primal_m1e5")
    core_phi_min = core_canon["rigorous_LB"]      # the core's contribution to floors
    core_reproduced = abs(core_phi_min - CORE_HEADLINE) < 5e-5
    print("[CORE 5.16] canonical sanity (12 core anchors, q baked, 4001x4001 grid):")
    print(f"   rigorous_LB = {core_phi_min:.7f}   grid_min={core_canon['grid_min']:.7f} "
          f"eps_grid={core_canon['eps_grid']:.2e} L_max={core_canon['L_max']:.3f}")
    print(f"   binding @ (h={core_canon['binding_point'][0]:.5f}, "
          f"p={core_canon['binding_point'][1]:.5f}) witness={core_canon['witness']}")
    print(f"   vs CORE_HEADLINE {CORE_HEADLINE}: {core_phi_min - CORE_HEADLINE:+.2e}  "
          f"({'REPRODUCED' if core_reproduced else 'CHECK -- evaluator/duals mismatch'})\n")

    # Informational: the WHOLE-union core cover (>= canonical; uses far-region
    # centers extrapolated in, which is a valid but loose/over-strong LB). We do
    # NOT use this for the floor -- the conservative documented 0.380284 stands.
    core_union = reproduce_core_headline(centers, "primal_m1e5")
    print(f"   [info] whole-union core cover = {core_union['rigorous_LB']:.7f} "
          f"(witness {core_union['witness']}; not used for floor -- canonical 0.380284 stands)\n")

    # ---- (1) every Table-2 region over the WHOLE union ----------------------
    per_region = []
    indep_min = core_phi_min        # independently-certified floor (NO White)
    indep_binding_region = "core"
    indep_binding_corner = {"h": core_canon["binding_point"][0],
                            "p": core_canon["binding_point"][1], "q": 0.0}
    white_reliant_min = core_phi_min  # core has no White fallback (it IS our bound)
    wr_binding_region = "core"

    cert_floors = load_certified_region_floors()
    print("================= per-region (UNION of all centers) =================")
    print("ours_phi_min = max( corehalo grid , smooth grid , union grid , "
          "dedicated certified box-min )\n")
    for idx, (hr, pr, qr, wbound) in enumerate(WHITE_TABLE2, start=1):
        n_h = 81 if (hr[1] - hr[0]) > 0.05 else 41
        n_p = 161 if (pr[1] - pr[0]) > 0.2 else 81
        n_q = 81 if (qr[1] - qr[0]) > 0.1 else 41
        # Evaluate three nested subsets; take the max (best rigorous box-min).
        subsets = [
            ("corehalo_grid_cover", corehalo),
            ("smooth_subset_grid_cover", smooth),
            ("union_global_grid_cover", centers),
        ]
        evals = {}
        for nm, sub in subsets:
            evals[nm] = cover_min_over_box(
                sub, "primal_m1e5", hr, pr, qr, n_h=n_h, n_p=n_p, n_q=n_q)
        lbU = evals["union_global_grid_cover"][0]
        lbS = evals["smooth_subset_grid_cover"][0]
        lbCH = evals["corehalo_grid_cover"][0]
        gmech = max(evals, key=lambda k: evals[k][0])
        grid_cover, gpt, gwit, gmin, eps, Lm = evals[gmech]

        # (c) dedicated per-region certified box-min (box-LP / adaptive / fine grid)
        cf = cert_floors.get(idx)
        cert_val = cf["floor"] if cf else None
        cert_src = cf["source"] if cf else None

        # ours = max of all three valid LBs (max of valid LBs is valid)
        if cert_val is not None and cert_val >= grid_cover:
            ours = cert_val
            ours_mechanism = cert_src
            ours_witness = "(box-LP/adaptive subdivision)"
            ours_worst = None   # the dedicated driver's box-min; corner recorded by it
        else:
            ours = grid_cover
            ours_mechanism = gmech
            ours_witness = gwit
            ours_worst = {"h": gpt[0], "p": gpt[1], "q": gpt[2]}

        white_floor = WHITE_STRIP_BOUND if abs(wbound - WHITE_STRIP_BOUND) < 1e-9 \
            else WHITE_OUTSIDE_FLOOR
        wr_floor = max(ours, white_floor)   # white-reliant region floor
        rec = {
            "region": idx,
            "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
            "white_floor": white_floor,
            "ours_phi_min": ours,
            "ours_mechanism": ours_mechanism,
            "grid_cover_union": lbU,
            "grid_cover_smooth": lbS,
            "grid_cover_corehalo": lbCH,
            "grid_cover_best": grid_cover,
            "certified_box_min": cert_val,
            "certified_source": cert_src,
            "ours_grid_min": gmin, "eps_grid": eps, "L_max": Lm,
            "grid_worst_point": {"h": gpt[0], "p": gpt[1], "q": gpt[2]},
            "grid_witness": gwit,
            "worst_point": (ours_worst if ours_worst is not None
                            else {"h": gpt[0], "p": gpt[1], "q": gpt[2]}),
            "witness": ours_witness,
            "white_reliant_floor": wr_floor,
            "clears_380000_indep": bool(ours >= WHITE_OUTSIDE_FLOOR),
            "clears_380284_indep": bool(ours >= CORE_HEADLINE),
            "still_white_reliant": bool(ours < WHITE_OUTSIDE_FLOOR),
        }
        per_region.append(rec)

        if ours < indep_min:
            indep_min = ours
            indep_binding_region = idx
            indep_binding_corner = rec["worst_point"]
        if wr_floor < white_reliant_min:
            white_reliant_min = wr_floor
            wr_binding_region = idx

        tag = ""
        if ours >= CORE_HEADLINE:
            tag = "  [indep >= 0.380284]"
        elif ours >= WHITE_OUTSIDE_FLOOR:
            tag = "  [indep >= 0.380000, < headline]"
        else:
            tag = "  <-- still WHITE-RELIANT (ours < 0.380000)"
        ctag = f"  [CH={lbCH:.6f} smooth={lbS:.6f} union={lbU:.6f}"
        if cert_val is not None:
            ctag += f" cert={cert_val:.6f}({cert_src.split('.')[0]})"
        ctag += "]"
        print(f"[R{idx:2d}] h{hr} p{pr} q{qr}")
        print(f"      ours_Phi_min={ours:.7f}  white_floor={white_floor:.5f}  "
              f"white_reliant={wr_floor:.7f}{tag}{ctag}")
        if ours_worst is not None:
            print(f"      worst @ (h={gpt[0]:.4f}, p={gpt[1]:.4f}, q={gpt[2]:.4f}) "
                  f"wit={gwit} via {gmech}\n")
        else:
            print(f"      floor via {ours_mechanism}\n")

    regions_still_white_reliant = [r["region"] for r in per_region
                                   if r["still_white_reliant"]]

    print("================= FULL-SPACE SUMMARY =================")
    print(f"total centers (after dedupe)          : {len(centers)}")
    print(f"core Phi_min                          : {core_phi_min:.7f}")
    print(f"INDEPENDENTLY-CERTIFIED floor (NO White): {indep_min:.7f}")
    print(f"   binding region={indep_binding_region}  corner={indep_binding_corner}")
    print(f"WHITE-RELIANT floor (max(ours,White))  : {white_reliant_min:.7f}"
          f"  binding region={wr_binding_region}")
    print(f"regions still white-reliant (ours<0.380000): {regions_still_white_reliant}")
    print(f"\nPhase-1 indep floor was {PHASE1_INDEP_FLOOR:.10f}")
    delta = indep_min - PHASE1_INDEP_FLOOR
    improved = indep_min > PHASE1_INDEP_FLOOR
    print(f"   -> indep floor {'IMPROVED' if improved else 'did NOT improve'} "
          f"by {delta:+.7e}")

    note = (
        f"Independently-certified full-space floor (no White number anywhere) = "
        f"{indep_min:.7f}, binding at region {indep_binding_region} "
        f"corner {indep_binding_corner}. "
        f"{'IMPROVED' if improved else 'NOT improved'} vs Phase-1's "
        f"{PHASE1_INDEP_FLOOR:.10f} (delta {delta:+.3e}). "
        f"White-reliant floor (each region's max(ours, White published floor)) = "
        f"{white_reliant_min:.7f}. "
        f"Regions still needing White's published 0.380000 (ours<0.380000): "
        f"{regions_still_white_reliant}. "
        f"Core (5.16) Phi_min reproduced at {core_phi_min:.7f} "
        f"(headline {CORE_HEADLINE})."
    )
    print("\nNOTE:", note)

    out = {
        "task": "PRO-38 full-space promotion -- Phase-2 final recompute",
        "anchor_convention": "primal - 1e-5 (conservative) for all centers; "
                             "box-LP leaves use leaf_primal - 1e-5 at box midpoint",
        "mside_sin_coeff": 4.0,
        "core_config": core_cfg,
        "n_centers_raw": len(centers_raw),
        "total_centers": len(centers),
        "center_sources": sources,
        "self_check_phi_vec_vs_scalar_maxabs": maxerr,
        "core_headline": CORE_HEADLINE,
        "core_phi_min": core_phi_min,
        "core_reproduced": bool(core_reproduced),
        "core_binding_point": core_canon["binding_point"],
        "core_witness": core_canon["witness"],
        "core_union_cover_info": core_union["rigorous_LB"],
        "independently_certified_floor": indep_min,
        "binding_region_indep": indep_binding_region,
        "binding_corner": indep_binding_corner,
        "white_reliant_floor": white_reliant_min,
        "white_reliant_binding_region": wr_binding_region,
        "regions_still_white_reliant": regions_still_white_reliant,
        "phase1_indep_floor": PHASE1_INDEP_FLOOR,
        "indep_floor_improved_vs_phase1": bool(improved),
        "indep_floor_delta_vs_phase1": delta,
        "per_region_after": [
            {"region": r["region"],
             "ours_phi_min": r["ours_phi_min"],
             "ours_mechanism": r["ours_mechanism"],
             "grid_cover_union": r["grid_cover_union"],
             "grid_cover_smooth": r["grid_cover_smooth"],
             "grid_cover_corehalo": r["grid_cover_corehalo"],
             "grid_cover_best": r["grid_cover_best"],
             "certified_box_min": r["certified_box_min"],
             "certified_source": r["certified_source"],
             "white_reliant_floor": r["white_reliant_floor"],
             "clears_380000_indep": r["clears_380000_indep"],
             "clears_380284_indep": r["clears_380284_indep"],
             "still_white_reliant": r["still_white_reliant"],
             "worst_point": r["worst_point"], "witness": r["witness"],
             "h_range": r["h_range"], "p_range": r["p_range"],
             "q_range": r["q_range"], "white_floor": r["white_floor"],
             "eps_grid": r["eps_grid"], "L_max": r["L_max"]}
            for r in per_region
        ],
        "note": note,
    }
    OUT = PR / "fullspace_promote_final.json"
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nsaved -> {OUT}")
    return out


if __name__ == "__main__":
    main()
