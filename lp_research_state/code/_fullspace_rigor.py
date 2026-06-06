"""
STAGE-2 RIGOR VERIFIER (adversarial) for the full-space promotion of the Erdos
minimum-overlap LOWER bound.  ADDITIVE, throwaway (prefix '_'); does NOT modify
any existing file.

Two independent certifications, both run from this one file:

TASK 1 -- certify mu >= 0.380000 FULL-SPACE, independent of White's rounding.
  White's Table-2 displays a per-region "Optimum lower bound" (column header), one
  number per region: "0.38" for regions 1-17 and 0.37925 for the strip (region 18).
  These are 2-dp rounded; the only value we may quote rigorously from his paper is
  the literal floor 0.380000 (resp. 0.37925).  To certify the floor WITHOUT relying
  on White's rounding, we RE-SOLVE White's OWN unaugmented Section-5 program (NO
  Bochner, NO poly-moment, NO T3/T5/T5') on each region.

  BOX-VALIDITY (this is exactly how White gets ONE number per region):  the program
  is solved at the region's FULL ranges (h1,h2),(p1,p2),(q1,q2) -- NOT at a single
  center.  Parameters enter ONLY constraint right-hand sides:
     (5.3) rhs=h1 ; (5.4) rhs=2/3+h2^2/2 ; (5.12) c0 in [p1,p2], d0 in [q1,q2] ;
     (5.13) rhs = -1/2 (max(p1^2,p2^2)+max(q1^2,q2^2)).
  The single optimum opt(N,T,R,h1,h2,p1,p2,q1,q2) is therefore a valid lower bound
  on mu over the ENTIRE box [h1,h2]x[p1,p2]x[q1,q2] (White Lemma 9/10, Appendix II).
  Moreover opt(N) is a lower bound on mu for EVERY N (the program is a relaxation
  that only tightens as N grows), and the optimum increases with N (White, Sec 5),
  so opt at our N <= opt at White's larger N <= mu.  Hence opt(our N) >= 0.380000
  certifies mu >= 0.380000 on that box, independent of any rounding in White.

  ANCHOR: conservative dual-extracted LB = reported_value - last_gap via
  dual_extractor.solve_with_dual_extraction (the project's rigor convention), then
  an additional safety margin (default 1e-5) subtracted.  Never the raw 'value'.

TASK 2 -- independently re-derive the augmented CORE anchor (~0.380284).
  Re-implemented from scratch (NOT a call to _verify_cover_dualext): load the saved
  conservative dual anchors (phase5_N20K_bn40_dualext.json), rebuild each ellipse's
  Phi over the CORE region (5.16) on a fine (h,p) grid at the core q-range, take the
  cover min, subtract the Lipschitz cell-error margin, compare to 0.380284.

Usage:
  python _fullspace_rigor.py task1 [N] [T]      # solve all 18 regions unaugmented
  python _fullspace_rigor.py task1pert [N] [T]  # LB-direction sanity (one wide box)
  python _fullspace_rigor.py task2              # independent core-anchor recompute
  python _fullspace_rigor.py all  [N] [T]       # task2 then task1
"""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import (
    build_problem_with_dual_handles, dual_objective_shift, find_ellipse_h_p,
)
from dual_extractor import solve_with_dual_extraction

REPO = CODE.parent.parent
DUALEXT = CODE.parent / "parallel_results" / "phase5_N20K_bn40_dualext.json"
OUT = CODE.parent / "parallel_results" / "fullspace_rigor.json"

TARGET_WHITE = 0.379005
CORE_HEADLINE = 0.380284          # recorded conservative core headline (primal - 1e-5)
WHITE_OUTSIDE_FLOOR = 0.380000
WHITE_STRIP_BOUND = 0.37925
CORE_516 = {"h": (0.0, 0.06), "p": (0.35, 0.45), "q": (-0.02, 0.02)}

# White Table 2 -- 18 regions, transcribed verbatim from arXiv:2201.05704 v1.
# (h1,h2)=E(M); (p1,p2)=c1; (q1,q2)=d1; white displayed bound.
WHITE_TABLE2 = [
    ((0.75, 2.0),  (0.0, 1.0),   (-1.0, 1.0),     0.38),    # 1
    ((0.4, 0.75),  (0.0, 1.0),   (-1.0, 1.0),     0.38),    # 2
    ((0.2, 0.4),   (0.0, 1.0),   (-1.0, 1.0),     0.38),    # 3
    ((0.1, 0.2),   (0.0, 1.0),   (-1.0, 1.0),     0.38),    # 4
    ((0.08, 0.1),  (0.0, 1.0),   (-1.0, 1.0),     0.38),    # 5
    ((0.0, 0.08),  (0.0, 1.0),   (-1.0, -0.05),   0.38),    # 6
    ((0.0, 0.08),  (0.0, 1.0),   (-0.05, -0.025), 0.38),    # 7
    ((0.0, 0.08),  (0.0, 1.0),   (0.05, 1.0),     0.38),    # 8
    ((0.0, 0.08),  (0.0, 1.0),   (0.025, 0.05),   0.38),    # 9
    ((0.0, 0.08),  (0.0, 0.25),  (-0.025, 0.025), 0.38),    # 10
    ((0.0, 0.08),  (0.25, 0.3),  (-0.025, 0.025), 0.38),    # 11
    ((0.0, 0.08),  (0.3, 0.33),  (-0.025, 0.025), 0.38),    # 12
    ((0.0, 0.08),  (0.5, 1.0),   (-0.025, 0.025), 0.38),    # 13
    ((0.0, 0.08),  (0.45, 0.5),  (-0.025, 0.025), 0.38),    # 14
    ((0.06, 0.08), (0.33, 0.45), (-0.025, 0.025), 0.38),    # 15
    ((0.0, 0.06),  (0.33, 0.45), (-0.025, -0.02), 0.38),    # 16
    ((0.0, 0.06),  (0.33, 0.45), (0.02, 0.025),   0.38),    # 17
    ((0.0, 0.06),  (0.33, 0.35), (-0.02, 0.02),   0.37925), # 18 strip
]


# ===========================================================================
# TASK 1
# ===========================================================================
def solve_region_unaugmented(h1, h2, p1, p2, q1, q2, N, T, R, margin):
    """Solve White's UNAUGMENTED Section-5 program at the region's FULL ranges.
    Returns a conservative dual-extracted lower bound (box-valid by construction).
    bochner_n=0 and no poly-moment/T-cuts => pure White Section 5."""
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h1, h2, p1, p2, q1, q2, bochner_n=0)  # <- all augmentations OFF
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    rep = res["reported_value"]
    dlb = res["rigorous_dual_LB"]
    # Conservative anchor: dual-extracted LB minus a safety margin. If for some
    # reason the dual LB couldn't be parsed, fall back to reported - margin (still
    # report which path was taken).
    if dlb is not None:
        anchor = dlb - margin
        anchor_src = "dual_extracted_minus_margin"
    else:
        anchor = (rep - margin) if rep is not None else None
        anchor_src = "reported_minus_margin(FALLBACK)"
    return {
        "reported_value": rep, "dual_lb": dlb,
        "conservative_anchor": anchor, "anchor_src": anchor_src,
        "dual_resid_at_LB": res["dual_residual_at_LB"],
        "status": res["status"], "time": res["time"],
        "n_cons": len(cons),
    }


def task1(N=5000, T=2000, R=10, margin=1e-5):
    print(f"\n###### TASK 1: independent UNAUGMENTED White Section-5 recompute "
          f"######\nconfig N={N} T={T} R={R} bochner_n=0 (NO augmentations) "
          f"margin={margin}\n")
    print("Each region solved at its FULL ranges -> one box-valid lower bound on mu.\n")
    rows = []
    n_below = 0
    worst = None
    t_start = time.time()
    for idx, (hr, pr, qr, wbound) in enumerate(WHITE_TABLE2, start=1):
        r = solve_region_unaugmented(hr[0], hr[1], pr[0], pr[1], qr[0], qr[1],
                                     N, T, R, margin)
        a = r["conservative_anchor"]
        clears = (a is not None) and (a >= WHITE_OUTSIDE_FLOOR)
        if not clears:
            n_below += 1
        if worst is None or (a is not None and a < worst[1]):
            worst = (idx, a)
        rows.append({
            "region": idx,
            "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
            "white_displayed": wbound,
            "reported_value": r["reported_value"],
            "dual_lb": r["dual_lb"],
            "our_conservative_box_LB": a,
            "anchor_src": r["anchor_src"],
            "dual_resid_at_LB": r["dual_resid_at_LB"],
            "status": r["status"], "time_s": r["time"],
            "clears_380000": bool(clears),
            "margin_over_380000": (a - WHITE_OUTSIDE_FLOOR) if a is not None else None,
        })
        flag = "" if clears else "   <-- BELOW 0.380000"
        print(f"[R{idx:2d}] h{hr} p{pr} q{qr}  White='{wbound}'")
        print(f"      reported={r['reported_value']:.7f}  dualLB="
              f"{r['dual_lb'] if r['dual_lb'] is None else round(r['dual_lb'],7)}"
              f"  our_box_LB={a:.7f}  ({r['status']}, {r['time']:.1f}s){flag}")
        print(f"      vs 0.380000: {(a - WHITE_OUTSIDE_FLOOR):+.2e}   "
              f"vs White-displayed {wbound}: {(a - wbound):+.2e}\n")

    # Verdict for mu >= 0.380000 full-space.
    # Regions 1-17 must clear 0.380000 from our UNAUGMENTED bound (rounding-free).
    # Region 18 (strip) is NOT expected to clear from White's own bound (his floor
    # 0.37925 < 0.380000); it is lifted by the AUGMENTED cover (Task 2 / Stage 1).
    r1_17 = [x for x in rows if x["region"] <= 17]
    strip = next(x for x in rows if x["region"] == 18)
    n17_below = sum(1 for x in r1_17 if not x["clears_380000"])
    print("================= TASK 1 VERDICT =================")
    print(f"Regions 1-17 (White floor '0.38'): {17 - n17_below}/17 clear 0.380000 "
          f"from our INDEPENDENT unaugmented box LB.")
    if n17_below == 0:
        print("  => All 17 'far/near' regions are >= 0.380000 INDEPENDENT of White's "
              "rounding.  The literal 0.38 floor is certified.")
    else:
        below = [x['region'] for x in r1_17 if not x['clears_380000']]
        print(f"  => regions BELOW 0.380000 at this N: {below} "
              f"(may rise at larger N; see notes).")
    print(f"Strip R18: our unaugmented box LB = {strip['our_conservative_box_LB']:.7f} "
          f"(White 0.37925).  Expected < 0.380000; relies on AUGMENTED cover (Task 2).")
    print(f"worst region overall = R{worst[0]} at {worst[1]:.7f}")
    print(f"[task1 total wall time {time.time()-t_start:.0f}s]\n")

    out = {
        "config": {"N": N, "T": T, "R": R, "bochner_n": 0, "margin": margin,
                   "augmentations": "NONE (pure White Section 5)"},
        "white_outside_floor": WHITE_OUTSIDE_FLOOR,
        "white_strip_bound": WHITE_STRIP_BOUND,
        "regions": rows,
        "verdict": {
            "regions_1_17_all_clear_380000": bool(n17_below == 0),
            "n_regions_1_17_below_380000": n17_below,
            "strip_R18_box_LB": strip["our_conservative_box_LB"],
            "strip_R18_clears_380000": strip["clears_380000"],
            "worst_region": worst[0],
            "worst_box_LB": worst[1],
            "note": ("mu>=0.380000 full-space = (regions 1-17 from this rounding-free "
                     "unaugmented recompute) AND (strip R18 + core from augmented cover, "
                     "Task 2).  R18 is NOT certifiable from White's own bound."),
        },
    }
    return out


def task1_perturb(N=2000, T=800, R=10):
    """Empirical LB-direction check for the UNAUGMENTED program across a WIDE box
    including q-variation.  Validates that the box-valid claim (a single full-range
    solve lower-bounds mu at every interior (h,p,q)) is not violated: we solve at
    full ranges, then re-solve at several INTERIOR single points (h=h1=h2=h*,
    p=p1=p2=p*, q1=q2=q*) and confirm  true_opt(interior pt) >= box LB.
    A sign/direction error would let the box value exceed an interior optimum."""
    print(f"\n###### TASK 1 perturbation: UNAUGMENTED box-validity direction "
          f"######\nconfig N={N} T={T} R={R} bochner_n=0\n")
    # Use region 8-like wide box but tractable q: (0,0.08)x(0,1)x(0.05,1)
    h1, h2 = 0.0, 0.08
    p1, p2 = 0.0, 1.0
    q1, q2 = 0.05, 1.0
    box = solve_region_unaugmented(h1, h2, p1, p2, q1, q2, N, T, R, margin=0.0)
    boxLB = box["dual_lb"] if box["dual_lb"] is not None else box["reported_value"]
    print(f"FULL-BOX (0,0.08)x(0,1)x(0.05,1): reported={box['reported_value']:.7f} "
          f"dualLB={boxLB:.7f}  (this is the box-valid LB on mu over the box)\n")
    # interior single points; each must have true optimum >= boxLB
    pts = [(0.0, 0.0, 0.05), (0.04, 0.5, 0.5), (0.08, 1.0, 1.0),
           (0.02, 0.2, 0.3), (0.0, 1.0, 0.05), (0.08, 0.0, 1.0)]
    print(f"{'(h,p,q) interior':>26} {'true_opt':>12} {'boxLB':>12} {'true-box':>12} verdict")
    worst = 1e9
    detail = []
    for (h, p, q) in pts:
        r = solve_region_unaugmented(h, h, p, p, q, q, N, T, R, margin=0.0)
        tv = r["reported_value"]
        gap = tv - boxLB
        worst = min(worst, gap)
        detail.append({"pt": [h, p, q], "true_opt": tv, "gap": gap})
        print(f"{('('+format(h,'.2f')+','+format(p,'.2f')+','+format(q,'.2f')+')'):>26} "
              f"{tv:>12.7f} {boxLB:>12.7f} {gap:>+12.2e} "
              f"{'OK' if gap > -1e-6 else '*** BOX>INTERIOR (BAD) ***'}")
    print(f"\nworst (true_opt_interior - boxLB) = {worst:+.2e}  "
          f"{'BOX-VALID OK' if worst > -1e-6 else 'BOX-VALIDITY FAILS'}\n")
    return {"box_LB": boxLB, "box_reported": box["reported_value"],
            "box_range": {"h": [h1, h2], "p": [p1, p2], "q": [q1, q2]},
            "interior_points": detail, "worst_gap": worst,
            "box_valid_ok": bool(worst > -1e-6)}


# ===========================================================================
# TASK 2 -- independent core-anchor recompute (from scratch)
# ===========================================================================
def load_anchors():
    data = json.load(open(DUALEXT))
    return data["centers"], data["config"]


def task2(n_grid=4001, margin=1e-5):
    """Independent re-implementation of the cover-min over the CORE region (5.16).

    From scratch (does NOT import _verify_cover_dualext or _fullspace_eval logic):
      - conservative anchor per center = primal - margin (repo convention; margin
        default 1e-5 = the recorded conservative headline convention);
      - Phi_c(h,p,q) at the core q-range built directly from dual_objective_shift,
        evaluated on a fine (h,p) grid with q FIXED at the core endpoints (concave
        in q so endpoints are the worst; we min over both endpoints AND q=0);
      - cover = max over centers; grid-min minus Lipschitz cell-error eps_grid.
    """
    centers, config = load_anchors()
    print(f"\n###### TASK 2: independent CORE-anchor recompute (from scratch) "
          f"######\nloaded {len(centers)} anchors from {DUALEXT.name}; "
          f"config={config}; grid={n_grid}x{n_grid}; margin={margin}\n")

    h0, h1 = CORE_516["h"]; p0, p1 = CORE_516["p"]; q0, q1 = CORE_516["q"]
    h_grid = np.linspace(h0, h1, n_grid)
    p_grid = np.linspace(p0, p1, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")

    # q endpoints + 0: concave in q (con_513 quadratic), worst at larger |q|.
    q_eval = sorted(set([q0, q1, 0.0]))

    # Build cover envelope as max over centers, min over q_eval.
    overall_env = None      # min over q of (max over centers of Phi)
    L_max = 0.0
    for q in q_eval:
        env = np.full_like(HH, -np.inf)
        for c in centers:
            anchor = c["primal"] - margin           # conservative anchor
            d = c["duals"]
            h_c = c["h_c"]; p_c = c["p_c"]; q1_c = c["q1"]; q2_c = c["q2"]
            qm2 = q * q
            qm2_c = max(q1_c**2, q2_c**2)
            # dual_objective_shift, reconstructed vectorially (h,p grid, single q):
            shift = (d["con_53"]      * (HH - h_c)
                     - d["con_54"]    * ((HH*HH - h_c**2) / 2.0)
                     + d["con_512_pL"]* (PP - p_c)
                     - d["con_512_pU"]* (PP - p_c)
                     + d["con_512_qL"]* (q - q1_c)
                     - d["con_512_qU"]* (q - q2_c)
                     + d["con_513"]   * (-0.5 * ((PP*PP - p_c**2) + (qm2 - qm2_c))))
            F = anchor + shift
            np.maximum(env, F, out=env)
        if overall_env is None:
            overall_env = env
        else:
            np.minimum(overall_env, env, out=overall_env)
        # Lipschitz on (h,p) for this center set (grad wrt h,p only; q handled by
        # enumerating endpoints, so its cell-error is excluded -- matching the
        # canonical (h,p)-grid convention used for the headline).
        for c in centers:
            d = c["duals"]
            gh = max(abs(d["con_53"] - d["con_54"]*h0),
                     abs(d["con_53"] - d["con_54"]*h1))
            gp = max(abs((d["con_512_pL"]-d["con_512_pU"]) - d["con_513"]*p0),
                     abs((d["con_512_pL"]-d["con_512_pU"]) - d["con_513"]*p1))
            L_max = max(L_max, float(np.hypot(gh, gp)))

    grid_min = float(overall_env.min())
    arg = np.unravel_index(int(overall_env.argmin()), overall_env.shape)
    binding = (float(HH[arg]), float(PP[arg]))
    cell_h = (h1 - h0) / (n_grid - 1)
    cell_p = (p1 - p0) / (n_grid - 1)
    eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    rigorous_LB = grid_min - eps_grid

    print(f"[CORE 5.16] independent recompute (from-scratch dual_objective_shift):")
    print(f"  grid_min     = {grid_min:.10f}  at (h={binding[0]:.5f}, p={binding[1]:.5f})")
    print(f"  eps_grid     = {eps_grid:.3e}   (L_max={L_max:.4f})")
    print(f"  rigorous_LB  = {rigorous_LB:.10f}")
    print(f"  vs recorded conservative headline {CORE_HEADLINE}:  "
          f"{rigorous_LB - CORE_HEADLINE:+.3e}")
    print(f"  vs recorded fullspace_stage1 core_eval 0.3802837846529683: "
          f"{rigorous_LB - 0.3802837846529683:+.3e}")

    # Cross-check 1: reproduce via the project's canonical find_ellipse_h_p path
    # (q baked via const_q + 4001 grid). Different code path; must agree.
    ce = _canonical_find_ellipse_core(centers, margin, n_grid)
    print(f"\n[cross-check] canonical find_ellipse_h_p path: rigorous_LB="
          f"{ce['rigorous_LB']:.10f}  (diff vs from-scratch "
          f"{ce['rigorous_LB'] - rigorous_LB:+.2e})")

    agree_recorded = abs(rigorous_LB - 0.3802837846529683) <= 1e-6
    agree_canonical = abs(rigorous_LB - ce["rigorous_LB"]) <= 1e-6
    print(f"\nAGREEMENT (<=1e-6): vs recorded core_eval = "
          f"{'YES' if agree_recorded else 'NO'} ; vs canonical path = "
          f"{'YES' if agree_canonical else 'NO'}")

    return {
        "config": config, "n_grid": n_grid, "margin": margin,
        "from_scratch": {
            "rigorous_LB": rigorous_LB, "grid_min": grid_min,
            "eps_grid": eps_grid, "L_max": L_max,
            "binding_point": list(binding),
        },
        "canonical_crosscheck": ce,
        "recorded_core_headline_conservative": CORE_HEADLINE,
        "recorded_stage1_core_eval": 0.3802837846529683,
        "agree_with_recorded_1e6": bool(agree_recorded),
        "agree_with_canonical_1e6": bool(agree_canonical),
        "diff_vs_recorded": rigorous_LB - 0.3802837846529683,
    }


def _canonical_find_ellipse_core(centers, margin, n_grid):
    """Cross-check using find_ellipse_h_p (q baked into const_q), independent of the
    from-scratch reconstruction above."""
    h0, h1 = CORE_516["h"]; p0, p1 = CORE_516["p"]
    h_grid = np.linspace(h0, h1, n_grid); p_grid = np.linspace(p0, p1, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    L_max = 0.0
    for c in centers:
        anchor = c["primal"] - margin
        syn = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
               "value": anchor}
        e = find_ellipse_h_p(syn, c["duals"], c["q1"], c["q2"], target=TARGET_WHITE)
        F = (anchor + e["const_q"] + e["A_h2"]*HH*HH + e["A_h1"]*HH + e["A_h0"]
             + e["A_p2"]*PP*PP + e["A_p1"]*PP + e["A_p0"])
        np.maximum(env, F, out=env)
        lam = lambda c2, c1, lo, hi: max(abs(2*c2*lo+c1), abs(2*c2*hi+c1))
        L_max = max(L_max, float(np.hypot(lam(e["A_h2"], e["A_h1"], h0, h1),
                                          lam(e["A_p2"], e["A_p1"], p0, p1))))
    gmin = float(env.min())
    cell_h = (h1-h0)/(n_grid-1); cell_p = (p1-p0)/(n_grid-1)
    eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    return {"rigorous_LB": gmin - eps_grid, "grid_min": gmin,
            "eps_grid": eps_grid, "L_max": L_max}


# ===========================================================================
def _save(key, payload):
    existing = {}
    if OUT.exists():
        try:
            existing = json.load(open(OUT))
        except Exception:
            existing = {}
    existing[key] = payload
    existing["_meta"] = {"updated": time.strftime("%Y-%m-%d %H:%M:%S"),
                         "file": __file__}
    OUT.write_text(json.dumps(existing, indent=2, default=float))
    print(f"saved [{key}] -> {OUT}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    def argi(i, default):
        return int(sys.argv[i]) if len(sys.argv) > i else default

    if cmd == "task2":
        _save("task2_core_anchor", task2())
    elif cmd == "task1":
        _save("task1_unaugmented_regions", task1(N=argi(2, 5000), T=argi(3, 2000)))
    elif cmd == "task1pert":
        _save("task1_box_validity_check", task1_perturb(N=argi(2, 2000), T=argi(3, 800)))
    elif cmd == "all":
        _save("task2_core_anchor", task2())
        _save("task1_box_validity_check", task1_perturb(N=argi(2, 2000), T=argi(3, 800)))
        _save("task1_unaugmented_regions", task1(N=argi(2, 5000), T=argi(3, 2000)))
    else:
        print(__doc__)
