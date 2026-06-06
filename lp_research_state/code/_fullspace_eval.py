"""
STAGE 1 — Full-space promotion of the Erdos minimum-overlap LOWER bound.

Evaluates our Bochner+poly-moment augmented dual cover Phi over White's FULL
(E(M), c1, d1) parameter space, i.e. over the CORE residual region (5.16) AND
each of White's Table-2 "outside" regions.

KEY MATH (verified, from project memory + path_b_analytical.py):
  For a single feasible dual point extracted at a center c (the 12 anchors in
  phase5_N20K_bn40_dualext.json), the dual objective is a GLOBALLY VALID lower
  bound on mu, varying with (h,p,q)=(E(M),c1,d1) only through the closed-form
  quadratic dual-objective shift:

      Phi_c(h, p, q) = anchor_c + shift_c(h, p, q1=q2=q)

  where shift_c is path_b_analytical.dual_objective_shift (sign convention there
  is validated by the project's perturbation tests; +lambda*Drhs for '>=' form,
  -lambda*Drhs for '<=' form). The cover value at a parameter point is

      Cover(h, p, q) = max_c Phi_c(h, p, q)

  (max of valid LBs is a valid LB). To lower-bound the cover over a BOX we grid
  finely and add a rigorous Lipschitz cell-error term eps_grid = L_max * half-diag
  (same convention as cde_evaluate.py / _verify_cover_dualext.py), where L_max is
  the max over centers of the gradient magnitude of the (concave-quadratic) Phi_c
  on the box. The min over a max-of-concaves is NOT necessarily at a corner, so
  the grid+Lipschitz bound is the rigorous box-min lower bound.

Q-DIMENSION SEMANTICS: a region's q-range (q1,q2) is the range of the SINGLE true
value d1*; the LB must hold at that single point. We evaluate Phi with q1=q2=q
(single point) and grid q across the region's [q1,q2]. The con_513 quadratic term
contributes -0.5*lambda_513*(q^2 - q_c^2): Phi decays as |q| grows (concave in q),
so the worst q is at whichever endpoint has the larger |q|.

NO EXPENSIVE SDP SOLVES here. Pure evaluation of saved duals.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import dual_objective_shift  # validated quadratic shift

REPO = CODE.parent.parent
DUALEXT = CODE.parent / "parallel_results" / "phase5_N20K_bn40_dualext.json"

TARGET_WHITE = 0.379005          # White's published Theorem-1 bound
CORE_HEADLINE = 0.380284         # our conservative core headline (primal - 1e-5)
WHITE_OUTSIDE_FLOOR = 0.380000   # literal floor of White's "0.38" rounded entries
WHITE_STRIP_BOUND = 0.37925      # White's strip #18 bound (the true global min)

# ----------------------------------------------------------------------------
# White's covering, transcribed EXACTLY from arXiv:2201.05704 v1, Tables 2 & 3.
# Variable map: White (h1,h2)=E(M) range; (p1,p2)=c1 range; (q1,q2)=d1 range.
# ----------------------------------------------------------------------------
CORE_516 = {  # eq (5.16); covered by Table-3 ellipses (our augmented core)
    "h": (0.0, 0.06), "p": (0.35, 0.45), "q": (-0.02, 0.02),
}

# Table 2: each (h_range, p_range, q_range, white_bound). 18 regions.
WHITE_TABLE2 = [
    ((0.75, 2.0),  (0.0, 1.0),   (-1.0, 1.0),     0.38),   # 1
    ((0.4, 0.75),  (0.0, 1.0),   (-1.0, 1.0),     0.38),   # 2
    ((0.2, 0.4),   (0.0, 1.0),   (-1.0, 1.0),     0.38),   # 3
    ((0.1, 0.2),   (0.0, 1.0),   (-1.0, 1.0),     0.38),   # 4
    ((0.08, 0.1),  (0.0, 1.0),   (-1.0, 1.0),     0.38),   # 5
    ((0.0, 0.08),  (0.0, 1.0),   (-1.0, -0.05),   0.38),   # 6
    ((0.0, 0.08),  (0.0, 1.0),   (-0.05, -0.025), 0.38),   # 7
    ((0.0, 0.08),  (0.0, 1.0),   (0.05, 1.0),     0.38),   # 8
    ((0.0, 0.08),  (0.0, 1.0),   (0.025, 0.05),   0.38),   # 9
    ((0.0, 0.08),  (0.0, 0.25),  (-0.025, 0.025), 0.38),   # 10
    ((0.0, 0.08),  (0.25, 0.3),  (-0.025, 0.025), 0.38),   # 11
    ((0.0, 0.08),  (0.3, 0.33),  (-0.025, 0.025), 0.38),   # 12
    ((0.0, 0.08),  (0.5, 1.0),   (-0.025, 0.025), 0.38),   # 13
    ((0.0, 0.08),  (0.45, 0.5),  (-0.025, 0.025), 0.38),   # 14
    ((0.06, 0.08), (0.33, 0.45), (-0.025, 0.025), 0.38),   # 15
    ((0.0, 0.06),  (0.33, 0.45), (-0.025, -0.02), 0.38),   # 16
    ((0.0, 0.06),  (0.33, 0.45), (0.02, 0.025),   0.38),   # 17
    ((0.0, 0.06),  (0.33, 0.35), (-0.02, 0.02),   0.37925),# 18  (the true global min)
]


def load_centers():
    data = json.load(open(DUALEXT))
    centers = []
    for c in data["centers"]:
        centers.append({
            "label": c["label"],
            "h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
            "primal": c["primal"], "dual_lb": c["dual_lb"], "duals": c["duals"],
        })
    return centers, data["config"]


def anchor_value(c, mode):
    """Conservative anchor for center c.
    mode='primal_m1e5' : primal - 1e-5  (conservative, matches CORE_HEADLINE)
    mode='primal_m1e6' : primal - 1e-6  (matched convention)
    """
    if mode == "primal_m1e5":
        return c["primal"] - 1e-5
    if mode == "primal_m1e6":
        return c["primal"] - 1e-6
    raise ValueError(mode)


def phi_center(c, anchor, h, p, q):
    """Phi_c(h,p,q) = anchor + dual_objective_shift with q1=q2=q (single point)."""
    center = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"]}
    return anchor + dual_objective_shift(h, p, q, q, center, c["duals"])


def phi_center_grid(c, anchor, HH, PP, q):
    """Vectorized Phi_c over an (h,p) meshgrid at a fixed single q value.

    dual_objective_shift is affine in the dual values and quadratic in (h,p),
    so we reconstruct it vectorially to match dual_objective_shift exactly.
    Cross-checked against the scalar phi_center in __main__.
    """
    d = c["duals"]
    h_c = c["h_c"]; p_c = c["p_c"]; q1_c = c["q1"]; q2_c = c["q2"]
    qm2 = q * q                      # max(q^2, q^2)
    qm2_c = max(q1_c**2, q2_c**2)
    # h terms
    Drhs_53 = HH - h_c
    Drhs_54 = (HH * HH - h_c**2) / 2.0
    # p terms
    Drhs_pL = PP - p_c
    Drhs_pU = PP - p_c
    # q terms (single point q)
    Drhs_qL = q - q1_c
    Drhs_qU = q - q2_c
    Drhs_513 = -0.5 * ((PP * PP - p_c**2) + (qm2 - qm2_c))
    shift = (d["con_53"] * Drhs_53
             - d["con_54"] * Drhs_54
             + d["con_512_pL"] * Drhs_pL
             - d["con_512_pU"] * Drhs_pU
             + d["con_512_qL"] * Drhs_qL
             - d["con_512_qU"] * Drhs_qU
             + d["con_513"] * Drhs_513)
    return anchor + shift


def cover_min_over_box(centers, anchor_mode, h_range, p_range, q_range,
                       n_h=121, n_p=121, n_q=41):
    """Rigorous lower bound on min over the box of Cover(h,p,q)=max_c Phi_c.

    Grid (h,p,q) finely; cover = max over centers; box-min = grid-min minus a
    Lipschitz cell-error term. Returns (box_min_lb, worst_point, worst_witness,
    grid_min_raw, eps_grid).
    """
    h0, h1 = h_range; p0, p1 = p_range; q0, q1 = q_range
    h_grid = np.linspace(h0, h1, n_h)
    p_grid = np.linspace(p0, p1, n_p)
    # q grid: include 0 and both endpoints; concave in q so endpoints are worst,
    # but grid the interior too for safety / to locate the worst point.
    if q0 == q1:
        q_grid = np.array([q0])
    else:
        q_grid = np.unique(np.concatenate([
            np.linspace(q0, q1, n_q), [0.0] if (q0 <= 0.0 <= q1) else []]))
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")

    best_anchor = {c["label"]: anchor_value(c, anchor_mode) for c in centers}

    overall_min = np.inf
    worst_pt = None
    worst_wit = None
    # Lipschitz constant: max over centers & box-corners of |grad_{h,p,q} Phi_c|.
    L_max = 0.0
    for c in centers:
        d = c["duals"]
        # grad_h = lambda53 - lambda54*h ; grad_p = (pL-pU) - lambda513*p ;
        # grad_q = (qL - qU) - lambda513*q  (from -0.5*lambda513*(q^2))
        def gmax(lin_const, quad_coeff, lo, hi):
            return max(abs(lin_const + quad_coeff * lo), abs(lin_const + quad_coeff * hi))
        gh = gmax(d["con_53"], -d["con_54"], h0, h1)
        gp = gmax(d["con_512_pL"] - d["con_512_pU"], -d["con_513"], p0, p1)
        gq = gmax(d["con_512_qL"] - d["con_512_qU"], -d["con_513"], q0, q1)
        L_max = max(L_max, float(np.sqrt(gh*gh + gp*gp + gq*gq)))

    for q in q_grid:
        env = np.full_like(HH, -np.inf)
        wit = np.empty(HH.shape, dtype=object)
        for c in centers:
            F = phi_center_grid(c, best_anchor[c["label"]], HH, PP, q)
            mask = F > env
            env[mask] = F[mask]
            wit[mask] = c["label"]
        qmin = float(env.min())
        if qmin < overall_min:
            overall_min = qmin
            arg = np.unravel_index(int(env.argmin()), env.shape)
            worst_pt = (float(HH[arg]), float(PP[arg]), float(q))
            worst_wit = str(wit[arg])

    # cell sizes
    cell_h = (h1 - h0) / (n_h - 1) if n_h > 1 else 0.0
    cell_p = (p1 - p0) / (n_p - 1) if n_p > 1 else 0.0
    cell_q = (q1 - q0) / (len(q_grid) - 1) if len(q_grid) > 1 else 0.0
    half_diag = 0.5 * float(np.sqrt(cell_h**2 + cell_p**2 + cell_q**2))
    eps_grid = L_max * half_diag
    return overall_min - eps_grid, worst_pt, worst_wit, overall_min, eps_grid, L_max


def reproduce_core_headline(centers, anchor_mode):
    """Reproduce the project's CORE (5.16) headline EXACTLY as _verify_cover_dualext.py
    does: q fixed at the core range (baked into const_q via find_ellipse_h_p), a
    4001x4001 (h,p) grid, and the (h,p)-only Lipschitz cell-error eps_grid. This is
    the project's canonical convention; reproducing it validates the evaluator and
    the saved duals before we trust Phi on the outside regions.
    """
    from path_b_analytical import find_ellipse_h_p
    H_BOX = CORE_516["h"]; P_BOX = CORE_516["p"]
    h_grid = np.linspace(*H_BOX, 4001); p_grid = np.linspace(*P_BOX, 4001)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf); wit = np.empty(HH.shape, dtype=object)
    L_max = 0.0
    for c in centers:
        anchor = anchor_value(c, anchor_mode)
        syn = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
               "value": anchor}
        e = find_ellipse_h_p(syn, c["duals"], c["q1"], c["q2"], target=TARGET_WHITE)
        F = (anchor + e["const_q"] + e["A_h2"]*HH*HH + e["A_h1"]*HH + e["A_h0"]
             + e["A_p2"]*PP*PP + e["A_p1"]*PP + e["A_p0"])
        mask = F > env; env[mask] = F[mask]; wit[mask] = c["label"]
        lam = lambda c2, c1, lo, hi: max(abs(2*c2*lo+c1), abs(2*c2*hi+c1))
        L_max = max(L_max, float(np.hypot(lam(e["A_h2"], e["A_h1"], *H_BOX),
                                          lam(e["A_p2"], e["A_p1"], *P_BOX))))
    gmin = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    cell_h = (H_BOX[1]-H_BOX[0])/4000; cell_p = (P_BOX[1]-P_BOX[0])/4000
    eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    return {"rigorous_LB": gmin - eps_grid, "grid_min": gmin, "eps_grid": eps_grid,
            "L_max": L_max, "witness": str(wit[arg]),
            "binding_point": [float(HH[arg]), float(PP[arg])]}


def main():
    centers, config = load_centers()
    print(f"loaded {len(centers)} centers from {DUALEXT.name}; config={config}\n")

    # ---- self-check: vectorized phi matches scalar phi at random points
    rng = np.random.default_rng(0)
    maxerr = 0.0
    for _ in range(200):
        c = centers[rng.integers(len(centers))]
        h = rng.uniform(0, 0.08); p = rng.uniform(0, 1); q = rng.uniform(-0.2, 0.2)
        a = anchor_value(c, "primal_m1e5")
        s1 = phi_center(c, a, h, p, q)
        s2 = float(phi_center_grid(c, a, np.array([[h]]), np.array([[p]]), q)[0, 0])
        maxerr = max(maxerr, abs(s1 - s2))
    print(f"[self-check] vectorized-vs-scalar Phi max abs diff = {maxerr:.2e} "
          f"({'OK' if maxerr < 1e-12 else 'FAIL'})\n")

    out = {"config": config, "white": {
        "target": TARGET_WHITE, "strip_bound": WHITE_STRIP_BOUND,
        "outside_floor_literal": WHITE_OUTSIDE_FLOOR},
        "core_headline": CORE_HEADLINE, "regions": []}

    for anchor_mode, tag in [("primal_m1e5", "CONSERVATIVE primal-1e-5")]:
        print(f"================= anchor mode: {tag} =================\n")

        # --- (0) Reproduce the CORE (5.16) headline (canonical convention) ---
        core = reproduce_core_headline(centers, anchor_mode)
        print(f"[CORE 5.16] canonical reproduction (q baked, 4001x4001 grid):")
        print(f"             rigorous_LB = {core['rigorous_LB']:.7f}")
        print(f"             grid_min={core['grid_min']:.7f} "
              f"eps_grid={core['eps_grid']:.2e} L_max={core['L_max']:.3f}")
        print(f"             binding @ (h={core['binding_point'][0]:.5f}, "
              f"p={core['binding_point'][1]:.5f}) witness={core['witness']}")
        print(f"             vs conservative core headline {CORE_HEADLINE}: "
              f"{core['rigorous_LB'] - CORE_HEADLINE:+.2e}  "
              f"({'REPRODUCED' if abs(core['rigorous_LB']-CORE_HEADLINE) < 5e-5 else 'CHECK'})\n")
        out["core_eval"] = {"phi_min_rigorous": core["rigorous_LB"],
                            "grid_min": core["grid_min"], "eps_grid": core["eps_grid"],
                            "L_max": core["L_max"], "binding_point": core["binding_point"],
                            "witness": core["witness"], "anchor_mode": anchor_mode}

        # --- (1) Each Table-2 outside region ---
        # CORRECT LOGIC: at every parameter point we have TWO valid lower bounds
        # on the true dual optimum: (a) White's published per-region bound, and
        # (b) our augmented Phi. The best certified LB at each point is the MAX
        # of the two. White's number is constant over the region, so the region's
        # certified floor is  max(min_box Phi_ours,  White_floor).
        #   - For White's "0.38" rows we use the literal rigorous floor 0.380000.
        #   - For the strip (row 18) White's exact value 0.37925 is the floor.
        # A region's certified floor >= 0.380000 ALWAYS for the "0.38" rows (White
        # alone). A region is a GATE for promoting the FULL bound to the core
        # headline (0.380284) iff its certified floor < 0.380284, i.e. iff our
        # Phi fails to reach 0.380284 there (White only guarantees 0.380000).
        region_results = []
        gate = {"certified_floor": np.inf}            # worst region for 0.380284
        strip_check = None
        for idx, (hr, pr, qr, wbound) in enumerate(WHITE_TABLE2, start=1):
            # adapt grid density to box size (wide boxes => more points)
            n_h = 81 if (hr[1]-hr[0]) > 0.05 else 41
            n_p = 161 if (pr[1]-pr[0]) > 0.2 else 81
            n_q = 81 if (qr[1]-qr[0]) > 0.1 else 41
            lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
                centers, anchor_mode, hr, pr, qr, n_h=n_h, n_p=n_p, n_q=n_q)
            ours = lb
            # White's rigorous floor for this region:
            white_floor = WHITE_STRIP_BOUND if abs(wbound - WHITE_STRIP_BOUND) < 1e-9 \
                else WHITE_OUTSIDE_FLOOR
            certified = max(ours, white_floor)        # best of two valid LBs
            rr = {
                "region": idx,
                "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
                "white_bound_displayed": wbound,
                "white_floor_rigorous": white_floor,
                "ours_phi_min": ours, "ours_grid_min": gmin, "eps_grid": eps,
                "L_max": Lm,
                "ours_worst_point": pt, "ours_witness": wit,
                "certified_floor": certified,          # max(ours, White_floor)
                "ours_reaches_core_headline": bool(ours >= CORE_HEADLINE),
                "certified_clears_380000": bool(certified >= WHITE_OUTSIDE_FLOOR),
                "certified_clears_core_headline": bool(certified >= CORE_HEADLINE),
                "gate_for_core_headline": bool(certified < CORE_HEADLINE),
            }
            region_results.append(rr)
            if idx == 18:
                strip_check = rr
            flag = ""
            if certified < CORE_HEADLINE:
                flag = "  <-- GATE for 0.380284 (ours below it; White floor=%.5f)" % white_floor
            if certified < gate["certified_floor"]:
                gate = {"region": idx, "certified_floor": certified,
                        "white_floor": white_floor, "ours_phi_min": ours,
                        "worst_point": pt, "witness": wit,
                        "h_range": list(hr), "p_range": list(pr),
                        "q_range": list(qr)}
            print(f"[R{idx:2d}] h{hr} p{pr} q{qr}")
            print(f"      White_floor={white_floor:.5f}  ours_Phi_min={ours:.6f}"
                  f"  certified=max={certified:.6f}{flag}")
            print(f"      ours worst @ (h={pt[0]:.4f}, p={pt[1]:.4f}, "
                  f"q={pt[2]:.4f}) wit={wit}\n")

        out["regions"] = region_results

        # --- full-space summary ---
        # (A) Is mu >= 0.380000 already established (literal White floor)?
        min_certified_380 = min(r["certified_floor"] for r in region_results)
        core_floor = out["core_eval"]["phi_min_rigorous"]
        fullspace_380 = min(min_certified_380, core_floor)
        # (B) Gate for promoting the FULL bound to the core headline 0.380284:
        gate_regions = [r["region"] for r in region_results
                        if r["gate_for_core_headline"]]

        print("================= FULL-SPACE SUMMARY =================")
        print(f"Core (5.16): our augmented cover >= {core_floor:.7f} "
              f"(conservative headline {CORE_HEADLINE})")
        print(f"Strip (R18): White=0.37925, ours_Phi_min={strip_check['ours_phi_min']:.6f}"
              f" -> certified {strip_check['certified_floor']:.6f}"
              f" ({'LIFTED above 0.380284' if strip_check['certified_clears_core_headline'] else 'NOT lifted to 0.380284'})")
        print(f"\n(A) FULL-SPACE mu LB taking White's '0.38' literally (floor 0.380000):")
        print(f"    min over ALL regions of max(ours, White_floor) = {min_certified_380:.7f}")
        print(f"    combined with core {core_floor:.6f}  ->  mu >= {fullspace_380:.7f}")
        print(f"    {'>>> mu >= 0.380000 ESTABLISHED' if fullspace_380 >= WHITE_OUTSIDE_FLOOR - 1e-12 else 'NOT yet 0.380000 (check core grid resolution)'}")
        print(f"\n(B) Promotion to core headline {CORE_HEADLINE}:")
        print(f"    GATE regions (ours fails to reach {CORE_HEADLINE}; White only 0.380000): "
              f"{gate_regions}")
        print(f"    worst certified floor among Table-2 = {gate['certified_floor']:.7f}"
              f" at region {gate['region']}")
        out["fullspace"] = {
            "min_certified_table2_floor": min_certified_380,
            "core_phi_min": core_floor,
            "fullspace_lb_literal_white_0p38": fullspace_380,
            "established_0p380000": bool(fullspace_380 >= WHITE_OUTSIDE_FLOOR - 1e-12),
            "gate_regions_for_core_headline": gate_regions,
            "worst_certified_floor_table2": gate["certified_floor"],
            "worst_region": gate["region"],
            "strip_R18_lifted": bool(strip_check["certified_clears_core_headline"]),
        }

    OUT = REPO / "lp_research_state" / "parallel_results" / "fullspace_stage1.json"
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nsaved -> {OUT}")


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
