"""
PRO-38 — Independently lift GATE region R6 via fresh augmented Phi-centers.

R6 = White Table-2 region 6: h in [0,0.08], c1=p in [0,1], d1=q in [-1,-0.05].
width_class = WIDE (q-width 0.95, p-width 1.0). Current best independent Phi floor
(core+halo+stage2-leaf union, from fullspace_recon.json) = 0.357546, worst at
(h=0.016, p=0.325, q=-0.1925). Shortfall to 0.380284 ~= 0.0227.

We solve FRESH single-point augmented dual-feasible centers, dual-extract exactly
as _verify_cover_dualext.py does (7 handles con_53/con_54/con_512_p{L,U}/con_512_q{L,U}
/con_513; poly-moment cuts added separately via build_even_moment_nonneg_constraints).
The path_b builder already bakes the CORRECTED mside coeff = 4.0 (line 106-107), so
every solve uses the corrected program (validity rule 3). Conservative anchor =
primal - 1e-5.

Centers are appended to parallel_results/fullspace_promote_R6.json in the SAME
center format as the halo file. After each solve we re-evaluate R6's box-min via
cover_min_over_box(core_centers + fresh_centers, ...).

WIDE-region budget rule: place ONE probe at the worst corner first; measure how much
of the box it lifts >= 0.380000 and ESTIMATE centers needed. Do NOT brute-force.
"""
from __future__ import annotations
import json, sys, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
import _fullspace_eval as FE

PR = CODE.parent / "parallel_results"
OUT = PR / "fullspace_promote_R6.json"

R6 = dict(h_range=(0.0, 0.08), p_range=(0.0, 1.0), q_range=(-1.0, -0.05))
TARGET = 0.380284
FLOOR = 0.380000
NEED = ("con_53", "con_54", "con_512_pL", "con_512_pU",
        "con_512_qL", "con_512_qU", "con_513")


def solve_center(label, h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    """Solve a single-point augmented center; dual-extract; return center dict in the
    halo/dualext format. Mirrors _verify_cover_dualext.py extraction exactly."""
    t0 = time.time()
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
    cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in NEED}
    primal = res["reported_value"]
    dual_lb = res["rigorous_dual_LB"]
    gap = (primal - dual_lb) if (primal is not None and dual_lb is not None) else None
    c = {
        "label": label, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "primal": primal, "dual_lb": dual_lb, "gap": gap,
        "dual_resid": res["dual_residual_at_LB"], "status": res["status"],
        "duals": duals, "time": time.time() - t0,
        "config": {"N": N, "T": T, "R": R, "bochner_n": bn, "pm_k_max": pm_k},
    }
    return c


def eval_R6(fresh_centers):
    """ours_phi_min over R6 = cover_min_over_box(core + fresh)."""
    core, _ = FE.load_centers()
    cset = core + fresh_centers
    lb, pt, wit, gm, eps, Lm = FE.cover_min_over_box(
        cset, "primal_m1e5", R6["h_range"], R6["p_range"], R6["q_range"],
        n_h=81, n_p=161, n_q=81)
    return dict(box_min_lb=lb, worst=pt, witness=wit, grid_min=gm,
                eps_grid=eps, L_max=Lm, n_centers=len(cset))


def fraction_lifted(fresh_centers, thresh, n_h=81, n_p=161, n_q=81):
    """Fraction of a fine R6 grid where Cover(core+fresh) >= thresh (no Lipschitz;
    pure coverage diagnostic to estimate centers needed)."""
    core, _ = FE.load_centers()
    cset = core + fresh_centers
    h0, h1 = R6["h_range"]; p0, p1 = R6["p_range"]; q0, q1 = R6["q_range"]
    hg = np.linspace(h0, h1, n_h); pg = np.linspace(p0, p1, n_p)
    qg = np.linspace(q0, q1, n_q)
    HH, PP = np.meshgrid(hg, pg, indexing="ij")
    anchors = {c["label"]: FE.anchor_value(c, "primal_m1e5") for c in cset}
    total = HH.size * len(qg)
    ok = 0
    worst = np.inf; worst_pt = None
    for q in qg:
        env = np.full_like(HH, -np.inf)
        for c in cset:
            F = FE.phi_center_grid(c, anchors[c["label"]], HH, PP, q)
            np.maximum(env, F, out=env)
        ok += int((env >= thresh).sum())
        m = float(env.min())
        if m < worst:
            worst = m
            a = np.unravel_index(int(env.argmin()), env.shape)
            worst_pt = (float(HH[a]), float(PP[a]), float(q))
    return ok / total, worst, worst_pt


if __name__ == "__main__":
    pass
