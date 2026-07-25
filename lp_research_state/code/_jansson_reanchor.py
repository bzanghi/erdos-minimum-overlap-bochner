"""
_jansson_reanchor.py — certify AND re-anchor the 12 core (5.16) envelope centers.

THE PROBLEM WITH THE CURRENT ANCHORS
------------------------------------
The core floor 0.3802838 is  min over the (h,p) box of  max_c Phi_c, with

    Phi_c(theta) = anchor_c + shift_c(theta),     anchor_c := primal_c - 1e-5

where `primal_c` is CLARABEL's reported objective at center c and `shift_c` is
White's exact dual-objective transport.  Two things are wrong with that anchor:

 1. It is not certified.  `primal_c` is a solver-reported value and
    `dual_extractor.py`'s "rigorous_dual_LB" subtracts the printed gap without
    correcting for dual infeasibility (its own docstring concedes a margin
    "can be absorbed" but none ever is).  The 1e-5 haircut is a convention, not
    a theorem.
 2. It is loose.  The Jansson interval-arithmetic certificate `p_lo` at these
    centers comes out ABOVE `primal_c - 1e-5` — at row4 by +4.2e-5, at row1 by
    +7.0e-4.  Every bit of that is thrown away.

THE FIX
-------
Anchor each center at its OWN Jansson `p_lo`, using the duals read from the
SAME solve.  Validity, in full:

  Let z be the conic dual from the solve at center theta_c, and
      p_lo = -b(theta_c)^T z  +  pen_zs  -  pen_Dx
  the Jansson bound (pen_zs <= 0 from the cone-distance term, pen_Dx >= 0 from
  the dual defect D = c + A^T z).  Then for any primal-feasible x,
      c^T x  >=  -b(theta)^T z + pen_zs(theta) - pen_Dx .
  Moving theta within the core box changes ONLY the right-hand side b of the
  four parameter-dependent constraints (5.3), (5.4), (5.12), (5.13) — the
  matrix A and the objective c are untouched — so:
      * -b(theta)^T z  =  -b(theta_c)^T z + shift_c(theta),  exactly, with
        shift_c the same transport `path_b_analytical.dual_objective_shift`
        already uses;
      * pen_Dx = sum_i |D_i| * xbar_i is theta-INDEPENDENT: D = c + A^T z does
        not involve b, and the xbar_i are model box bounds (Omega <= 1,
        w,v in [0,1], |c_k|,|d_k| <= 2/pi, fixed tail caps);
      * pen_zs = sum_j min(0, lambda_min(z_j)) * sbar_j is theta-INDEPENDENT
        WHEN IT VANISHES, i.e. when z lies exactly in the dual cone so every
        lambda_min(z_j) >= 0.  This script ASSERTS pen_zs == 0 at every center
        and refuses to emit a re-anchored center otherwise.  (In the runs
        recorded in docs/RND_WHITESPACE/L2_PROD.json it is identically 0.)

  Hence  Phi_c(theta) := p_lo_c + shift_c(theta)  <=  c^T x  for every feasible
  x at parameters theta — a rigorous lower bound on mu conditional on theta,
  which is exactly what the cover needs.

So this is not a re-tuning: it replaces an uncertified convention by a
certified quantity, and the certified quantity happens to be larger.

OUTPUT
------
`parallel_results/jansson_core12_reanchored.json` in the same shape as
`phase5_N20K_bn40_dualext.json` (label / h_c / p_c / q1 / q2 / duals) plus
`p_lo`, so the existing evaluator can consume it with anchor = p_lo.

    ../../.venv/bin/python _jansson_reanchor.py --all
    ../../.venv/bin/python _jansson_reanchor.py --evaluate
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pathlib

warnings.filterwarnings("ignore")

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

REPO = CODE.parent.parent
DUALEXT = CODE.parent / "parallel_results" / "phase5_N20K_bn40_dualext.json"
OUT_JSON = CODE.parent / "parallel_results" / "jansson_core12_reanchored.json"
# overridable so a re-run at a different N writes its own file (see --out)

DUAL_KEYS = ["con_53", "con_54", "con_512_pL", "con_512_pU",
             "con_512_qL", "con_512_qU", "con_513"]

H_BOX = (0.0, 0.06)
P_BOX = (0.35, 0.45)


def load_anchors():
    d = json.loads(DUALEXT.read_text())
    return d["centers"], d["config"]


def load_db():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"meta": {
        "charter": "Jansson p_lo re-anchoring of the 12 core (5.16) envelope centers",
        "validity": "Phi_c(theta) = p_lo_c + shift_c(theta); requires penalty_zs == 0",
        "old_anchor_convention": "primal - 1e-5 (uncertified)",
        "new_anchor": "p_lo (Jansson interval-arithmetic certificate)",
    }, "centers": []}


def save_db(db):
    OUT_JSON.write_text(json.dumps(db, indent=2, default=float))


def solve_certify_and_read_duals(N, T, R, h_c, p_c, q1, q2, bn, pm, slack_infl):
    """Run the Jansson certificate AND read the matched model-constraint duals
    from the SAME solve.  Returns (p_lo, duals, diag)."""
    import cvxpy as cp
    from _jansson_verify import (
        build_center_problem, canonical_x_bounds, split_cone_blocks,
        cone_lambda_min_lower, slack_size_upper, _iv,
    )
    import mpmath
    from mpmath import iv
    import scipy.sparse as sp

    prob, H, _pm_tb = build_center_problem(N, T, R, h_c, p_c, q1, q2, bn, pm)
    # solver_opts={} (not the default None) so the chain can be inverted below:
    # clarabel_conif.invert reads inverse_data.solver_options unconditionally.
    data, chain, inv = prob.get_problem_data(cp.CLARABEL, solver_opts={})
    A = data["A"]; b = np.asarray(data["b"]); c = np.asarray(data["c"])
    dims = data["dims"]

    t0 = time.time()
    sol = chain.solve_via_data(prob, data)
    solve_s = time.time() - t0
    x = np.asarray(sol.x); z = np.asarray(sol.z); s = np.asarray(sol.s)

    # populate .dual_value on the model constraints from THIS solve
    prob.unpack_results(sol, chain, inv)
    duals = {}
    for k in DUAL_KEYS:
        if k in H and H[k].dual_value is not None:
            dv = H[k].dual_value
            duals[k] = float(np.asarray(dv).reshape(-1)[0]) if np.ndim(dv) else float(dv)
        else:
            duals[k] = 0.0

    # ---- Jansson penalty terms (interval arithmetic) ----
    xbar, xflags, _ = canonical_x_bounds(prob, H, N, T, R, x_solution=x)
    Acsc = sp.csc_matrix(A)
    z_iv = [_iv(float(zr)) for zr in z]
    pen_Dx = iv.mpf(0)
    for i in range(A.shape[1]):
        s0, e0 = Acsc.indptr[i], Acsc.indptr[i + 1]
        acc = _iv(float(c[i]))
        for r_, v_ in zip(Acsc.indices[s0:e0].tolist(), Acsc.data[s0:e0].tolist()):
            acc = acc + _iv(float(v_)) * z_iv[r_]
        pen_Dx = pen_Dx + abs(acc) * _iv(float(xbar[i]))
    pen_Dx_hi = pen_Dx.b

    z_blocks = split_cone_blocks(z, dims)
    s_blocks = split_cone_blocks(s, dims)
    pen_zs = iv.mpf(0)
    worst_d = None
    for (zk, zp, zn), (sk, sp_, sn) in zip(z_blocks, s_blocks):
        if zk == "zero":
            continue
        d_j, _ = cone_lambda_min_lower(zk, zp, zn)
        sbar_j, _ = slack_size_upper(sk, sp_, sn, infl=slack_infl)
        d_neg = min(mpmath.mpf(0), mpmath.mpf(d_j))
        worst_d = float(d_j) if worst_d is None else min(worst_d, float(d_j))
        pen_zs = pen_zs + _iv(float(d_neg)) * _iv(float(sbar_j))
    pen_zs_lo = pen_zs.a

    neg_bz = _iv(0)
    for i in range(len(b)):
        neg_bz = neg_bz - _iv(float(b[i])) * z_iv[i]
    p_lo_iv = neg_bz + pen_zs - pen_Dx
    p_lo = float(p_lo_iv.a)

    diag = {
        "obj_val": float(sol.obj_val),
        "obj_val_dual": float(sol.obj_val_dual),
        "neg_bz_lo": float(neg_bz.a),
        "penalty_Dx_upper": float(pen_Dx_hi),
        "penalty_zs_lower": float(pen_zs_lo),
        "penalty_zs_is_zero": bool(float(pen_zs_lo) == 0.0),
        "worst_cone_lambda_min_lower": worst_d,
        "r_prim": float(sol.r_prim), "r_dual": float(sol.r_dual),
        "status": str(sol.status),
        "unmapped_cols": xflags.get("unmapped_cols"),
        "solve_s": solve_s,
    }
    return p_lo, duals, diag


def run_one(center, N, T, R, bn, pm, slack_infl):
    db = load_db()
    cfg = dict(N=N, T=T, R=R, bochner_n=bn, pm_k_max=pm)
    # key the skip on the CONFIG too — otherwise a re-run at a different N is
    # silently skipped and the old N's numbers are reported as the new N's
    if any(r.get("label") == center["label"] and r.get("ok")
           and r.get("config") == cfg for r in db["centers"]):
        print(f"[skip] {center['label']} already done at {cfg}", flush=True)
        return
    lbl = center["label"]
    print(f"\n### {lbl}  h_c={center['h_c']} p_c={center['p_c']} "
          f"N={N} T={T} bn={bn} pm={pm} ###", flush=True)
    t0 = time.time()
    try:
        p_lo, duals, diag = solve_certify_and_read_duals(
            N, T, R, center["h_c"], center["p_c"], center["q1"], center["q2"],
            bn, pm, slack_infl)
        old_anchor = center["primal"] - 1e-5
        rec = {
            "label": lbl, "ok": True, "config": cfg,
            "h_c": center["h_c"], "p_c": center["p_c"],
            "q1": center["q1"], "q2": center["q2"],
            "p_lo": p_lo,
            "old_anchor": old_anchor,
            "gain_vs_old_anchor": p_lo - old_anchor,
            "duals": duals,
            "old_duals": center["duals"],
            "diag": diag,
            "wall_s": time.time() - t0,
        }
        flag = "" if diag["penalty_zs_is_zero"] else "  *** pen_zs != 0: NOT USABLE ***"
        print(f"[{lbl}] p_lo={p_lo:.12f}  old_anchor={old_anchor:.12f}  "
              f"gain={p_lo - old_anchor:+.3e}{flag}", flush=True)
    except Exception as e:
        rec = {"label": lbl, "ok": False, "config": cfg,
               "error": f"{type(e).__name__}: {e}",
               "tb": traceback.format_exc(), "wall_s": time.time() - t0}
        print(f"[{lbl}] ERROR: {e}", flush=True)
    db = load_db()
    db["centers"].append(rec)
    save_db(db)


def envelope_floor(centers, anchors, n_grid=4001):
    """min over the core (h,p) box of max_c (anchor_c + shift_c), with the
    rigorous Lipschitz cell-error term.  Same convention as
    _fullspace_eval.reproduce_core_headline."""
    from path_b_analytical import find_ellipse_h_p
    h_grid = np.linspace(*H_BOX, n_grid)
    p_grid = np.linspace(*P_BOX, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    wit = np.empty(HH.shape, dtype=object)
    L_max = 0.0
    for c, a in zip(centers, anchors):
        syn = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
               "value": a}
        e = find_ellipse_h_p(syn, c["duals"], c["q1"], c["q2"], target=0.379005)
        F = (a + e["const_q"] + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        m = F > env
        env[m] = F[m]; wit[m] = c["label"]
        lam = lambda c2, c1, lo, hi: max(abs(2 * c2 * lo + c1), abs(2 * c2 * hi + c1))
        L_max = max(L_max, float(np.hypot(lam(e["A_h2"], e["A_h1"], *H_BOX),
                                          lam(e["A_p2"], e["A_p1"], *P_BOX))))
    gmin = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    cell_h = (H_BOX[1] - H_BOX[0]) / (n_grid - 1)
    cell_p = (P_BOX[1] - P_BOX[0]) / (n_grid - 1)
    eps_grid = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    return {"rigorous_LB": gmin - eps_grid, "grid_min": gmin,
            "eps_grid": eps_grid, "L_max": L_max, "witness": str(wit[arg]),
            "binding_point": [float(HH[arg]), float(PP[arg])]}


def _box_min(centers, anchors, h_box, p_box, n_grid):
    """(grid_min - eps, grid_min, eps, L_max, witness, argmin) over one (h,p) box.

    Identical convention to _fullspace_eval.reproduce_core_headline, but with the
    box as a parameter so it can be subdivided.  L_max is recomputed PER BOX,
    which is the whole point: a smaller box has both a smaller half-diagonal and
    (usually) a smaller local gradient bound.
    """
    from path_b_analytical import find_ellipse_h_p
    h_grid = np.linspace(*h_box, n_grid)
    p_grid = np.linspace(*p_box, n_grid)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    wit = np.empty(HH.shape, dtype=object)
    L_max = 0.0
    for c, a in zip(centers, anchors):
        syn = {"h_c": c["h_c"], "p_c": c["p_c"], "q1": c["q1"], "q2": c["q2"],
               "value": a}
        e = find_ellipse_h_p(syn, c["duals"], c["q1"], c["q2"], target=0.379005)
        F = (a + e["const_q"] + e["A_h2"] * HH * HH + e["A_h1"] * HH + e["A_h0"]
             + e["A_p2"] * PP * PP + e["A_p1"] * PP + e["A_p0"])
        m = F > env
        env[m] = F[m]; wit[m] = c["label"]
        lam = lambda c2, c1, lo, hi: max(abs(2 * c2 * lo + c1), abs(2 * c2 * hi + c1))
        L_max = max(L_max, float(np.hypot(lam(e["A_h2"], e["A_h1"], *h_box),
                                          lam(e["A_p2"], e["A_p1"], *p_box))))
    gmin = float(env.min())
    arg = np.unravel_index(int(env.argmin()), env.shape)
    cell_h = (h_box[1] - h_box[0]) / (n_grid - 1)
    cell_p = (p_box[1] - p_box[0]) / (n_grid - 1)
    eps = L_max * 0.5 * float(np.hypot(cell_h, cell_p))
    return (gmin - eps, gmin, eps, L_max, str(wit[arg]),
            (float(HH[arg]), float(PP[arg])))


def envelope_floor_adaptive(centers, anchors, target, n_grid=801,
                            max_depth=7, h_box=None, p_box=None, depth=0):
    """Rigorous box-min of the envelope over the core box, by subdivision.

    A single global grid pays `eps = L_max * half_diag` everywhere, and L_max is
    set by the steepest center anywhere on the box.  Splitting the longest axis
    shrinks both factors locally, so the certified floor rises toward the true
    grid minimum.  Recursion stops once a sub-box clears `target` — so the
    result is "at least target" when it clears, and the true tightest bound this
    method can give when it does not.

    Returns (floor, worst_box_report).
    """
    h_box = h_box or H_BOX
    p_box = p_box or P_BOX
    lb, gmin, eps, L, wit, pt = _box_min(centers, anchors, h_box, p_box, n_grid)
    if lb >= target or depth >= max_depth:
        return lb, {"h_box": h_box, "p_box": p_box, "lb": lb, "grid_min": gmin,
                    "eps": eps, "L_max": L, "witness": wit, "argmin": pt,
                    "depth": depth}
    if (h_box[1] - h_box[0]) / (H_BOX[1] - H_BOX[0]) >= \
       (p_box[1] - p_box[0]) / (P_BOX[1] - P_BOX[0]):
        mid = 0.5 * (h_box[0] + h_box[1])
        parts = [((h_box[0], mid), p_box), ((mid, h_box[1]), p_box)]
    else:
        mid = 0.5 * (p_box[0] + p_box[1])
        parts = [(h_box, (p_box[0], mid)), (h_box, (mid, p_box[1]))]
    worst, rep = np.inf, None
    for hb, pb in parts:
        l, r = envelope_floor_adaptive(centers, anchors, target, n_grid,
                                       max_depth, hb, pb, depth + 1)
        if l < worst:
            worst, rep = l, r
    return worst, rep


def emit_dualext(out_path=None):
    """Write a drop-in replacement for phase5_N20K_bn40_dualext.json carrying the
    certified anchors, so every downstream evaluator can be pointed at it via
    $LP_DUALEXT with no code change.

    `primal` is set to `p_lo + 1e-5` ON PURPOSE: the evaluators call
    anchor_value(., 'primal_m1e5') = primal - 1e-5, so this makes their anchor
    come out at exactly `p_lo`.  The honest value is also stored under `p_lo`,
    and `anchor_value(., 'p_lo')` reads it directly — prefer that where the call
    site can be changed.
    """
    out_path = Path(out_path or (CODE.parent / "parallel_results" /
                                 "phase5_N20K_bn40_dualext_reanchored.json"))
    old_centers, cfg = load_anchors()
    db = load_db()
    done = {r["label"]: r for r in db["centers"] if r.get("ok")}
    missing = [c["label"] for c in old_centers if c["label"] not in done]
    if missing:
        raise SystemExit(f"cannot emit: {len(missing)} centers uncertified: {missing}")
    bad = [l for l, r in done.items() if not r["diag"]["penalty_zs_is_zero"]]
    if bad:
        raise SystemExit(f"cannot emit: penalty_zs != 0 at {bad}")

    order = [c["label"] for c in old_centers] + \
            [l for l in done if l not in {c["label"] for c in old_centers}]
    centers = []
    for lbl in order:
        r = done[lbl]
        centers.append({
            "label": lbl, "h_c": r["h_c"], "p_c": r["p_c"],
            "q1": r["q1"], "q2": r["q2"],
            "p_lo": r["p_lo"],
            "primal": r["p_lo"] + 1e-5,   # so primal_m1e5 == p_lo; see docstring
            "dual_lb": r["p_lo"],
            "duals": r["duals"],
            "_anchor_source": "jansson p_lo, duals from the same solve",
        })
    out = {"config": cfg,
           "_note": ("Jansson-re-anchored core centers. 'primal' is synthetic "
                     "(= p_lo + 1e-5) so anchor_value(.,'primal_m1e5') yields the "
                     "certified p_lo; the real certificate is the 'p_lo' field."),
           "centers": centers}
    Path(out_path).write_text(json.dumps(out, indent=2, default=float))
    print(f"emitted {len(centers)} re-anchored centers -> {out_path}")
    print(f"use with:  LP_DUALEXT={out_path} ../../.venv/bin/python <evaluator>.py")
    return out_path


def evaluate():
    old_centers, cfg = load_anchors()
    db = load_db()
    done = {r["label"]: r for r in db["centers"] if r.get("ok")}

    missing = [c["label"] for c in old_centers if c["label"] not in done]
    unusable = [l for l, r in done.items() if not r["diag"]["penalty_zs_is_zero"]]
    print(f"centers certified: {len(done)}/{len(old_centers)}")
    if missing:
        print(f"  MISSING: {missing}")
    if unusable:
        print(f"  UNUSABLE (pen_zs != 0): {unusable}")

    base = envelope_floor(old_centers, [c["primal"] - 1e-5 for c in old_centers])
    print(f"\n[baseline]  old anchors (primal-1e-5), old duals")
    print(f"   rigorous_LB = {base['rigorous_LB']:.10f}  grid_min={base['grid_min']:.10f} "
          f"eps={base['eps_grid']:.2e} L_max={base['L_max']:.4f} "
          f"witness={base['witness']} @ {base['binding_point']}")

    if missing or unusable:
        print("\n[re-anchored] cannot evaluate until every center is certified "
              "and has penalty_zs == 0.")
        return None

    # the original 12, plus any FRESH center certified since (--new-center)
    original = [c["label"] for c in old_centers]
    extra = [l for l in done if l not in original]
    if extra:
        print(f"  plus {len(extra)} fresh center(s): {extra}")
    use = [done[l] for l in original] + [done[l] for l in extra]
    new_centers = [{"label": r["label"], "h_c": r["h_c"], "p_c": r["p_c"],
                    "q1": r["q1"], "q2": r["q2"], "duals": r["duals"]} for r in use]
    new_anchors = [r["p_lo"] for r in use]
    new = envelope_floor(new_centers, new_anchors)
    print(f"\n[re-anchored]  anchors = Jansson p_lo, duals from the SAME solve")
    print(f"   rigorous_LB = {new['rigorous_LB']:.10f}  grid_min={new['grid_min']:.10f} "
          f"eps={new['eps_grid']:.2e} L_max={new['L_max']:.4f} "
          f"witness={new['witness']} @ {new['binding_point']}")
    print(f"\n   CORE FLOOR {base['rigorous_LB']:.10f} -> {new['rigorous_LB']:.10f}"
          f"   ({new['rigorous_LB'] - base['rigorous_LB']:+.3e})")

    print(f"\n{'center':<24} {'old anchor':>14} {'p_lo':>16} {'gain':>12}")
    for r in use:
        oa = r["old_anchor"]
        if oa == float("-inf"):
            print(f"{r['label']:<24} {'(fresh)':>14} {r['p_lo']:>16.12f} {'--':>12}")
        else:
            print(f"{r['label']:<24} {oa:>14.9f} {r['p_lo']:>16.12f} "
                  f"{r['gain_vs_old_anchor']:>+12.3e}")

    db["evaluation"] = {"baseline": base, "reanchored": new}
    save_db(db)
    return new


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--center", type=str, default=None)
    ap.add_argument("--new-center", type=str, default=None,
                    help="LABEL,h_c,p_c[,q1,q2] — certify a FRESH center (e.g. at the "
                         "current binding point) and add it to the envelope")
    ap.add_argument("--evaluate", action="store_true")
    ap.add_argument("--emit-dualext", action="store_true",
                    help="write a drop-in dualext file with the certified anchors")
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=40)
    ap.add_argument("--pm_k_max", type=int, default=20)
    ap.add_argument("--slack_infl", type=float, default=4.0)
    ap.add_argument("--out", type=str, default=None,
                    help="output JSON (default jansson_core12_reanchored.json)")
    args = ap.parse_args()
    if args.out:
        OUT_JSON = pathlib.Path(args.out)

    if args.emit_dualext:
        emit_dualext()
    elif args.evaluate:
        evaluate()
    elif args.new_center:
        parts = args.new_center.split(",")
        if len(parts) not in (3, 5):
            raise SystemExit("--new-center wants LABEL,h_c,p_c[,q1,q2]")
        c = {"label": parts[0], "h_c": float(parts[1]), "p_c": float(parts[2]),
             "q1": float(parts[3]) if len(parts) == 5 else -0.02,
             "q2": float(parts[4]) if len(parts) == 5 else 0.02,
             # a fresh center has no prior anchor; -inf makes `gain` meaningless
             # rather than fake, and keeps it out of the baseline comparison
             "primal": float("-inf"), "duals": {}}
        run_one(c, args.N, args.T, args.R, args.bochner_n, args.pm_k_max,
                args.slack_infl)
    else:
        centers, cfg = load_anchors()
        todo = centers if args.all else [c for c in centers if c["label"] == args.center]
        if not todo:
            raise SystemExit(f"no such center; have {[c['label'] for c in centers]}")
        for c in todo:
            run_one(c, args.N, args.T, args.R, args.bochner_n, args.pm_k_max,
                    args.slack_infl)
        evaluate()
