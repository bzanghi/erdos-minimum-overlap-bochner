"""
STAGE 2 — Full-space promotion of the Erdos minimum-overlap LOWER bound.

GOAL: lift every one of White's 18 Table-2 "outside" regions to a certified
lower bound >= TARGET (= 0.380284, the conservative core headline), so that
combined with the core (5.16) cover (>= 0.3802838) the bound mu >= 0.380284 holds
UNCONDITIONALLY over the full (E(M), c1, d1) parameter space.

------------------------------------------------------------------------------
METHOD (rigorous box-LP subdivision -- White's own framework, augmented)
------------------------------------------------------------------------------
Stage 1 evaluated the EXISTING 12-center cover's quadratic dual-shift Phi over the
boxes and found 12 "gate" regions where Phi dips below TARGET. Investigation
(see FULLSPACE_PROMOTION_STAGE2.md) showed the dips occur at LARGE |d1| corners
(c1~0, |d1| up to 1) which are INFEASIBLE moment triples -- and that the
quadratic Phi from any single center necessarily DECAYS in q (the (5.13)
multiplier con_513 ~ 0.3 gives a -0.5*con_513*(q^2 - q_c^2) penalty), so no set
of centers can lift Phi at those wide-q corners. A center solved with a matching
wide q-range has primal ~0.374 (the (5.13) rhs -0.5*max(q^2) loosens the program),
structurally below TARGET. So the Phi-center-addition route CANNOT close the wide-q
gate regions.

The correct rigorous tool is the SAME one White uses, applied per sub-box and
augmented: solve the augmented LP over a sub-box's RANGES (h in [h1,h2], c1 in
[p1,p2], d1 in [q1,q2]). cvxpy `build_problem_with_dual_handles(N,T,R, h1,h2,
p1,p2, q1,q2, bochner_n, +poly_moment)` builds exactly this: con_53 uses h1,
con_54 uses h2, con_512_pL/pU pin c1 in [p1,p2], con_512_qL/qU pin d1 in [q1,q2],
con_513 uses max(p^2), max(q^2). Minimizing Omega over this program is a RIGOROUS
LOWER BOUND on mu for every admissible function whose moments lie in the sub-box
(weak duality + the augmentation constraints -- Bochner PSD, even-moment-nonneg --
are NECESSARY conditions satisfied by every 0<=f<=1).

Per sub-box outcome:
  * status infeasible  -> NO admissible function has moments in this sub-box
                          (the augmentation constraints are necessary conditions,
                          so relaxation-infeasible => truly empty). EXCLUDED.
  * optimal, val >= TARGET -> sub-box cleared (mu >= val >= TARGET there).
  * optimal, val <  TARGET -> subdivide further (the box-min over a wide range is
                          conservative -- it uses the loosest constraints over the
                          whole range; splitting tightens it). If still < TARGET at
                          max depth -> GENUINE residual gate (report precisely).

A region is CERTIFIED >= TARGET iff every leaf sub-box is infeasible or >= TARGET.
Because the true extremal function's moment triple is feasible and lies in exactly
one leaf, that leaf is the ">= TARGET" kind, so mu >= TARGET. This is fully rigorous
and needs NO Phi-extrapolation.

RIGOR NOTES
  * We use the augmented program (bochner_n + even-moment-nonneg poly-moment cuts).
    poly_moment.even_moment_tail_bound is the rigorous tail (j_part=200000 + analytic
    remainder 4k/(pi^2 j_part)) -- verified by reading poly_moment.py for this run.
  * The reported per-sub-box value is cvxpy's primal `prob.value`. For a clean
    rigorous margin we ALSO dual-extract (verbose CLARABEL parse) and certify the
    sub-box with the DUAL objective `dual_LB` (<= true LP opt always). We clear a
    sub-box only if dual_LB - dual_margin >= TARGET (dual_margin absorbs the tiny
    dual residual). The primal is used only to flag sub-boxes for subdivision.
  * Infeasibility is detected from cvxpy status; we additionally require the solve
    to be unambiguous (status contains 'infeasible' and not 'inaccurate' for the
    accept; 'infeasible_inaccurate' triggers a re-solve at higher fidelity).

CONFIG: light & parallel-safe by default (N=3000, T=1200, bn=30, pm_k=20). The
wide far regions clear with huge margin at this scale (feasible vals 0.39-0.53), so
light is sufficient; production escalation only if a leaf is borderline.

OUTPUT: lp_research_state/parallel_results/fullspace_stage2_centers.json
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction

REPO = CODE.parent.parent
OUT = CODE.parent / "parallel_results" / "fullspace_stage2_centers.json"

TARGET = 0.380284
CORE_HEADLINE = 0.380284

# White Table-2 outside regions (transcription identical to _fullspace_eval.py).
# (h_range, p_range, q_range, white_bound, region_id)
WHITE_TABLE2 = [
    ((0.75, 2.0),  (0.0, 1.0),   (-1.0, 1.0),     0.38,    1),
    ((0.4, 0.75),  (0.0, 1.0),   (-1.0, 1.0),     0.38,    2),
    ((0.2, 0.4),   (0.0, 1.0),   (-1.0, 1.0),     0.38,    3),
    ((0.1, 0.2),   (0.0, 1.0),   (-1.0, 1.0),     0.38,    4),
    ((0.08, 0.1),  (0.0, 1.0),   (-1.0, 1.0),     0.38,    5),
    ((0.0, 0.08),  (0.0, 1.0),   (-1.0, -0.05),   0.38,    6),
    ((0.0, 0.08),  (0.0, 1.0),   (-0.05, -0.025), 0.38,    7),
    ((0.0, 0.08),  (0.0, 1.0),   (0.05, 1.0),     0.38,    8),
    ((0.0, 0.08),  (0.0, 1.0),   (0.025, 0.05),   0.38,    9),
    ((0.0, 0.08),  (0.0, 0.25),  (-0.025, 0.025), 0.38,    10),
    ((0.0, 0.08),  (0.25, 0.3),  (-0.025, 0.025), 0.38,    11),
    ((0.0, 0.08),  (0.3, 0.33),  (-0.025, 0.025), 0.38,    12),
    ((0.0, 0.08),  (0.5, 1.0),   (-0.025, 0.025), 0.38,    13),
    ((0.0, 0.08),  (0.45, 0.5),  (-0.025, 0.025), 0.38,    14),
    ((0.06, 0.08), (0.33, 0.45), (-0.025, 0.025), 0.38,    15),
    ((0.0, 0.06),  (0.33, 0.45), (-0.025, -0.02), 0.38,    16),
    ((0.0, 0.06),  (0.33, 0.45), (0.02, 0.025),   0.38,    17),
    ((0.0, 0.06),  (0.33, 0.35), (-0.02, 0.02),   0.37925, 18),
]

# Stage-1 gate regions (the ones needing Stage-2 work). The other 6 already clear
# from the existing cover's Phi alone (Stage 1). We re-certify the gate ones here.
GATE_REGIONS = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 17}


def box_lp_solve(h1, h2, p1, p2, q1, q2, N, T, R, bn, pm_k, extract_dual=True):
    """Solve the augmented LP over the sub-box ranges. Returns dict with status,
    primal value, and (if extract_dual) the dual objective LB and 7 dual handles."""
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h1, h2, p1, p2, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        if extract_dual:
            res = solve_with_dual_extraction(prob)
            status = res["status"]
            primal = res["reported_value"]
            dual_LB = res["rigorous_dual_LB"]
            dual_resid = res["dual_residual_at_LB"]
            elapsed = res["time"]
        else:
            t0 = time.time()
            prob.solve(solver="CLARABEL", verbose=False)
            elapsed = time.time() - t0
            status = prob.status
            primal = float(prob.value) if (prob.value is not None and np.isfinite(prob.value)) else None
            dual_LB = None
            dual_resid = None
    except Exception as e:
        # CLARABEL can fail numerically on a degenerate sub-box. Treat as a
        # "solver_failed" outcome -> caller subdivides (or escalates). NOT an
        # infeasibility claim (those must come from a clean 'infeasible' status).
        return {
            "status": "solver_failed", "primal": None, "dual_LB": None,
            "dual_resid": None, "duals": {}, "time": 0.0,
            "box": [h1, h2, p1, p2, q1, q2], "error": f"{type(e).__name__}: {e}",
        }
    duals = {}
    for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
              "con_512_qL", "con_512_qU", "con_513"):
        dv = H[k].dual_value
        duals[k] = float(dv) if dv is not None else 0.0
    return {
        "status": status, "primal": primal, "dual_LB": dual_LB,
        "dual_resid": dual_resid, "duals": duals, "time": elapsed,
        "box": [h1, h2, p1, p2, q1, q2],
    }


def is_infeasible(status):
    """Rigorous infeasibility: accept only an unambiguous 'infeasible' status.
    'infeasible_inaccurate' is treated as ambiguous (caller re-solves / subdivides)."""
    if status is None:
        return False
    s = str(status).lower()
    return ("infeasible" in s) and ("inaccurate" not in s)


def is_ambiguous_infeasible(status):
    s = str(status).lower() if status is not None else ""
    return "infeasible" in s and "inaccurate" in s


def certify_region(region_id, hr, pr, qr, cfg, dual_margin, max_depth, log):
    """Adaptively subdivide a region until every leaf is infeasible or dual_LB>=TARGET.
    Returns (certified: bool, leaves: list, residual_gaps: list)."""
    leaves = []
    residual = []
    # work queue of boxes with a depth counter; split the longest *relevant* axis.
    # Seed subdivision: split q first for wide-q regions (the hard axis), then p.
    h0, h1 = hr; p0, p1 = pr; q0, q1 = qr
    # initial coarse seed depending on region shape
    wide_q = (q1 - q0) > 0.1
    wide_p = (p1 - p0) > 0.2
    wide_h = (h1 - h0) > 0.05
    # seed edges
    def seed_edges(lo, hi, n):
        return list(np.linspace(lo, hi, n + 1))
    nq = 6 if wide_q else 1
    np_ = 3 if wide_p else 1
    nh = 2 if wide_h else 1
    h_e = seed_edges(h0, h1, nh)
    p_e = seed_edges(p0, p1, np_)
    q_e = seed_edges(q0, q1, nq)
    queue = []
    for hi in range(len(h_e) - 1):
        for pi in range(len(p_e) - 1):
            for qi in range(len(q_e) - 1):
                queue.append((h_e[hi], h_e[hi+1], p_e[pi], p_e[pi+1],
                              q_e[qi], q_e[qi+1], 0))
    while queue:
        h1b, h2b, p1b, p2b, q1b, q2b, depth = queue.pop(0)
        # On deep / borderline boxes use higher fidelity.
        N, T, bn, pm_k = cfg["N"], cfg["T"], cfg["bn"], cfg["pm_k"]
        if depth >= cfg.get("escalate_depth", 3):
            N, T, bn, pm_k = cfg["N2"], cfg["T2"], cfg["bn2"], cfg["pm_k2"]
        r = box_lp_solve(h1b, h2b, p1b, p2b, q1b, q2b, N, T, cfg["R"], bn, pm_k,
                         extract_dual=True)
        cert_val = (r["dual_LB"] - dual_margin) if r["dual_LB"] is not None else None
        box_str = (f"R{region_id} d{depth} h[{h1b:.4f},{h2b:.4f}] "
                   f"p[{p1b:.3f},{p2b:.3f}] q[{q1b:.4f},{q2b:.4f}]")
        if str(r["status"]) == "solver_failed":
            # numerical failure: retry once at production fidelity, else subdivide
            r2 = box_lp_solve(h1b, h2b, p1b, p2b, q1b, q2b, cfg["N2"], cfg["T2"],
                              cfg["R"], cfg["bn2"], cfg["pm_k2"], extract_dual=True)
            if str(r2["status"]) != "solver_failed":
                r = r2
                cert_val = (r["dual_LB"] - dual_margin) if r["dual_LB"] is not None else None
            else:
                log(f"  {box_str}: SOLVER_FAILED (retry failed) -> subdivide")
                if depth < max_depth:
                    pass  # fall through to subdivision below (cert_val stays None)
                else:
                    leaves.append({**r, "verdict": "residual_gate",
                                   "cert_val": None, "depth": depth})
                    residual.append({**r, "cert_val": None, "depth": depth})
                    continue
        if is_infeasible(r["status"]):
            log(f"  {box_str}: INFEAS (excluded)")
            leaves.append({**r, "verdict": "infeasible", "depth": depth})
            continue
        if is_ambiguous_infeasible(r["status"]) and depth < max_depth:
            # ambiguous: re-solve at production fidelity once before deciding
            r2 = box_lp_solve(h1b, h2b, p1b, p2b, q1b, q2b, cfg["N2"], cfg["T2"],
                              cfg["R"], cfg["bn2"], cfg["pm_k2"], extract_dual=True)
            if is_infeasible(r2["status"]):
                log(f"  {box_str}: INFEAS@prod (excluded)")
                leaves.append({**r2, "verdict": "infeasible", "depth": depth})
                continue
            r = r2
            cert_val = (r["dual_LB"] - dual_margin) if r["dual_LB"] is not None else None
        if cert_val is not None and cert_val >= TARGET:
            log(f"  {box_str}: dualLB={r['dual_LB']:.6f} cert={cert_val:.6f} OK")
            leaves.append({**r, "verdict": "cleared", "cert_val": cert_val, "depth": depth})
            continue
        # below target -> subdivide if depth budget remains
        if depth < max_depth:
            # split the axis that is widest *and* relevant; prefer q, then p, then h
            lens = []
            if (q2b - q1b) > 1e-9:
                lens.append(("q", q2b - q1b))
            if (p2b - p1b) > 1e-9:
                lens.append(("p", p2b - p1b))
            if (h2b - h1b) > 1e-9:
                lens.append(("h", h2b - h1b))
            axis = max(lens, key=lambda t: t[1])[0]
            if axis == "q":
                qm = 0.5 * (q1b + q2b)
                queue.append((h1b, h2b, p1b, p2b, q1b, qm, depth + 1))
                queue.append((h1b, h2b, p1b, p2b, qm, q2b, depth + 1))
            elif axis == "p":
                pm = 0.5 * (p1b + p2b)
                queue.append((h1b, h2b, p1b, pm, q1b, q2b, depth + 1))
                queue.append((h1b, h2b, pm, p2b, q1b, q2b, depth + 1))
            else:
                hm = 0.5 * (h1b + h2b)
                queue.append((h1b, hm, p1b, p2b, q1b, q2b, depth + 1))
                queue.append((hm, h2b, p1b, p2b, q1b, q2b, depth + 1))
            cv = f"{cert_val:.6f}" if cert_val is not None else "None"
            log(f"  {box_str}: cert={cv} < TARGET -> split {axis}")
        else:
            cv = f"{cert_val:.6f}" if cert_val is not None else "None"
            log(f"  {box_str}: cert={cv} < TARGET at MAX DEPTH -> RESIDUAL GATE")
            leaves.append({**r, "verdict": "residual_gate",
                           "cert_val": cert_val, "depth": depth})
            residual.append({**r, "cert_val": cert_val, "depth": depth})
    certified = (len(residual) == 0)
    return certified, leaves, residual


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", type=str, default="",
                    help="comma-separated region ids to process (default: all gate regions)")
    ap.add_argument("--N", type=int, default=3000)
    ap.add_argument("--T", type=int, default=1200)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=30)
    ap.add_argument("--pm_k", type=int, default=20)
    # production escalation config (for deep/borderline leaves)
    ap.add_argument("--N2", type=int, default=10000)
    ap.add_argument("--T2", type=int, default=4000)
    ap.add_argument("--bn2", type=int, default=40)
    ap.add_argument("--pm_k2", type=int, default=20)
    ap.add_argument("--escalate_depth", type=int, default=3)
    ap.add_argument("--max_depth", type=int, default=6)
    ap.add_argument("--dual_margin", type=float, default=1e-5,
                    help="subtracted from dual_LB before certifying (absorbs dual residual)")
    args = ap.parse_args()

    if args.regions:
        want = set(int(x) for x in args.regions.split(","))
    else:
        want = set(GATE_REGIONS)

    cfg = {"N": args.N, "T": args.T, "R": args.R, "bn": args.bn, "pm_k": args.pm_k,
           "N2": args.N2, "T2": args.T2, "bn2": args.bn2, "pm_k2": args.pm_k2,
           "escalate_depth": args.escalate_depth}

    # load existing results (append/merge so multiple invocations accumulate)
    existing = {}
    if OUT.exists():
        try:
            prev = json.load(open(OUT))
            for r in prev.get("regions", []):
                existing[r["region"]] = r
        except Exception:
            existing = {}

    print(f"=== STAGE 2: rigorous box-LP subdivision (TARGET={TARGET}) ===")
    print(f"config light: N={args.N} T={args.T} bn={args.bn} pm_k={args.pm_k}; "
          f"escalate>=d{args.escalate_depth} to N={args.N2} bn={args.bn2}")
    print(f"regions to process: {sorted(want)}\n", flush=True)

    region_results = dict(existing)
    t_start = time.time()
    for hr, pr, qr, wb, rid in WHITE_TABLE2:
        if rid not in want:
            continue
        print(f"--- Region {rid}: h{hr} p{pr} q{qr} (White {wb}) ---", flush=True)
        logs = []
        def log(s, _logs=logs):
            print(s, flush=True); _logs.append(s)
        t0 = time.time()
        certified, leaves, residual = certify_region(
            rid, hr, pr, qr, cfg, args.dual_margin, args.max_depth, log)
        dt = time.time() - t0
        min_cleared = min((lf.get("cert_val", np.inf) for lf in leaves
                           if lf["verdict"] == "cleared"), default=None)
        region_results[rid] = {
            "region": rid, "h_range": list(hr), "p_range": list(pr),
            "q_range": list(qr), "white_bound": wb,
            "certified_ge_target": bool(certified),
            "target": TARGET,
            "n_leaves": len(leaves),
            "n_infeasible": sum(1 for lf in leaves if lf["verdict"] == "infeasible"),
            "n_cleared": sum(1 for lf in leaves if lf["verdict"] == "cleared"),
            "n_residual": len(residual),
            "min_cleared_cert": (float(min_cleared) if min_cleared is not None else None),
            "residual_gates": [{"box": rr["box"], "cert_val": rr["cert_val"],
                                "primal": rr["primal"], "status": rr["status"]}
                               for rr in residual],
            "leaves": leaves,
            "time_s": dt,
            "config": cfg,
        }
        status_str = "CERTIFIED >= TARGET" if certified else f"RESIDUAL ({len(residual)} gaps)"
        print(f"  => Region {rid}: {status_str}  "
              f"(min cleared cert={min_cleared}, {len(leaves)} leaves, {dt:.0f}s)\n",
              flush=True)
        # persist incrementally after each region
        save(region_results, args)

    save(region_results, args)
    print(f"\n=== Stage-2 done in {time.time()-t_start:.0f}s ===")
    allcert = all(region_results[r]["certified_ge_target"] for r in want
                  if r in region_results)
    print(f"all processed gate regions certified >= {TARGET}: {allcert}")
    for r in sorted(region_results):
        rr = region_results[r]
        print(f"  R{r}: {'OK' if rr['certified_ge_target'] else 'RESIDUAL'} "
              f"(inf={rr['n_infeasible']} cleared={rr['n_cleared']} "
              f"resid={rr['n_residual']} mincert={rr['min_cleared_cert']})")


def save(region_results, args):
    out = {
        "method": "rigorous box-LP subdivision (augmented White program)",
        "target": TARGET,
        "anchor": "dual_LB - dual_margin (rigorous dual objective)",
        "config": {"N": args.N, "T": args.T, "R": args.R, "bn": args.bn,
                   "pm_k": args.pm_k, "N2": args.N2, "T2": args.T2, "bn2": args.bn2,
                   "pm_k2": args.pm_k2, "dual_margin": args.dual_margin,
                   "max_depth": args.max_depth},
        "regions": [region_results[r] for r in sorted(region_results)],
    }
    OUT.write_text(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
