"""
PRO-38 R8 promotion — PROBE phase.

R8 box: h in [0,0.08], p(=c1) in [0,1], q(=d1) in [0.05,1.0]  (WIDE).
With the 12 core centers (q_c=0.02) the combined cover Phi-min over R8 is ~0.275,
worst at q=1.0 (con_513 q-decay -0.5*lambda_513*(q^2 - q_c^2) is catastrophic when
the strip's |q| reaches 1.0 while the center sits at q_c=0.02).

Strategy: place fresh AUGMENTED centers whose own q-range sits HIGH on the strip so
Phi has small q-decay over the local sub-strip. Each fresh center is independently
dual-feasible => adding it only raises the cover (validity automatic). Anchor is the
conservative dual-extracted LB minus 1e-5.

This PROBE solves a handful of centers at moderate config to map the landscape:
  - how high is the SDP primal as q_c climbs toward 1.0?
  - how large is lambda_513 (=> how fast does Phi decay off the center in q)?
  - how much of the q-axis can ONE center cover >= 0.380000 / 0.380284?
so we can ESTIMATE centers_needed before committing to a full grid (WIDE => don't
brute force; budget rule).

Uses path_b_analytical.build_problem_with_dual_handles (which hardcodes the corrected
mside_sin_coeff=4.0 on its 5.6/5.7 RHS) + poly_moment, dual-extracted exactly like
_verify_cover_dualext.py / _fullspace_stage2_halo_centers.py.
"""
from __future__ import annotations
import json, sys, time, warnings, argparse
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints
from dual_extractor import solve_with_dual_extraction
from _fullspace_eval import (load_centers, cover_min_over_box, phi_center,
                             CORE_HEADLINE, WHITE_OUTSIDE_FLOOR)

OUT = CODE.parent / "parallel_results" / "fullspace_promote_R8.json"
TARGET = CORE_HEADLINE  # 0.380284


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return res, duals


def make_center_record(label, h_c, p_c, q1, q2, res, duals, dt, config):
    anchor = res["rigorous_dual_LB"]
    return {
        "label": label, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "primal": res["reported_value"], "dual_lb": anchor,
        "anchor_primal_m1e5": (res["reported_value"] - 1e-5) if res["reported_value"] is not None else None,
        "dual_resid": res["dual_residual_at_LB"], "status": res["status"],
        "duals": duals, "time": dt, "config": config,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=5000)
    ap.add_argument("--T", type=int, default=2000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=30)
    ap.add_argument("--pm_k", type=int, default=20)
    # probe centers: list of (h_c, p_c, q1, q2). Single-point q (q1=q2) per task.
    ap.add_argument("--probes", type=str,
                    default="0.0,0.5,1.0,1.0;0.0,0.5,0.5,0.5;0.0,0.5,0.1,0.1")
    args = ap.parse_args()
    config = {"N": args.N, "T": args.T, "R": args.R, "bochner_n": args.bn, "pm_k_max": args.pm_k,
              "mside_sin_coeff": 4.0}

    probes = []
    for tok in args.probes.split(";"):
        parts = [float(x) for x in tok.split(",")]
        probes.append(tuple(parts))

    centers = []
    print(f"=== R8 PROBE: config N={args.N} T={args.T} bn={args.bn} pm_k={args.pm_k} ===\n", flush=True)
    for (h_c, p_c, q1, q2) in probes:
        label = f"R8probe_h{h_c}_c{p_c}_q{q1}_{q2}"
        t0 = time.time()
        res, duals = solve_center(h_c, p_c, q1, q2, args.N, args.T, args.R, args.bn, args.pm_k)
        dt = time.time() - t0
        rec = make_center_record(label, h_c, p_c, q1, q2, res, duals, dt, config)
        centers.append(rec)
        print(f"[{label}] status={res['status']} primal={res['reported_value']} "
              f"dualLB={res['rigorous_dual_LB']} dual_resid={res['dual_residual_at_LB']} "
              f"con_513={duals['con_513']:.4f} ({dt:.0f}s)", flush=True)
        # how far does this center alone reach >= floor along the q-axis (at its own p_c, h=0)?
        if res["rigorous_dual_LB"] is not None:
            anchor = res["reported_value"] - 1e-5  # conservative
            tmpc = {"h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2, "duals": duals}
            qs = np.linspace(0.05, 1.0, 96)
            phis = np.array([phi_center(tmpc, anchor, h_c, p_c, q) for q in qs])
            ge380 = qs[phis >= WHITE_OUTSIDE_FLOOR]
            ge_core = qs[phis >= TARGET]
            def rng(a):
                return (float(a.min()), float(a.max())) if a.size else None
            print(f"        alone @ (h={h_c},p={p_c}): q-range >=0.380000 = {rng(ge380)}; "
                  f">=0.380284 = {rng(ge_core)}; Phi@q=0.05={phis[0]:.5f} Phi@q=1.0={phis[-1]:.5f}",
                  flush=True)
        # save incrementally
        OUT.write_text(json.dumps({"region": 8, "phase": "probe", "target": TARGET,
                                   "config": config, "centers": centers}, indent=2, default=float))

    # combined cover over R8 with core + probes
    core, _ = load_centers()
    for c in centers:
        c["primal_for_anchor"] = c["primal"]
    combined = core + [{"label": c["label"], "h_c": c["h_c"], "p_c": c["p_c"],
                        "q1": c["q1"], "q2": c["q2"], "primal": c["primal"],
                        "duals": c["duals"]} for c in centers if c["primal"] is not None]
    R8 = ((0.0, 0.08), (0.0, 1.0), (0.05, 1.0))
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        combined, "primal_m1e5", *R8, n_h=41, n_p=161, n_q=81)
    print(f"\n=== R8 combined ({len(core)} core + {len(centers)} probe) ===")
    print(f"phi_min_lb={lb:.6f} grid_min={gmin:.6f} eps_grid={eps:.2e} L_max={Lm:.3f}")
    print(f"worst @ (h={pt[0]:.4f}, p={pt[1]:.4f}, q={pt[2]:.4f}) wit={wit}")
    print(f"clears 0.380000: {lb >= WHITE_OUTSIDE_FLOOR}; clears 0.380284: {lb >= TARGET}")


if __name__ == "__main__":
    main()
