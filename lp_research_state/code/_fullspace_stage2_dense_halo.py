"""
STAGE 2 (dense halo cover) — production N=20000 Phi-centers for the near-core halo.

Established facts:
  * Far regions R1-R5, R10 are CERTIFIED >= TARGET by box-LP subdivision
    (_fullspace_stage2_solve.py).
  * The near-core HALO (the c1~[0.31,0.45], small |d1| in [0.02,0.05] parts of
    R7, R9, R16, R17, plus the |d1| in [0.05,0.1] slivers of R6, R8) cannot be
    cleared by box-LP at feasible compute (range relaxation conservatism ~1.5e-3).
  * But SINGLE-POINT N=20000 floors in the halo are SAFELY above target:
    (h=0.003,c1=0.3915,d1=-0.025)=0.381140; (h=0.008,c1=0.3835,d1=-0.025)=0.380853;
    (h=0,c1=0.39,d1=-0.025)=0.381050. So the true floor ~0.3808-0.3811 >> 0.380284;
    a DENSE N=20000 Phi-cover (like the core's 12 centers) should reach it.

This script builds that dense cover: production N=20000 centers at h=0 across
c1 in [0.33,0.45] (extendable to [0.25,0.5]) at two q-ranges -- a NARROW one
[-0.025,0.025] for the inner halo (R16/R17 + d1<=0.025 of R7/R9) and a WIDER one
[-0.05,0.05] for the outer halo (d1 in [0.025,0.05] of R7/R9). It then evaluates the
COMBINED cover (existing 12 + these) over each halo region with the rigorous
grid+Lipschitz box-min (_fullspace_eval.cover_min_over_box), and CDE-iterates: add a
center at the current worst point until the region's box-min >= TARGET or it stalls.

Anchors: conservative dual-extracted (anchor = dual_LB - 1e-5).

NOTE on h: the halo's worst h is h=0 (concave in h; con_54>0). Centers at h=0 anchor
the worst E(M); Phi decays slightly toward h=0.08 (~ -7e-4) but the existing core
centers (at h up to 0.03) and the far-region box-LP cover larger h. We evaluate the
halo over the regions' actual h-ranges so the grid+Lipschitz captures any h-decay.
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
from _fullspace_eval import load_centers, cover_min_over_box

OUT = CODE.parent / "parallel_results" / "fullspace_stage2_halo_centers.json"
TARGET = 0.380284

# Halo regions and their boxes. For R7/R9 we evaluate only the c1 sub-range that
# box-LP cannot clear (the rest is certified separately by box-LP); the inner c1
# clears at light (we restrict to [0.31,0.45]).
HALO_EVAL = {
    16: ((0.0, 0.06), (0.33, 0.45), (-0.025, -0.02)),
    17: ((0.0, 0.06), (0.33, 0.45), (0.02, 0.025)),
    7:  ((0.0, 0.08), (0.31, 0.45), (-0.05, -0.025)),
    9:  ((0.0, 0.08), (0.31, 0.45), (0.025, 0.05)),
    # R6/R8 near-core slivers (|d1| in [0.05,0.1], c1 in [0.25,0.5]).
    61: ((0.0, 0.08), (0.25, 0.5), (-0.1, -0.05)),
    81: ((0.0, 0.08), (0.25, 0.5), (0.05, 0.1)),
}


def solve_center(h_c, p_c, q1, q2, N, T, R, bn, pm_k):
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bn)
    if pm_k > 0:
        pm_cons, _ = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    try:
        res = solve_with_dual_extraction(prob)
    except Exception as e:
        return None, None, f"{type(e).__name__}: {e}"
    duals = {k: (float(H[k].dual_value) if H[k].dual_value is not None else 0.0)
             for k in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                       "con_512_qL", "con_512_qU", "con_513")}
    return res, duals, None


def to_cover_center(label, h_c, p_c, q1, q2, anchor, duals):
    # primal := anchor so cover_min_over_box's 'primal_m1e5' uses anchor-1e-5.
    # We set anchor := dual_LB (then -1e-5 in eval gives dual_LB-1e-5, conservative).
    return {"label": label, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
            "primal": anchor, "dual_lb": anchor, "duals": duals}


def eval_region(combined, hr, pr, qr):
    n_h = 41
    n_p = 161 if (pr[1] - pr[0]) > 0.2 else 101
    n_q = 41
    return cover_min_over_box(combined, "primal_m1e5", hr, pr, qr,
                             n_h=n_h, n_p=n_p, n_q=n_q)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=40)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--margin", type=float, default=1e-5)
    ap.add_argument("--max_cde_iters", type=int, default=10)
    # initial dense seed grid (c1 values), per q-range band
    ap.add_argument("--c1_narrow", type=str, default="0.34,0.37,0.40,0.43")
    ap.add_argument("--c1_wide", type=str, default="0.34,0.39,0.44")
    args = ap.parse_args()

    existing, _ = load_centers()
    halo = []
    if OUT.exists():
        try:
            prev = json.load(open(OUT))
            for c in prev.get("centers", []):
                # normalize to cover-center schema
                halo.append(to_cover_center(c["label"], c["h_c"], c["p_c"],
                                            c["q1"], c["q2"], c["dual_lb"], c["duals"]))
        except Exception:
            halo = []
    have = {(round(c["h_c"], 5), round(c["p_c"], 5), round(c["q1"], 5), round(c["q2"], 5))
            for c in halo}
    raw_records = []  # for persistence with full primal/dual info

    def add_center(h_c, p_c, q1, q2, tagprefix):
        key = (round(h_c, 5), round(p_c, 5), round(q1, 5), round(q2, 5))
        if key in have:
            return None
        t0 = time.time()
        res, duals, err = solve_center(h_c, p_c, q1, q2, args.N, args.T, args.R, args.bn, args.pm_k)
        dt = time.time() - t0
        if res is None or res["rigorous_dual_LB"] is None:
            print(f"  [skip] center (h={h_c:.4f},c1={p_c:.4f},q[{q1},{q2}]): "
                  f"{(err or res['status'])} ({dt:.0f}s)", flush=True)
            return None
        anchor = res["rigorous_dual_LB"]
        lab = f"{tagprefix}_h{h_c:.3f}_c{p_c:.4f}_q{q1}_{q2}"
        cc = to_cover_center(lab, h_c, p_c, q1, q2, anchor, duals)
        halo.append(cc); have.add(key)
        raw_records.append({"label": lab, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
                            "primal": res["reported_value"], "dual_lb": anchor,
                            "dual_resid": res["dual_residual_at_LB"],
                            "status": res["status"], "duals": duals, "time": dt})
        print(f"  [add] center (h={h_c:.4f},c1={p_c:.4f},q[{q1},{q2}]): "
              f"primal={res['reported_value']:.6f} dualLB={anchor:.6f} "
              f"con513={duals['con_513']:.4f} ({dt:.0f}s)", flush=True)
        persist(raw_records, existing, halo, {})
        return cc

    print(f"=== Dense halo cover (N={args.N}, bn={args.bn}) ===\n", flush=True)

    # --- Phase A: seed dense grid ---
    print("Phase A: seed narrow-q [-0.025,0.025] grid (inner halo) ...", flush=True)
    for p in [float(x) for x in args.c1_narrow.split(",")]:
        add_center(0.0, p, -0.025, 0.025, "narrow")
    print("Phase A: seed wide-q [-0.05,0.05] grid (outer halo) ...", flush=True)
    for p in [float(x) for x in args.c1_wide.split(",")]:
        add_center(0.0, p, -0.05, 0.05, "wide")

    # --- Phase B: evaluate + CDE iterate per region ---
    def eval_all():
        out = {}
        for rid, (hr, pr, qr) in HALO_EVAL.items():
            combined = existing + halo
            lb, pt, wit, gmin, eps, Lm = eval_region(combined, hr, pr, qr)
            out[rid] = {"lb": lb, "pt": pt, "wit": wit, "eps": eps, "L": Lm,
                        "hr": hr, "pr": pr, "qr": qr}
        return out

    print("\nPhase B: evaluate combined cover over halo regions + CDE iterate\n", flush=True)
    region_eval = eval_all()
    for rid in sorted(region_eval):
        e = region_eval[rid]
        print(f"  [R{rid}] init box-min={e['lb']:.6f} "
              f"{'OK' if e['lb'] >= TARGET else 'BELOW'} worst@(h={e['pt'][0]:.4f},"
              f"c1={e['pt'][1]:.4f},d1={e['pt'][2]:.4f}) wit={e['wit']}", flush=True)

    for it in range(1, args.max_cde_iters + 1):
        # find the worst region below target
        worst_rid = None; worst_lb = np.inf
        for rid, e in region_eval.items():
            if e["lb"] < worst_lb:
                worst_lb = e["lb"]; worst_rid = rid
        if worst_lb >= TARGET:
            print(f"\n[CDE] all halo regions >= TARGET (min={worst_lb:.6f}). Done.", flush=True)
            break
        e = region_eval[worst_rid]
        h_w, c_w, d_w = e["pt"]
        # place a center at the worst point's (c1) at h=0, q-range matched to region band.
        qr = e["qr"]
        # choose center q-range: narrow if |d1|<=0.025 else wide covering the band
        if abs(qr[0]) <= 0.0251 and abs(qr[1]) <= 0.0251:
            cq1, cq2 = -0.025, 0.025
        elif max(abs(qr[0]), abs(qr[1])) <= 0.0501:
            cq1, cq2 = -0.05, 0.05
        else:
            cq1, cq2 = -0.1, 0.1
        print(f"\n[CDE iter {it}] worst R{worst_rid} lb={worst_lb:.6f} at c1={c_w:.4f}; "
              f"add center (h=0, c1={c_w:.4f}, q[{cq1},{cq2}])", flush=True)
        cc = add_center(0.0, round(c_w, 4), cq1, cq2, f"cde{it}")
        if cc is None:
            # center already exists or failed -> nudge c1 slightly to break stall
            cc = add_center(0.0, round(c_w + 0.005, 4), cq1, cq2, f"cde{it}b")
            if cc is None:
                print(f"[CDE iter {it}] could not add a new center; stalling.", flush=True)
                break
        region_eval = eval_all()
        for rid in sorted(region_eval):
            ev = region_eval[rid]
            print(f"  [R{rid}] box-min={ev['lb']:.6f} "
                  f"{'OK' if ev['lb'] >= TARGET else 'BELOW'} worst@(h={ev['pt'][0]:.4f},"
                  f"c1={ev['pt'][1]:.4f},d1={ev['pt'][2]:.4f}) wit={ev['wit']}", flush=True)

    persist(raw_records, existing, halo, region_eval)
    halo_min = min(e["lb"] for e in region_eval.values())
    halo_min_rid = min(region_eval, key=lambda r: region_eval[r]["lb"])
    print(f"\n=== Dense halo cover result ===")
    print(f"  halo min over regions = {halo_min:.6f} at R{halo_min_rid}")
    print(f"  all halo regions >= {TARGET}: {halo_min >= TARGET}")
    print(f"  n halo centers = {len(halo)}")


def persist(raw_records, existing, halo, region_eval):
    out = {
        "method": "dense production Phi-cover (halo) + combined box-min eval + CDE",
        "target": TARGET, "anchor": "dual_LB - 1e-5 (conservative)",
        "n_existing": len(existing), "n_halo": len(halo),
        "centers": raw_records,
        "region_eval": {str(rid): {"lb": region_eval[rid]["lb"],
                                   "worst_point": region_eval[rid]["pt"],
                                   "witness": region_eval[rid]["wit"],
                                   "eps_grid": region_eval[rid]["eps"],
                                   "h_range": list(region_eval[rid]["hr"]),
                                   "p_range": list(region_eval[rid]["pr"]),
                                   "q_range": list(region_eval[rid]["qr"])}
                        for rid in region_eval},
    }
    OUT.write_text(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
