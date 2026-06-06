"""
STAGE 1 (cont.) — Characterize the GATE regions to scope Stage 2.

For each gate region (where our existing 12-center augmented cover Phi fails to
reach the core headline 0.380284), determine:
  (1) how our Phi varies across the box (which sub-direction drives it down);
  (2) whether SUBDIVISION of the box into sub-boxes lets the EXISTING cover clear
      0.380284 over most of it, isolating the genuinely-hard sub-region;
  (3) for the two near-core q-strips (R16, R17) how fine a q-subdivision lifts them.

This is still pure EVALUATION of the saved duals — NO SDP solves. It tells us, per
region, the minimal Stage-2 work: a few sub-box re-solves vs many.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _fullspace_eval import (load_centers, anchor_value, phi_center_grid,
                             cover_min_over_box, CORE_HEADLINE,
                             WHITE_OUTSIDE_FLOOR, WHITE_TABLE2)

ANCHOR = "primal_m1e5"


def cover_grid_values(centers, hr, pr, q, n_h, n_p):
    """Return the cover (max over centers of Phi) on an (h,p) grid at fixed q."""
    h_grid = np.linspace(*hr, n_h); p_grid = np.linspace(*pr, n_p)
    HH, PP = np.meshgrid(h_grid, p_grid, indexing="ij")
    env = np.full_like(HH, -np.inf)
    for c in centers:
        a = anchor_value(c, ANCHOR)
        F = phi_center_grid(c, a, HH, PP, q)
        np.maximum(env, F, out=env)
    return HH, PP, env


def subdivide_clear_fraction(centers, hr, pr, qr, n_sub_h, n_sub_p, n_sub_q,
                             target=CORE_HEADLINE):
    """Partition the box into n_sub_h x n_sub_p x n_sub_q sub-boxes; for each,
    compute our cover's rigorous min (grid+Lipschitz). Return how many sub-boxes
    ALREADY clear `target` with the existing cover, and the worst sub-box.
    This shows whether the gate is a thin sliver (few hard sub-boxes) or pervasive.
    """
    h_edges = np.linspace(*hr, n_sub_h + 1)
    p_edges = np.linspace(*pr, n_sub_p + 1)
    q_edges = np.linspace(*qr, n_sub_q + 1) if qr[0] != qr[1] else np.array([qr[0], qr[1]])
    total = 0; cleared = 0
    worst = {"lb": np.inf}
    for i in range(len(h_edges) - 1):
        for j in range(len(p_edges) - 1):
            for k in range(max(1, len(q_edges) - 1)):
                hh = (h_edges[i], h_edges[i+1])
                pp = (p_edges[j], p_edges[j+1])
                qq = (q_edges[k], q_edges[k+1]) if len(q_edges) > 1 else qr
                lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
                    centers, ANCHOR, hh, pp, qq, n_h=21, n_p=21, n_q=11)
                total += 1
                if lb >= target:
                    cleared += 1
                elif lb < worst["lb"]:
                    worst = {"lb": lb, "h": hh, "p": pp, "q": qq,
                             "point": pt, "witness": wit}
    return {"total_subboxes": total, "cleared_subboxes": cleared,
            "fraction_cleared": cleared / total, "worst_subbox": worst}


def main():
    centers, config = load_centers()
    print(f"GATE ANALYSIS (anchor={ANCHOR}, target={CORE_HEADLINE})\n")

    out = {"target": CORE_HEADLINE, "anchor": ANCHOR, "gate_detail": []}

    # gate region indices from the main eval (recompute quickly):
    gate_idx = []
    for idx, (hr, pr, qr, wb) in enumerate(WHITE_TABLE2, start=1):
        n_h = 81 if (hr[1]-hr[0]) > 0.05 else 41
        n_p = 161 if (pr[1]-pr[0]) > 0.2 else 81
        n_q = 81 if (qr[1]-qr[0]) > 0.1 else 41
        lb, *_ = cover_min_over_box(centers, ANCHOR, hr, pr, qr,
                                    n_h=n_h, n_p=n_p, n_q=n_q)
        wfloor = 0.37925 if abs(wb-0.37925) < 1e-9 else WHITE_OUTSIDE_FLOOR
        if max(lb, wfloor) < CORE_HEADLINE:
            gate_idx.append(idx)
    print(f"Gate regions (certified < {CORE_HEADLINE}): {gate_idx}\n")

    # Classify each gate region by subdivision behaviour.
    for idx in gate_idx:
        hr, pr, qr, wb = WHITE_TABLE2[idx-1]
        # near-core q-strips (16,17) -> subdivide q finely; wide regions -> subdivide h,p
        wide_p = (pr[1]-pr[0]) > 0.2
        wide_h = (hr[1]-hr[0]) > 0.05
        wide_q = (qr[1]-qr[0]) > 0.1
        nsh = 4 if wide_h else 2
        nsp = 8 if wide_p else 2
        nsq = 8 if wide_q else (4 if (qr[1]-qr[0]) > 1e-6 else 1)
        res = subdivide_clear_fraction(centers, hr, pr, qr, nsh, nsp, nsq)
        w = res["worst_subbox"]
        print(f"[R{idx:2d}] h{hr} p{pr} q{qr}  (White={wb})")
        print(f"      subdiv {nsh}x{nsp}x{nsq}={res['total_subboxes']} sub-boxes: "
              f"{res['cleared_subboxes']} already clear {CORE_HEADLINE} "
              f"({100*res['fraction_cleared']:.0f}%)")
        print(f"      worst sub-box: Phi_min={w['lb']:.5f} at h{w['h']} p{w['p']} "
              f"q{w['q']}  (wit {w['witness']})\n")
        out["gate_detail"].append({
            "region": idx, "h_range": list(hr), "p_range": list(pr),
            "q_range": list(qr), "white_bound": wb,
            "subdiv": [nsh, nsp, nsq], **res,
        })

    OUT = CODE.parent / "parallel_results" / "fullspace_stage1_gate.json"
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"saved -> {OUT}")


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
