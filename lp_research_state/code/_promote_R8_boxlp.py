"""
PRO-38 R8 promotion via the PROVEN rigorous box-LP subdivision (_fullspace_stage2_solve).

R8 = h[0,0.08] p(=c1)[0,1] q(=d1)[0.05,1.0], a WIDE gate region. The stage2 box-LP
subdivision already certified R1-R5,R10 >= 0.380284 by this exact mechanism but never
processed R8. It is the cleanest rigorous tool here because:
  - it solves the augmented LP over each sub-BOX (range relaxation), so a cleared leaf's
    dual_LB - margin is a RIGOROUS LB over the whole sub-box (no Lipschitz/Phi-extension);
  - clean 'infeasible' leaves are vacuously covered (no admissible f with |d1| that large /
    c1 that large given sum_squares<=0.5 and |d|<=2/pi) -- the established convention
    (is_infeasible: unambiguous 'infeasible' only; solver_failed/ambiguous => re-solve).

This driver runs ONLY R8, writes the leaf tree + region verdict to fullspace_promote_R8.json,
and reports the independent floor = min over CLEARED leaves of cert_val (infeasible leaves
contribute +inf). If certified, OUR program (not White's 0.38) certifies R8 >= 0.380284.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _fullspace_stage2_solve import certify_region, TARGET

OUT = CODE.parent / "parallel_results" / "fullspace_promote_R8.json"
R8 = ((0.0, 0.08), (0.0, 1.0), (0.05, 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=3000)
    ap.add_argument("--T", type=int, default=1200)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bn", type=int, default=30)
    ap.add_argument("--pm_k", type=int, default=20)
    ap.add_argument("--N2", type=int, default=8000)
    ap.add_argument("--T2", type=int, default=3000)
    ap.add_argument("--bn2", type=int, default=40)
    ap.add_argument("--pm_k2", type=int, default=20)
    ap.add_argument("--escalate_depth", type=int, default=3)
    ap.add_argument("--max_depth", type=int, default=6)
    ap.add_argument("--dual_margin", type=float, default=1e-5)
    args = ap.parse_args()

    cfg = {"N": args.N, "T": args.T, "R": args.R, "bn": args.bn, "pm_k": args.pm_k,
           "N2": args.N2, "T2": args.T2, "bn2": args.bn2, "pm_k2": args.pm_k2,
           "escalate_depth": args.escalate_depth}
    hr, pr, qr = R8
    print(f"=== R8 box-LP subdivision (TARGET={TARGET}) ===")
    print(f"box h{hr} p{pr} q{qr}; light N={args.N} bn={args.bn}, escalate>=d{args.escalate_depth} "
          f"to N={args.N2} bn={args.bn2}; max_depth={args.max_depth}\n", flush=True)

    logs = []
    def log(s, _l=logs):
        print(s, flush=True); _l.append(s)

    t0 = time.time()
    certified, leaves, residual = certify_region(8, hr, pr, qr, cfg, args.dual_margin,
                                                 args.max_depth, log)
    dt = time.time() - t0

    cleared = [lf for lf in leaves if lf["verdict"] == "cleared"]
    infeas = [lf for lf in leaves if lf["verdict"] == "infeasible"]
    min_cleared = min((lf["cert_val"] for lf in cleared), default=None)
    out = {
        "region": 8, "method": "rigorous box-LP subdivision (augmented White program)",
        "anchor": "dual_LB - dual_margin (rigorous dual objective)",
        "target": TARGET, "h_range": list(hr), "p_range": list(pr), "q_range": list(qr),
        "config": cfg, "dual_margin": args.dual_margin, "max_depth": args.max_depth,
        "certified_ge_target": bool(certified),
        "n_leaves": len(leaves), "n_infeasible": len(infeas), "n_cleared": len(cleared),
        "n_residual": len(residual),
        "min_cleared_cert": (float(min_cleared) if min_cleared is not None else None),
        "independent_floor": (float(min_cleared) if min_cleared is not None else None),
        "residual_gates": [{"box": r["box"], "cert_val": r.get("cert_val"),
                            "status": r.get("status")} for r in residual],
        "leaves": leaves, "time_s": dt, "log": logs,
    }
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"\n=== R8 verdict: certified_ge_{TARGET} = {certified} ===")
    print(f"leaves={len(leaves)} infeasible={len(infeas)} cleared={len(cleared)} residual={len(residual)}")
    print(f"min cleared cert (independent floor) = {min_cleared}")
    print(f"time {dt:.0f}s; saved -> {OUT}")


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    main()
