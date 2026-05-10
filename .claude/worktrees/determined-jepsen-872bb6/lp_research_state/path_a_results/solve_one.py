"""Solve one (h, p, qm, qp) Bochner-augmented SDP and append result to JSONL."""
import sys, os, time, warnings, json, argparse
warnings.filterwarnings("ignore")
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp
import io, contextlib


def solve_one(N, T, R, h1, h2, p1, p2, qm, qp, bn):
    Omega, w, v, c, d, eps_v, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, qm, qp, bochner_n=bn)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    buf = io.StringIO()
    t0 = time.time()
    with contextlib.redirect_stdout(buf):
        res = solve_with_dual_extraction(prob)
    last = res['raw_iterations'][-1] if res['raw_iterations'] else None
    LB = (res['reported_value'] - last['gap']) if last else None
    return {"reported": res['reported_value'], "LB": LB,
            "status": res['status'], "time": time.time() - t0,
            "last_gap": last['gap'] if last else None,
            "dual_residual": last['dual_residual'] if last else None}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("--N", type=int, default=3000)
    p.add_argument("--T", type=int, default=1200)
    p.add_argument("--R", type=int, default=10)
    p.add_argument("--bn", type=int, default=30)
    p.add_argument("--h1", type=float, required=True)
    p.add_argument("--h2", type=float, required=True)
    p.add_argument("--p1", type=float, required=True)
    p.add_argument("--p2", type=float, required=True)
    p.add_argument("--qm", type=float, required=True)
    p.add_argument("--qp", type=float, required=True)
    p.add_argument("--label", default="")
    args = p.parse_args()
    r = solve_one(args.N, args.T, args.R, args.h1, args.h2, args.p1, args.p2, args.qm, args.qp, args.bn)
    rec = {"label": args.label, "N": args.N, "T": args.T, "R": args.R, "bn": args.bn,
           "h1": args.h1, "h2": args.h2, "p1": args.p1, "p2": args.p2, "qm": args.qm, "qp": args.qp,
           **r}
    with open(args.out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"DONE {args.label}: LB={r['LB']:.7f} t={r['time']:.1f}s status={r['status']}", flush=True)
