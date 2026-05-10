"""Parameterized driver: dual-extract a rigorous LB for one (row, N, T, R, bochner_n) config.

Usage:
    python3 run_dual_extract.py row6 3000 1200 10 30
    python3 run_dual_extract.py row1 3000 1200 10 30
    ...

Persists JSON to lp_research_state/dual_extract_<row>_N<N>_n<n>.json
"""
import sys, os, time, json, warnings
warnings.filterwarnings("ignore")

STATE_DIR = os.environ.get(
    "ERDOS_STATE_DIR",
    os.path.dirname(os.path.abspath(__file__)),
)
sys.path.insert(0, STATE_DIR + "/code")

from white_full_convex import build_problem, WHITE_TABLE3
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp

POINTS_BY_ROW = {label: (h, p, q1, q2) for (h, p, q1, q2, label) in WHITE_TABLE3}

def main():
    if len(sys.argv) < 6:
        print("usage: run_dual_extract.py <row> <N> <T> <R> <bochner_n>")
        sys.exit(2)
    row = sys.argv[1]
    N = int(sys.argv[2]); T = int(sys.argv[3]); R = int(sys.argv[4]); n_b = int(sys.argv[5])

    h, p, qm, qp = POINTS_BY_ROW[row]
    print(f"Building {row} (h={h}, p={p}, q∈[{qm},{qp}]) N={N} T={T} R={R} bochner_n={n_b} ...", flush=True)
    Omega, w, v, c, d, eps_v, dlt, cons = build_problem(
        N, T, R, h, h, p, p, qm, qp, bochner_n=n_b,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)

    # The most precise rigorous LB is: reported_value - last_gap
    # (CLARABEL's printed dual_obj column is only ~5 sig figs;
    #  prob.value and the gap column carry full ~10 sig fig precision).
    last_gap = None
    if res["raw_iterations"]:
        last_gap = res["raw_iterations"][-1]["gap"]

    rigorous_LB_value_minus_gap = None
    if res["reported_value"] is not None and last_gap is not None:
        rigorous_LB_value_minus_gap = res["reported_value"] - last_gap

    out = {
        "row": row, "N": N, "T": T, "R": R,
        "h": h, "p": p, "q_range": [qm, qp],
        "bochner_n": n_b, "use_T5p": False,
        "kind": "lp_run_bochner_dual",
        "value": res["reported_value"],
        "status": res["status"],
        "rigorous_dual_LB_low_precision": res["rigorous_dual_LB"],
        "rigorous_dual_LB": rigorous_LB_value_minus_gap,
        "last_gap": last_gap,
        "dual_residual_at_LB": res["dual_residual_at_LB"],
        "best_iter": res["best_iter"],
        "n_iters_total": res["n_iters_total"],
        "time": res["time"],
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    out_path = f"{STATE_DIR}/dual_extract_{row}_N{N}_n{n_b}.json"
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=2)

    print(f"  reported_value         = {out['value']:.10g}")
    print(f"  rigorous_dual_LB (val-gap) = {rigorous_LB_value_minus_gap:.10g}"
          if rigorous_LB_value_minus_gap is not None else "  rigorous_dual_LB = None")
    print(f"  last_gap               = {last_gap:.3e}" if last_gap else "")
    print(f"  status                 = {out['status']}")
    print(f"  iters                  = {out['n_iters_total']}")
    print(f"  time                   = {out['time']:.1f}s")
    print(f"  wrote {out_path}")

if __name__ == "__main__":
    main()
