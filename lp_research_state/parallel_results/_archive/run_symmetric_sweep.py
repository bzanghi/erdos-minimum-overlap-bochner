"""Conditional bound assuming f* is even: sweep all 7 White Table-3 rows.

Solves the SDP with assume_even=True (d=0, dlt=0, v=w) at all 7 ellipse
centers. Reports MIN over rows = conditional lower bound on µ.

NOTE: For even f we have ∫ x f(x) dx = 0 so h_1 = 0. Rows with h>0 are
infeasible under the even assumption. We record those as +∞ and take MIN
only over feasible rows (effectively the h=0, q=0 axis).
"""
import sys, os, json, time, traceback, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
from white_full_convex import solve_full_program, WHITE_TABLE3
import numpy as np

OUT = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/symmetric_conditional.json"

# Choose bochner_n. With assume_even, only c-block matters in the Bochner PSD
# (since d=0); the matrix size is the same but the problem is much sparser.
# Try n=30 first; fall back if memory/time issues.
BOCHNER_LEVELS = [30, 25, 20]

# Solver tries
SOLVER = "CLARABEL"
R = 10
N, T = 10000, 4000

results = []
for (h, p, qm, qp, label) in WHITE_TABLE3:
    print(f"\n========= {label}: h={h:.3f} p={p:.4f} q∈[{qm:.2f},{qp:.2f}] =========", flush=True)
    PARAMS = dict(h1=h, h2=h, p1=p, p2=p, q1=qm, q2=qp, R=R, solver=SOLVER)
    chain = []
    found = None
    for bn in BOCHNER_LEVELS:
        for (NN, TT) in [(N, T), (5000, 2000), (3000, 1500)]:
            print(f"  trying  N={NN} T={TT} bochner_n={bn} ...", flush=True)
            t0 = time.time()
            try:
                res = solve_full_program(N=NN, T=TT, bochner_n=bn,
                                          assume_even=True, **PARAMS)
                el = time.time() - t0
                v = res["value"] if res["value"] is not None else None
                st = str(res["status"])
                print(f"    status={st}  Omega*={v}  t={el:.1f}s", flush=True)
                chain.append({"N": NN, "T": TT, "bochner_n": bn, "status": st,
                              "value": float(v) if v is not None else None, "time": el})
                if st == "infeasible":
                    # Even f forces h_1 = 0; row is infeasible, treat as +∞.
                    found = {
                        "row": label, "h": h, "p": p, "q1": qm, "q2": qp,
                        "N": NN, "T": TT, "R": R, "bochner_n": bn,
                        "value": float("inf"), "status": "infeasible_under_even",
                        "time": el, "chain": chain,
                    }
                    break
                if st in ("optimal", "optimal_inaccurate") and v is not None:
                    # numerical sanity: d should be ~0
                    dmax = float(np.max(np.abs(res["d"]))) if res["d"] is not None else None
                    vwmax = (float(np.max(np.abs(res["v"] - res["w"])))
                              if (res["v"] is not None and res["w"] is not None) else None)
                    found = {
                        "row": label, "h": h, "p": p, "q1": qm, "q2": qp,
                        "N": NN, "T": TT, "R": R, "bochner_n": bn,
                        "value": float(v), "status": st, "time": el,
                        "max_abs_d": dmax, "max_abs_v_minus_w": vwmax,
                        "chain": chain,
                    }
                    break
            except Exception as e:
                el = time.time() - t0
                print(f"    EXC: {e}  t={el:.1f}s", flush=True)
                chain.append({"N": NN, "T": TT, "bochner_n": bn,
                              "status": f"exception: {e}", "value": None, "time": el})
        if found is not None:
            break
    if found is None:
        found = {"row": label, "h": h, "p": p, "q1": qm, "q2": qp,
                 "value": None, "status": "all_failed", "chain": chain}
    # JSON does not allow Infinity; substitute string for storage.
    if found.get("value") == float("inf"):
        found_for_json = {**found, "value": "infeasible"}
    else:
        found_for_json = found
    results.append(found)
    # Save partial after each row
    res_for_json = []
    for r in results:
        if r.get("value") == float("inf"):
            res_for_json.append({**r, "value": "infeasible"})
        else:
            res_for_json.append(r)
    with open(OUT, "w") as f:
        json.dump({
            "description": ("CONDITIONAL bound on Erdős minimum overlap constant µ "
                            "ASSUMING the optimal f* is even (f(x)=f(-x)). "
                            "Modified SDP: d_k=0, dlt=0, v_j=w_j for all k,j."),
            "interpretation": ("If conjecture (f* even) holds: µ ≥ MIN over rows. "
                                "If MIN > 0.379005, either µ > White's bound OR f* is not even."),
            "white_full": 0.379005,
            "bochner_non_conditional": 0.3799,
            "rows": res_for_json,
        }, f, indent=2)
    print(f"  -> {label}: Omega*={found['value']}", flush=True)

# Aggregate
vals = [r["value"] for r in results if r["value"] is not None and r["value"] != float("inf")]
mn = min(vals) if vals else None
infeasible_rows = [r["row"] for r in results if r.get("status") == "infeasible_under_even"]
print(f"  infeasible-under-even rows: {infeasible_rows}")
print("\n========= SUMMARY =========")
for r in results:
    print(f"  {r['row']}: Omega*={r['value']}  (status={r.get('status')}, "
          f"bochner_n={r.get('bochner_n')}, N={r.get('N')}, T={r.get('T')})")
print(f"  MIN over rows (CONDITIONAL bound assuming f* even) = {mn}")
print(f"  White's full bound: 0.379005")
print(f"  Bochner non-conditional bound: 0.3799")

# Final write with min
res_for_json = []
for r in results:
    if r.get("value") == float("inf"):
        res_for_json.append({**r, "value": "infeasible"})
    else:
        res_for_json.append(r)
with open(OUT, "w") as f:
    json.dump({
        "description": ("CONDITIONAL bound on Erdős minimum overlap constant µ "
                        "ASSUMING the optimal f* is even (f(x)=f(-x)). "
                        "Modified SDP: d_k=0, dlt=0, v_j=w_j for all k,j."),
        "interpretation": ("If conjecture (f* even) holds: µ ≥ MIN over rows. "
                            "If MIN > 0.379005, either µ > White's bound OR f* is not even."),
        "note_h1": ("For even f, ∫ x f(x) dx = 0 = h_1, so rows with h>0 are "
                     "infeasible under the even assumption (recorded as 'infeasible'). "
                     "The MIN is taken over feasible rows only."),
        "white_full": 0.379005,
        "bochner_non_conditional": 0.3799,
        "min_over_rows_conditional": mn,
        "infeasible_rows_under_even": infeasible_rows,
        "rows": res_for_json,
    }, f, indent=2)
print(f"\nWROTE {OUT}")
