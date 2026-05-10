import sys, os, json, time, traceback
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
import warnings; warnings.filterwarnings("ignore")
from white_full_convex import solve_full_program

OUT = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/row7_bochner.json"

PARAMS = dict(h1=0.030, h2=0.030, p1=0.375, p2=0.375, q1=-0.02, q2=0.02,
              bochner_n=20, R=10, solver="CLARABEL")

fallback_chain = []
result = None
for N, T in [(10000, 4000), (5000, 2000), (3000, 1500)]:
    print(f"=== Trying N={N}, T={T} ===", flush=True)
    t0 = time.time()
    try:
        res = solve_full_program(N=N, T=T, **PARAMS)
        elapsed = time.time() - t0
        fallback_chain.append({"N": N, "T": T, "status": str(res["status"]),
                               "value": float(res["value"]) if res["value"] is not None else None,
                               "time": elapsed})
        print(f"   status={res['status']}  value={res['value']}  time={elapsed:.1f}s", flush=True)
        if res["status"] in ("optimal", "optimal_inaccurate") and res["value"] is not None:
            result = {
                "row": "row7", "N": N, "T": T, "R": 10, "bochner_n": 20,
                "value": float(res["value"]), "status": str(res["status"]),
                "time": elapsed, "fallback_chain": fallback_chain,
            }
            break
        else:
            print("   non-optimal, falling back", flush=True)
    except Exception as e:
        elapsed = time.time() - t0
        tb = traceback.format_exc()
        print(f"   exception: {e}\n{tb}", flush=True)
        fallback_chain.append({"N": N, "T": T, "status": f"exception: {e}", "value": None, "time": elapsed})

if result is None:
    result = {"row": "row7", "N": None, "T": None, "R": 10, "bochner_n": 20,
              "value": None, "status": "all_failed", "time": None,
              "fallback_chain": fallback_chain}

with open(OUT, "w") as f:
    json.dump(result, f, indent=2)
print("WROTE", OUT, flush=True)
print(json.dumps(result, indent=2))
