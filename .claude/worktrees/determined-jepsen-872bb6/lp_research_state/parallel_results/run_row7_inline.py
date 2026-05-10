import sys, os, json, time, traceback
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
import warnings; warnings.filterwarnings("ignore")
from white_full_convex import solve_full_program
OUT = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/row7_bochner.json"
PARAMS = dict(h1=0.030, h2=0.030, p1=0.375, p2=0.375, q1=-0.02, q2=0.02, bochner_n=20, R=10, solver="CLARABEL")
N = int(sys.argv[1]); T = int(sys.argv[2])
print(f"=== N={N}, T={T} ===", flush=True)
t0 = time.time()
try:
    res = solve_full_program(N=N, T=T, **PARAMS)
    elapsed = time.time() - t0
    payload = {"N": N, "T": T, "status": str(res["status"]), "value": float(res["value"]) if res["value"] is not None else None, "time": elapsed}
except Exception as e:
    elapsed = time.time() - t0
    payload = {"N": N, "T": T, "status": f"exception:{type(e).__name__}", "value": None, "time": elapsed, "err": str(e)}
    traceback.print_exc()
print(json.dumps(payload, indent=2), flush=True)
log_path = OUT + ".attempts"
arr = []
if os.path.exists(log_path):
    with open(log_path) as f: arr = json.load(f)
arr.append(payload)
with open(log_path, "w") as f: json.dump(arr, f, indent=2)
