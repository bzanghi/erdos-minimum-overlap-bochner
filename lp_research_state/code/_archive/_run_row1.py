import sys, json, time, traceback
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
import warnings; warnings.filterwarnings("ignore")
from white_full_convex import solve_full_program

N = int(sys.argv[1]); T = int(sys.argv[2]); R = int(sys.argv[3]); nb = int(sys.argv[4])
out_path = sys.argv[5]
t0 = time.time()
print(f"[start] N={N} T={T} R={R} nb={nb}", flush=True)
try:
    res = solve_full_program(N=N, T=T, R=R, h1=0.015, h2=0.015,
        p1=0.381, p2=0.381, q1=-0.02, q2=0.02, bochner_n=nb, solver="CLARABEL")
    elapsed = time.time() - t0
    out = {
        "N": N, "T": T, "R": R, "bochner_n": nb,
        "value": float(res["value"]) if res["value"] is not None else None,
        "status": res["status"],
        "time": elapsed,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"[done] {json.dumps(out)}", flush=True)
except Exception:
    traceback.print_exc()
    sys.exit(1)
