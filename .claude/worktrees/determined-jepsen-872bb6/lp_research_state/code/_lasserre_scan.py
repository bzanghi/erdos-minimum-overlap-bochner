"""Lasserre (T_max, T_loc) scan at row 4 N=2000, T=200 — sequential, writes per-config JSON."""
import warnings; warnings.filterwarnings('ignore')
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from white_full_convex import solve_full_program

if len(sys.argv) < 5:
    print("usage: <Tmax> <Tloc> <N> <T> [out]")
    sys.exit(1)

Tmax = int(sys.argv[1])
Tloc = int(sys.argv[2])
N = int(sys.argv[3])
T = int(sys.argv[4])
R = 10
out = sys.argv[5] if len(sys.argv) > 5 else f"/tmp/lasscan_T{Tmax}_L{Tloc}_N{N}_T{T}.json"

H, P, QM, QP = 0.004, 0.3875, -0.02, 0.02
t0 = time.time()
r = solve_full_program(N, T, R, H, H, P, P, QM, QP, bochner_n=20,
                       lasserre_T_max=Tmax, lasserre_T_loc=Tloc)
elapsed = time.time() - t0
res = {"value": float(r["value"]) if r["value"] is not None else None,
       "status": r["status"], "time": elapsed,
       "T_max": Tmax, "T_loc": Tloc, "N": N, "T": T}
with open(out, "w") as f:
    json.dump(res, f, indent=2)
print(json.dumps(res))
