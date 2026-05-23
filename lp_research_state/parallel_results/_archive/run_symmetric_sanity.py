"""Sanity test of the assume_even=True modification at small scale."""
import sys, time, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
from white_full_convex import solve_full_program

# small sanity: row6 (h=0, p=0.381, q in [-0.02, 0.02])
PARAMS = dict(h1=0.0, h2=0.0, p1=0.381, p2=0.381, q1=-0.02, q2=0.02, R=8, solver="CLARABEL")

print("=== sanity: row6 small-scale, NON-EVEN, no Bochner ===")
t0 = time.time()
r1 = solve_full_program(N=2000, T=1000, bochner_n=0, assume_even=False, **PARAMS)
print(f"  status={r1['status']}  Omega*={r1['value']}  t={time.time()-t0:.1f}s")

print("=== sanity: row6 small-scale, EVEN, no Bochner ===")
t0 = time.time()
r2 = solve_full_program(N=2000, T=1000, bochner_n=0, assume_even=True, **PARAMS)
print(f"  status={r2['status']}  Omega*={r2['value']}  t={time.time()-t0:.1f}s")

print("=== sanity: row6 small-scale, EVEN, bochner_n=10 ===")
t0 = time.time()
r3 = solve_full_program(N=2000, T=1000, bochner_n=10, assume_even=True, **PARAMS)
print(f"  status={r3['status']}  Omega*={r3['value']}  t={time.time()-t0:.1f}s")

# Verify d ~ 0
import numpy as np
if r2["d"] is not None:
    print(f"  EVEN(no boch): max|d|={np.max(np.abs(r2['d'])):.2e}, max|v-w|={np.max(np.abs(r2['v']-r2['w'])):.2e}")
if r3["d"] is not None:
    print(f"  EVEN(boch=10): max|d|={np.max(np.abs(r3['d'])):.2e}, max|v-w|={np.max(np.abs(r3['v']-r3['w'])):.2e}")
