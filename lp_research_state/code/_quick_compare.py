"""
Quick comparison: bochner-only vs Lasserre-2 (existing) vs Lasserre-2 (independent).
Small problem to fit in budget.
"""
from __future__ import annotations
import time, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np
import cvxpy as cp
import importlib.util as _ilu
from pathlib import Path

HERE = Path("/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")
def _load(name, fname):
    spec = _ilu.spec_from_file_location(name, HERE / fname)
    mod = _ilu.module_from_spec(spec); spec.loader.exec_module(mod); return mod

wfc = _load("wfc", "white_full_convex.py")
las_e = _load("las_e", "lasserre.py")
las_i = _load("las_i", "lasserre_independent.py")

mode = sys.argv[1] if len(sys.argv) > 1 else "all"
N    = int(sys.argv[2]) if len(sys.argv) > 2 else 200
T    = int(sys.argv[3]) if len(sys.argv) > 3 else 40
R    = int(sys.argv[4]) if len(sys.argv) > 4 else 5
bn   = int(sys.argv[5]) if len(sys.argv) > 5 else 8
Tmx  = int(sys.argv[6]) if len(sys.argv) > 6 else 4
Tlc  = int(sys.argv[7]) if len(sys.argv) > 7 else 4
h    = float(sys.argv[8]) if len(sys.argv) > 8 else 0.004
p    = float(sys.argv[9]) if len(sys.argv) > 9 else 0.3875
qm   = float(sys.argv[10]) if len(sys.argv) > 10 else -0.02
qp   = float(sys.argv[11]) if len(sys.argv) > 11 else 0.02

print(f"N={N} T={T} R={R} bochner_n={bn} T_max={Tmx} T_loc={Tlc}")
print(f"h={h} p={p} q in [{qm},{qp}]   mode={mode}")

def run(label, add_las=None):
    Omega, w, v, c, d, eps, dlt, cons = wfc.build_problem(
        N=N, T=T, R=R, h1=h, h2=h, p1=p, p2=p, q1=qm, q2=qp,
        cell_mode="exact", bochner_n=bn, lasserre_T_max=0,
    )
    if add_las is not None:
        add_las(cons, c, d, T_max=Tmx, T_loc=Tlc)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver="CLARABEL")
    print(f"  {label:25s}  Omega* = {prob.value:.10f}   ({prob.status}, {time.time()-t0:.1f}s)")
    return prob.value

if mode in ("all", "b"):
    val_b = run("bochner-only", None)
if mode in ("all", "e"):
    val_e = run("bochner+existingLasserre", las_e.add_lasserre2_constraint)
if mode in ("all", "i"):
    val_i = run("bochner+independentLasserre", las_i.add_lasserre2_constraint_indep)
