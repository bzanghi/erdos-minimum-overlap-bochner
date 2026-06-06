"""Check both runs at the requested configuration with verbose output to
inspect actual Clarabel convergence."""
import warnings; warnings.filterwarnings("ignore")
import time, sys
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

mode = sys.argv[1]
N=2000; T=200; R=10; bn=20; Tmx=10; Tlc=10
h=0.004; p=0.3875; qm=-0.02; qp=0.02

Omega, w, v, c, d, eps, dlt, cons = wfc.build_problem(
    N=N, T=T, R=R, h1=h, h2=h, p1=p, p2=p, q1=qm, q2=qp,
    cell_mode="exact", bochner_n=bn, lasserre_T_max=0,
)
if mode == 'i':
    las_i.add_lasserre2_constraint_indep(cons, c, d, T_max=Tmx, T_loc=Tlc)
elif mode == 'e':
    las_e.add_lasserre2_constraint(cons, c, d, T_max=Tmx, T_loc=Tlc)

prob = cp.Problem(cp.Minimize(Omega), cons)
t0 = time.time()
prob.solve(solver="CLARABEL", verbose=True,
           eps_abs=1e-9, eps_rel=1e-9,
           max_iter=300)
print(f"\n>>> mode={mode}: Omega* = {prob.value:.10f}  status = {prob.status}  time={time.time()-t0:.1f}s")
