"""
Compare independent Lasserre-2 implementation vs the existing lasserre.py.
Solve the row-4 SDP (h=0.004, p=0.3875, q in [-0.02,0.02]) at the requested
parameters, both ways, and print Omega* to 10 digits.
"""
from __future__ import annotations
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import cvxpy as cp

# Import the *base* program builder (no Lasserre attached).
import importlib.util as _ilu
from pathlib import Path
HERE = Path("/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code")

def _load(name, fname):
    spec = _ilu.spec_from_file_location(name, HERE / fname)
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

wfc = _load("wfc", "white_full_convex.py")
las_existing = _load("las_existing", "lasserre.py")
las_indep = _load("las_indep", "lasserre_independent.py")


def solve_with(mode: str, T_max: int, T_loc: int,
               N=400, T=80, R=10, bochner_n=15,
               h=0.004, p=0.3875, qm=-0.02, qp=0.02):
    """mode in {'bochner_only', 'existing', 'independent'}"""
    Omega, w, v, c, d, eps, dlt, cons = wfc.build_problem(
        N=N, T=T, R=R, h1=h, h2=h, p1=p, p2=p, q1=qm, q2=qp,
        cell_mode="exact", bochner_n=bochner_n, lasserre_T_max=0,
    )
    if mode == "existing":
        las_existing.add_lasserre2_constraint(cons, c, d,
                                              T_max=T_max, T_loc=T_loc)
    elif mode == "independent":
        las_indep.add_lasserre2_constraint_indep(cons, c, d,
                                                 T_max=T_max, T_loc=T_loc)
    elif mode == "bochner_only":
        pass
    else:
        raise ValueError(mode)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver="CLARABEL")
    return prob.value, prob.status, time.time() - t0


if __name__ == "__main__":
    print("Configuration: row 4 (h=0.004, p=0.3875, q in [-0.02,0.02])")
    print("              N=400, T=80, R=10, bochner_n=15, T_max=10  (smaller for time)")
    print()

    Tmx = 6
    Tlc = 6
    # 1. Bochner-only (sanity)
    val_b, st_b, tt_b = solve_with("bochner_only", T_max=Tmx, T_loc=Tlc)
    print(f"  Bochner-only:           Omega* = {val_b:.10f}   ({st_b}, {tt_b:.1f}s)")

    # 2. Existing lasserre.py
    val_e, st_e, tt_e = solve_with("existing", T_max=Tmx, T_loc=Tlc)
    print(f"  Bochner + existing-Las: Omega* = {val_e:.10f}   ({st_e}, {tt_e:.1f}s)")

    # 3. Independent re-implementation
    val_i, st_i, tt_i = solve_with("independent", T_max=Tmx, T_loc=Tlc)
    print(f"  Bochner + indep-Las   : Omega* = {val_i:.10f}   ({st_i}, {tt_i:.1f}s)")

    print()
    print(f"  Lasserre contribution (existing) : {val_e - val_b:+.6e}")
    print(f"  Lasserre contribution (indep)    : {val_i - val_b:+.6e}")
    print(f"  Disagreement existing vs indep   : {abs(val_e - val_i):.3e}")
