import sys, time, json; sys.path.insert(0,'.')
import numpy as np, cvxpy as cp
from path_b_analytical import build_problem_with_dual_handles

def feas(h,p,q, N=4000, T=1600, R=10, bn=20):
    Omega, cons, H = build_problem_with_dual_handles(N,T,R, h,h, p,p, q,q, bochner_n=bn)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0=time.time()
    try:
        prob.solve(solver='CLARABEL', verbose=False)
    except Exception as e:
        return 'error:'+str(e)[:40], float('nan'), time.time()-t0
    v = prob.value if prob.value is not None else float('nan')
    return prob.status, float(v), time.time()-t0

# Feasibility map. Box: h[0,0.08], p[0,1], q[0.05,1].
# Probe q-range at fixed p (near binding strip p~0.37, h=0) and p near 0,
# plus the worst point. Focus: where does q become infeasible?
pts = []
# (A) q-scan at p=0.37, h=0  (binding strip; find feasible/infeasible boundary)
for q in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.70, 1.0]:
    pts.append((0.0, 0.37, q, 'q-scan p0.37 h0'))
# (B) q-scan at p=0.0, h=0
for q in [0.05, 0.20, 0.40, 0.60, 1.0]:
    pts.append((0.0, 0.0, q, 'q-scan p0 h0'))
# (C) p-scan at q=0.05, h=0 (the low-q feasible edge -- the binding strip)
for p in [0.0, 0.15, 0.25, 0.30, 0.33, 0.42, 0.50, 0.70, 1.0]:
    pts.append((0.0, p, 0.05, 'p-scan q0.05 h0'))
# (D) check h=0.08 top edge near binding
for q in [0.05, 0.20, 0.40]:
    pts.append((0.08, 0.37, q, 'q-scan p0.37 h0.08'))

out = []
for (h,p,q,tag) in pts:
    st,val,dt = feas(h,p,q)
    rec = {'h':h,'p':p,'q':q,'tag':tag,'status':st,'val':val,'t':round(dt,1)}
    out.append(rec)
    feasible = st in ('optimal','optimal_inaccurate')
    print(f"  ({h:.2f},{p:.2f},{q:.2f}) {tag:22s} -> {st:20s} val={val:.5f} {'FEAS' if feasible else 'INFEAS'} ({dt:.0f}s)", flush=True)

json.dump(out, open('_verify_R8_feasmap.json','w'), indent=1)
print('DONE feasmap, wrote _verify_R8_feasmap.json', flush=True)
