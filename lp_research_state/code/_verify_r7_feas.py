import sys; sys.path.insert(0,'.')
import time, json
import cvxpy as cp
from path_b_analytical import build_problem_with_dual_handles

PTS = [
    # (h, p, q) label
    (0.014, 0.3862, -0.025, "BINDING_covermin"),   # the cover grid-min point
    (0.003, 0.3875, -0.05,  "candidate_worst"),     # candidate's reported worst pt
    (0.0,   0.39,   -0.025, "h0_p039_qm025"),
    (0.0,   0.0,    -0.025, "h0_p0_qm025"),         # corner low p
    (0.0,   1.0,    -0.025, "h0_p1_qm025"),         # corner high p
    (0.08,  0.0,    -0.05,  "h08_p0_qm05"),          # corner
    (0.08,  1.0,    -0.05,  "h08_p1_qm05"),          # corner
    (0.04,  0.5,    -0.0375,"mid_box"),
    (0.0,   0.3875, -0.025, "h0_p03875_qm025"),
    (0.0,   0.3875, -0.05,  "h0_p03875_qm05"),
]
N,T,R,bn = 4000,1600,10,20
out=[]
for h,p,q,lab in PTS:
    t0=time.time()
    try:
        Omega,cons,handles = build_problem_with_dual_handles(N,T,R,h,h,p,p,q,q,bochner_n=bn)
        prob = cp.Problem(cp.Minimize(Omega), cons)
        prob.solve(solver='CLARABEL')
        st=prob.status; val=prob.value
    except Exception as e:
        st=f"ERR:{type(e).__name__}:{str(e)[:60]}"; val=None
    dt=time.time()-t0
    rec={"h":h,"p":p,"q":q,"label":lab,"status":st,"value":(None if val is None else float(val)),"sec":round(dt,1)}
    out.append(rec)
    print(json.dumps(rec), flush=True)
json.dump(out, open("../parallel_results/_r7_feas_map.json","w"), indent=2)
print("DONE")
