"""
PRO-49 — how HIGH does centro reach (centro-only, to avoid OOM-ing the co-running
L2 solve with the heavy real form). Each bn in its own child process; reports peak
RSS, wall, PSD side, row4 value + rigorous dual LB. We compare against the real
form's MEASURED footprint at lower bn (sym_centro_result.json) + the established
~3.2x memory factor to state the equal-budget max-bn claim honestly.
"""
from __future__ import annotations
import sys, os, json, time, subprocess, textwrap

CODE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
OUT = os.path.normpath(os.path.join(CODE, "..", "..", "docs", "NEW_APPROACHES", "sym_centro_result.json"))
H, P, QM, QP = 0.004, 0.3875, -0.02, 0.02

CHILD = textwrap.dedent(r"""
    import sys, os, json, time, resource
    sys.path.insert(0, "__CODE__")
    import warnings; warnings.filterwarnings("ignore")
    import cvxpy as cp
    form, N, T, R, bn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    H,P,QM,QP = __H__, __P__, __QM__, __QP__
    if form=="real":
        from white_full_convex import build_problem as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    else:
        from white_full_convex_centro import build_problem_centro as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    prob = cp.Problem(cp.Minimize(O), cons)
    try:
        data,_,_ = prob.get_problem_data(cp.CLARABEL); psd=sorted(list(data["dims"].psd))
    except Exception as e: psd=["ERR:%s"%e]
    t0=time.time()
    try:
        from dual_extractor import solve_with_dual_extraction
        res = solve_with_dual_extraction(prob)
        status=res["status"]; value=res["reported_value"]; rig=res["rigorous_dual_LB"]
    except MemoryError:
        status="OOM"; value=None; rig=None
    except Exception as e:
        status="EXC:%s"%type(e).__name__; value=None; rig=None
    wall=time.time()-t0
    pk=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    pk_mb=pk/(1024.0*1024.0) if sys.platform=="darwin" else pk/1024.0
    print(json.dumps({"form":form,"bn":bn,"psd_blocks":psd,"wall":wall,"peak_mb":pk_mb,"status":status,"value":value,"rigorous_dual_LB":rig}))
""").replace("__CODE__",CODE).replace("__H__",repr(H)).replace("__P__",repr(P)).replace("__QM__",repr(QM)).replace("__QP__",repr(QP))


def run(form, bn, N, T, R, timeout=1800):
    try:
        out = subprocess.run([PY,"-c",CHILD,form,str(N),str(T),str(R),str(bn)],
                             capture_output=True,text=True,timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"form":form,"bn":bn,"status":"TIMEOUT","wall":timeout,"peak_mb":None,"value":None}
    line=None
    for ln in out.stdout.splitlines():
        if ln.strip().startswith("{"): line=ln.strip()
    if line is None:
        return {"form":form,"bn":bn,"status":"NOJSON","stderr_tail":out.stderr[-400:],"value":None,"peak_mb":None,"wall":None}
    return json.loads(line)


if __name__ == "__main__":
    N,T,R = 5000,2000,10
    bns = [60, 80]
    if len(sys.argv)>1: bns=[int(x) for x in sys.argv[1].split(",")]
    if len(sys.argv)>2: N,T,R=int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])
    form = sys.argv[5] if len(sys.argv)>5 else "centro"

    data = json.load(open(OUT)) if os.path.exists(OUT) else {}
    data.setdefault("highbn_measurements", [])
    for bn in bns:
        print(f"[{time.strftime('%H:%M:%S')}] form={form} bn={bn} (N={N},T={T},R={R}) ...", flush=True)
        rec = run(form, bn, N, T, R)
        print(f"    -> status={rec.get('status')} wall={rec.get('wall')!r} peak_mb={rec.get('peak_mb')!r} "
              f"psd={rec.get('psd_blocks')} value={rec.get('value')!r} rigLB={rec.get('rigorous_dual_LB')!r}", flush=True)
        data["highbn_measurements"].append(rec)
        json.dump(data, open(OUT,"w"), indent=2)
    print(f"Saved -> {OUT}")
