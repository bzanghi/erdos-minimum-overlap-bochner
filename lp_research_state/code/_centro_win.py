"""
PRO-49 centrosymmetric split — WIN MEASUREMENT (real form vs centro).

For each bochner_n, a child process builds row4 at fixed (N,T,R), records the
PSD cone CLARABEL receives, solves with CLARABEL, and reports wall time + peak
RSS (resource.getrusage; macOS ru_maxrss is BYTES). We then report:
  - time_factor  = real_wall / centro_wall   (>1 => centro faster)
  - mem_factor   = real_peak / centro_peak   (>1 => centro lighter)
  - max bochner_n each form reaches before OOM/timeout at fixed memory
  - row4 rigorous dual LB at the highest feasible bn (via dual_extractor)

Moderate scale (N<=5000) to coexist with the co-running production L2 solve.
"""
from __future__ import annotations
import sys, os, json, time, subprocess, textwrap

CODE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
OUT = os.path.join(CODE, "..", "..", "docs", "NEW_APPROACHES", "sym_centro_result.json")
OUT = os.path.normpath(OUT)

# row4 = binding row
H, P, QM, QP = 0.004, 0.3875, -0.02, 0.02

CHILD = textwrap.dedent(r"""
    import sys, os, json, time, resource
    sys.path.insert(0, "__CODE__")
    import warnings; warnings.filterwarnings("ignore")
    import cvxpy as cp
    form, N, T, R, bn, want_lb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]=="1"
    H,P,QM,QP = __H__, __P__, __QM__, __QP__
    if form == "real":
        from white_full_convex import build_problem as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    else:
        from white_full_convex_centro import build_problem_centro as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    prob = cp.Problem(cp.Minimize(O), cons)
    try:
        data,_,_ = prob.get_problem_data(cp.CLARABEL)
        psd = sorted(list(data["dims"].psd)); Ashape = list(data["A"].shape); Annz = int(data["A"].nnz)
    except Exception as e:
        psd = ["ERR:%s"%e]; Ashape=None; Annz=None
    t0 = time.time()
    rig_lb = None
    try:
        if want_lb:
            from dual_extractor import solve_with_dual_extraction
            res = solve_with_dual_extraction(prob)
            status = res["status"]; value = res["reported_value"]; rig_lb = res["rigorous_dual_LB"]
        else:
            prob.solve(solver="CLARABEL")
            status = prob.status; value = (float(prob.value) if prob.value is not None else None)
    except MemoryError:
        status = "OOM"; value=None
    except Exception as e:
        status = "EXC:%s"%type(e).__name__; value=None
    wall = time.time()-t0
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mb = peak_kb/(1024.0*1024.0) if sys.platform=="darwin" else peak_kb/1024.0
    print(json.dumps({"form":form,"N":N,"T":T,"R":R,"bn":bn,"psd_blocks":psd,
                      "A_shape":Ashape,"A_nnz":Annz,"wall":wall,"peak_mb":peak_mb,
                      "status":status,"value":value,"rigorous_dual_LB":rig_lb}))
""")
CHILD = (CHILD.replace("__CODE__", CODE).replace("__H__", repr(H))
              .replace("__P__", repr(P)).replace("__QM__", repr(QM))
              .replace("__QP__", repr(QP)))


def run_child(form, N, T, R, bn, want_lb=False, timeout=1200):
    try:
        out = subprocess.run([PY, "-c", CHILD, form, str(N), str(T), str(R), str(bn),
                              "1" if want_lb else "0"],
                             capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"form": form, "N": N, "T": T, "R": R, "bn": bn, "status": "TIMEOUT",
                "wall": timeout, "peak_mb": None, "value": None, "psd_blocks": None}
    line = None
    for ln in out.stdout.splitlines():
        ln = ln.strip()
        if ln.startswith("{"):
            line = ln
    if line is None:
        return {"form": form, "N": N, "T": T, "R": R, "bn": bn, "status": "NOJSON",
                "stderr_tail": out.stderr[-500:], "wall": None, "peak_mb": None, "value": None}
    return json.loads(line)


def load():
    if os.path.exists(OUT):
        try:
            with open(OUT) as f: return json.load(f)
        except Exception: return {}
    return {}


def save(d):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f: json.dump(d, f, indent=2)


if __name__ == "__main__":
    N, T, R = 5000, 2000, 10
    bns = [10, 20, 30, 40, 50, 60]
    want_lb_at = None  # bn at which to also extract rigorous LB; default = max successful
    if len(sys.argv) > 1: bns = [int(x) for x in sys.argv[1].split(",")]
    if len(sys.argv) > 2: N, T, R = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    if len(sys.argv) > 5: want_lb_at = int(sys.argv[5])

    data = load()
    data.setdefault("win_measurements", [])
    data["win_config"] = {"N": N, "T": T, "R": R, "row4": [H, P, QM, QP], "scale_note": "moderate (N<=5000) to coexist with co-running L2 solve"}

    for bn in bns:
        pair = {"bochner_n": bn}
        do_lb = (want_lb_at == bn)
        for form in ("real", "centro"):
            print(f"[{time.strftime('%H:%M:%S')}] form={form} bn={bn} (N={N},T={T},R={R}) lb={do_lb} ...", flush=True)
            rec = run_child(form, N, T, R, bn, want_lb=do_lb)
            pair[form] = rec
            print(f"    -> status={rec.get('status')} wall={rec.get('wall')!r} "
                  f"peak_mb={rec.get('peak_mb')!r} psd={rec.get('psd_blocks')} "
                  f"value={rec.get('value')!r} rigLB={rec.get('rigorous_dual_LB')!r}", flush=True)
        try:
            wr, wc = pair["real"].get("wall"), pair["centro"].get("wall")
            mr, mc = pair["real"].get("peak_mb"), pair["centro"].get("peak_mb")
            pair["time_factor_real_over_centro"] = (wr / wc) if (wr and wc) else None
            pair["mem_factor_real_over_centro"] = (mr / mc) if (mr and mc) else None
            print(f"    => time_factor(real/centro)={pair['time_factor_real_over_centro']}, "
                  f"mem_factor(real/centro)={pair['mem_factor_real_over_centro']}", flush=True)
        except Exception:
            pass
        data["win_measurements"].append(pair)
        save(data)
    print(f"\nSaved -> {OUT}")
