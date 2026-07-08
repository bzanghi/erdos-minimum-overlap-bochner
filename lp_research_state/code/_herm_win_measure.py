"""
Approach ③, step 3 — measure the cost win (or its absence).

Solve binding row4 at increasing bochner_n with BOTH the real-embedding and the
complex-Hermitian forms. Record:
  * canonicalized PSD block side-lengths CLARABEL actually sees (the decisive
    number — if these are equal, there is NO solver-level win);
  * wall-time (compile + solve);
  * peak RSS (via resource.getrusage of a child process — measured per-solve);
  * max bochner_n each form reaches before OOM/slowdown.

Writes to docs/NEW_APPROACHES/sym_reduction_result.json.
"""
from __future__ import annotations
import sys, os, json, time, subprocess, textwrap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
JSON_PATH = os.path.join(REPO, "docs", "NEW_APPROACHES", "sym_reduction_result.json")
PY = sys.executable
CODE = os.path.dirname(os.path.abspath(__file__))

H, P, QM, QP = 0.004, 0.3875, -0.02, 0.02


def load_json():
    if os.path.exists(JSON_PATH):
        with open(JSON_PATH) as f:
            return json.load(f)
    return {"win_measurements": []}


def save_json(data):
    os.makedirs(os.path.dirname(JSON_PATH), exist_ok=True)
    with open(JSON_PATH, "w") as f:
        json.dump(data, f, indent=2)


# Child script: solve ONE (form, N, T, R, bn) and print a JSON line with
# peak RSS (resource), wall time, status, value, and PSD block dims.
# (Use __PLACEHOLDER__ markers instead of str.format to avoid JSON-brace clashes.)
CHILD = textwrap.dedent(r"""
    import sys, os, json, time, resource
    sys.path.insert(0, "__CODE__")
    import warnings; warnings.filterwarnings("ignore")
    import cvxpy as cp
    form, N, T, R, bn = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    H,P,QM,QP = __H__, __P__, __QM__, __QP__
    if form == "real":
        from white_full_convex import build_problem as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    else:
        from white_full_convex_hermitian import build_problem_hermitian as B
        O,*_,cons = B(N,T,R,H,H,P,P,QM,QP,bochner_n=bn)
    prob = cp.Problem(cp.Minimize(O), cons)
    # canonical PSD dims
    try:
        data,_,_ = prob.get_problem_data(cp.CLARABEL)
        psd = list(data["dims"].psd); Ashape = list(data["A"].shape); Annz = int(data["A"].nnz)
    except Exception as e:
        psd = ["ERR:%s"%e]; Ashape=None; Annz=None
    t0 = time.time()
    try:
        prob.solve(solver="CLARABEL")
        status = prob.status; value = (float(prob.value) if prob.value is not None else None)
    except Exception as e:
        status = "EXC:%s"%type(e).__name__; value=None
    wall = time.time()-t0
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mb = peak_kb/1024.0 if sys.platform=="darwin" else peak_kb/1024.0  # darwin: bytes? handled below
    # macOS ru_maxrss is in BYTES; linux in KB. Detect:
    if sys.platform == "darwin":
        peak_mb = peak_kb/(1024.0*1024.0)
    print(json.dumps({"form":form,"N":N,"T":T,"R":R,"bn":bn,"psd_blocks":psd,
                      "A_shape":Ashape,"A_nnz":Annz,"wall":wall,"peak_mb":peak_mb,
                      "status":status,"value":value}))
""")
CHILD = (CHILD.replace("__CODE__", CODE).replace("__H__", repr(H))
              .replace("__P__", repr(P)).replace("__QM__", repr(QM))
              .replace("__QP__", repr(QP)))


def run_child(form, N, T, R, bn, timeout=900):
    try:
        out = subprocess.run([PY, "-c", CHILD, form, str(N), str(T), str(R), str(bn)],
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
                "stderr_tail": out.stderr[-400:], "wall": None, "peak_mb": None, "value": None}
    return json.loads(line)


if __name__ == "__main__":
    # Fixed row4 program size; sweep bochner_n upward for both forms.
    N, T, R = 5000, 2000, 10
    bns = [10, 15, 20, 25, 30, 40]
    if len(sys.argv) > 1:
        bns = [int(x) for x in sys.argv[1].split(",")]
    if len(sys.argv) > 2:
        N, T, R = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

    data = load_json()
    data.setdefault("win_measurements", [])
    data["win_config"] = {"N": N, "T": T, "R": R, "row4": [H, P, QM, QP]}

    for bn in bns:
        pair = {"bochner_n": bn}
        for form in ("real", "hermitian"):
            print(f"[{time.strftime('%H:%M:%S')}] solving form={form} bn={bn} (N={N},T={T},R={R}) ...", flush=True)
            rec = run_child(form, N, T, R, bn)
            pair[form] = rec
            print(f"    -> status={rec.get('status')} wall={rec.get('wall')!r} "
                  f"peak_mb={rec.get('peak_mb')!r} psd_blocks={rec.get('psd_blocks')} "
                  f"value={rec.get('value')!r}", flush=True)
        # factors
        try:
            wr, wh = pair["real"]["wall"], pair["hermitian"]["wall"]
            mr, mh = pair["real"]["peak_mb"], pair["hermitian"]["peak_mb"]
            pair["time_factor_real_over_herm"] = (wr / wh) if (wr and wh) else None
            pair["mem_factor_real_over_herm"] = (mr / mh) if (mr and mh) else None
        except Exception:
            pass
        data["win_measurements"].append(pair)
        save_json(data)
        print(f"    time_factor(real/herm)={pair.get('time_factor_real_over_herm')!r} "
              f"mem_factor={pair.get('mem_factor_real_over_herm')!r}", flush=True)
