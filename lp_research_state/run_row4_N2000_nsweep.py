"""Bochner n-sweep on row4 at N=2000: try n in [25, 30] in order.

Persist each result to experiments_done.json IMMEDIATELY after solving so
that a bash-kill mid-batch retains the partial progress (mirrors
run_baselines.py / run_row4_n50.py pattern).

Logic:
  - n=25 always attempted.
  - n=30 only attempted if n=25 finished in <30s (so we don't blow the
    bash budget on something the prior data point already showed is slow).
"""
import json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE_DIR = HERE / "code"
STATE_DIR = HERE
sys.path.insert(0, str(CODE_DIR))
import white_full_convex as w

POINTS_BY_ROW = {
    "row4": (0.004, 0.3875, -0.02, 0.02),
}

DONE = STATE_DIR / "experiments_done.json"
LOG = STATE_DIR / "cron_log.txt"

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def already_have(done, N, row, bochner_n):
    for r in done["results"]:
        if (r.get("kind") == "lp_run_bochner"
            and r.get("N") == N
            and r.get("row") == row
            and r.get("bochner_n") == bochner_n
            and r.get("use_T5p") == False):
            return True
    return False

def run_one(N, T, R, row, bochner_n):
    done = json.loads(DONE.read_text())
    if already_have(done, N, row, bochner_n):
        log(f"  skip {row} N={N} bochner_n={bochner_n} (already have)")
        # find existing time for the gating decision
        for r in done["results"]:
            if (r.get("kind") == "lp_run_bochner" and r.get("N") == N
                and r.get("row") == row and r.get("bochner_n") == bochner_n):
                return r.get("time", 0.0)
        return 0.0
    h, p, qm, qp = POINTS_BY_ROW[row]
    t0 = time.time()
    res = w.solve_full_program(N, T, R, h, h, p, p, qm, qp,
                                use_T5p=False, bochner_n=bochner_n)
    elapsed = time.time() - t0
    rec = {
        "kind": "lp_run_bochner",
        "N": N, "T": T, "R": R,
        "row": row, "h": h, "p": p, "q_range": [qm, qp],
        "use_T5p": False,
        "bochner_n": bochner_n,
        "value": float(res["value"]) if res["value"] is not None else None,
        "status": res["status"],
        "time": elapsed,
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    # re-load to avoid clobbering concurrent writes
    done = json.loads(DONE.read_text())
    done["results"].append(rec)
    DONE.write_text(json.dumps(done, indent=2))
    val_str = f"{rec['value']:.7f}" if rec['value'] is not None else "None"
    log(f"  {row} N={N} bochner_n={bochner_n}: Ω*={val_str} ({rec['status']}, {elapsed:.1f}s) [persisted; total={len(done['results'])}]")
    return elapsed

def main():
    N, T, R, row = 2000, 800, 10, "row4"
    log(f"run_row4_N2000_nsweep.py starting (N={N}, T={T}, R={R}, row={row})")
    # n=25 first
    t25 = run_one(N, T, R, row, 25)
    # gate n=30 on n=25 time
    if t25 < 30.0:
        log(f"  n=25 finished in {t25:.1f}s (<30s); attempting n=30")
        run_one(N, T, R, row, 30)
    else:
        log(f"  n=25 took {t25:.1f}s (>=30s); skipping n=30 to stay in budget")

if __name__ == "__main__":
    main()
