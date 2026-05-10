"""Bypass: run Bochner_n=50 LP on row4 at N=3000, persisting result immediately.

Workaround for the bash-kill issue documented in queue _notes.bash_kill_workaround.
Mirrors run_baselines.py pattern but for a single Bochner LP.
"""
import json, sys, time
from pathlib import Path

# Detect mount path
import os
HERE = Path(__file__).resolve().parent
CODE_DIR = HERE / "code"
STATE_DIR = HERE
sys.path.insert(0, str(CODE_DIR))
import white_full_convex as w

POINTS_BY_ROW = {
    "row1": (0.015, 0.381,  -0.02, 0.02),
    "row4": (0.004, 0.3875, -0.02, 0.02),
    "row6": (0.000, 0.381,  -0.02, 0.02),
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

def main():
    N, T, R, row, bochner_n = 3000, 1200, 10, "row4", 50
    done = json.loads(DONE.read_text())
    log(f"run_row4_n50.py starting; existing total = {len(done['results'])}")
    if already_have(done, N, row, bochner_n):
        log(f"  skip {row} bochner_n={bochner_n} (already have)")
        return
    h, p, qm, qp = POINTS_BY_ROW[row]
    t0 = time.time()
    res = w.solve_full_program(N, T, R, h, h, p, p, qm, qp,
                                use_T5p=False, bochner_n=bochner_n)
    rec = {
        "kind": "lp_run_bochner",
        "N": N, "T": T, "R": R,
        "row": row, "h": h, "p": p, "q_range": [qm, qp],
        "use_T5p": False,
        "bochner_n": bochner_n,
        "value": float(res["value"]) if res["value"] is not None else None,
        "status": res["status"],
        "time": time.time() - t0,
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    done["results"].append(rec)
    DONE.write_text(json.dumps(done, indent=2))
    log(f"  {row} bochner_n={bochner_n}: Ω*={rec['value']:.7f} ({rec['status']}, {rec['time']:.1f}s) [persisted; total={len(done['results'])}]")

if __name__ == "__main__":
    main()
