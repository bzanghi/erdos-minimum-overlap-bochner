"""Bypass: run baseline LPs (no Bochner, no T5p) for the 7 rows at N=2000,
saving each one to experiments_done.json IMMEDIATELY after solving.

This avoids the issue where cron_runner saves only after all LPs finish, so a
bash-induced kill loses all progress.
"""
import json, os, sys, time
from pathlib import Path

CODE_DIR = Path("/sessions/trusting-serene-gauss/mnt/Erdos/lp_research_state/code")
STATE_DIR = Path("/sessions/trusting-serene-gauss/mnt/Erdos/lp_research_state")
sys.path.insert(0, str(CODE_DIR))
import white_full_convex as w

POINTS_BY_ROW = {
    "row1": (0.015, 0.381,  -0.02, 0.02),
    "row2": (0.015, 0.385,  -0.02, 0.02),
    "row3": (0.020, 0.375,  -0.02, 0.02),
    "row4": (0.004, 0.3875, -0.02, 0.02),
    "row5": (0.000, 0.4,    -0.02, 0.02),
    "row6": (0.000, 0.381,  -0.02, 0.02),
    "row7": (0.030, 0.375,  -0.02, 0.02),
}

DONE = STATE_DIR / "experiments_done.json"
LOG = STATE_DIR / "cron_log.txt"

def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def already_have(done_results, N, row, use_T5p, bochner_n=0):
    for r in done_results:
        if (r.get("kind") == "lp_run"
            and r.get("N") == N
            and r.get("row") == row
            and r.get("use_T5p") == use_T5p
            and r.get("bochner_n", 0) == bochner_n):
            return True
    return False

def main():
    N, T, R = 2000, 800, 10
    done = json.loads(DONE.read_text())
    log(f"run_baselines.py starting; existing total = {len(done['results'])}")
    new_count = 0
    target = ["row1", "row2", "row3", "row4", "row6", "row7"]  # row5 already exists
    # We only need use_T5p=False for the baseline MIN computation.
    for row in target:
        if already_have(done["results"], N, row, False, 0):
            log(f"  skip {row} (already have)")
            continue
        h, p, qm, qp = POINTS_BY_ROW[row]
        t0 = time.time()
        res = w.solve_full_program(N, T, R, h, h, p, p, qm, qp, use_T5p=False)
        rec = {
            "kind": "lp_run",
            "N": N, "T": T, "R": R,
            "row": row, "h": h, "p": p, "q_range": [qm, qp],
            "use_T5p": False,
            "value": float(res["value"]) if res["value"] is not None else None,
            "status": res["status"],
            "time": time.time() - t0,
            "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        done["results"].append(rec)
        DONE.write_text(json.dumps(done, indent=2))
        new_count += 1
        log(f"  {row}: Ω*={rec['value']:.7f} ({rec['status']}, {rec['time']:.1f}s) [persisted; total={len(done['results'])}]")
    log(f"run_baselines.py done; added {new_count} new results")

if __name__ == "__main__":
    main()
