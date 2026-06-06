"""Run row 4 (h=0.004, p=0.3875, q∈[-0.02,0.02]) at N=2000, T=800, R=10 with
four configs: baseline, +Bochner_n=20, +Bochner_n=20 +M-Schur_n=10, +M-Schur_n=20.
Writes JSON to /Erdos/lp_research_state/parallel_results/mside_schur_results.json.
"""
import warnings; warnings.filterwarnings("ignore")
import sys, json, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from white_full_convex import solve_full_program

OUT = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/mside_schur_results.json"

PARAMS = dict(
    N=2000, T=800, R=10,
    h1=0.004, h2=0.004,
    p1=0.3875, p2=0.3875,
    q1=-0.02, q2=0.02,
    cell_mode="exact", solver="CLARABEL", verbose=False,
)

CONFIGS = [
    ("baseline",                   dict()),
    ("bochner_n20",                dict(bochner_n=20)),
    ("bochner_n20_mschur_n10",     dict(bochner_n=20, mside_bochner_schur_n=10)),
    ("bochner_n20_mschur_n20",     dict(bochner_n=20, mside_bochner_schur_n=20)),
]

def run_one(label, cfg):
    print(f"\n=== {label}  {cfg} ===", flush=True)
    t0 = time.time()
    res = solve_full_program(**PARAMS, **cfg)
    return {
        "label": label, "config": cfg,
        "value": float(res["value"]) if res["value"] is not None else None,
        "status": res["status"],
        "time_sec": round(time.time() - t0, 2),
    }

if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else None
    out = {}
    if os.path.exists(OUT):
        try:
            with open(OUT) as fh:
                out = json.load(fh)
        except Exception:
            out = {}
    for label, cfg in CONFIGS:
        if only and label != only:
            continue
        rec = run_one(label, cfg)
        out[label] = rec
        print(f"  -> value={rec['value']} status={rec['status']} t={rec['time_sec']}s", flush=True)
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        with open(OUT, "w") as fh:
            json.dump(out, fh, indent=2)
    print("\nFinal results:")
    print(json.dumps(out, indent=2))
