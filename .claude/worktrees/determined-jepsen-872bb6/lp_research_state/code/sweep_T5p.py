"""
Sweep over White's Table 3 points with and without T5'.
Persist results to disk so we can resume.
"""
import warnings, time, json, sys, os
warnings.filterwarnings("ignore")
from white_full_convex import solve_full_program

OUT = "T5p_sweep.json"

POINTS = [
    (0.015, 0.381,  -0.02, 0.02, "row1"),
    (0.015, 0.385,  -0.02, 0.02, "row2"),
    (0.020, 0.375,  -0.02, 0.02, "row3"),
    (0.004, 0.3875, -0.02, 0.02, "row4"),
    (0.000, 0.4,    -0.02, 0.02, "row5"),
    (0.000, 0.381,  -0.02, 0.02, "row6"),
    (0.030, 0.375,  -0.02, 0.02, "row7"),
]


def load():
    if os.path.exists(OUT):
        return json.load(open(OUT))
    return {}


def save(d):
    json.dump(d, open(OUT, "w"), indent=2)


def run_one(N, T, R, point, label, use_T5p, db):
    h, p, qm, qp, row = point
    key = f"N{N}_T{T}_R{R}_{row}_{label}"
    if key in db:
        return
    t0 = time.time()
    res = solve_full_program(N, T, R, h, h, p, p, qm, qp, use_T5p=use_T5p)
    db[key] = {
        "N": N, "T": T, "R": R, "row": row,
        "h": h, "p": p, "q_range": [qm, qp],
        "use_T5p": use_T5p,
        "value": float(res["value"]) if res["value"] is not None else None,
        "status": res["status"],
        "time": time.time() - t0,
    }
    save(db)
    return db[key]


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 800
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    rows = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    db = load()
    print(f"Sweep: N={N} T={T} R={R}")
    for pt in POINTS:
        if rows is not None and pt[4] not in rows:
            continue
        for use_T5p, label in [(False, "base"), (True, "T5p")]:
            r = run_one(N, T, R, pt, label, use_T5p, db)
            if r is None:
                continue
            print(f"  {pt[4]} {label:4s}: Ω*={r['value']:.7f}  ({r['status']}, {r['time']:.1f}s)")
    # Print summary across rows
    print()
    print(f"--- Summary at N={N} T={T} R={R} ---")
    for pt in POINTS:
        kb = f"N{N}_T{T}_R{R}_{pt[4]}_base"
        kt = f"N{N}_T{T}_R{R}_{pt[4]}_T5p"
        if kb in db and kt in db:
            b = db[kb]["value"]; t = db[kt]["value"]
            print(f"  {pt[4]}: base={b:.7f}  +T5p={t:.7f}  Δ={t-b:+.7f}")
    bs = [db[f"N{N}_T{T}_R{R}_{pt[4]}_base"]["value"] for pt in POINTS if f"N{N}_T{T}_R{R}_{pt[4]}_base" in db]
    ts = [db[f"N{N}_T{T}_R{R}_{pt[4]}_T5p"]["value"] for pt in POINTS if f"N{N}_T{T}_R{R}_{pt[4]}_T5p" in db]
    if bs and ts:
        print(f"  MIN over rows:  base={min(bs):.7f}  +T5p={min(ts):.7f}  Δ={min(ts)-min(bs):+.7f}")
