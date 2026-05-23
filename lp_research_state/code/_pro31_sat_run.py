"""PRO-31 extension driver: certified M(n) via SAT, pushing past the ILP wall.

Validates against Haugland (n=10,12,15,16,17,18) then extends n=19,20,21,...
upward until a per-n time wall. Persists incrementally (merges, no clobber) to
lp_research_state/data/discrete_M_values.json under key "extension_sat".
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from discrete_M_sat import solve_Mn, max_overlap  # noqa: E402

DATA = CODE.parent / "data" / "discrete_M_values.json"
TOGETHER_UB = 0.380871

KNOWN = {10: 5, 12: 5, 13: 6, 14: 6, 15: 6, 16: 7, 17: 7, 18: 8}

PER_N_BUDGET = 1800.0          # cap per n (seconds)
TOTAL_BUDGET = float(sys.argv[1]) if len(sys.argv) > 1 else 9000.0


def load():
    if DATA.exists():
        return json.loads(DATA.read_text())
    return {"validation": [], "extension": []}


def save(doc):
    DATA.write_text(json.dumps(doc, indent=2))


def main():
    doc = load()
    doc.setdefault("extension_sat", [])
    done_ns = {r["n"] for r in doc["extension_sat"] if r.get("optimal")}

    # Validation phase (fast; confirms encoding incl. n=18 the ILP missed)
    print("=== VALIDATION (SAT) ===", flush=True)
    val = []
    for n in [10, 12, 15, 16, 17, 18]:
        r = solve_Mn(n, time_budget_sec=120, verbose=False)
        exp = KNOWN[n]
        ok = r["M"] == exp and r["optimal"]
        print(f"  n={n}: M={r['M']} exp={exp} opt={r['optimal']} "
              f"t={r['wall_time_sec']:.2f}s {'OK' if ok else 'FAIL'}", flush=True)
        if not ok:
            print("  ABORT: validation mismatch — SAT encoding bug.", flush=True)
            return
        val.append({"n": n, "M": r["M"], "ratio": r["ratio"],
                    "optimal": r["optimal"], "time": r["wall_time_sec"]})
    doc["validation_sat"] = val
    save(doc)

    # Extension phase. M is non-decreasing-ish; bracket tightly using last M.
    print("\n=== EXTENSION (SAT) ===", flush=True)
    t_ext = time.time()
    prev_M = 8  # M(18)
    n = 19
    while True:
        if n in done_ns:
            # use recorded M as the new prev
            rec = next(r for r in doc["extension_sat"] if r["n"] == n)
            prev_M = rec["M"]
            n += 1
            continue
        elapsed_total = time.time() - t_ext
        if elapsed_total > TOTAL_BUDGET:
            print(f"Total budget {TOTAL_BUDGET}s exhausted at n={n}.", flush=True)
            break
        budget = min(PER_N_BUDGET, TOTAL_BUDGET - elapsed_total)
        # bracket: M(n) ≥ prev_M is NOT guaranteed (ratios fluctuate), but
        # M jumps by at most ~1 between consecutive n in practice. Be safe:
        lo = max(0, prev_M - 1)
        hi = prev_M + 2
        print(f"\n--- n={n} (bracket [{lo},{hi}], budget {budget:.0f}s) ---",
              flush=True)
        r = solve_Mn(n, M_lo=lo, M_hi=hi, time_budget_sec=budget, verbose=True)

        beats = (r["ratio"] is not None and r["ratio"] < TOGETHER_UB)
        rec = {
            "n": n, "N": 2 * n, "M": r["M"], "ratio": r["ratio"],
            "optimal": r["optimal"],
            "sat_cert_at_M": r["sat_certificate_at_M"],
            "unsat_cert_at_M_minus_1": r["unsat_certificate_at_M_minus_1"],
            "verified_overlap": r["verified_overlap"],
            "A_star": r["A_star"],
            "wall_time_sec": r["wall_time_sec"],
            "beats_together_ub": beats,
            "solver": r["solver"],
            "n_sat_calls": len(r["calls"]),
        }
        # replace any prior partial entry for this n
        doc["extension_sat"] = [x for x in doc["extension_sat"] if x["n"] != n]
        doc["extension_sat"].append(rec)
        doc["extension_sat"].sort(key=lambda x: x["n"])
        save(doc)

        if r["M"] is None:
            print(f"  n={n}: no witness within budget — WALL. Stopping.",
                  flush=True)
            break
        print(f"  -> M({n})={r['M']} ratio={r['ratio']:.6f} "
              f"optimal={r['optimal']} t={r['wall_time_sec']:.1f}s "
              f"beats_together={beats}", flush=True)
        if not r["optimal"]:
            print(f"  n={n}: optimality NOT proven within budget "
                  f"(have UB M≤{r['M']}). WALL. Stopping.", flush=True)
            break

        prev_M = r["M"]
        n += 1

    print(f"\nSaved to {DATA}", flush=True)
    # summary
    print("\n=== EXTENSION SUMMARY ===", flush=True)
    print(f"{'n':>4} {'M':>4} {'M/n':>10} {'opt':>5} {'time(s)':>9} {'<UB?':>5}",
          flush=True)
    for r in doc["extension_sat"]:
        print(f"{r['n']:>4} {str(r['M']):>4} "
              f"{(r['ratio'] if r['ratio'] else 0):>10.6f} "
              f"{str(r['optimal']):>5} {r['wall_time_sec']:>9.1f} "
              f"{str(r['beats_together_ub']):>5}", flush=True)


if __name__ == "__main__":
    main()
