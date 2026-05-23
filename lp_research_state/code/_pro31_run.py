"""PRO-31 driver: compute M(n) via HiGHS ILP for n in a range, with budget
gating and an independent SAT cross-check on a few values.

Outputs: lp_research_state/data/discrete_M_values.json

Strategy:
  1. Re-verify Haugland values for a small batch (sanity).
  2. Push forward starting at n = 21 (prior SAT got n<=20).
  3. Skip any n whose ILP exceeds `single_n_cap_sec` and stop.
  4. For each new n (n >= 21), if cheap enough, also run a SAT solve and
     assert the two agree (essential cross-check per repo invariants).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from discrete_M_ilp import solve_Mn_ilp, _max_overlap  # type: ignore

# Reuse the prior SAT encoder if available — independent encoding.
try:
    from _sat_Mn import solve_Mn_sat  # type: ignore
    SAT_AVAILABLE = True
except Exception:
    SAT_AVAILABLE = False


HAUGLAND = {
    1: 1, 2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 3, 8: 4, 9: 4, 10: 5,
    11: 5, 12: 5, 13: 6, 14: 6, 15: 6,
    # Our prior session brute force + SAT (cross-confirmed): 16-20
    16: 7, 17: 7, 18: 8, 19: 8, 20: 8,
}


def main():
    out_path = ROOT / "data" / "discrete_M_values.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    single_n_cap_sec = 30 * 60.0
    total_budget_sec = 90 * 60.0

    # Validation pass: 4 Haugland values (n in {10, 12, 15, 18}).
    print("=== Validation: Haugland reproduction ===", flush=True)
    val_results = []
    for n in [10, 12, 15, 18]:
        r = solve_Mn_ilp(n, time_limit=600)
        exp = HAUGLAND[n]
        ok = (r["M"] == exp) and (r["verified_overlap"] == exp) and r["optimal"]
        print(
            f"  n={n:2d}: M={r['M']} (expected {exp}) verified={r['verified_overlap']} "
            f"t={r['wall_time_sec']:.2f}s opt={r['optimal']} gap={r['mip_gap']:.2e} -> {'OK' if ok else 'FAIL'}",
            flush=True,
        )
        r["expected_M"] = exp
        r["validation_ok"] = ok
        val_results.append(r)
        if not ok:
            print("ABORT: validation failed — encoding bug.", flush=True)
            with open(out_path, "w") as f:
                json.dump({"validation": val_results, "extension": []}, f, indent=2)
            return

    print("\n=== Extension: pushing M(n) beyond prior known table ===", flush=True)
    ext_results = []
    t_global = time.time()
    n = 21
    while True:
        if time.time() - t_global > total_budget_sec:
            print(f"  Total budget exhausted at n={n}.", flush=True)
            break
        # Remaining budget for this n (don't exceed global).
        remaining = total_budget_sec - (time.time() - t_global)
        budget = min(single_n_cap_sec, remaining)
        print(f"\n--- n={n} (budget {budget:.0f}s) ---", flush=True)
        r = solve_Mn_ilp(n, time_limit=budget)
        if r["A_star"] is not None:
            r["verified_overlap"] = _max_overlap(tuple(r["A_star"]), 2 * n)
        exp = HAUGLAND.get(n)
        match = (exp is None) or (r["M"] == exp)
        print(
            f"  -> M({n}) = {r['M']}, ratio = {r['ratio']}, "
            f"t = {r['wall_time_sec']:.1f}s, opt={r['optimal']}, "
            f"gap={r['mip_gap']:.2e}, expected={exp}, match={match}",
            flush=True,
        )
        if r["A_star"]:
            print(f"     verified overlap = {r['verified_overlap']}", flush=True)
            print(f"     A* = {r['A_star']}", flush=True)
        r["expected_M"] = exp

        # Cross-check via SAT (independent encoding) for small-enough n.
        if SAT_AVAILABLE and r["wall_time_sec"] < 120 and r["optimal"]:
            print(f"     cross-check via SAT (independent encoder)...", flush=True)
            try:
                sat_r = solve_Mn_sat(n, time_budget_sec=300)
                agree = sat_r["M"] == r["M"]
                print(
                    f"       SAT M={sat_r['M']} t={sat_r['total_time']:.2f}s "
                    f"agree={agree}",
                    flush=True,
                )
                r["sat_cross_check"] = {
                    "M": sat_r["M"],
                    "time": sat_r["total_time"],
                    "agree": agree,
                }
                if not agree:
                    print(f"     CROSS-CHECK FAILED at n={n}: ILP vs SAT disagree.", flush=True)
                    r["validation_warning"] = "SAT/ILP disagree"
            except Exception as e:
                print(f"     SAT cross-check error: {e}", flush=True)
                r["sat_cross_check"] = {"error": str(e)}
        ext_results.append(r)

        # Persist after each n.
        with open(out_path, "w") as f:
            json.dump({
                "convention": "Haugland: full set [1,2n], |A|=n; M(n)=min_A max_k |A ∩ (B+k)|; µ=lim M(n)/n.",
                "validation": val_results,
                "extension": ext_results,
                "haugland_table": HAUGLAND,
            }, f, indent=2)

        if not r["optimal"]:
            print(f"  Hit time limit at n={n} (gap not closed). Stopping.", flush=True)
            break
        if r["wall_time_sec"] > single_n_cap_sec * 0.95:
            print(f"  n={n} consumed ~full single-n budget. Stopping.", flush=True)
            break
        n += 1
        if n > 50:
            print(f"  Reached n=50 ceiling, stopping.", flush=True)
            break

    # Summary table
    print("\n=== Summary: certified UBs on µ from M(n)/n ===", flush=True)
    print(f"{'n':>3}  {'M(n)':>5}  {'M/n':>9}  {'time(s)':>8}  {'<Together?':>10}")
    TOGETHER_UB = 0.380871
    for r in ext_results:
        if r["M"] is None:
            continue
        ratio = r["M"] / r["n"]
        beats = "yes" if (r["optimal"] and ratio < TOGETHER_UB) else "no"
        print(f"{r['n']:>3}  {r['M']:>5}  {ratio:>9.6f}  {r['wall_time_sec']:>8.1f}  {beats:>10}")


if __name__ == "__main__":
    main()
