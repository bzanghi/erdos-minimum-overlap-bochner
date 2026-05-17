"""Sweep all candidate constraint families at fast-eval scale.

Measure ΔΩ for each family vs baseline. Save results to JSON.
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from fast_eval import baseline_solve, solve_with_extra, FAST_CONFIG
from dsl import FAMILIES


def main():
    print("=" * 70)
    print("Phase 2: sweep candidate families at fast-eval scale")
    print("=" * 70)
    print(f"Config: {FAST_CONFIG}")
    print()

    t0 = time.time()
    base, base_status = baseline_solve()
    base_dt = time.time() - t0
    print(f"Baseline: Omega = {base:.8f}, status = {base_status}, time = {base_dt:.2f}s")
    print()

    results = {"baseline": {"Omega": base, "status": base_status, "time_s": base_dt},
               "families": {}}

    print(f"{'family':<28s} {'Omega':>12s} {'ΔΩ':>12s} {'time':>7s} {'status':>15s}")
    print("-" * 80)

    for name, factory in FAMILIES.items():
        constraint_fn = factory()
        t0 = time.time()
        val, status = solve_with_extra(constraint_fn)
        dt = time.time() - t0
        if val is None:
            delta = None
            print(f"{name:<28s} {'FAIL':>12s} {'-':>12s} {dt:>6.2f}s {status:>15s}")
        else:
            delta = val - base
            print(f"{name:<28s} {val:>12.8f} {delta:>+12.4e} {dt:>6.2f}s {status:>15s}")
        results["families"][name] = {
            "Omega": val,
            "delta_Omega": delta,
            "status": status,
            "time_s": dt,
        }

    # Top candidates (positive delta = tightens the SDP)
    print()
    print("=" * 70)
    print("Top candidates (sorted by ΔΩ):")
    print("=" * 70)
    valid = [(n, r) for n, r in results["families"].items() if r["delta_Omega"] is not None]
    valid.sort(key=lambda x: x[1]["delta_Omega"] or -1, reverse=True)
    for name, r in valid[:10]:
        print(f"  {name:<28s} ΔΩ = {r['delta_Omega']:+.4e}  status={r['status']}")

    out = Path(__file__).parent.parent.parent / "data" / "ai_constraint_candidates.json"
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
