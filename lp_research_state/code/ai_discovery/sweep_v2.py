"""Phase 2 sweep v2: larger N, include k=1 (existing T5p) as calibration.

Tests:
1. baseline at N=500
2. T5p_k=1..10 (each k tested individually)
3. F4_sumcos with various theta vectors (parameter search)
4. Two cell-envelope high-freq extensions
5. Combinations: T5p_k1 + T5p_k_other (does it help to STACK)
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from fast_eval import baseline_solve, solve_with_extra
from dsl import family_T5pk, family_T5p_sumcos, family_cell_env_high_freq, family_fejer


# Larger fast-eval config to amplify signals
CFG = {
    "N": 500,
    "T": 200,
    "R": 8,
    "h1": 0.004,
    "h2": 0.004,
    "p1": 0.3875,
    "p2": 0.3875,
    "q1": -0.02,
    "q2": 0.02,
    "bochner_n": 10,
}


def measure(name, cfn, base, results=None, out_path=None):
    t0 = time.time()
    val, status = solve_with_extra(cfn, **CFG)
    dt = time.time() - t0
    if val is None:
        print(f"  {name:<35s} FAIL ({status}) in {dt:.1f}s", flush=True)
        r = None
    else:
        delta = val - base
        flag = "***" if delta > 1e-4 else ("** " if delta > 1e-5 else "   ")
        print(f"  {name:<35s} Ω={val:.7f} ΔΩ={delta:+9.3e} {flag} {dt:>5.1f}s {status}", flush=True)
        r = {"Omega": val, "delta": delta, "status": status, "time_s": dt}
    # Incremental save after every measurement
    if results is not None and out_path is not None:
        results["families"][name] = r
        out_path.write_text(json.dumps(results, indent=2, default=str))
    return r


def main():
    print("=" * 80)
    print(f"Phase 2 sweep v2 — config: N={CFG['N']}, T={CFG['T']}, R={CFG['R']}, bn={CFG['bochner_n']}")
    print("=" * 80)
    t0 = time.time()
    base, base_status = baseline_solve(**CFG)
    dt = time.time() - t0
    print(f"Baseline: Ω = {base:.8f}, status = {base_status}, time = {dt:.2f}s")
    print()

    out = Path(__file__).parent.parent.parent / "data" / "ai_constraint_sweep_v2.json"
    results = {"config": CFG, "baseline": {"Omega": base, "status": base_status},
               "families": {}}

    print("Series A: T5pk for k = 1..15 (individual)", flush=True)
    for k in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15]:
        measure(f"T5pk_k{k}", family_T5pk(k), base, results, out)

    print(flush=True)
    print("Series B: T5p_sumcos with varied theta", flush=True)
    series_B = [
        ("uniform_5",       np.ones(5)),
        ("uniform_10",      np.ones(10)),
        ("decay_1/k_5",     np.array([1/k for k in range(1, 6)])),
        ("decay_1/k_10",    np.array([1/k for k in range(1, 11)])),
        ("decay_1/k2_10",   np.array([1/k**2 for k in range(1, 11)])),
        ("peak_k3",     np.array([0, 0, 1, 0, 0])),
        ("peak_k5",     np.array([0, 0, 0, 0, 1])),
        ("peak_k7",     np.array([0]*6 + [1])),
        ("sparse_k357",  np.array([0, 0, 1, 0, 1, 0, 1])),
    ]
    for name, theta in series_B:
        measure(f"sumcos_{name}", family_T5p_sumcos(theta), base, results, out)

    print(flush=True)
    print("Series C: Cell-envelope at higher freq", flush=True)
    for m_high in [2, 4, 8, 12]:
        measure(f"cell_env_high_m{m_high}", family_cell_env_high_freq(m_high), base, results, out)

    print(flush=True)
    print("Series D: Fejér kernel of various degrees", flush=True)
    for n in [3, 5, 10, 20, 30]:
        measure(f"fejer_n{n}", family_fejer(n), base, results, out)

    # ===== Sorted summary =====
    print()
    print("=" * 80)
    print("Top 15 by ΔΩ:")
    print("=" * 80)
    valid = [(n, r) for n, r in results["families"].items() if r.get("delta") is not None]
    valid.sort(key=lambda x: x[1]["delta"], reverse=True)
    for name, r in valid[:15]:
        flag = "***" if r["delta"] > 1e-4 else ("**" if r["delta"] > 1e-5 else "  ")
        print(f"  {flag} {name:<35s} ΔΩ = {r['delta']:+.3e}")

    print(f"\nFinal results saved to {out}")


if __name__ == "__main__":
    main()
