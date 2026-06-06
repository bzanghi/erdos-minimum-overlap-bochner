"""Test additivity of T5pk constraints.

If baseline + T5pk_k1 gives ΔΩ_1, and baseline + T5pk_k1 + T5pk_k3 gives
ΔΩ_{1+3}, then the marginal gain from k=3 is ΔΩ_{1+3} - ΔΩ_1.

This determines whether T5pk_k=k for k > 1 is a genuinely new lever
or already subsumed by the existing T5p_k=1.
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from fast_eval import baseline_solve, solve_with_extra
from dsl import family_T5pk


CFG = {
    "N": 500, "T": 200, "R": 8, "h1": 0.004, "h2": 0.004,
    "p1": 0.3875, "p2": 0.3875, "q1": -0.02, "q2": 0.02, "bochner_n": 10,
}


def combined_constraint(*k_list):
    """Constraint fn that adds T5p_k for each k in k_list."""
    fns = [family_T5pk(k) for k in k_list]
    def cfn(Omega, w, v, c, d, eps, dlt, cfg):
        out = []
        for fn in fns:
            out.extend(fn(Omega, w, v, c, d, eps, dlt, cfg))
        return out
    return cfn


def measure(name, cfn, base):
    t0 = time.time()
    val, status = solve_with_extra(cfn, **CFG)
    dt = time.time() - t0
    if val is None:
        print(f"  {name:<28s} FAIL ({status})", flush=True)
        return None
    delta = val - base
    flag = "***" if delta > 1e-4 else ("** " if delta > 1e-5 else "   ")
    print(f"  {name:<28s} Ω={val:.7f} ΔΩ={delta:+9.3e} {flag} {dt:>5.1f}s", flush=True)
    return delta


def main():
    base, base_status = baseline_solve(**CFG)
    print(f"Baseline (N=500): Ω = {base:.8f}, status = {base_status}\n", flush=True)

    print("=== Individual constraints ===", flush=True)
    d1 = measure("k=1 only", combined_constraint(1), base)
    d3 = measure("k=3 only", combined_constraint(3), base)
    d5 = measure("k=5 only", combined_constraint(5), base)
    d7 = measure("k=7 only", combined_constraint(7), base)
    d9 = measure("k=9 only", combined_constraint(9), base)

    print("\n=== Pairs ===", flush=True)
    d13 = measure("k=1+3", combined_constraint(1, 3), base)
    d15 = measure("k=1+5", combined_constraint(1, 5), base)
    d35 = measure("k=3+5", combined_constraint(3, 5), base)

    print("\n=== Triples ===", flush=True)
    d135 = measure("k=1+3+5", combined_constraint(1, 3, 5), base)
    d137 = measure("k=1+3+7", combined_constraint(1, 3, 7), base)
    d1357 = measure("k=1+3+5+7", combined_constraint(1, 3, 5, 7), base)

    print("\n=== Full odd sweep ===", flush=True)
    d_all = measure("k=1,3,5,7,9", combined_constraint(1, 3, 5, 7, 9), base)
    d_all_odd = measure("k=1,3,5,7,9,11", combined_constraint(1, 3, 5, 7, 9, 11), base)
    d_all_more = measure("k=1,3,5,7,9,11,13,15", combined_constraint(1, 3, 5, 7, 9, 11, 13, 15), base)

    # ==== Marginals ====
    print("\n=== Marginal analysis ===", flush=True)
    if d1 and d13:
        marg3_given1 = d13 - d1
        print(f"  Marginal gain from k=3 given k=1:     ΔΔΩ = {marg3_given1:+.4e}")
    if d1 and d3 and d13:
        sup_additive = d1 + d3
        print(f"  Independence prediction (d1 + d3):     {sup_additive:+.4e}")
        print(f"  Actual d13:                            {d13:+.4e}")
        print(f"  Subadditivity (d13 / (d1+d3)):         {d13/sup_additive:.3f}")
    if d_all and d1:
        marg_all_given1 = d_all - d1
        print(f"\n  Marginal gain from all odd 3-9 given k=1: ΔΔΩ = {marg_all_given1:+.4e}")
    if d_all_more and d1:
        marg_more = d_all_more - d1
        print(f"  Marginal gain from all odd 3-15 given k=1: ΔΔΩ = {marg_more:+.4e}")

    out = {
        "baseline": base, "config": CFG,
        "individual": {"k1": d1, "k3": d3, "k5": d5, "k7": d7, "k9": d9},
        "pairs": {"k1_3": d13, "k1_5": d15, "k3_5": d35},
        "triples": {"k1_3_5": d135, "k1_3_7": d137, "k1_3_5_7": d1357},
        "all_odd": {"k1_3_5_7_9": d_all, "k1_3_5_7_9_11": d_all_odd, "k1_3_5_7_9_11_13_15": d_all_more},
    }
    out_path = Path(__file__).parent.parent.parent / "data" / "ai_t5pk_additivity.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
