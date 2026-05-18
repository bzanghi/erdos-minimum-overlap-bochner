"""PRO-11 full check: smoke (N=200, T=80, bochner=4) + medium (N=2000, T=800, bochner=10)."""
from __future__ import annotations
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
from sdpa_runner import solve_with_sdpa_gmp
import cvxpy as cp


def run(N, T, R, bochner_n, sdpa_timeout, label=""):
    print(f"\n=== {label}  N={N} T={T} R={R} bochner_n={bochner_n} ===", flush=True)
    # CLARABEL first (cheap)
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bochner_n
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    cl = solve_with_dual_extraction(prob)
    print(f"CLARABEL: reported={cl['reported_value']!r}", flush=True)
    print(f"CLARABEL: rigorous LB={cl['rigorous_dual_LB']!r}  status={cl['status']}  time={cl['time']:.2f}s", flush=True)

    # SDPA-GMP
    Omega2, w2, v2, c2, d2, eps2, dlt2, cons2 = build_problem(
        N, T, R, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=bochner_n
    )
    prob2 = cp.Problem(cp.Minimize(Omega2), cons2)
    print(f"SDPA-GMP: serializing + solving (timeout={sdpa_timeout}s) ...", flush=True)
    sd = solve_with_sdpa_gmp(prob2, timeout_sec=sdpa_timeout)
    print(f"SDPA-GMP: phase={sd.get('phase')}  primal={sd.get('primal_obj')!r}", flush=True)
    print(f"SDPA-GMP: dual={sd.get('dual_obj')!r}  gap={sd.get('duality_gap')!r}  digits={sd.get('precision_digits')!r}", flush=True)
    print(f"SDPA-GMP: m={sd.get('m')} blocks={len(sd.get('block_structure') or ())}  "
          f"serialize={sd.get('serialize_sec'):.2f}s  solve={sd.get('solve_sec')}s  total={sd.get('runtime_sec'):.2f}s", flush=True)

    if (cl.get("reported_value") is not None
        and sd.get("primal_obj") is not None
        and sd.get("dual_obj") is not None):
        diff_pri = sd["primal_obj"] - cl["reported_value"]
        diff_dual_vs_rep = sd["dual_obj"] - cl["reported_value"]
        if cl.get("rigorous_dual_LB") is not None:
            diff_dual_vs_lb = sd["dual_obj"] - cl["rigorous_dual_LB"]
        else:
            diff_dual_vs_lb = None
        print(f"\n  >> SDPA primal - CLARABEL reported = {diff_pri:+.3e}", flush=True)
        print(f"  >> SDPA dual   - CLARABEL reported = {diff_dual_vs_rep:+.3e}", flush=True)
        if diff_dual_vs_lb is not None:
            print(f"  >> SDPA dual   - CLARABEL rig LB  = {diff_dual_vs_lb:+.3e}", flush=True)

    return {
        "label": label, "N": N, "T": T, "R": R, "bochner_n": bochner_n,
        "clarabel_reported": cl.get("reported_value"),
        "clarabel_rigorous_LB": cl.get("rigorous_dual_LB"),
        "clarabel_status": cl.get("status"),
        "clarabel_time": cl.get("time"),
        "sdpa_phase": sd.get("phase"),
        "sdpa_primal_obj": sd.get("primal_obj"),
        "sdpa_dual_obj": sd.get("dual_obj"),
        "sdpa_duality_gap": sd.get("duality_gap"),
        "sdpa_rigorous_LB": sd.get("rigorous_dual_LB"),
        "sdpa_precision_digits": sd.get("precision_digits"),
        "sdpa_iterations": sd.get("iterations"),
        "sdpa_runtime_sec": sd.get("runtime_sec"),
        "sdpa_serialize_sec": sd.get("serialize_sec"),
        "sdpa_solve_sec": sd.get("solve_sec"),
        "sdpa_m": sd.get("m"),
        "sdpa_n_blocks": len(sd.get("block_structure") or ()),
    }


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    results = []
    results.append(run(200, 80, 10, bochner_n=4, sdpa_timeout=600, label="SMOKE"))
    # Medium scale (the brief: N=2000, T=800, bochner_n=10)
    results.append(run(2000, 800, 10, bochner_n=10, sdpa_timeout=10800, label="MEDIUM"))

    out_path = os.path.join(os.path.dirname(__file__),
                            "..", "parallel_results", "pro11_sdpa_s_serializer.json")
    out_path = os.path.abspath(out_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nResults saved to: {out_path}")
