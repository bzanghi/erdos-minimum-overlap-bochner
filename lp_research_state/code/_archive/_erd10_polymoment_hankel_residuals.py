"""ERD-10: extract poly_moment and Hankel-PSD duals at row 4 N=3000, derive
per-family residual bounds, and add to F3 full-stack picture.

Residual bound derivations:

(A) poly_moment Hausdorff (k=2..14):
    Constraint: m_k_truncated ≥ -tail_bound_k(T)  [LINEAR; one per even k]
    The exact constraint is m_k_exact ≥ 0; relaxation slack = tail_bound_k.
    Residual on Ω from tightening: μ_k · tail_bound_k per constraint.
    Total: Σ_k μ_k · tail_bound_k.

(B) Hankel-PSD: M_n(f)_Hankel ≥ 0 via truncated moments m_var ± ε_k.
    Relaxation: m_var[k] ≈ m_truncated[k] ± tail_bound_k (slack vars).
    Residual: ||Z_H||_2 · sum of tail_bound_k absorbed.

We use the existing path_b_with_polymoment infrastructure but extract all
constraint duals (not just the path-B handles).
"""
from __future__ import annotations
import sys, time, json, warnings
import cvxpy as cp
import numpy as np

CODE = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code"
sys.path.insert(0, CODE)
warnings.simplefilter("ignore")

from path_b_analytical import build_problem_with_dual_handles
from poly_moment import build_even_moment_nonneg_constraints, build_even_hankel_psd, even_moment_tail_bound


def main():
    # Match Phase 5 production except smaller N for speed
    N, T, R = 3000, 1200, 10
    bochner_n = 20
    pm_k_max = 14
    hankel_n = 6
    h_c, p_c = 0.004, 0.3875
    q1, q2 = -0.02, 0.02

    print(f"=== ERD-10: poly_moment + Hankel residual extraction ===")
    print(f"Config: N={N}, T={T}, R={R}, bn={bochner_n}, pm_k_max={pm_k_max}, hankel_n={hankel_n}")
    print(f"Row 4 (h={h_c}, p={p_c}, q=[{q1},{q2}])\n")

    # Build SDP
    Omega, cons_base, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bochner_n
    )

    # Add poly_moment constraints; capture handle for dual extraction
    pm_cons, pm_tails = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k_max)

    # Add Hankel-PSD constraints
    hk_cons, m_var, hk_tails = build_even_hankel_psd(H["c"], H["d"], T, n_hankel=hankel_n)

    all_cons = list(cons_base) + list(pm_cons) + list(hk_cons)
    print(f"Total constraints: base={len(cons_base)}, pm={len(pm_cons)}, hankel={len(hk_cons)}, total={len(all_cons)}")

    prob = cp.Problem(cp.Minimize(Omega), all_cons)
    t0 = time.time()
    prob.solve(solver="CLARABEL", verbose=False)
    dt = time.time() - t0
    print(f"\nSolve done: Ω = {prob.value:.7f}, status = {prob.status}, time = {dt:.1f}s")

    # ===== Extract poly_moment dual multipliers =====
    print(f"\n--- (A) poly_moment dual multipliers ---")
    print(f"{'k':>3s} {'tail_bound':>14s} {'μ_k (dual)':>14s} {'μ_k · tail':>14s}")
    pm_residual_total = 0.0
    pm_per_k = {}
    for i, ck in enumerate(pm_cons):
        # Recover k from pm_tails dict (ordered by k from k_min to k_max)
        # Actually pm_cons is in order, and pm_tails is dict {k: tb}. The k's
        # iterate over even k ∈ {2, 4, ..., k_max}
        k_list = sorted(pm_tails.keys())
        k = k_list[i] if i < len(k_list) else None
        if k is None: continue
        tb = pm_tails[k]
        dv = ck.dual_value
        if dv is None: dv = 0.0
        if not np.isscalar(dv):
            dv = float(np.asarray(dv).ravel()[0]) if np.asarray(dv).size > 0 else 0.0
        else:
            dv = float(dv)
        residual_k = abs(dv) * tb
        pm_residual_total += residual_k
        pm_per_k[k] = {"tail_bound": tb, "mu": dv, "residual": residual_k}
        print(f"{k:>3d} {tb:>14.4e} {dv:>+14.4e} {residual_k:>14.4e}")
    print(f"\nTotal poly_moment residual = Σ_k |μ_k|·tail_k = {pm_residual_total:.4e}")

    # ===== Extract Hankel-PSD dual =====
    print(f"\n--- (B) Hankel-PSD dual multiplier ---")
    # The Hankel-PSD constraint is the LAST one in hk_cons typically (PSD constraint
    # on the moment Hankel matrix). Other constraints are slack-variable bounds.
    hankel_residual_total = 0.0
    for i, ck in enumerate(hk_cons):
        dv = ck.dual_value
        s = str(ck)[:80]
        if dv is None:
            print(f"  hk[{i}] dual=None: {s}")
            continue
        if np.isscalar(dv):
            print(f"  hk[{i}] dual={float(dv):+.4e}: {s}")
        else:
            arr = np.asarray(dv)
            if arr.ndim == 2:
                # PSD constraint matrix dual
                Z = arr
                spec = float(np.max(np.abs(np.linalg.eigvalsh(Z))))
                trace = float(np.trace(Z))
                # Hankel truncation residual: ||Z||_2 * Σ tail_bounds
                tail_sum = sum(hk_tails.values()) if isinstance(hk_tails, dict) else 0.0
                hankel_residual_total += spec * tail_sum
                print(f"  hk[{i}] PSD-block: ||Z||_2 = {spec:.4e}, tr(Z) = {trace:.4e}")
                print(f"         hankel residual bound = ||Z||·Σtail = {spec:.4e}·{tail_sum:.4e} = {spec*tail_sum:.4e}")
            else:
                # Bound on slack variable (scalar)
                arr_flat = arr.ravel()
                print(f"  hk[{i}] vec-dual: shape={arr.shape}, |max|={float(np.abs(arr_flat).max()):.4e}")
    print(f"\nTotal Hankel-PSD residual estimate = {hankel_residual_total:.4e}")

    # ===== Summary + F3 update =====
    print(f"\n{'='*70}")
    print("SUMMARY: full-stack per-family residual at row 4 N=3000")
    print(f"{'='*70}")
    print(f"  Cell-envelope cos+sin (Step E baseline at N=3000): 6.2e-4 (estimated)")
    print(f"  Bochner-PSD truncation (F3, row7 N=3000):           2.2e-4")
    print(f"  poly_moment k=2..{pm_k_max} (THIS RUN):              {pm_residual_total:.2e}")
    print(f"  Hankel-PSD n={hankel_n} (THIS RUN):                  {hankel_residual_total:.2e}")
    sum_naive = 6.2e-4 + 2.2e-4 + pm_residual_total + hankel_residual_total
    print(f"  Naive sum:                                            {sum_naive:.2e}")
    LB = 0.3801279
    UB = 0.380871
    C_naive = LB + sum_naive
    print(f"  C_total (naive sum)                                 = {C_naive:.7f}")
    print(f"  Compare to μ_UB = {UB}")
    print(f"  C_total - UB = {C_naive - UB:+.4e}")
    if C_naive < UB:
        print(f"  ✓ Full-stack naive-sum theorem NON-VACUOUS")
    else:
        print(f"  ✗ Full-stack naive-sum theorem VACUOUS")

    out = {
        "config": dict(N=N, T=T, R=R, bochner_n=bochner_n, pm_k_max=pm_k_max, hankel_n=hankel_n),
        "row": "row4", "h_c": h_c, "p_c": p_c,
        "Omega": float(prob.value), "status": prob.status,
        "poly_moment": {"per_k": pm_per_k, "total_residual": pm_residual_total},
        "hankel_psd": {"total_residual": hankel_residual_total},
        "summary": {
            "cell_env_estimate_N3000": 6.2e-4,
            "bochner_residual_N3000": 2.2e-4,
            "polymom_residual": pm_residual_total,
            "hankel_residual": hankel_residual_total,
            "naive_sum": sum_naive,
            "C_naive": C_naive,
            "UB": UB,
            "margin": UB - C_naive,
            "non_vacuous": C_naive < UB,
        },
    }
    json.dump(out, open("/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/data/full_saturation_residuals.json", 'w'),
              indent=2, default=str)
    print(f"\nWrote lp_research_state/data/full_saturation_residuals.json")


if __name__ == "__main__":
    main()
