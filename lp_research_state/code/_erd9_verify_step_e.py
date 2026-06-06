"""ERD-9 verification: solve relaxed SDP → extract Σm·λ → predict C_explicit →
solve exact SDP → assert Ω_exact ≤ C_explicit (Step E theorem).

This is a falsification test of LEVER_I_PRIME_FINAL.md Theorem 3.
"""
from __future__ import annotations
import sys, time, json, warnings
import numpy as np
import cvxpy as cp

CODE = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code"
sys.path.insert(0, CODE)
warnings.simplefilter("ignore")

from white_full_convex import build_problem
from white_full_convex_exact import build_problem_exact


def extract_cell_envelope_multipliers(cons, R: int):
    """Cell-envelope cosine: cons[8..8+2R-1]; sine: cons[8+2R..8+6R-1]."""
    cos_start = 8
    cos_end = cos_start + 2 * R
    sin_start = cos_end
    sin_end = sin_start + 4 * R

    def _scalar(ci):
        dv = ci.dual_value
        if dv is None: return None
        if np.isscalar(dv): return float(dv)
        arr = np.asarray(dv).ravel()
        return float(arr[0]) if arr.size == 1 else None

    lam = [_scalar(cons[i]) or 0.0 for i in range(cos_start, cos_end)]
    sigma_pairs = []
    for m in range(1, 2 * R + 1):
        s1 = _scalar(cons[sin_start + 2 * (m - 1)]) or 0.0
        s2 = _scalar(cons[sin_start + 2 * (m - 1) + 1]) or 0.0
        sigma_pairs.append((s1, s2))
    return lam, sigma_pairs


def verify(N, T, R, bochner_n, h1, h2, p1, p2, q1, q2, label):
    print(f"\n{'='*70}")
    print(f"ERD-9 verification: {label}")
    print(f"  config: N={N}, T={T}, R={R}, bn={bochner_n}")
    print(f"{'='*70}")

    # Step 1: solve RELAXED, extract multipliers
    print("\n[1] Solve RELAXED (cell-min)...")
    t0 = time.time()
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, q1, q2, bochner_n=bochner_n
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    omega_rel = prob.value
    print(f"  Ω_relaxed = {omega_rel:.7f}, status = {prob.status}, time = {time.time()-t0:.1f}s")

    lam, sigma_pairs = extract_cell_envelope_multipliers(cons, R)
    lam_arr = np.array([abs(x) for x in lam])
    sig_arr = np.array([abs(p[0]) + abs(p[1]) for p in sigma_pairs])
    ms = np.arange(1, 2 * R + 1)
    Sml = float((ms * lam_arr).sum())
    Sm3l = float((ms ** 3 * lam_arr).sum())
    Sms = float((ms * sig_arr).sum())
    Sm3s = float((ms ** 3 * sig_arr).sum())
    print(f"  Σ m·λ (cos) = {Sml:.4f}")
    print(f"  Σ m·σ (sin) = {Sms:.4f}")

    # Step 2: predict C_explicit via corrected residual formula
    Omega_est = max(omega_rel, 0.25)  # for Ω in residual formula
    factorA = np.pi / (2 * N)
    factorB = np.pi ** 2 * Omega_est / (3 * N ** 3)
    residual_cos = factorA * Sml + factorB * Sm3l
    residual_sin = factorA * Sms + factorB * Sm3s
    residual_total = residual_cos + residual_sin
    C_explicit = omega_rel + residual_total

    print(f"\n[2] Predicted C_explicit via corrected formula:")
    print(f"  residual_cos = {residual_cos:.4e}")
    print(f"  residual_sin = {residual_sin:.4e}")
    print(f"  total        = {residual_total:.4e}")
    print(f"  C_explicit   = {omega_rel:.7f} + {residual_total:.4e} = {C_explicit:.7f}")

    # Step 3: solve EXACT, measure Ω_exact
    print(f"\n[3] Solve EXACT (cell integral)...")
    t0 = time.time()
    Omega, *_ , cons = build_problem_exact(
        N, T, R, h1, h2, p1, p2, q1, q2, bochner_n=bochner_n
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    omega_exact = prob.value
    print(f"  Ω_exact = {omega_exact:.7f}, status = {prob.status}, time = {time.time()-t0:.1f}s")

    # Step 4: verify Ω_exact ≤ C_explicit (Step E theorem)
    margin = C_explicit - omega_exact
    delta_measured = omega_exact - omega_rel
    print(f"\n[4] Step E theorem check:")
    print(f"  Ω_exact - Ω_relaxed = {delta_measured:+.4e}  (measured residual)")
    print(f"  Predicted residual  = {residual_total:.4e}")
    print(f"  C_explicit - Ω_exact = {margin:+.4e}")
    if margin >= -1e-7:
        print(f"  ✓ THEOREM HOLDS (Ω_exact ≤ C_explicit by {margin:+.2e})")
    else:
        print(f"  ✗ THEOREM VIOLATED (Ω_exact > C_explicit by {-margin:.2e})")

    return dict(
        label=label, N=N, T=T, R=R, bochner_n=bochner_n,
        omega_relaxed=omega_rel,
        Sml=Sml, Sms=Sms,
        residual_total=residual_total,
        C_explicit=C_explicit,
        omega_exact=omega_exact,
        delta_measured=delta_measured,
        margin=margin,
        theorem_holds=margin >= -1e-7,
    )


def main():
    h1, h2, p1, p2, q1, q2 = 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02
    results = []
    for cfg in [
        dict(N=200, T=100, R=5, bochner_n=8),
        dict(N=500, T=200, R=8, bochner_n=10),
        dict(N=1000, T=500, R=10, bochner_n=15),
    ]:
        try:
            r = verify(cfg["N"], cfg["T"], cfg["R"], cfg["bochner_n"],
                       h1, h2, p1, p2, q1, q2,
                       label=f"row4 N={cfg['N']}")
            results.append(r)
        except Exception as e:
            print(f"  FAILED: {e}")

    out = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/data/exact_integral_verification.json"
    json.dump(results, open(out, 'w'), indent=2, default=str)
    print(f"\nWrote {out}")
    print("\n" + "=" * 70)
    print(f"SUMMARY: theorem held in {sum(1 for r in results if r['theorem_holds'])}/{len(results)} configs")


if __name__ == "__main__":
    main()
