"""Test genuinely-novel constraint families NOT subsumed by existing stack.

Candidates:
A. Entropy bound: ∫f log f ≥ -log 2 (Jensen on convex x log x).
   Validity: f is a prob density on [-2,2] with ∫f=1, |support|≤2;
   the uniform distribution on [-1,1] minimizes ∫f log f at -log 2.
   But MIN-overlap optimizer wants spread-out f (high entropy),
   so this constraint is AWAY from binding. UNLIKELY to give ΔΩ.

B. UPPER entropy bound: ∫f log f ≤ -log 2 + ε  (forces f to be
   spread out). NOT VALID — f can have lower entropy easily.
   SKIP.

C. (1-f) cell-envelope: apply (W.1)-style cell-min relaxation to
   (1-f). Coefficients flip sign. Validity: standard.

D. Bochner of (f * f): the autocorrelation M = f*f also is a density
   on [-2,2] with ∫M=1. Apply Bochner to M's Fourier coefficients.
   M̂(k) = 2|f̂(k)|² — QUADRATIC in c, d. Bochner is PSD of Toeplitz
   matrix of these — encoding may need lifting.

E. Higher-order Hausdorff (k > 14): extension of poly_moment.

Each tested in isolation and combined with the existing T5p stack.
"""
from __future__ import annotations
import sys, time, json
from pathlib import Path
import numpy as np
import cvxpy as cp

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from fast_eval import baseline_solve, solve_with_extra
from dsl import family_T5pk


CFG = {
    "N": 500, "T": 200, "R": 8, "h1": 0.004, "h2": 0.004,
    "p1": 0.3875, "p2": 0.3875, "q1": -0.02, "q2": 0.02, "bochner_n": 10,
}


# ============================================================================
# Family A: Entropy lower bound ∫f log f ≥ -log 2
# ============================================================================
def family_entropy_lower(threshold: float = -float(np.log(2))):
    """∫f log f ≥ threshold.

    Discretely: L·Σ(w_j log w_j + v_j log v_j) ≥ threshold.

    cvxpy entropy: cp.entr(x) = -x log x. So our LHS = -L·Σ(entr(w) + entr(v)).
    Constraint: -L·Σ(entr(w)+entr(v)) ≥ threshold  ⇔  L·Σ(entr(w)+entr(v)) ≤ -threshold = log 2.
    """
    def cfn(Omega, w, v, c, d, eps, dlt, cfg):
        L = 2.0 / cfg["N"]
        return [L * (cp.sum(cp.entr(w)) + cp.sum(cp.entr(v))) <= -threshold]
    return cfn


# ============================================================================
# Family C: Cell-envelope on (1-f)
# ============================================================================
def family_one_minus_f_cell_env():
    """Apply the (W.1) cell-envelope cosine constraint to (1-f).

    With f → 1-f: ∫(1-f)·cos(πmx/2) gives the analog of (am, bm) but for
    (1-f). Since (1-f) ≥ 0 and bounded by 1, the cell-envelope cosine
    constraint applies with f → 1-f.

    For (1-f): the SDP variables transform as w_j → 1 - w_j, v_j → 1 - v_j.
    f̂(0) = 1/2 → 1 - 1/2 = 1/2 (since ∫f=1, ∫(1-f) over [-2,2] = 4-1 = 3,
    and (1-f)̂(0) = (1/4)·3 = 3/4). Other Fourier coefs flip sign.

    Hmm — this complicates because the support of 1-f is different from f.
    Let me try a simpler form: just transcribe the cell-envelope constraint
    with sign flips. This MAY NOT BE VALID without careful derivation.
    """
    def cfn(Omega, w, v, c, d, eps, dlt, cfg):
        # Skip for now: needs careful derivation
        return []
    return cfn


# ============================================================================
# Family E: Higher poly_moment via direct encoding
# ============================================================================
def family_higher_hausdorff(k_extra: int = 16):
    """Add ∫x^k f(x) dx ≥ tail_bound for k=k_extra (must be even).

    Discretely on cells centered at jL - L/2: x_j = (j - 1/2)L for j=1..N
    (positive side), x_j = -(j - 1/2)L for negative side.

    ∫x^k f(x) dx ≈ L·Σ x_j^k · w_j + L·Σ (-x_j)^k · v_j (positive + negative)
                 = L·Σ x_j^k (w_j + v_j) for even k.

    For k even and ≥ 2: ∫x^k f ≥ 0 always (since x^k ≥ 0 and f ≥ 0).
    The Hausdorff moments are nonneg.

    To make it nontrivially binding: replace ≥ 0 with ≥ some Fourier-tail
    LOWER bound from poly_moment.py. For now, just use ≥ 0 as a baseline.
    """
    def cfn(Omega, w, v, c, d, eps, dlt, cfg):
        N = cfg["N"]
        L = 2.0 / N
        j = np.arange(1, N+1)
        x = (j - 0.5) * L
        # Even k → both sides same:
        coefs = x ** k_extra
        # Constraint: L · sum(coefs · (w + v)) ≥ 0
        # This is vacuous (sum of nonneg quantities is nonneg). But to give
        # tightening, link to Fourier coefficients via x^k = ... cosines.
        # The simplest binding form: ∫x^k f ≤ ∫x^k · Ω·1[-1,1] = Ω · 2/(k+1)
        # (since on |x|≤1, max f = Ω, integrated x^k = 2/(k+1) for even k).
        # Actually this needs to be tighter. Skip for now.
        return []
    return cfn


# ============================================================================
# Run
# ============================================================================
def measure(name, cfn, base):
    t0 = time.time()
    val, status = solve_with_extra(cfn, **CFG)
    dt = time.time() - t0
    if val is None:
        print(f"  {name:<35s} FAIL ({status}) {dt:.1f}s", flush=True)
        return None
    delta = val - base
    flag = "***" if delta > 1e-4 else ("** " if delta > 1e-5 else "   ")
    print(f"  {name:<35s} Ω={val:.7f} ΔΩ={delta:+9.3e} {flag} {dt:>5.1f}s", flush=True)
    return delta


def main():
    base, status = baseline_solve(**CFG)
    print(f"Baseline (N=500): Ω = {base:.8f}, status = {status}\n", flush=True)

    print("=== Test novel families ===", flush=True)
    measure("entropy_lower(-log2)", family_entropy_lower(-float(np.log(2))), base)
    measure("entropy_lower(-0.5)", family_entropy_lower(-0.5), base)
    measure("entropy_lower(-0.2)", family_entropy_lower(-0.2), base)
    measure("entropy_lower(-0.1)", family_entropy_lower(-0.1), base)
    measure("entropy_lower(-0.05)", family_entropy_lower(-0.05), base)

    print("\n=== Combine entropy with T5p_k1 (the existing lever) ===", flush=True)
    def combined(threshold):
        ent_fn = family_entropy_lower(threshold)
        t5p_fn = family_T5pk(1)
        def cfn(Omega, w, v, c, d, eps, dlt, cfg):
            return ent_fn(Omega, w, v, c, d, eps, dlt, cfg) + t5p_fn(Omega, w, v, c, d, eps, dlt, cfg)
        return cfn
    measure("entropy(-log2) + T5p_k1", combined(-float(np.log(2))), base)
    measure("entropy(-0.5) + T5p_k1", combined(-0.5), base)
    measure("entropy(-0.2) + T5p_k1", combined(-0.2), base)
    measure("entropy(-0.1) + T5p_k1", combined(-0.1), base)
    measure("entropy(-0.05) + T5p_k1", combined(-0.05), base)

    print("\n=== Just T5p_k1 baseline for comparison ===", flush=True)
    measure("T5p_k1 alone", family_T5pk(1), base)


if __name__ == "__main__":
    main()
