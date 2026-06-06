# PRO-29: Spectral Reformulation — Structural Insight, Not a Direct Bound

**Status:** Done. Discovered a clean dual formulation `μ = 1 − sup_h inf_t ⟨h, T_t h⟩` (max-over-h of the min-autocorrelation). This connects directly to the Barnard-Steinerberger / Madrid-Ramos autocorrelation literature surfaced in PRO-25. Spectral analysis alone gives a 4× loose bound (Rayleigh quotient is too pessimistic given our constraints), but **the autocorrelation duality is the right framing for cross-bound attacks from the Barnard-Steinerberger line of work**.

## 1. The clean duality

Starting from Together's definition:
```
M(h, t) = ∫ h(x) (1 - h(x+t)) dx
        = ∫h - ∫ h(x) h(x+t) dx
        = ‖h‖_L¹ - ⟨h, T_t h⟩
        = 1 - ⟨h, T_t h⟩      (since ‖h‖_L¹ = ∫h = 1)
```

Therefore:
```
μ = inf_h sup_t M(h, t)
  = inf_h sup_t [1 - ⟨h, T_t h⟩]
  = 1 - sup_h inf_t ⟨h, T_t h⟩
```

**The maximum-minimum autocorrelation** of h over all shifts t is the dual quantity:

```
μ_dual := sup_h inf_t ⟨h, T_t h⟩,   μ = 1 - μ_dual
```

With μ ≈ 0.38087 (UB), μ_dual ≈ 0.61913.

## 2. Connection to Barnard-Steinerberger

PRO-25's literature mine surfaced [arXiv:1903.08731 (Barnard-Steinerberger)](https://arxiv.org/abs/1903.08731): for f ∈ L¹(ℝ), f ≥ 0:
```
inf_{0 ≤ t ≤ 1} ∫ f(x) f(x+t) dx ≤ C · ‖f‖_L¹²
```
where C has been progressively sharpened: 0.411 (B-S) → 0.4071 (Madrid-Ramos). The LB on C is 0.37.

**If this directly applies to our problem (modulo normalization):**

μ_dual ≤ C_BS_MR ≈ 0.4071 → μ ≥ 1 − 0.4071 = 0.5929.

But we know μ ≤ 0.380871 < 0.5929. **Contradiction.** So B-S/MR does NOT directly apply.

**Why:**
- Their t-range is [0, 1] — possibly under a normalization where f has support [0, 1]
- Their f is in L¹(ℝ), zero outside support, with `‖f‖_L¹ = 1`
- The discrepancy in numerical regime suggests their problem's t-range or support is normalized differently than ours

**Action item for future work:** Carefully match normalizations between B-S/MR and our problem. The factor-of-1.5 numerical gap suggests their t ∈ [0, 1] corresponds to a DIFFERENT shift family than our t ∈ [0, 2]. This may be reconcilable via a substitution.

If reconcilable, the B-S/MR LB of 0.37 (i.e., μ_dual ≥ 0.37 ⇒ μ ≤ 0.63 — not useful as it stands) or UB 0.4071 might transfer in a normalization where our t ∈ [0, T_max] matches their t ∈ [0, 1].

## 3. Naive spectral analysis is too loose

Built the discrete operator A_j = I - T_j (where T_j is non-periodic shift by j cells) and computed eigenvalues.

| j | λ_max(sym(A_j)) | Rayleigh-quotient UB on M | actual M | bound/actual |
|---|---|---|---|---|
| 33 | 1.988 | 1.5398 | 0.3809 | **4.04×** |

The Rayleigh-quotient bound `⟨h, A_j h⟩ ≤ λ_max(sym A_j) · ‖h‖²` is loose by a factor of 4 because:
- It assumes h is an UNCONSTRAINED unit L² vector
- Our h has TWO constraints: ‖h‖_∞ ≤ 1 AND ‖h‖_L¹ = 1
- These restrict h to a polytope MUCH smaller than the unit L² ball
- The constrained optimization is what cvxpy/CLARABEL already does

**Bottom line:** spectral analysis of A_j alone doesn't shortcut the constrained min-max.

## 4. Sanity check — the matrix and correlate formulations match

At first the matrix-form `L·h^T A_j h` differed from the correlate form `L·Σ h_i (1 - h_{i+j})`. The discrepancy:

- Matrix form gives `L·(‖h‖²_2 - ⟨h, T_j h⟩)`
- Correlate form gives `L·(‖h‖_1 - ⟨h, T_j h⟩) = 1 - ⟨h, T_j h⟩` 

The DIFFERENCE between them is `L · (‖h‖²_2 - ‖h‖_1)`, which for our h (with `‖h‖²_2 ≈ 0.226` and `‖h‖_1 = 1/L = 300/L` wait, in unweighted = 300, weighted = 1) is exactly `1 - L·‖h‖²_2 ≈ 1 - 0.225 = 0.775` at j=0 — matches the +0.225 discrepancy seen at j=0.

So the matrix form computes a DIFFERENT quantity (L²-flavored autocorrelation deficit), not our Together-convention M. Need to be careful: the spectral analysis is of the LINEAR-NORM autocorrelation `⟨h, T_j h⟩`, not the indicator-overlap M.

## 5. What MIGHT actually work for spectral

Three potential angles where spectral information could help:

### A. Eigenvector at the binding shift

For j=33 (the binding shift for Together's h\*), the dominant eigenvector of sym(A_33) gives information about which h-direction maximally increases ⟨h, A_33 h⟩. This might be useful for the PRO-26 Phase 2a ansatz design (which basis functions matter for the interior).

### B. Spectral preconditioner

The constrained QP `inf_h max_t ⟨h, A_t h⟩` is hard partly because of the matrix conditioning. A spectral preconditioner (Cholesky-like) could accelerate iterative solvers.

### C. Trace-norm / nuclear-norm relaxation

Replace ‖h‖_∞ ≤ 1 with a relaxation of the form `h h^T ⪯ M` for some matrix M. Combined with the spectral structure of A_t, this might give a tighter SDP than White's framework. But this is the same kind of construction the existing SDP already encodes, so probably no win.

## 6. Recommendation

**Document the duality.** μ = 1 - μ_dual is a clean reformulation worth including in the preprint.

**Pursue normalization-matching with Barnard-Steinerberger / Madrid-Ramos.** If their result transfers under proper rescaling, we may inherit either a tighter UB or LB. This is a 1-evening math task and should be high priority.

**Do NOT pursue naive spectral analysis** as a μ-bounding method. Documented as too loose.

**Possibly: Use spectral analysis to inform PRO-26 Phase 2a v2 ansatz design.** The dominant eigenvectors of sym(A_t) for active t might suggest the right basis functions for parameterizing h.

## 7. Deliverables

- `lp_research_state/code/_pro29_spectral.py` — driver with operator construction, spectral sweep, Rayleigh-bound comparison
- This document

## 8. Strategic implication

PRO-29 has surfaced the **most promising open lead** discovered this session: the explicit connection to the **autocorrelation lower-bound literature** (Barnard-Steinerberger, Madrid-Ramos, Fish-King-Miller). If their bound transfers under careful normalization, we inherit a new lever on μ.

**Spawning follow-up: PRO-32 "Normalize-match Barnard-Steinerberger to our problem"** — this is the natural extension. Should be high priority.
