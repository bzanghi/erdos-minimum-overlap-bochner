# PRO-26 Phase 1: Rechnitzer Technique Transfer to μ — Feasibility Analysis

**Verdict:** 🟢 **Partial transfer: UB-side likely feasible; LB-side blocked by structural mismatch.** Recommend prototyping a UB-side ansatz to potentially tighten Together's 0.380871 to ~30 digits, but expect the LB will still need SDP infrastructure (PRO-11).

## 1. Rechnitzer's pipeline (decoded from arXiv:2602.07292)

Rechnitzer's problem: `ν₂² = inf_f ‖f∗f‖²₂` over `f ∈ L¹(-1/2, 1/2)` with `∫f = 1`, `f ≥ 0`. **No pointwise upper bound on f.** Final result: 128 verified digits.

The pipeline is **four ingredients**:

### A. White's decomposition (eq 8 of preprint)

`C(f) = ‖F∗F‖²₂ = (1/2) Σ_k |F̂(k)|⁴ + (8/π⁴) Σ_m |Σ_k (-1)^k F̂(k)/(2k - (2m+1))|⁴`

Truncation of f̂ to |k| ≤ T converts the variational problem into finite-dim optimization on T coefficients.

### B. Asymptotic ansatz (Section 2)

From the empirical structure of White's near-optimum, derive:
- **Real-space ansatz** (eq 16, 23): `f(x) ≈ (2/π) / √(1 - 4x²)` — **inverse-square-root singularity at ±1/2**.
- **Fourier ansatz** (eq 10): `(-1)^k f̂(k) = (1/√k) Σ_{j=0}^{P-1} a_j / k^j`
- **Generalized real-space ansatz** (eq 26): `f(x) = Σ a_j · C(j, 1/2) · (1-4x²)^{j-1/2}`, with `Σa_j = 1`

This collapses ~4096 free coefficients into ~32 parameters (a_0, ..., a_{31}).

### C. Rigorous evaluation via ball-arithmetic + Bessel asymptotics (Section 3)

Fourier coefficients of (1-4x²)^{j-1/2} are Bessel functions: `F̂(k) = (1/2) Σ_j a_j · J(j, πk/2) · j!·(4/πk)^j`

Define `𝒥(p, k) := J(p, πk/2) · p! · (4/πk)^p`. Then `C(a⃗) = 1/2 + Σ_k [Σ_j a_j 𝒥(j, πk/2)]^4`.

The large-k Bessel asymptotic has period-4 sign behavior (eq 32-35). Tail of the infinite k-sum is bounded by Hurwitz zeta values. **C++ + FLINT ball-arithmetic** to track floating-point errors rigorously.

### D. Hölder-Plancherel LB (Section 4)

For the LB, define dual G(x) periodically with `Ĝ(0) = 0`. Plancherel gives `1 = 2 Σ F̂(k) · Ĝ(k)` (eq 42). Hölder:

`1/16 ≤ (Σ|F̂(k)|⁴)(Σ|Ĝ(k)|^{4/3})³`

Rearranged: `ν₂² ≥ 1/2 + (1/2)(Σ|Ĝ|^{4/3})^{-3}` (eq 44).

For Hölder equality, need `Ĝ(k) ∝ F̂(k)³`. The function `(F∗F∗F)(x)` is constant on (-1/2, 1/2) at the OPTIMUM — so use a curve-fit `b₁√(1-4x²) + b₂(1-4x²)^{3/2}` of the empirical F∗F∗F for a near-tight G.

## 2. Our problem μ — structural comparison

| Property | Rechnitzer ν₂² | Our μ |
|---|---|---|
| Objective | `inf ‖f∗f‖²₂` (smooth L² norm) | `inf_h sup_t M(t)` (min-max, non-smooth) |
| Function class | `f ∈ L¹(-1/2, 1/2), f ≥ 0, ∫f=1` | `h:[0,2]→[0,1], ∫h=1` (pointwise bounded) |
| Optimal function geometry | Inverse-square-root singularities at endpoints | **Bang-bang: h=0 on 28%, h=1 on 10%, interior 62% (PRO-23)** |
| Fourier asymptotics | `(-1)^k f̂(k) ∼ a/√k` | Likely `1/k` decay (indicator-like jumps) |
| Active set | Smooth — F∗F is at sup over all of [-1/2, 1/2] | Discrete: 437 shifts t active to 1e-9 (PRO-23) |
| LB approach | Hölder via Plancherel on L² | Min-max duality (γ-weighted KKT) |

## 3. What transfers, what doesn't

### ✅ Transfers cleanly

**The general METHODOLOGY:**
- Ansatz + finite parameter vector
- BFGS / Newton-Raphson in mpmath
- Rigorous floating-point via ball arithmetic
- Asymptotic-series acceleration of slow-convergent sums

This is universal numerical machinery and applies to ANY high-precision optimization.

### 🟡 Partially transfers (UB side)

**Rechnitzer's real-space ansatz `(1-4x²)^{j-1/2}`** uses inverse-square-root basis functions. For h with bang-bang structure, a natural analog is:

```
h(x) = 1 · 1_{A_+}(x) + 0 · 1_{A_0}(x) + Σ_j a_j φ_j(x) · 1_{A_int}(x)
```

where (A_+, A_0, A_int) are upper-active / lower-active / interior sets, and φ_j are basis functions on A_int.

**Crucial wrinkle:** A_+, A_0 are themselves unknowns. Either:
- Fix them from Together's h* (gives near-tight UB but not improvable beyond his global-search precision)
- Parameterize their boundaries (introduces additional unknowns — boundary positions)

**If we fix sets and only optimize the interior:** the problem becomes finite-dim LP-style. mpmath at 50 digits is feasible. **Expected payoff:** UB tighter than 0.380871, perhaps 10-30 digits, but no better than Together's underlying h* allows.

**If we parameterize set boundaries:** the problem becomes mixed continuous-discrete, much harder. This is where Together's SA + grad descent stops at n=600.

### ❌ Does NOT transfer

**The Hölder-Plancherel LB** (Section 4 of Rechnitzer): specific to L² objectives. For min-max `sup_t M(t)`, the analogous step would be:

`μ = inf_h sup_t M(t) = inf_h sup_{γ ∈ Δ} Σγ_t M(t)`

By Fenchel-Sion: `μ = sup_γ inf_h Σγ_t M(t)`. The inner problem is linear in h. **This IS the KKT functional equation route (LEVER_FUNCTIONAL_EQUATION)**, which we already explored in PRO-23 and got blocked because no clean γ achieves equality.

So Rechnitzer's LB trick doesn't give us anything new.

**The Fourier ansatz `(-1)^k h̃(k) = a/√k`** is wrong for μ. Together's h* has bang-bang structure → Fourier coefficients decay as 1/k (jump-induced), not 1/√k (singularity-induced).

## 4. Recommended Phase 2 (if we proceed)

### Phase 2a — UB side (probability of success ~50%)

1. **Take Together's h* as a starting point.** Discrete with 600 cells.
2. **Identify A_+, A_0, A_int** from the cell values.
3. **Fit interior cells** to a smooth ansatz using a Chebyshev or modified polynomial basis on A_int.
4. **mpmath BFGS** at 50-digit precision on the interior ansatz coefficients, holding A_+ and A_0 fixed.
5. **Evaluate sup_t M(h̃; t)** at 50 digits via FFT-accelerated correlation.
6. **Iterate** until M decreases below 0.3808703.

**Estimated effort:** 1-2 focused sessions of coding + 1 of debugging. The key risk: the bang-bang structure means the gradient w.r.t. interior coefficients has a `(1-h)·h` factor that vanishes near A_+ and A_0 — could lead to flat directions.

### Phase 2b — Boundary refinement (riskier)

Parameterize A_+, A_0 by smooth functions (curve fits). Optimize together with interior coefficients.

This is essentially **a high-precision version of Together's pipeline.** Their global SA + grad descent stops at n=600 due to combinatorial complexity; mpmath BFGS at much higher resolution could plausibly continue.

### What we punt on

**LB side improvement.** Rechnitzer's Hölder trick doesn't transfer; the existing KKT-functional-equation route is the analog, and it's blocked. To tighten μ_LB further, we still need either:
- SDPA-GMP serializer (PRO-11) for ultra-precise SDP solves
- Or a fundamentally new constraint family (open research)

## 5. Honest probability estimate

| Outcome | P |
|---|---|
| Phase 2a succeeds → μ_UB tightened to e.g. 0.38083 with 20+ digits | 30% |
| Phase 2a yields modest tightening (e.g. 0.380869, 8 digits) | 35% |
| Phase 2a stuck at Together's local min (no improvement) | 30% |
| Phase 2b unblocks Phase 2a stuck case | 10% (conditional on 2a stuck) |

Expected value of Phase 2a: substantial. Worth pursuing.

## 6. Compare to PRO-11 (SDPA-GMP serializer)

| Metric | PRO-26 Phase 2a | PRO-11 serializer |
|---|---|---|
| Effort | 1-2 sessions | 1-2 sessions |
| Risk | h's bang-bang gradient issues | SDPA-S format quirks |
| Output | High-precision μ_UB | High-precision μ_LB (via Ω*) |
| Side affected | UB | LB |
| Synergy | If both succeed, bracket halves on BOTH sides | Same |

**Conclusion:** PRO-11 and PRO-26 are **complementary, not competing**. Different bracket sides. Both worth doing in sequence — PRO-26 Phase 2a next (faster path to a publishable result), then PRO-11 (more robust verification machinery).

## 7. Phase 1 deliverables (this document)

- Decoded Rechnitzer's full pipeline (4 ingredients: White decomposition, ansatz, ball-arithmetic evaluation, Hölder LB)
- Structural comparison ν₂² vs μ
- UB-side transfer: 🟢 feasible with caveats
- LB-side transfer: ❌ blocked (Hölder is L²-specific)
- Phase 2 plan with success-probability estimates

## 8. Recommendation

**Proceed with Phase 2a** (UB-side ansatz refinement). If μ_UB tightens below 0.3808 with ≥10 verified digits, this is a publishable result on its own. If it gets stuck at Together's level, document the negative result and refocus on PRO-11.

**Do not abandon PRO-11** — needed for high-precision LB.

**Updated strategic order:**
1. 🥇 **PRO-26 Phase 2a** — UB-side ansatz prototype (new top priority)
2. 🥈 PRO-11 serializer
3. 🥉 PRO-5 preprint v2 with all this session's findings
4. PRO-8 White email v3
