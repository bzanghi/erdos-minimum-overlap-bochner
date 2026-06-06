# PRO-27: Lean/Mathlib Lemma Inventory for the Erdős μ Problem

**Status:** Done. Mathlib has the basic Fourier/Plancherel/Hölder primitives we'd need, but **lacks the specialized infrastructure** for our problem class — no autocorrelation extremal lemmas, no Beurling-Selberg majorants, no SDP duality, no Bochner-Herglotz PSD-trigonometric-moment theorem, no Lasserre hierarchy. Formalizing PRO-6's saturation theorem requires substantial new development; **Rechnitzer's UB pipeline (PRO-26 Phase 2a) is a more tractable formalization target**.

## Top-line verdict

**Tractable now (~weeks):** basic Fourier setup on `AddCircle 2`, Parseval, Hölder, Cauchy-Schwarz, convolution-Fourier duality, Riemann-Lebesgue.

**Buildable on Mathlib (~3-6 months):** cell-envelope `α_m^-(j) = cos(πmL/2)` lemma, per-cell Taylor bound `δ_m(j) = O((πm)²L³/24)`. These are direct Taylor-expansion calculations in scope of `Mathlib.Analysis.Calculus.Taylor`.

**Not buildable today (12+ months gap):**
- SDP primal-dual / KKT framework — Mathlib has zero
- PSD moment matrices `M_n(f) ⪰ 0` Bochner-style — `Matrix.PosSemidef` exists but the Toeplitz Fourier-coefficient PSD theorem doesn't
- Bochner-Herglotz / trigonometric moment problem — absent
- Beurling-Selberg / Vaaler extremal majorants — absent
- Lasserre / SOS certificates — absent

**Math-side gaps (math not formalization):** the uniform shadow-price bounds `|ξ| ≤ Ξ, τ ≤ T, ν_3 ≤ V` are unproven conjectures (LEVER_I_PRIME_THEOREM.md §2.2). Formalization can't get ahead of the math.

## 18 candidate Mathlib lemmas

| # | Mathlib name | Application | Applies? |
|---|---|---|---|
| 1 | `MeasureTheory.Lp.fourierTransformₗᵢ` | L² isometry, central for Plancherel chain | ✅ |
| 2 | `MeasureTheory.Lp.norm_fourier_eq` | `‖𝓕f‖₂ = ‖f‖₂` baseline | ✅ |
| 3 | `MeasureTheory.Lp.inner_fourier_eq` | Parseval `⟨f,g⟩ = ⟨f̂,ĝ⟩` | ✅ |
| 4 | `hasSum_sq_fourierCoeff` | Parseval on AddCircle — our natural domain | ✅ **most useful** |
| 5 | `orthonormal_fourier` | Foundational orthonormality | ✅ |
| 6 | `fourierBasis` | Hilbert basis from Fourier monomials | ✅ |
| 7 | `hasSum_fourier_series_L2` | L² convergence of Fourier series | ✅ |
| 8 | `Real.fourier_mul_convolution_eq` | `𝓕(f∗g) = 𝓕f · 𝓕g` — converts White's `‖F∗F‖²₂` to `Σ|F̂|⁴` | ✅ **central** |
| 9 | `MeasureTheory.convolution_assoc` | Needed for `F∗F∗F` in Rechnitzer §4 | ✅ |
| 10 | `MeasureTheory.Integrable.fourierInv_fourier_eq` | Fourier inversion | Maybe |
| 11 | `Real.tsum_eq_tsum_fourier_of_rpow_decay` | Poisson summation | Maybe |
| 12 | `ENNReal.lintegral_mul_le_Lp_mul_Lq` | Hölder — but blocks at math level for min-max μ | ✅ but limited |
| 13 | `ENNReal.lintegral_Lp_add_le` | Minkowski | ✅ (trivial) |
| 14 | `inner_mul_le_norm_mul_norm` | Cauchy-Schwarz | ✅ |
| 15 | `orthonormal.tsum_inner_products_le` | Bessel inequality, for tail bounds on Fourier coefs | ✅ |
| 16 | `Real.tendsto_integral_cocompact_fourierIntegral` | Riemann-Lebesgue qualitative truncation | ✅ |
| 17 | `Finset.addEnergy` | Discrete `Σ\|1̂_A\|⁴` energy — no Fourier identity proven | 🟡 partial |
| 18 | `SchwartzMap.fourier_convolution` | Smooth-truncation tail args | Maybe |

## Strongest hit: `MeasureTheory.Lp.fourierTransformₗᵢ` + `Real.fourier_mul_convolution_eq`

Sufficient to formalize White's decomposition:
```
‖F∗F‖²₂ = ‖𝓕(F∗F)‖²₂ = ‖(𝓕F)²‖²₂ = Σ|F̂(k)|⁴
```
This is **mechanically reachable from current Mathlib**. The wrinkle: `fourierTransformₗᵢ` is on ℝⁿ (continuous Fourier transform), while our `h : [0, 2] → [0, 1]` lives on `AddCircle 2` (discrete-index Fourier coefficients). Translation via `Real.fourierCoeff_tsum_comp_add` (Poisson summation, #11).

## What's missing for our actual saturation theorem

PRO-6's theorem requires:
1. ❌ Convex programming framework (Mathlib has zero SDP/primal-dual machinery)
2. ❌ Bochner-style "PSD iff finite Hilbert representation" on Toeplitz matrices
3. ❌ Beurling-Selberg / Vaaler extremal majorants
4. ❌ Lasserre hierarchy / SOS certificates

Estimated effort to bring even one of these to publish-quality formalization: 12+ months per item. The Singer-Sidon arxiv preprint (2605.03274) confirms this — they explicitly list "cyclic Fourier API: Parseval, convolution, L⁴ energy" as future work.

## Recommendation

**If we pursue Lean formalization (PRO-7):**

✅ **Target: Rechnitzer's UB pipeline** (1-2 month effort, ~1-2k lines). Bessel-function bounds and ball-arithmetic-style interval bounds are present in Mathlib (`Mathlib.Analysis.SpecialFunctions.Bessel`, `Mathlib.Topology.Algebra.Order.IntermediateValue`).

❌ **Avoid: PRO-6's saturation theorem.** The SDP/PSD/Bochner infrastructure gaps make this 12+ months of foundational work before the specific theorem becomes addressable.

✅ **Standalone mini-result feasible:** the cell-envelope Taylor lemmas `δ_m(j) ≤ (πm)²L³/24` (LEVER_I_PRIME_THEOREM.md §3.1). Pure calculus, fully in scope. Roughly 1-2k lines.

## Sources

- [Mathlib.Analysis.Fourier.FourierTransform](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/Fourier/FourierTransform.html)
- [Mathlib.Analysis.Fourier.AddCircle](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/Fourier/AddCircle.html)
- [Mathlib.Analysis.Convolution](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Analysis/Convolution.html)
- [Mathlib.MeasureTheory.Integral.MeanInequalities](https://leanprover-community.github.io/mathlib4_docs/Mathlib/MeasureTheory/Integral/MeanInequalities.html)
- [Mathlib.Combinatorics.Additive.Energy](https://leanprover-community.github.io/mathlib4_docs/Mathlib/Combinatorics/Additive/Energy.html)
- [arXiv:2605.03274 — Formalizing Singer Sidon Constructions in Lean 4](https://arxiv.org/pdf/2605.03274) (confirms Fourier API gap)
- [Floris van Doorn — BonnAnalysis Plancherel blueprint](https://florisvandoorn.com/BonnAnalysis/blueprint/chap-plancherel.html)
