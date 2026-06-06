# PRO-32: B-S/M-R Autocorrelation Bounds Don't Transfer to μ

**Status:** Done. **❌ No clean transfer.** Three independent structural obstructions identified. The numerical paradox (B-S/M-R's ~0.41 vs our μ ≈ 0.38) is resolved.

## 1. The setup

Per PRO-29's duality:
```
μ = inf_h sup_t M(h, t) = 1 - sup_h inf_t ⟨h, T_t h⟩
```
where μ_dual = `sup_h inf_t ⟨h, T_t h⟩` is the **max-over-h of the minimum autocorrelation**.

Barnard-Steinerberger / Madrid-Ramos prove:
- B-S (arXiv:1903.08731), Theorem 2: `min_{t∈[0,1]} ∫ f(x)f(x+t) dx ≤ 0.411 · ‖f‖_L¹²`, sharp constant in (0.37, 0.411]
- Fish-King-Miller (arXiv:2001.02326): perturbation conditions on extremizers; same setup
- Madrid-Ramos (arXiv:2003.06962): Theorem 1.3 — `min_{t∈[-1/2,1/2]} ∫f(x)f(x+t)dx ≤ 0.829604 · ‖f‖_L¹·‖f‖_L²` (mixed norm; different scale); Theorem 1.7 — strict-but-non-quantitative improvement on B-S's 0.411 for compactly-supported nonneg L¹

## 2. Setup comparison

| | B-S/FKM | M-R | Our μ |
|---|---|---|---|
| Function class | nonneg `f ∈ L¹(ℝ)` | same | `h:[0,2]→[0,1]` |
| Support | unrestricted | unrestricted | compact in `[0, 2]` |
| Upper bound | **none** | **none** | **`h ≤ 1` pointwise** |
| t-range | `[0, 1]` | `[-1/2, 1/2]` | `t ∈ ℝ` (effectively `[-2, 2]`) |
| L¹ norm | `‖f‖_L¹` arbitrary | same | `‖h‖_L¹ = 1` (fixed) |
| Best published constant | 0.411 (B-S) | 0.829604 (mixed norm); strict improvement on 0.411 (no number) | — |

## 3. Substitution attempt

Tried `g(x) = (1/c) h(x/c)` to map our problem onto theirs:
- `‖g‖_L¹ = ‖h‖_L¹ = 1` ✓
- `∫g(x)g(x+t')dx = (1/c) ∫h(y)h(y+t'/c)dy`
- B-S applied: `min_{t'∈[0,1]} (1/c)∫h(y)h(y+t'/c)dy ≤ 0.411`
- Substituting `s = t'/c`: `min_{s∈[0,1/c]}∫h·h(·+s)dy ≤ 0.411·c`

To cover our `s ∈ [0, 2]` we need `1/c ≥ 2`, so `c ≤ 1/2`. Then `0.411·c ≤ 0.2055`.

This would give `inf_t ⟨h, T_t h⟩ ≤ 0.205`, so `μ ≥ 1 − 0.205 = 0.795`. **Contradicts** our verified `μ ≤ 0.380871`.

## 4. The bug in the substitution

When `g(x) = (1/c) h(x/c)` with `c = 1/2`:
- Support of `g` is `[0, 2c] = [0, 1]` ✓ (compact)
- But `g(x) = 2·h(2x)` reaches `2·1 = 2` ≫ 1

**g violates the `g ≤ 1` constraint by a factor of 2.** B-S/FKM/M-R don't impose this constraint; their extremizers are SPIKE-LIKE (large concentrated bumps). Our h is bounded above by 1 — a fundamentally different function class.

## 5. Three obstructions ranked by severity

1. **(decisive) Missing L∞ bound.** Our `h ≤ 1` is not scale-invariant. Any substitution that fits their `t ∈ [0,1]` into our `t ∈ [-2, 2]` crashes against the ceiling.

2. **(directional) t-range goes the wrong way.** B-S's `min over t ∈ [0,1]` is a min over a **smaller** set than our `inf over t ∈ [-2, 2]`. Shrinking the min-set gives a *larger* min (`min over subset ≥ min over superset`), which is **not** an upper bound on the larger-set min.

3. **(structural) Support constraint.** Their f lives on all of ℝ; ours is compactly supported in `[0, 2]`. Even fixing this via specialization `f = h · 1_{[0,2]}` only gives the t-range `[0,1]` issue (obstruction 2).

## 6. Conclusion

❌ **B-S, FKM, and M-R results do not directly transfer to our μ problem.** The bracket `μ ∈ [0.3803027, 0.380871]` remains state-of-the-art; no immediate cross-bound is available from this literature.

**Resolved paradox:** the apparent claim `μ ≥ 0.5929` from naive substitution was wrong because B-S's 0.411 constant applies to a strictly weaker function class (no L∞ bound), and the t-range direction is incompatible.

## 7. What COULD still apply

Madrid-Ramos's class `𝓛_g(I)` (nonneg L¹ compactly supported, dominated by g outside I) is the closest formal analogue. Their Theorem 1.7 proves a strict-but-non-quantitative improvement on 0.411 for this class. **Quantifying their improvement under the additional `h ≤ 1` constraint** would require re-deriving the constant in their proof — this is a research problem, not a simple substitution.

Re-deriving would be a substantial mathematical undertaking (~weeks of math research), and **the proof technique might not give a tighter bound than our existing SDP work**. Their bound 0.411 corresponds to μ ≥ 0.589, which is much weaker than our μ ≥ 0.3803027. **Even if we could make their approach work, it would likely give a weaker LB.**

## 8. Strategic implications

- **Stop pursuing B-S/M-R transfer for new bounds.** Document as resolved-negative.
- The autocorrelation duality `μ = 1 − μ_dual` is still **useful for preprint exposition** — it positions our work in the autocorrelation literature even if no constant transfers.
- **Cite Tao's `optimizationproblems` C_{1b}** as the authoritative name for our constant.

## 9. Deliverables

- This document
- (No code; pure math analysis)
