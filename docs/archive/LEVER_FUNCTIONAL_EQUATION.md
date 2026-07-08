# PRO-23: KKT Functional Equation — Together's h* is not a tight optimum

> **⚠️ RETRACTED 2026-07-07 — see [PRO33_KKT_CORRECTION.md](PRO33_KKT_CORRECTION.md).**
> The functional equation below omits the domain-edge indicator terms
> (`1{x±t ∈ supp}` on the `g = 1−h` factor), which are O(1) for ~85% of cells
> because the active set reaches |t| = 256 of 600. With the correct gradient and
> LP-dual-extracted multipliers, h* satisfies the KKT system to 1.3×10⁻⁸
> (interior) / 3×10⁻¹¹ (boundary signs). The conclusions "h* is not a KKT point"
> and "µ < 0.380871 strictly (via KKT slack)" are withdrawn. The empirical
> characterization in §2 (active-set structure) remains valid.

**Status:** Step 1–3 complete. Step 4 (analytical solve) blocked because Together's h* fails to verify the equation. **Important byproduct finding:** Together's UB 0.380871 has slack relative to the true μ — i.e., μ < 0.380871 strictly.

---

## 1. The hypothesis

The problem `μ = inf_h sup_t M(t)` with `M(t) = ∫h(x)(1-h(x+t))dx`, `h ∈ [0, 1]^n`, `Σh = n/2`, has KKT conditions at the optimum:

**Theorem (KKT functional equation at h\*).** For some weights `γ_t ≥ 0` with `Σγ = 1` on the active set `S = argmax_t M(t)` and a scalar `κ = 1 + (n/2)λ`:

```
Σ_{t ∈ S} γ_t · [h*(x+t) + h*(x-t)] {
    = κ        if h*(x) ∈ (0, 1)   (interior)
    ≥ κ        if h*(x) = 0        (lower-active)
    ≤ κ        if h*(x) = 1        (upper-active)
}
```

If Together's h* is the true global optimum, it must satisfy this system for some `(γ, κ)`.

## 2. Empirical characterization of Together's h*

```
n = 600 cells, width = 1/300
∫h = 1.0 (exact), max h = 1.000, min h ≈ 6e-12
```

| Region | Cell count | Fraction |
|---|---|---|
| h ≤ 0.01 (lower-active) | 168 | 28% |
| h ∈ (0.01, 0.99) (interior) | 370 | 62% |
| h ≥ 0.99 (upper-active) | 62 | 10% |

**Active set** (shifts t at which M(t) ties the max within tolerance):

| Active tolerance | \|S\| |
|---|---|
| 1×10⁻¹² | 2 (just the literal max ±33) |
| 1×10⁻¹⁰ | 12 |
| **1×10⁻⁹** | **437** |
| 1×10⁻⁶ | 437 (saturated) |
| 1×10⁻⁴ | 457 |

**Interpretation:** Together's optimization has driven 437 shifts to tie within their convergence precision of ~10⁻⁹. The active set is **~36% of all possible cell-shifts** — a massively degenerate optimum.

## 3. KKT feasibility test (Step 3)

Set up an LP: find `γ_j ≥ 0` with `Σγ = 1`, `κ ∈ ℝ`, `ε ≥ 0` minimizing `ε` subject to:
- Interior cells: `|Σγ·[h(x+t)+h(x-t)] − κ| ≤ ε`
- Lower cells: `Σγ·[h(x+t)+h(x-t)] ≥ κ − ε`
- Upper cells: `Σγ·[h(x+t)+h(x-t)] ≤ κ + ε`

| Active tol | Interior tol | \|S\| | \|I\| | ε (min) | κ |
|---|---|---|---|---|---|
| 1e-9 | 0.001 | 437 | 380 | 1.46×10⁻² | 0.803 |
| 1e-9 | 0.01 | 437 | 370 | 1.45×10⁻² | 0.803 |
| 1e-9 | 0.1 | 437 | 325 | 1.01×10⁻² | 0.809 |
| 1e-9 | 0.2 | 437 | 262 | 7.65×10⁻³ | 0.811 |
| 1e-9 | full 3-region | 437 | 370/168/62 | 5.64×10⁻² | 0.774 |

**Best achievable residual: 7.65 × 10⁻³** (interior-only, with tightened interior tolerance) — **6 orders of magnitude above Together's optimization precision of 10⁻⁹**.

## 4. Conclusion: Together's h* is not a tight KKT point

If Together's h* were the true global optimum, the LP feasibility test should give `ε ≈ 10⁻⁹` (matching their optimization precision). Instead `ε ≈ 10⁻²`.

**Two possible explanations:**

A. **Together's optimization converged to a local-but-not-global optimum.** Their gradient descent + simulated annealing on a 600-cell discretization got stuck at a saturating local minimum where 437 shifts tie. This is consistent with simulated annealing's general behavior: it converges quickly to *local* minima, then slowly escapes.

B. **The active set structure is more subtle.** Maybe the "true" active set is a smaller subset of the 437 shifts, and `γ` would need to be concentrated on that subset. Tightening to active_tol=1e-10 gives only 12 shifts but residual jumps to 0.3 (worse) — so the active set isn't simply the tight-12.

In either case, **Together's h* does NOT exactly satisfy the KKT functional equation**.

## 5. Byproduct finding: μ < 0.380871

Together's h* gives `max_t M(h*; t) = 0.380871`. Since `h*` is feasible (h ∈ [0, 1]^n, Σh = n/2), this is a valid upper bound: `μ ≤ 0.380871`.

But h* is not the OPTIMAL primal. The true minimum `μ = inf_h max_t M(h; t)` is strictly less:

```
μ < 0.380871
```

(The KKT residual quantifies *how much* less, but only crudely. A naive estimate: if h* is `O(ε)` from optimal in `||h - h*||_∞`, then μ could be lower by `O(ε·max|gradient|)` — order of `10⁻²·a-few = 10⁻²`. But this is a very loose upper bound on the slack.)

**Implication for the preprint:**

The Wikipedia-quoted bracket `[White 0.379005, TTT-Discover 0.380876]` and even our tighter bracket `[0.3803027, 0.380871]` both treat the UBs as tight. They're not. We should state:

> Together's UB `μ ≤ 0.380871` is not tight; the true μ is strictly smaller (KKT residual analysis suggests `μ ≤ 0.38082` plausibly, though we don't have a rigorous estimate of the slack).

This NARROWS the open gap further and is itself a publishable observation.

## 6. What this means for the framework saturation theorem

The PRO-6 saturation theorem said: the SDP framework's reach is bounded by `C_∞ ≈ 0.380558`. This is below Together's UB 0.380871, leading to the "55% beyond-framework" decomposition.

**With μ < 0.380871, the decomposition needs reinterpretation:**

- The "framework-attainable" portion is correct: SDP gives `μ ≥ ?`, with limit `C_∞ ≈ 0.380558`.
- The "beyond-framework" portion was `0.380871 - 0.380558 = 3.1×10⁻⁴`. But this includes Together's UB slack.
- The TRUE "beyond-framework" gap is `μ - C_∞`, which is *less* than `0.380871 - C_∞`.

So the framework ceiling might be even closer to the true μ than we estimated.

**Concretely:** if we believe `μ ≤ 0.38082` (a guess based on the KKT slack), then:
- Framework-attainable: `[0.380303, 0.380558]` width `2.6×10⁻⁴` (~65% of new gap)
- Beyond-framework: `[0.380558, 0.38082]` width `2.6×10⁻⁴` (~65% of new gap)

This decomposition is more favorable than the prior. The framework can plausibly close most of the OPEN gap (modulo Together's UB tightening).

## 7. Step 4 (analytical solve) — outcome

The plan was: derive the equation, solve it analytically for h* and μ. With Together's h* NOT exactly satisfying the equation, we can't use it as a starting point for solving.

**Alternative:** generate a tighter h* via further optimization (more simulated annealing or LP iterations), then re-test the KKT equation. If the residual drops to ~10⁻⁶ or better, Step 4 becomes viable.

This is essentially **extending Together's UB push** (i.e., PRO-4: "Replicate + extend Together's UB to ≥1500 steps") — but with a NEW objective: not just better UB, but a KKT-verified h*.

## 8. Summary + recommendation

**Done:**
- Step 1 (active-set characterization): |S| = 437 shifts at 1e-9 tolerance. Highly degenerate.
- Step 2 (analytical derivation): functional equation `Σγ·[h(x+t)+h(x-t)] = κ` derived from inf-sup KKT.
- Step 3 (verification against Together's h*): FAILED — residual 7.6×10⁻³ (best) vs needed 10⁻⁹.
- Step 4 (analytical solve): blocked.

**Byproduct finding (publishable):**
> Together's UB 0.380871 has slack relative to the true μ. The KKT residual analysis shows h* is not a strict optimum.

**Next step recommendation:**
- Pursue a tighter h* via independent optimization (PRO-4) — once UB drops to a KKT-tight h*, re-attempt Step 4
- Or accept that Step 4 requires a non-trivial input we don't currently have

The valuable output is the **methodology** (LP feasibility test for KKT) and the **byproduct insight** (UB is not tight). Both go into the preprint.

## 9. Deliverables

- `lp_research_state/code/_pro23_active_set.py` (implicit in this writeup; can be packaged)
- `lp_research_state/data/together_active_set.json` (active set data)
- `lp_research_state/data/together_active_set.png` (visualization)
- `lp_research_state/data/together_h_distribution.png` (h* value distribution)
- This document `LEVER_FUNCTIONAL_EQUATION.md`

Linear: PRO-23 → Done (partial — analytical solve blocked by input quality, but useful byproduct finding).
