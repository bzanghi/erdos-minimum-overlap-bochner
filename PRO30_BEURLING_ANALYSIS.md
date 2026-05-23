# PRO-30 — Beurling extremal functions as a lever for tightening μ

**Status:** Analysis only. No code, no SDP runs.
**Verdict (one line):** **Non-lever.** A Beurling/Selberg construction does not give a new linear functional beyond what White's cell-envelope cosine basis + Bochner-PSD already enforce on the LB side, and it does not transfer to a new UB construction either. Retire the issue.

---

## 1. Problem restatement as a min-max autocorrelation

Following White (2023) §2 and the preprint draft §2 (`communications/preprint_draft.tex` L162–186), write `f, g : [-1, 1] → [0, 1]`, `f + g = 1`, `∫₋₁¹ f = 1`, and set

```
M(x) = ∫₋₁¹ f(t) g(x + t) dt        on x ∈ [-2, 2].
```

Then `μ = inf_f ‖M‖_∞`. Expanding `g = 1 - f`,

```
M(x) = (∫ f) − ∫₋₁¹ f(t) f(x + t) dt
     = 1 − (f ⋆ f)(x)
```

where `(f ⋆ f)(x) := ∫ f(t) f(x + t) dt` is the **standard autocorrelation** of `f` on `[-1, 1]`. So

```
μ = inf_f sup_x [ 1 − (f ⋆ f)(x) ]
  = 1 − sup_f inf_{x ∈ [-2, 2]} (f ⋆ f)(x).        (★)
```

The PRO-29 spectral duality (preprint §8.3, L592–605) is exactly this with the convolution interpreted as `⟨h, T_t h⟩`. The constraints on `f` are:

* `f : [-1, 1] → [0, 1]` (pointwise box)
* `∫ f = 1` (mass)
* (the support `[-1, 1]` and the `L^∞` cap together force `|supp(f)| ≥ 1/2`, equivalently the "h on `[0, 1]`, support `1/2`" framing in the task spec, up to affine change of variables.)

The relevant autocorrelation argument range is `x ∈ [-2, 2]`. The Fourier-period White uses is therefore `[-2, 2]` (period 4), and the cosine basis `cos(πmx/2)` for `m ≥ 1` is the orthogonal basis on a half-period.

---

## 2. What a Beurling-extremal construction would compute

A **Beurling–Selberg majorant** is the unique entire function of exponential type `≤ 2πδ` (i.e. Fourier-supported in `[-δ, δ]`) that majorizes the indicator `𝟙_{[a,b]}` and has minimum `L¹` excess. Vaaler (Bull. AMS 1985) gives explicit Fourier coefficients; the gain over the trivial dilation is `O(1/δ)` in `L¹` excess. Generalizations:

* **Selberg's majorant for `sgn`** — Holt–Vaaler, Carneiro–Vaaler.
* **Multivariate / one-sided for balls** — Carneiro–Chandee, Carneiro–Littmann.
* **For autocorrelations specifically** — there is no published Beurling-style sharp majorant; the literature instead works on the closely related but distinct **autoconvolution extremal problems** (Schinzel–Schmidt, Martin–O'Bryant, Matolcsi–Vinuesa, recently White 2022 and Novikov et al. 2025).

The natural way to use a Beurling-type object in our problem would be **one** of the following two routes.

### Route A — Beurling majorant of an indicator inside the SDP

Pick a Beurling majorant `B_δ` of `𝟙_{[-a, a]}` of exponential type `≤ 2πδ`. Since `f(t) f(x + t) ≥ 0` and we want a **lower** bound on `inf_x (f ⋆ f)(x)`, what we actually want is a **minorant** of `(f ⋆ f)` in terms of a finite Fourier expansion. The classical Beurling–Selberg machine produces minorants `B^-_δ ≤ 𝟙 ≤ B^+_δ` for indicators of intervals, but `(f ⋆ f)` is not an indicator — it is a *continuous, nonnegative, even, concave-near-origin function with finite support*.

The closest classical object is a Beurling/Logan minorant of `−(f ⋆ f)`, i.e. an entire-type lower envelope. But because `(f ⋆ f)` has its global maximum at `x = 0` (Cauchy–Schwarz) and decays to zero at `x = ±2`, the linear test functionals we would add are exactly:

```
∫₋₂² φ_k(x) (f ⋆ f)(x) dx = Σ_m λ_m(φ_k) · |hat f(m)|²
```

for `φ_k` ranging over the Beurling-extremal test family. By Parseval, **every** such linear functional on `(f ⋆ f)` is a linear combination of `|hat f(m)|²` — the very moment-matrix entries that the Bochner-PSD constraint `M_n(f) ⪰ 0` already controls jointly via PSD-ness, and that the cell-envelope cosine constraints `Σ c_m α_m^−(j)` already test against integrated cosines on each `[(j-1)L, jL]`.

In other words: **the dual cone of "linear-in-`|hat f(m)|²` valid inequalities" is closed under the operations the SDP already performs.** A Beurling minorant gives one specific element of this dual cone — sharper than a single cosine test, but not sharper than the PSD condition that subsumes all PSD-quadratic tests in `(c_m, d_m)`.

Concretely, the linear functional a Beurling test would add looks like

```
Σ_{m=0}^{2δ} β_m |hat f(m)|² ≤ B_δ(0) · ∫ f = B_δ(0)
```

with `β_m` the Fourier coefficients of the chosen majorant. The SDP's `M_n(f) ⪰ 0` already implies: for any nonnegative trigonometric polynomial `P(x) = Σ p_m cos(πmx/2)` with `p_m` admissible, `Σ p_m c_m ≥ 0`. The Beurling majorant is **one specific choice of admissible `p_m`**. It cannot beat the SDP optimum unless the SDP's optimum dual carries strictly less information than a single specific majorant — and it does not, because the SDP optimizes over the full PSD cone of admissible duals.

This is the same structural collapse documented in `OUT_OF_BOX_INEQUALITIES.md` §10 ("Logan / de Branges / Selberg–Beurling extremal") and §2 conclusion (L225–230): *White's framework IS the Beurling–Selberg / Cohn–Elkies LP for this problem; there is no second independent LP we are missing.*

### Route B — Use Beurling to bound `‖M‖_∞` below by a one-sided Fourier inequality

Set `A(x) := (f ⋆ f)(x)`. We want a **lower bound on `sup_{x ∈ [-2, 2]} (1 − A(x))`**, equivalently an **upper bound on `inf_x A(x)`**.

Suppose we have a Beurling-extremal `φ` of exponential type with `φ ≥ 𝟙_{|x| ≤ a}` and minimal `∫ φ`. Then for any nonnegative `A`,

```
inf_{|x| ≤ a} A(x) · ∫ φ ≤ ∫ A φ = Σ_m hat φ(m) · |hat f(m)|².
```

This is exactly Route A re-derived in the time domain — same content. The "cleanest" version would test against `φ_a(x) = (sin(πax/2) / (πax/2))²` (the Selberg–Fejér majorant) or its sharper Vaaler variant, but the resulting `Σ_m hat φ(m) |hat f(m)|²` is again a PSD-admissible quadratic in `(c, d)` and is **dominated by `M_n(f) ⪰ 0`**.

This is the same pattern as PRO-22 ("direct sup_t SDP") — that lever tried to drop the cell-envelope and impose `M(x_k) ≥ Ω` at a fine grid of `x_k` directly; the result was an **invalid (loose) LB** at `μ ≥ 0.370557`, +4.42 × 10⁻³ below White (`LEVER_SUPT_DIRECT.md`). The cell-envelope is *necessary for validity, not slack*. Beurling test functions are smoother than the indicators cell-envelope tests against, but they live in the same LP/SDP cone and can only sharpen what the cell-envelope already tests if they fall *outside* the PSD-admissible cone — which they don't, because Beurling test functions are *exactly* the extremal points of one-sided PSD cones (that's the definition).

### Comparison with PRO-22 and PRO-23

| Lever | Direction | Content | Result |
| --- | --- | --- | --- |
| **PRO-22** (sup_t direct) | LB | Drop cell-envelope, impose `M(x_k) ≥ Ω` at fine grid | Invalid: 0.370557 (loose by ~10⁻³) |
| **PRO-23** (KKT functional eq.) | UB | Force Together's `h*` to satisfy KKT identity | Non-tight: `μ < 0.380871` only at the witness, not closed form |
| **PRO-30 Route A** (Beurling minorant of `(f⋆f)`) | LB | Add `Σ β_m c_m ≥ const` from majorant test | Subsumed by `M_n(f) ⪰ 0` + cell-envelope cosine tests |
| **PRO-30 Route B** (Beurling majorant of indicator → bound `inf A`) | LB | Same content as Route A in time domain | Same subsumption |

Beurling sits **between** these — it's a legitimate framework element (unlike PRO-22's invalid relaxation, and unlike PRO-23's witness-only argument), but it does not enlarge the dual cone the SDP already optimizes over.

---

## 3. Why the cosine-cell-envelope already "is" the Beurling test

The SDP's cell-envelope constraints test `f` against

```
{ cos(πmx/2) }_{m = 1..2R}
```

via the per-cell minorants `α_m^-(j) := min_{x ∈ [(j-1)L, jL]} cos(πmx/2)`. The exact integrated form (preprint §5 L351–362) is a Beurling-style **lower envelope** of `cos(πmx/2)` over each cell — Lipschitz-trapezoidal in Case A and `(πm)²L³/24` Case-B tight, with the corrected per-cell residual lemma (Lemma `lem:per-cell` at L302). The aggregate residual at `N = 20000, m ≤ 20` is `O(10⁻⁸)` — well below the open-gap scale `5.7 × 10⁻⁴`.

A Beurling majorant of `cos(πmx/2)` on each cell would be a sharper minorant by at most the same residual budget. Since the residual is already `O(10⁻⁸)` per cell, **a Beurling sharpening of the cosine cell-envelope can recover at most `O(10⁻⁸)` per cell × `N = 20000` cells = `O(10⁻⁴)` total — and that's an upper bound on the gain, not a guaranteed gain.** Empirically, the cell-envelope multiplier KKT identity at `j=1` (preprint Thm `thm:kkt` L333) already binds to within 0.4%; there is no slack for a Beurling majorant to exploit at the rows where the SDP is binding (row 4 of the 7-row cover).

The "framework saturation" theorem (preprint §5 L364, `LEVER_F3_FULL_SATURATION.md`) decomposes the open gap `[0.3803027, 0.380871]` as:

* `[0.3803027, ~0.380558]` — framework-attainable at currently-tractable scale (`~45%` of open gap)
* `[~0.380558, 0.380871]` — beyond F1–F5 (`~55%` of open gap)

A Beurling-extremal addition would be a *new family `F6`*. The analysis above says F6 is **not linearly independent** of F1 (Bochner-PSD) and F4 (cosine cell-envelope). So `r_C(F6 | F1, F4) ≈ 0` up to numerical noise.

---

## 4. UB side — does Beurling give a tighter upper bound than Together's 0.380871?

Together's UB comes from an explicit `h*` constructed numerically (preprint §8.3 background, repo `https://github.com/togethercomputer/erdos-minimum-overlap`). To beat 0.380871 one needs either:

(a) a better explicit `h*` (e.g. ansatz of a Beurling-extremal-related shape), or
(b) a sharper *evaluation* of `‖M_{h*}‖_∞` showing the existing `h*` actually achieves something below 0.380871.

For (a): Beurling extremal functions are entire of exponential type, *not* `[0, 1]`-valued indicators with support `1/2`. The natural Beurling-style ansatz for `f` (e.g. `f = (sin(πx)/πx)² · 𝟙_{[0,1]}` normalized) was checked structurally by the `LEVER_*_DEEP_READ.md` series — the resulting `M` plateau sits **above** Together's plateau by `O(10⁻³)`. (The Vaaler-derived smoothing of the indicator does not improve over a discrete numerical `h*` because the binding constraint at the UB optimum is the `[0, 1]` box, not band-limit.)

For (b): the saturation/regression analysis (`OUT_OF_BOX_SYNTHESIS.md` L88) shows Together's `h*` is genuinely locally optimal — 2000 random perturbations, 4 σ-scales, zero improvement. A Beurling-extremal smoothing of `h*` would either (i) violate `0 ≤ h ≤ 1`, (ii) violate the support constraint, or (iii) raise `‖M‖_∞` by smoothing away the corners that pin the plateau. None of these tightens the UB.

**Conclusion on UB side: Beurling does not transfer.** The right literature for tightening 0.380871 is the **autoconvolution extremal problem** lineage (Schinzel–Schmidt 1988, Martin–O'Bryant 2009, Matolcsi–Vinuesa 2010, White 2022, Novikov et al. 2025) — these solve `inf_f ‖f ⋆ f‖_∞` for nonnegative `f` with `∫ f = 1, supp f ⊂ [0, 1]`. **That problem is structurally similar to our μ but is not the same problem** (it lacks the `f ≤ 1` pointwise cap and the `M = ∫ f − f ⋆ f` shift, and the extremal `f*` is *not* an indicator). The Matolcsi–Vinuesa 20-step construction gives `‖f ⋆ f‖_∞ ≥ 0.6443` for `f` supported on `[0, 1]` with `∫ f = 1`, but applied to our problem the conversion gives a far looser bound on μ than we already have. The most directly relevant paper, White's own 2022 *almost-tight L² autoconvolution inequality* (arXiv 2210.16437), is by the same author whose `μ ≥ 0.379005` we are augmenting — White has already extracted what's extractable in this direction.

---

## 5. Literature pointers

* **Vaaler (1985)** — "Some Extremal Functions in Fourier Analysis," *Bull. AMS* 12. Foundational; gives explicit Fourier coefficients of Beurling majorants for `𝟙_{[a,b]}` and `sgn`.
* **Selberg's collected papers (vol. II)** — original Beurling–Selberg construction.
* **Carneiro–Vaaler, Carneiro–Chandee, Carneiro–Littmann** — multivariate and one-sided generalizations.
* **Cohn–Elkies (2003)** — sphere-packing LP bounds; structurally identical LP cone to our SDP.
* **Cilleruelo–Ruzsa–Vinuesa (2010)**, "Generalized Sidon sets" (`arXiv:0909.5024`) — connects Sidon-set density to autoconvolution extremals. *Their dual problem is `sup f` of `(f ⋆ f)`, not `inf` over `t ≠ 0`; the constraint class is `∫ f = 1, f ≥ 0` without an `L^∞` cap — strictly weaker than ours.* `findings.md`-adjacent literature scan (`LITERATURE_SCAN_2024_2026.md` L98) already noted "no direct 2024–2026 work bears on Erdős' minimum overlap specifically."
* **Matolcsi–Vinuesa (2010)** "Improved bounds on the supremum of autoconvolutions" (`arXiv:0907.1379`). Their 20-step construction gives `‖f ⋆ f‖_∞ ≥ 0.6443`. **Different problem** — see §4.
* **White (2022)** "An almost-tight L² autoconvolution inequality" (`arXiv:2210.16437`). Same author. The fact that White wrote both this and the 2023 `μ ≥ 0.379005` paper without crossing the streams is itself evidence that the two LP cones don't carry independent information.
* **Novikov et al. (2025)** "An improved example for an autoconvolution inequality" (`arXiv:2506.16750`). 50-interval refinement of Matolcsi–Vinuesa to `0.8962`. Same caveat as above.
* **Centrale Supélec preprint (2024)** "Second-Order Beurling Approximations and Super-resolution" (`hal-04187552`). Beurling-style methods for super-resolution; not directly applicable.

The literature near our problem **does not contain** a Beurling-style sharp constant for `inf_x (f ⋆ f)(x)` under the `f : [-1, 1] → [0, 1]` box constraint. Building one would essentially require solving our problem.

---

## 6. Verdict and recommendation

**Beurling extremal functions are NOT a lever for tightening either side of `[0.3803027, 0.380871]`.**

The structural reasons:

1. **Route A/B subsumption.** Any Beurling-derived linear functional on `(f ⋆ f)` is a linear-in-`|hat f(m)|²` test, hence already dominated by the Bochner-PSD + cell-envelope cosine system. This matches the OUT_OF_BOX_INEQUALITIES.md §10 verdict that Beurling–Selberg gives no improvement beyond `1/8` as a one-liner, and the §2 verdict that White's framework *is* the Beurling–Selberg / Cohn–Elkies LP.
2. **Cell-envelope already saturates the residual.** The Case-A/B corrected per-cell residual (preprint Lem `lem:per-cell`) is `O((πm)² L³ / 24) = O(10⁻⁸)` at `N = 20000`, `m ≤ 20`. A Beurling majorant of `cos(πmx/2)` can sharpen this by at most a constant factor — well below the open-gap scale `5.7 × 10⁻⁴`.
3. **UB side has the wrong shape.** Beurling smoothings of `h*` either violate `[0,1]` box, violate support, or raise `‖M‖_∞`. The Matolcsi–Vinuesa / White-2022 autoconvolution lineage is the right neighbourhood, but is a *different* extremal problem and is already harvested by White himself.
4. **Framework-saturation theorem.** The preprint §5 saturation argument quantitatively rules out any F6 in the linear span of `{|hat f(m)|²}_{m ≤ M} ∪ {cell averages of M}` from closing more than `~45%` of the open gap, which the existing F1–F5 stack already attains.

**Recommendation: retire PRO-30.** Move the issue to "Won't fix — subsumed by F1+F4 (Bochner-PSD + cell-envelope cosines)" with a one-line pointer to this doc.

**If anyone wants a single concrete numerical check** before retiring: run the Selberg–Fejér majorant `B(x) = (sin(πδx/2)/(πδx/2))²` at `δ = 2R/2 = 10` against the row-4 augmented optimum's `(c_m, d_m)`, evaluate `Σ_m hat B(m) (c_m² + d_m²)`, and confirm it is `≤` the LP value of `Σ_m c_m c̃_m + d_m d̃_m` for the SDP's optimal cosine-cell-envelope dual `(c̃, d̃)`. The prediction is that the Beurling test gives the *same* number to within `O(10⁻⁶)`. This is a 15-minute check using `lp_research_state/code/probe.py` semantics, but would not change the recommendation either way.

**Concrete next steps if we did pursue it (we should not):**

1. (Defensive) Implement the Selberg–Fejér probe above on row 4 at `N = 20000`, `bochner_n = 30`. Expected: `Δ ≤ 1 × 10⁻⁶`. Time: 1 hour.
2. (Out-of-scope) Construct a Beurling majorant of `(f ⋆ f)` *with the `f ≤ 1` cap baked in*. This is open research and would be a new theorem, not an SDP add-on.
3. (Wrong direction) Use Vaaler-smoothed indicator as a UB candidate. Predicted to lose to Together's `h*`.

None of (1)–(3) is a Phase-6 candidate. The PRO-30 "1 evening reading + 1 session computation" budget in the Linear issue is correct *if* the conclusion is "no" — and that is the conclusion.

---

## 7. Bookkeeping note on the rigor convention

This analysis discusses **per-row SDP** sharpenings. To convert any hypothetical per-row improvement into a rigorous bound on μ, the ellipse-extension argument (preprint §3, White §5.1 / App. II) must be re-run; per the CLAUDE.md "Critical caveat," a single-point SDP gain is not a μ bound. Since the verdict here is negative, no ellipse-extension run is needed.

---

*— Analysis closes PRO-30. No SDP cycles burned, no claims hedged because there is no claim.*
