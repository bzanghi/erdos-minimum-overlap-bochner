# Lever H Deep Read — arXiv:2210.16437 transferability to min-overlap

**Date:** 2026-05-10
**Author:** investigative session, post-CDE Phase 5 (µ ≥ 0.3801279)
**Hypothesis under test:** White's "almost-tight" uniqueness/eigenvalue technique
for L² autoconvolution (arXiv:2210.16437) transfers to the Erdős minimum-overlap
problem and closes (or significantly tightens) the rigorous lower bound on `µ`.
**Sources consulted:** arXiv:2210.16437 (HTML/ar5iv rendering — PDF binary fetch
failed, ar5iv HTML succeeded); arXiv:2508.02803 (follow-up LB construction);
arXiv:2506.16750 (follow-up LB construction); `lp_research_state/code/white_full_convex.py`
(our min-overlap encoding).

---

## 1. TL;DR

**Verdict: TRANSFER FAILS.** The technique in arXiv:2210.16437 is *not* an
eigenvalue / second-variation analysis as the cross-domain note hypothesized.
It is a Fourier-truncation convex program (same genus as White's Acta Arith.
program that our SDP already implements), plus a uniqueness argument that
relies critically on the **strict convexity of the L² (= L⁴-of-Fourier-coefficients)
norm**. The min-overlap functional `sup_t (f⋆f)(t)` is an L∞ norm in `t`, which
is **not strictly convex**, so the uniqueness machinery does not transfer. The
"almost-tight" 0.0014% bound is a *numerical* convergence-rate statement
(`R⁻¹/⁶`) of the same Fourier-truncated convex program — analogous to our
own SDP's convergence to its continuum value, with the same intrinsic
truncation error. No new analytical bound, no compact-operator eigenvalue
extraction, no smooth-proxy passage that survives `p → ∞`. **Expected gain
on `µ`: 0.** Recommendation: do not pursue this lever; redirect to the
non-sharpness / duality-gap diagnostic (top-2a in OUT_OF_BOX_CROSS_DOMAIN).

---

## 2. White's L² autoconvolution technique (arXiv:2210.16437)

**Theorem 1.** `0.574635728 ≤ µ₂² ≤ 0.574643711`, where
`µ₂ = inf{‖f*f‖₂ : f : [-½, ½] → ℝ, ∫f = 1}`. The ratio of upper to lower bound
is `~1.0000139`, i.e. the value is pinned to **0.0014%**.

**Proof skeleton (sections 2–4).**

- **§2 (qualitative).** Direct method in calculus of variations: weak-*
  compactness in `L¹ × L∞` gives existence of a minimizer `f^◇` (Prop. 4).
  Uniqueness is then obtained from **strict convexity of the L² norm of
  `f*f`** — equivalently, strict convexity of `‖f̂‖₄⁴` via Plancherel. The
  proof is a Minkowski-inequality equality-case argument: two minimizers
  must be linearly dependent, and the normalization `∫f = ∫g = 1` then forces
  `f = g`. Lemma 5 says one can drop the `L^∞` upper bound without changing
  the infimum. Lemma 6 says `f^◇ * f^◇ * f^◇` is constant on `(-½, ½)` (an
  Euler–Lagrange / Lagrange-multiplier consequence).

- **§3 (Plancherel identity, Lemma 7).** The KEY identity is the closed-form
  expansion
  ```
  ‖f*f‖₂² = ½ + Σ_{m≥1} f̂(m)⁴ + (16/π⁴) Σ_{m odd} (rational in f̂(k))⁴
  ```
  This converts the L² objective into a **polynomial in Fourier coefficients**.
  Lemma 8 gives Hölder-type lower bounds used to turn a numerical optimum
  into a certified lower bound.

- **§4 (quantitative, Prop. 9).** Define a finite-Fourier-truncated convex
  program `𝒪(T, R)` with `T` Fourier modes and an `R`-parameter discretization
  of the constraint. This program is a **quadratically-constrained linear program**
  (per the paper) — it is *not* an SDP and does *not* use Bochner-PSD or
  positive-trigonometric-polynomial constraints. Prop. 9 proves the rate
  `|𝒪(R, √R) − µ₂²| < 10·R⁻¹/⁶`. The "almost-tight" 0.0014% comes from
  picking `R` large enough to make `10·R⁻¹/⁶` small in absolute terms.

**Where uniqueness comes from (the critical point).** Strict convexity of
the L² norm (equivalently the L⁴-of-Fourier-coefficients norm). NOT from
eigenvalue analysis. NOT from a second-variation operator. NOT from a
spectral gap. The proof never constructs a compact operator and never
diagonalizes anything. Our cross-domain note's hypothesis about "principal
eigenfunction of a compact operator" was wrong.

---

## 3. The transfer question — explicit comparison

**Functional 1 (White's L²).** `F₂(f) = ‖f*f‖₂² = ∫(f*f)(t)² dt`. By
Plancherel, `F₂(f) = ∫|f̂(y)|⁴ dy` — a smooth, strictly convex,
**polynomial-in-Fourier-coefficients** functional (Lemma 7). Second
variation: well-defined and positive-definite at the minimizer (this is
what uniqueness via Minkowski actually encodes).

**Functional 2 (min-overlap).** `F_∞(h) = sup_t (h⋆h)(t)`. This is an L∞
functional in `t`. Its second variation does not exist in the usual sense —
`F_∞` is *not* Fréchet-differentiable; its subdifferential at the optimum
is a *convex set of measures* supported on the (possibly multi-point)
argmax set `{t : (h⋆h)(t) = F_∞(h)}`. This is a "sup-type" extremal problem,
not a smooth one.

**The first natural attempt: replace `sup` by `L^p` and let `p → ∞`.** Define
`F_p(h) = (∫|(h⋆h)(t)|^p dt)^{1/p}` for `p < ∞`. For each finite `p`, `F_p`
is smooth and strictly convex (for `p > 1`), and has a White-style Fourier
expansion. For `p = 2`, this *is* the L² autoconvolution functional, just
on the autocorrelation instead of the autoconvolution — and that yields
the L² autoconvolution constant `µ₂`, not `µ`.

The problem: **`F_p → F_∞` non-uniformly in `f`**, and the rate is governed
by the *width* `w(h) = meas{t : (h⋆h)(t) ≥ F_∞(h) − ε}` of the argmax set,
which depends on `h`. The triangle inequality gives, on the support of
`(h⋆h)` of total mass `M`:
```
F_p(h) ≤ M^{1/p} · F_∞(h),       F_p(h) ≥ w(h)^{1/p} · (F_∞(h) − ε).
```
So `F_p` agrees with `F_∞` only up to a multiplicative `M^{1/p}` slack. For
`p = 100`, `M^{1/p} ≈ M^{0.01} ≈ 1 + 0.01·log M`. With `M ≈ 4` (support of
the autocorrelation of an indicator on `[0,2]`), this is `~1.014` — i.e.
1.4% slack. To get within `10⁻⁴` of `µ` we would need `p ≥ 14,000` or so,
at which point the polynomial-in-Fourier expansion of `F_p` has degree
`2p ≥ 28,000` in `ĥ` and is utterly intractable.

Equivalently: White's Lemma 7 (the `f̂(k)⁴` expansion) is the special case
`p = 2`, which already taxes our SDP scale. The expansion at `p = 100`
generates `ĥ(k)^{200}` terms; Lasserre relaxations of such constraints have
already been shown (in our own `lasserre.py` retraction) to require tail
bounds that kill the gain at currently-tractable scale.

**The second natural attempt: directly mimic uniqueness from strict convexity
of `F_∞`.** The L∞ norm in `t` is **not strictly convex**. Specifically,
`F_∞(½h + ½g) = ½ F_∞(h) + ½ F_∞(g)` whenever `h, g` achieve their suprema
at *distinct* points (the typical case). Minkowski's equality case for
`L∞` requires `h⋆h` and `g⋆g` to be *proportional on the joint argmax
set*, a vastly weaker condition than linear dependence. The Minkowski
argument from Prop. 4 collapses.

This is not a fixable obstruction — it is intrinsic to the L∞ functional.

**The third natural attempt: use Lemma 6's Euler–Lagrange constancy
property.** White's Lemma 6: `f^◇ * f^◇ * f^◇` is constant on `(-½, ½)`.
The analogous statement for min-overlap would be: at the optimum `h^*`,
the *autocorrelation* `(h^* ⋆ h^*)(t)` is constant on the *argmax set*,
which is trivially true (it's literally the level set of the sup). The
analogue gives NO additional information — it is the obvious first-order
condition, not a strong identity. White's identity has analytical content
because the L² problem has a *smooth, polynomial-quartic* Lagrangian; the
L∞ problem's Lagrangian is a measure supported on a level set.

**Where does this leave the analogy?** The L² convex program of §4 is
**already the analogue our SDP implements** for the min-overlap functional
(via Bochner-PSD relaxation of `sup_t (h⋆h)(t) ≤ c` into trigonometric-
polynomial-positivity constraints). White's quantitative result is
`|𝒪(R, √R) − µ₂²| < 10·R⁻¹/⁶`; ours is the empirical 5.4×10⁻⁴ ceiling
that the diagnostic chain has already mapped. The "almost-tight" is a
matter of `R` being large, not of any technique beyond Fourier truncation.

---

## 4. Verdict — TRANSFER FAILS

**Mathematical obstruction.** The hypothesized "second-variation eigenvalue
technique" is not what arXiv:2210.16437 actually does. The actual technique
is:
  - Fourier-truncated convex program (we already have this — it's our SDP).
  - Strict convexity of the L² (= L⁴-of-Fourier) norm for uniqueness.

Min-overlap functional `sup_t (h⋆h)(t)` is L∞ in `t`, which is **not
strictly convex**, killing the uniqueness argument. Smoothing to `F_p`
recovers strict convexity but at exponentially-growing polynomial degree
in the Fourier expansion (`ĥ^{2p}`), pushing the SDP/LP into a regime
already proven intractable by our own retracted Lasserre work. There is
no "free" smoothing that gives a tighter bound on `µ` than what the
present SDP yields.

**Quantitative estimate of expected gain.** Strictly zero, modulo
re-deriving via the smoothed `F_p` what our SDP already computes for
`p = ∞`. White's own paper notes (paraphrased): "the L∞ method of Cloninger
& Steinerberger is computational, limited by a nonconvex optimization
program" — i.e. White himself flags L∞ as the harder regime and does not
claim the L² technique transfers. The 2025 follow-ups (arXiv:2508.02803,
arXiv:2506.16750) are *constructive lower bounds* via step functions —
same flavor as Together's min-overlap UB. They do not transfer to the
analytical-LB side.

**Sanity check from our own retraction.** `lasserre.py` was withdrawn
precisely because truncating `(f²)̂(m)` without a tail bound is not
rigorous, and the natural Fejér–Riesz tail bound kills the gain at
tractable scale (`communications/lasserre_tail_bound.md`). The L² method
of White avoids this only because its native objective is a polynomial
of degree 4 in `f̂(k)` — finitely many terms with explicit, computable
tails. Min-overlap's native objective is `sup_t`, which has no polynomial
expansion at all in the Fourier basis without the same kind of relaxation
the SDP already employs.

---

## 5. Recommended action

1. **Do not pursue this lever further.** Redirect the email to Ethan White
   away from "does the L² technique transfer" and toward two more
   productive questions:
   - *(a)* Does White believe the Bochner-PSD + ellipse-extension family
     has an intrinsic duality gap to `µ`? (his §5.1 in the Acta Arith.
     paper is the most likely source of insight on this).
   - *(b)* Is there a Lemma-7-style closed-form Fourier identity for
     `sup_t (h⋆h)(t)` via the *test-function-against-`δ_t`* dual that
     our SDP doesn't currently encode? (Unlikely, but White is the
     person who'd know.)

2. **Pivot to Lever 2a (Cohn–de Laat–Salmon non-sharpness).** This is now
   the highest-EV remaining lever: rather than push the LB *up*, prove
   the rigorous *ceiling* of the Bochner+ellipse family is strictly
   below `µ` (or, possibly, equal to `µ` — both outcomes are publishable).
   Concrete next experiment: adapt the discrete-reduction dual construction
   of arXiv:2206.09876 to our setting and check whether the dual SDP
   admits a witness exceeding `0.3801279` at the SDP scale we already
   solve at.

3. **Document the Lever H closure in the research note.** Update
   `erdos_lower_bound_research_note.md` (or `findings.md`) with a
   one-line entry: *"Lever H (transfer of arXiv:2210.16437 L² uniqueness
   technique to min-overlap) fails — L∞-not-strictly-convex obstruction.
   See LEVER_H_DEEP_READ.md."*

4. **(Optional) Direct experimental falsification.** Implement the `F_p`
   smoothed-proxy SDP for `p = 4, 8, 16` and verify empirically that
   the LB it yields *degrades* relative to the SDP's `sup_t` encoding,
   not improves. Expected outcome: smoothing strictly weakens the LB.
   If the empirical result contradicts this, revisit.

---

## Sources

- White, *Canad. Math. Bull.*, arXiv:2210.16437 (L² autoconvolution).
  Theorem 1, Propositions 4 & 9, Lemmas 5–8. Accessed via
  https://ar5iv.labs.arxiv.org/html/2210.16437 ; PDF binary fetch failed
  on arxiv.org/pdf/2210.16437 (saved-to-disk-only, no text).
- arXiv:2508.02803 (2025), arXiv:2506.16750 (2025) — constructive LB via
  step functions; classified as upper-bound-analogue work (LB on a
  *ratio*, not on `µ` itself).
- `lp_research_state/code/white_full_convex.py` — our SDP encoding, same
  Fourier-truncation genus as White's `𝒪(T, R)`.
- `OUT_OF_BOX_CROSS_DOMAIN.md` — prior hypothesis under test (now falsified).
- `communications/lasserre_tail_bound.md` — sanity check on the
  "polynomial-in-Fourier" route's intrinsic limits at our scale.
