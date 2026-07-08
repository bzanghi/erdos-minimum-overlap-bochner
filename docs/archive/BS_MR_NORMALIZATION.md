# BS / MR Autocorrelation Inequalities — Normalization Map to µ (PRO-32)

**Question.** Do the Barnard–Steinerberger (BS, arXiv:1903.08731) and
Madrid–Ramos (MR, arXiv:2003.06962) extremal *autocorrelation* inequalities
transfer to a useful bound on the Erdős minimum-overlap constant µ?

**Verdict: NO — they do not transfer.** The naive map gives the spurious
"µ ≥ 0.59" precisely because it (i) confuses the cross-correlation
`f ⋆ (1−f)` with the autocorrelation `f ⋆ f`, and (ii) replaces the Erdős
*sup over all real shifts* with BS's *inf over a unit window* `t ∈ [0,1]`.
Both substitutions are individually fatal; together they manufacture a
contradiction. There is no clean substitution (`x → x/L`, support length,
`‖f‖₁` rescaling) that repairs them, because the two functionals optimize
genuinely different objects. **No number from BS/MR bears on the bracket
`[0.379544, 0.380871]`.** Details and numerics below.

---

## 1. Exact theorem statements (verbatim from the papers)

### Barnard–Steinerberger, *Three Convolution Inequalities on the Real Line with Connections to Additive Combinatorics* (arXiv:1903.08731)

Abstract gives three inequalities. The relevant one is the **minimum
autocorrelation** bound:

> For `f ∈ L¹(ℝ)`,
> `min_{0 ≤ t ≤ 1} ∫_ℝ f(x) f(x+t) dx ≤ 0.42 ‖f‖_{L¹}²`,
> where the constant 1/2 is trivial, **0.42 cannot be replaced by 0.37**.

- **Type:** AUTOcorrelation — `f` against `f` (not `1−f`).
- **Support:** none — `f ∈ L¹(ℝ)`, full real line.
- **t-range:** the unit window `t ∈ [0, 1]` (dimensionless).
- **Normalization:** RHS is `‖f‖₁²` (scale-invariant: both sides scale
  like `‖f‖₁²` under amplitude scaling, but the LHS is NOT scale-invariant
  under the *spatial* dilation `f(x) → f(x/λ)` — see §3).
- "0.42 cannot be replaced by 0.37" means: there exists `f ≥ 0` whose
  `min_{t∈[0,1]} (f⋆f)(t)` is `≥ 0.37 ‖f‖₁²`, so the optimal constant
  `C₄ ∈ [0.37, 0.42]`. The 0.37 example is compactly supported (§4.2 of MR).

### Madrid–Ramos, *On Optimal Autocorrelation Inequalities on the Real Line* (arXiv:2003.06962)

MR sharpen BS. Defining `C₄` as the smallest constant with
`min_{t∈[0,1]} ∫ f(x)f(t+x) dx ≤ C₄ ‖f‖₁²` (Thm 1.7):

- BS's explicit upper constant is `C₄ ≤ 1/(2(1+θ₀))` with
  `θ₀ = −inf_x sin(x)/x = 0.217…`, i.e. `C₄ ≤ 0.41075…`.
- **MR Theorem 1.7:** `C₄` is *strictly smaller* than `1/(2(1+θ₀))`
  (qualitative — they explicitly state they "cannot quantify" the
  improvement; *"a new method altogether is needed"*).
- **MR §4.2 / Conjecture 4.2:** the lower bound `C₄ ≥ 0.37` comes from the
  compactly supported Chebyshev-weighted example
  `f(x) = 1_{[−1/4,1/4]}/√(1−4x²) − 1_{[−1/2,1/2]}/(4√(1−4x²))`,
  for which `(f⋆f)(t) ≥ 4/π` on `t∈[0,1]` and `‖f‖₁ ≤ 1.439`. MR's
  Conjecture 4.2 ties the *compactly-supported* version to
  `f ∈ L¹([−1/2,1/2])` (**support length 1**), and Prop. 4.3 relates it to
  `f ∈ L¹([−1,1])` with a factor 4.
- MR also improve the *mean* autocorrelation bound (Thm 1.1: `0.91 → 0.864`)
  and the *min-over-`[−1/2,1/2]`* mixed-norm bound (Thm 1.3: `0.8296`).
  None of these is the functional in our problem either.

So the BS/MR object is: **the largest possible MINIMUM of the
AUTOcorrelation over a UNIT window, for a nonnegative L¹ function**,
`C₄ = sup_{f≥0} inf_{t∈[0,1]} (f⋆f)(t) / ‖f‖₁² ∈ [0.37, 0.41075]`.

---

## 2. Our functional, and the naive (broken) reduction

Erdős/White/Haugland (see CLAUDE.md, `together_loader.py`):

> `µ = inf_f sup_{x∈ℝ} ∫ f(t)(1 − f(x+t)) dt`,  `f:[−1,1]→[0,1]`, `∫f = 1`,
> f zero-extended.

Write the **cross-correlation** `M(x) := ∫ f(t)f(x+t) dt = (f ⋆ f)(x)`
(here it coincides with the autocorrelation because the cross term is
`f` vs `f`; the `1−f` only enters through the constant). The naive
algebra:

```
∫ f(t)(1 − f(x+t)) dt = ∫ f(t) dt − ∫ f(t) f(x+t) dt = 1 − A(x),
   where A(x) = (f ⋆ f)(x) is the autocorrelation, ‖f‖₁ = 1.
⇒ sup_x [1 − A(x)] = 1 − inf_x A(x)
⇒ µ = inf_f sup_x[1−A(x)] = 1 − sup_f inf_x A(x).
```

BS says `sup_{f} inf_{t∈[0,1]} A(t) = C₄ ≤ 0.41`, so the naive reader
concludes `sup_f inf_x A ≤ C₄` and `µ ≥ 1 − C₄ ≈ 0.59`. **This contradicts
`µ ≤ 0.380871` and is wrong.** The two errors:

### Error A — sup is over ALL real shifts, not a unit window
The Erdős sup ranges over **all** `x ∈ ℝ`. As `|x| → 2` (support length
of `f⋆f`), `A(x) → 0`, so `inf_{x∈ℝ} A(x) = 0` and the naive identity gives
`sup_x (1 − A(x)) = 1`, a vacuous bound. BS's inf is restricted to the
**unit window** `t ∈ [0,1]`, where `A` is bounded *away from 0*. These are
different operations on the same `A`: `inf_{[0,1]} A ≫ inf_{ℝ} A = 0`.

### Error B — the Erdős sup is NOT `1 − inf A`; it picks an interior shift
The honest Erdős/Haugland functional is the **cross-correlation of `f`
with `1−f`**, `(f ⋆ (1−f))(x) = ∫ f(t)(1−f(x+t)) dt`. Its argmax is at an
*interior* shift, NOT at the tail where `A` vanishes. So `sup_x` of the
overlap is decided by where `f` and `1−f` overlap maximally, a balance
condition that has nothing to do with `inf A`. (`A` minimized = `f`
shifted off itself; the overlap maximized = `f` and the *complement*
`1−f` aligned.) The identity `µ = 1 − sup_f inf_x A` is FALSE; the correct
statement is `µ = inf_f sup_x (f ⋆ (1−f))(x)`.

In short: BS/MR study `sup_f inf_t (f ⋆ f)`; Erdős is
`inf_f sup_x (f ⋆ (1−f))`. Same letters, different problem.

---

## 3. Normalization algebra — is there any clean map?

Attempted substitution to align the settings:

| Feature | BS/MR | Erdős/White |
|---|---|---|
| object | autocorr `f⋆f` | cross-corr `f⋆(1−f)` |
| extremal type | `sup_f inf_t` (max-min) | `inf_f sup_x` (min-max) |
| shift range | unit window `t∈[0,1]` | all real `x∈ℝ` |
| support of `f` | `ℝ` (or `[−1/2,1/2]`, len 1) | `[−1,1]`, len 2 |
| pointwise bound | `f ≥ 0` only | `0 ≤ f ≤ 1` (box) |
| normalization | `‖f‖₁²` on RHS | `∫f = 1` |

- **Spatial dilation.** Map our `f` on `[−1,1]` to support length 1 via
  `g(x) = 2 f(2x)` on `[−1/2,1/2]` keeps `∫g = 1`. Under `x→x/λ`,
  `(f⋆f)(t)` rescales but the *window* `[0,1]` does not, so a dilation that
  fixes the support breaks the window and vice-versa. No single `λ` aligns
  both the support AND the t-window with BS's normalization.
- **The `0 ≤ f ≤ 1` box constraint** is absent in BS/MR. Their extremal
  `C₄`-functions are *unbounded* (Chebyshev singularities `1/√(1−4x²)`).
  Our `f` is hard-capped at 1, so even the autocorrelation `sup_f inf A`
  would be a *different* (smaller) constant in our class.
- **Direction.** Even granting Error A/B away, BS gives an UPPER bound on
  `inf_t A`, which through `µ = 1 − sup inf A` would give a LOWER bound on
  µ. But that very identity is false (Error B). MR's *strict* improvement
  is qualitative (no number), so even the direction carries no usable
  numeric.

**Conclusion of the algebra: no clean substitution exists.** The mismatch is
structural (max-min vs min-max, auto vs cross, window vs full line, box vs
unbounded), not a units bug.

---

## 4. Numerical sanity check (`lp_research_state/code/_bsmr_check.py`)

Run on Together's near-optimal `f*` (`together_f_star.json`, n=600 step
function on `[0,2]`, `∫f = 1`, `‖f‖₁ = 1`) and on test densities.

```
M_T = sup_k ∫ h(x)(1−h(x+k)) dx   = 0.380870   (Erdős/Together UB on µ)
A(0) = ∫ h²                        = 0.774675
inf_t A(t) over t ∈ [0,1] (BS)     = 0.173075
inf_t A(t) over ALL shifts         = 0.000000   (→ 0 at the tail |t|→2)

naive identity  sup_k(1 − A(k))    = 1.000000   attained at lag ≈ −1.997 (tail!)
1 − inf_{t∈[0,1]} A                = 0.826925
TRUE cross-correlation argmax      at physical shift ≈ 0.11, value 0.380870
```

Readings:

1. The identity `overlap(k) = 1 − A(k)` holds **exactly** pointwise (it's
   just `∫f − ∫f·f(·+k)`). But `sup_k (1 − A(k)) = 1`, attained at the
   **tail** `k ≈ −2` where `A → 0` — this is Error A made concrete.
   The actual Erdős sup (0.38087) is attained at an **interior** shift
   `k ≈ 0.11` by the cross-correlation `f ⋆ (1−f)` — Error B made concrete.
2. `inf_{t∈[0,1]} A = 0.173` on this near-optimal `f*` is already **well
   below** BS's `C₄ ≥ 0.37` envelope, confirming `f*` is nowhere near a
   BS extremizer — the two problems pull `f` in different directions.
3. Plugging `inf_{[0,1]}A = 0.173` into the (false) identity gives
   `1 − 0.173 = 0.827`, nowhere near µ ≈ 0.38 — direct proof the map is
   broken.
4. Uniform control: `h ≡ 1/2` on `[0,2]` gives `M_T = 0.5` and triangular
   `A` (peak 0.5 at lag 0, → 0 at lag 2), consistent with the above.

The BS `C₄ ≥ 0.37` reference example, recomputed, is a **signed** function
(`f⋆f` taken in absolute value in MR Conj. 4.2); it is *not* admissible in
our `0 ≤ f ≤ 1` class at all, reinforcing that the extremal regimes don't
overlap. (Our reproduction yields a negative-valued `f`, so its raw
`min f⋆f` is negative; MR use `|f⋆f|`. Either way it lies outside our
feasible set.)

**The algebra and the numbers agree: the substitution is not self-consistent
for our functional.** A mismatch would have shown up as `1 − inf_{[0,1]}A ≈ µ`;
instead it is `0.827 ≠ 0.381`.

---

## 5. What would have to change for any transfer

A BS/MR-style bound could touch µ only if one had an extremal theory for the
**min-max cross-correlation** `inf_f sup_x (f ⋆ (1−f))` over the *box class*
`0 ≤ f ≤ 1, ∫f = 1, supp f ⊆ [−1,1]`, with the sup over the **full** real
line. That is exactly White's program — i.e. the thing we are already
solving. BS/MR's max-min *autocorrelation* over unbounded `L¹` functions is a
sibling extremal problem with the same Fourier-analytic flavor but a
different objective, feasible set, and extremizer. Their constants
(0.37–0.42, 0.864, 0.8296, etc.) do not bound µ in either direction.

**Bottom line for the bracket `[0.379544, 0.380871]`: BS/MR contribute
nothing.** This is the expected clean negative. The earlier "done" mark on
PRO-32 with no committed artifact is now superseded by this note.

---

### Provenance
- BS abstract fetched from arXiv:1903.08731 (verbatim, §1 above).
- MR text fetched from arXiv:2003.06962v2 (PDF parsed; Thm 1.1, 1.3, 1.7,
  §4.2, Conj. 4.2, Prop. 4.3 quoted/paraphrased above).
- Numerics: `lp_research_state/code/_bsmr_check.py` on
  `lp_research_state/data/together_f_star.json` via
  `together_loader.compute_overlap_from_f` and direct `np.correlate`.
