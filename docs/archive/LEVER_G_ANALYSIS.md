# Lever G analysis — (f, g) variables with f + g = 1

**Date:** 2026-05-10
**Question:** Does rewriting White's SDP in Together's `(f, g)` variables (with `f + g = 1` enforced as a linear constraint) admit a strictly tighter convex relaxation than White's f-only formulation?

**Verdict (one line, defended in §4):** **Provably equivalent at the convex-hull level for the *encoding* of the constraint set `{(f, 1−f) : 0 ≤ f ≤ 1}`. But the two formulations minimize over *different objective functionals*; the (f, g)/Haugland objective `M_T` admits a different — and in principle differently relaxable — SDP family than White's `Ω`. Lever G is therefore not a no-op, but its potential gain is a separate, well-defined research question on the M-side, not a structural fix to White's f-side. Recommendation: do not prototype as drop-in tightening; instead, treat as an open M-side encoding question — concrete prototype spec in §5.**

---

## 1. The two formulations side by side

### 1.1 Together / Haugland (the "(f, g)" formulation)

From Together's `README.md`:

> Let `C_5` be the largest constant satisfying
> `sup_{x ∈ [−2, 2]} ∫_{−1}^{1} f(t) g(x + t) dt ≥ C_5`
> for all non-negative `f, g: [−1, 1] → [0, 1]` with `f + g = 1` on `[−1, 1]` and `∫_R f = 1`,
> where `f, g` are extended by zero outside `[−1, 1]`.

Substituting `g = 1 − f` (forced by the linear constraint `f + g = 1` on the support), this is **identical to Haugland's (2016) formulation**:

```
C_5  =  inf_h  sup_{k ∈ R}  ∫_R  h(x) (1 − h(x + k)) dx,
        h : [0, 2] → [0, 1],  ∫_0^2 h = 1,  h zero-extended off [0, 2].
```

(Together's README cites Haugland for this equivalence verbatim.) The translation between Together's "f on [−1, 1]" and Haugland's "h on [0, 2]" is a unit translation `h(x) = f(x − 1)`. Both have **support length 2** and **shift range [−2, 2]**.

Call this the **M-functional**:

```
  M_T(h) := sup_{k ∈ R}  ∫_R  h(x) (1 − h(x + k)) dx.       (*)
```

### 1.2 White (the "f-only" formulation)

From `white_full_convex.py` and `lp_research_state/code/_fourier_convention_notes.md`:

- `f : [−2, 2] → [0, 1]` is the density variable, periodically extended with period 4.
- `∫_{−2}^{2} f = 1` (line 142).
- White's objective: `Ω = sup_t M_W(t)` where

```
  M_W(t) := ∫_{−2}^{2} f(x) f(x + t) dx       (autocorrelation, period-4 on [−2, 2]).
```

White's SDP minimizes `Ω` subject to a Fourier-analytic upper bound on `sup_t M_W(t)`, the moment / parameter constraints `(h_1, h_2, p, q)` from the §5.1 ellipse covering, Bochner PSD on `f` and `1 − f`, and optional augmentations.

Call this the **Ω-functional**.

### 1.3 The domain mismatch is real and load-bearing

White uses `f` on `[−2, 2]` (length 4), Haugland/Together use `h` on `[0, 2]` (length 2). The bridge identity (White 2023, §3, and also Haugland 2016):

```
  µ  =  inf_f sup_t M_W(t)  =  inf_h sup_k M_T(h)(k).
```

Both `inf`s equal the same constant `µ`. But at *fixed* admissible inputs the two functionals give *different numbers* — the diagnostic in `_fourier_convention_notes.md` §7.2 measures this directly: at Together's `h*`,

```
  M_T(h*)             = 0.380871   (Together's published value)
  Ω(f_even(h*))       = 0.387337   (White's autocorrelation, even embedding)
```

A gap of `~7 × 10⁻³`. This is not a discretization artifact; it is the **structural gap between two related-but-not-equal objective functionals**.

---

## 2. Relating `M_T(h)` to `Ω_W(f)`

This is the algebraic core of the question. Two relations to keep straight.

### 2.1 The naïve "(1 − f)" expansion

Suppose we work on Haugland's domain `[0, 2]` directly (i.e. `f : [0, 2] → [0, 1]`, zero-extended off `[0, 2]`). Then

```
  ∫_R f(t) (1 − f(t + k)) dt
    =  ∫_R f(t) dt  −  ∫_R f(t) f(t + k) dt
    =  1  −  R_f(k),
```

where `R_f(k) := ∫_R f(t) f(t + k) dt` is the autocorrelation of `f` extended by zero.

Crucially: because `f` is zero-extended, the integrand `f(t) f(t + k)` has compact support of length at most `2 − |k|` for `|k| ≤ 2` and is identically zero for `|k| > 2`. So `R_f(0) = ∫f² ≤ 1` and `R_f(k) ≤ 1 − k/2 · (something)` (just by support-length / Cauchy–Schwarz).

Therefore

```
  M_T(h)  =  sup_k  [ 1 − R_h(k) ]  =  1 − inf_k R_h(k).
```

**This is a striking and clean identity.** Together / Haugland is *minimizing the supremum of `1 − R_h(k)`* — equivalently, **maximizing the infimum of `R_h(k)`** over shifts `k`.

But: the inf is achieved at `|k|` near the boundary (typically `k → 2⁻` where the support windows barely overlap), at which point `R_h(k) → 0`. So the inf over all real `k` of `R_h(k)` is `0` (attained in the limit `|k| > 2`). That makes the supremum of `1 − R_h(k)` equal to `1`, which is obviously not Haugland's `µ ≈ 0.38`.

**Subtlety I missed in the brief.** Haugland's `M_T` is `sup_k`, not `sup_t over the whole line` of the integrand — but it's also not *just* over integer cell-shifts. Look at notebook cell 3 of Together's analysis (cited in `together_loader.py:32–34`):

> "Discretely the max is attained on the grid of integer cell shifts."

So in the **discrete** (piecewise-constant) setting the max is over `k ∈ {−(n−1)L, …, (n−1)L}`, the relevant set. Equivalently, in the continuous setting, the sup is over `k ∈ [−2, 2]` since for `|k| > 2` the supports don't overlap and the integrand vanishes identically — `∫h(1 − h(·+k)) = ∫h · 1 − 0 = 1`, which is *outside* the supremum search range because the integrand of `h(x)(1 − h(x+k))` requires `h(x+k)` to be defined where `h(x)` is also nonzero (otherwise we get the trivial value `∫h = 1`).

Concretely Haugland defines:

```
  M_T(h) = sup_k  Φ(k),  where  Φ(k) := ∫_R h(t) (1 − h(t + k)) dt.
```

For `|k| > 2`, `Φ(k) = ∫_R h(t) · 1 dt − 0 = 1`. For `|k| ≤ 2`,

```
  Φ(k) = (∫_R h) − R_h(k) = 1 − R_h(k) ∈ [1 − ∫h², 1].
```

So `sup_k Φ(k) = 1` regardless of `h` — **the literal sup is degenerate at `|k| ≥ 2`**. Haugland must mean the supremum *restricted to* `|k| ≤ 2 − support_width` or similar. Looking at the actual definition (Haugland 2016, the original combinatorial Erdős problem): the supremum is over `k` in the *interior of the overlap region*, where both `h(t)` and `h(t+k)` are "live" — equivalently, the **max** over integer cell-shifts in the discretized version.

The clean restatement: Haugland's `µ` is

```
  µ  =  inf_h  max_{k ∈ Z, |k| < n}  L · Σ_i h_i (1 − h_{i+k}),
```

where `L = 2/n`, `h ∈ [0, 1]^n`, `Σh = n/2`. The sup over real `k` would degenerate, so the discretization to integer cell-shifts (with zero-padding) is the right object — and indeed it's what Together implements. **In the continuous version, the right reading is `sup_k Φ(k)` over the open interval `(−2, 2)`, where `Φ(k) → 1` only in the limit `|k| → 2⁻`.**

This `Φ`-discretization is what `together_loader.py:431–458 compute_overlap_from_f` implements: `L · max(np.correlate(h, 1 − h, mode="full"))`.

So the correct relation:

```
  M_T(h)  =  sup_{k ∈ (−2, 2)}  [ 1 − R_h(k) ]   (zero-extended autocorrelation; h on [0, 2]).
```

`R_h(k)` is **maximized at `k = 0`** (Cauchy–Schwarz / autocorrelation peak), so `inf_{k ∈ (−2,2)} R_h(k)` is achieved away from `k = 0`. White's `Ω` is **`R_h(0) = ∫h²`** in this picture — the autocorrelation peak (or actually the supremum of autocorrelation, which is achieved at `t = 0` for nonnegative `h`).

**So: White's Ω-functional looks at the PEAK of `R_h`, while Together's M-functional looks at `1 minus` the INFIMUM of `R_h` (over a bounded shift range). These are different statistics of the same autocorrelation function.**

### 2.2 The minimum-overlap identity (Haugland)

Haugland (2016) proves that the infima coincide:

```
  inf_f sup_t R_f(t)  =  inf_h sup_k [ 1 − R_h(k) ]  =  µ.
```

This is a non-trivial theorem — it uses the duality between maximizing the autocorrelation peak (subject to a mass constraint) and minimizing the "off-peak depth" (subject to symmetric constraints). The proof is essentially that for the *optimal* `f` and the *optimal* `h`, one can construct the other and they achieve the same value.

But for **non-optimal** inputs, the two functionals disagree. The disagreement at Together's `h*` (which is M-optimal but not Ω-optimal) is the `0.38734 − 0.38087 = 0.00647` gap reported in §1.3.

---

## 3. The SDP relaxation comparison

This is the actual question Lever G poses.

### 3.1 White's Ω-functional SDP relaxation

White's program (lines 99–339 of `white_full_convex.py`) relaxes `Ω = sup_t M_W(t)` via:

1. **Fourier truncation** + **cell-Fourier inequalities** (lines 176–190). These encode

   ```
   sup_t M_W(t)  ≤  (Ω/2) · 2  =  Ω
   ```

   tested against trig polynomials of bounded degree, using cell-kernel cosine envelopes `[α^−_{jm}, α^+_{jm}]`.

2. **Bochner-PSD on `f` and `1 − f`** (lines 233–258). These encode `f ≥ 0` and `f ≤ 1` in the dual (Fourier) sense via Hermitian Toeplitz PSD constraints on the moment matrices `M_n(f̂)` and `M_n((1−f)̂)`.

3. **M-side Bochner-PSD** (`mside_bochner.py`, lines 260–292). This encodes `M_W(t) ≥ 0` for all `t` as a Hermitian Toeplitz PSD on the Fourier coefficients of `M_W`, with the bilinear term `−4|f̂(m)|²` SOC-relaxed.

The objective `Ω` is bounded from above by all of these constraints; the SDP value is a rigorous lower bound on `µ`.

### 3.2 The (hypothetical) M_T-functional SDP relaxation

A Lever-G prototype would minimize a different objective: `M_T` rather than `Ω`. With `g = 1 − f` substituted, the M-functional is

```
  M_T(f)  =  sup_k  ∫ f(t) (1 − f(t + k)) dt
        =  sup_k  [ ∫f  −  R_f(k) ]
        =  ∫f  −  inf_k R_f(k)        (over the allowed shift range).
```

Since `∫f = 1` is a hard constraint, **minimizing `M_T` is equivalent to MAXIMIZING `inf_k R_f(k)`** — i.e., maximizing the *flatness* of the autocorrelation across shifts, not minimizing its peak.

This is **fundamentally a different SDP**:

- White's SDP: `min Ω` such that `R_f(t) ≤ Ω/2` for all `t` (PSD-style upper bound on autocorrelation peak).
- Lever-G SDP: `min (1 − μ)` such that `R_f(k) ≥ μ` for all `k ∈ K` (PSD-style *lower* bound on autocorrelation infimum over restricted shift set `K`).

The relaxations are **not the same family**:

- The peak-control (White) uses Bochner-PSD on `M_W ≥ 0` PLUS a *scalar upper bound* `M̂(0) = Ω/2`. This is a *one-sided* constraint on the function `M_W`.
- The off-peak-depth control (Lever G) requires that some shifted version of `1 − M_W` (specifically, `1 − R_f(k)` for `k ∈ (−2, 2)`) has its supremum ≤ Ω. This is *also* a one-sided constraint but on a *different* derived function.

### 3.3 Are the two SDPs's convex hulls comparable?

**Claim.** The feasible set in `f`-space is the *same*: `{f : f ∈ [0, 1], ∫f = 1, f-domain = [0, 2] (or [−2, 2] after embedding)}`. Both SDPs over-approximate this set via the same Bochner / moment / Fourier-truncation machinery.

What differs is the **objective**:

```
  White:  Ω(f) = max(R_f, taken over all shifts; supported on autocorr peak)
  G   :  Ω(f) = max(1 − R_f, taken over restricted shifts; supported off-peak)
```

Each is a *convex* functional of the moment matrix (both are suprema of linear-in-moments expressions, when the autocorrelation is parametrized via its Fourier coefficients).

**Now: a SDP relaxation of `min Φ(f)` over a convex set `F` is determined by:**

1. The over-approximation `F̃ ⊇ F` of the feasible set (same in both SDPs).
2. The structure of the objective `Φ`.

Both `Ω` and `M_T` are suprema of linear-in-`f̂` functions (after the Bochner-PSD encoding makes `f̂` accessible). The standard SDP relaxation of a sup is a *single Toeplitz-PSD inequality* that majorizes the function pointwise. **For the same `F̃`, the SDP optima of `min Ω` and `min M_T` differ only insofar as the objectives differ as functions on `F̃`.**

**Both objectives are convex in `f̂`** (after Bochner encoding), and **both equal `µ` at the optimum** by Haugland's theorem. So:

```
  inf_{F̃} Ω  ≤  inf_F Ω  =  µ
  inf_{F̃} M_T  ≤  inf_F M_T  =  µ
```

The question is whether `inf_{F̃} M_T` could be **strictly larger** than `inf_{F̃} Ω` (giving a tighter LB on `µ`).

This is **not automatic**. In general, given two convex functions `φ_1, φ_2` with `inf_F φ_1 = inf_F φ_2`, their relaxed minima on `F̃ ⊋ F` can satisfy *either* direction of inequality, or coincide. **The relaxation can favor one objective over the other depending on which is "more curved" or which has a "shallower" relaxation gap on `F̃ \ F`.**

### 3.4 A specific structural argument: why the Ω-SDP might be **looser** than an M_T-SDP

White's Ω-SDP has the well-documented weakness:

- **Together's `h*` is a much better M_T-point than it is an Ω-point** (M_T(h*) = 0.3809 vs Ω(h*) = 0.3873).
- This suggests there exist `f` in White's *relaxed* feasible set `F̃` with `Ω(f) < 0.3795` (the current rigorous LB) but with `M_T(f) > 0.3795`.
- Such an `f` is **infeasible for the µ ≥ 0.3795 bound under M_T** but **feasible under Ω**. Equivalently: the M_T-objective is MORE sensitive than Ω to the slack of `F̃` over `F`.

This is suggestive but **not a proof**. The actual question is: at the *minimum* of the relaxed problem (different problems with different objectives have different optima), do we get the same value or different ones?

**Empirically (per `_fourier_convention_notes.md` §7.2):** the SDP-optimal `f̃` (White's program at row 4 with Bochner-PSD) has `M_T(f̃) ≈ 0.66`, far higher than `Ω(f̃) ≈ 0.38`. So at White's optimum, M_T is *not* a tight relaxation either — the SDP optimizer just doesn't care about M_T.

If we *re-ran* the SDP with `min M_T` instead, the new optimum would be a different `f̃'`, and we cannot a priori predict its `M_T(f̃')` from current data.

### 3.5 The decisive structural fact

There is one decisive structural fact that **bounds Lever G's potential gain**:

**Both objectives, on the same relaxed feasible set `F̃`, are minorized by the same dual lower bound.** Specifically, take the dual certificate for White's SDP that produces `Ω* ≥ 0.379544`. This certificate is built from:
- a SOS-like nonneg trig polynomial `p(t)` with `M̂_W(t) ≥ 0` pointwise,
- combined with constraints that force `(Ω/2) p(0) − Σ_m p̂(m) M̂_W(m) ≤ 0`.

The dual certificate for M_T would involve a different family of nonneg trig polynomials `q(k)` that test `R_f(k) ≥ µ` over `k ∈ K`. These are *different* polynomial cones — and **either could in principle be larger or smaller** than the other.

**However:** the Bochner-PSD cone on `f` and `1 − f` is **the same** in both SDPs, and these are the *binding* constraints in the current best SDP per the diagnostic (CLAUDE.md). The cell-kernel cosine envelopes (lines 176–190) are also the same. The objective change alone, with the *same* relaxed feasible set, cannot increase the SDP optimum *beyond the value that the dual cone can certify*.

**Concretely:** if Lever G's SDP has optimum `µ_G* > 0.379544 (White's current)`, then the dual cone of Lever G's SDP must certify this value. Since the *primal* feasible-set constraints are identical between the two formulations (same Bochner-PSD, same cell-kernel, same moments), **any tighter dual certificate Lever G enables must come from the M_T-objective's specific algebraic structure** — specifically, from a *negative-curvature* or *concavity* exploitation in the `−R_f(k)` term that the Ω-objective cannot see.

### 3.6 Lever G in the (f, g)-with-`f+g=1`-as-an-explicit-variable form

The brief's framing is "add `g` as an explicit variable with `f + g = 1` as a linear constraint." This is a **trivial reformulation** at the LP/SDP level: `g` is determined by `f` via the linear constraint, so adding `g` as a variable + the constraint adds zero information (no new feasible solutions, no new infeasible ones, no new dual variables that aren't already implicit).

**Unless:** we use `g` in the *objective* in a way that exploits structure not visible in `f` alone. The natural such use is:

```
  M_T(f, g) = sup_k ∫ f(t) g(t + k) dt
```

which, when `g = 1 − f`, becomes Haugland's M-functional. This is the only nontrivial use of the (f, g) reparameterization — and it amounts to **changing the objective from Ω to M_T**, as analyzed in §3.2–§3.5.

So Lever G as stated is equivalent to: **switch the SDP's objective from Ω to M_T, while keeping the same feasible-set encoding.**

---

## 4. Verdict and reasoning

**Verdict: could be tighter in principle, but the gain is bounded by — and aligned with — the M-side Bochner-PSD machinery already in the codebase. The (f, g)-with-`f+g=1` framing offers no new convex-hull tightening beyond what M-side Bochner achieves.**

### 4.1 Justification

1. **`g = 1 − f` adds no new variables in any informative sense.** The constraint `f + g = 1` is the *definition* of `g`. Adding `g` to the SDP variable list with this constraint is a no-op for the feasible set.

2. **The interesting object is the objective `M_T`** rather than `Ω`. Changing the objective alone, with the same Bochner-PSD / cell-kernel / moment constraints, **does not change the relaxed feasible set `F̃`** — it changes which point in `F̃` minimizes.

3. **`M_T(f) = ∫f − inf_k R_f(k) = 1 − inf_k R_f(k)`** (per §2). So `min_F̃ M_T` is equivalent to `1 − max_F̃ inf_k R_f(k)`. The inner `inf_k R_f(k)` is concave in `f̂` (it's the infimum of linear functions of `f̂` — actually wait, `R_f(k)` is *quadratic* in `f` and hence in `f̂` via `R̂_f(m) = |f̂(m)|²`). The inf of concave functions is concave, so `inf_k R_f(k)` is concave in `f̂`. **Concave-in-`f̂` maximization over a convex set is the right shape for a tractable SDP** — it's the same shape as the current M-side Bochner approach.

4. **The codebase already has the M-side Bochner-PSD machinery** (`mside_bochner.py`, `mside_bochner_schur.py`, `mside_via_lasserre.py`), which encodes `M_W ≥ 0` (where `M_W` is White's autocorrelation-derived M-function) as a relaxed Toeplitz-PSD constraint. **This is structurally the same family** as what Lever G would need for `M_T`-objective relaxation — both are "PSD-test a derived autocorrelation function." The Lasserre level-2 variant `mside_via_lasserre.py` is the tightest known relaxation; if there is a Lever-G gain, it lives in that family.

5. **The diagnostic chain has tested M-side Bochner at various orders**, and findings.md / CLAUDE.md show it does not break past 0.379544 at the currently-tractable problem size. **This implicitly bounds Lever G's gain too**: whatever tightening Lever G offers is structurally upper-bounded by what an exactly-solved M-side Bochner-PSD at high order achieves.

6. **The brief's specific claim ("the (f, g) rewrite with explicit `f + g = 1`") is a no-op at the SDP level.** The substantive question is "could a different objective (M_T instead of Ω) give a tighter SDP," and the answer to *that* is "in principle, by exactly the amount that M-side Bochner-PSD captures, which has already been tested."

### 4.2 Verdict refined

- **Provably equivalent (Lever G as literal "add `g` as a variable with `f + g = 1`"): YES.** This adds zero new information.
- **Provably equivalent (Lever G as "switch the objective to `M_T`"): NO** — these are different SDPs. But they share the same feasible-set encoding, so the gain is bounded by the M-side Bochner / Lasserre family already explored.

**Final verdict: PROVABLY EQUIVALENT** for the literal reading. The interesting reading ("change objective to `M_T`") is **not new** — it is the M-side Bochner approach, already in the codebase as `mside_bochner.py` and its Lasserre lift `mside_via_lasserre.py`. **No new prototype warranted as a Lever-G investigation.** If a research lever exists in this direction, it is *not* "use (f, g) variables" — it is "find a higher-order or different M-side Bochner encoding," which is a known research question, not a Lever-G one.

---

## 5. Prototype spec (negative recommendation)

**Do not prototype Lever G as `(f, g)`-with-`f+g=1` variables.** It is provably equivalent to White's f-only formulation for the *constraint encoding*, and the only non-trivial reading (different objective) is already covered by the M-side Bochner machinery in `mside_bochner.py`.

If future research time goes here, the right test would be:

1. Add a `cp.Variable` `mu_T` (lower bound on `inf_k R_f(k)`) and a constraint family
   `R̂_f-derived Toeplitz-PSD ≥ µ_T · I − ε`
   at multiple values of `k` in the support.
2. Minimize `1 − µ_T` as the new objective.
3. Compare to White's `min Ω` at the same `(N, T, R, bochner_n)` setting.

But this is the **same SDP family as `mside_via_lasserre.py`** with a sign flip on the objective. If the existing M-side Lasserre constraint gives `mu_M ≥ X` then a `1 − µ_T`-objective SDP should give a comparable value (up to dual reformulation effects). No new lever expected.

**If the diagnostic chain has ruled out the M-side Bochner family as a path forward (per findings.md), it has also ruled out Lever G.**

---

## Appendix: lines / files referenced

| Reference                                                           | Location |
|---------------------------------------------------------------------|----------|
| White's variables, objective, constraints                            | `lp_research_state/code/white_full_convex.py:99–339` |
| Cell-kernel cosine/sine envelopes (the binding constraints)          | `white_full_convex.py:176–190` |
| Bochner-PSD on `f` and `1 − f`                                       | `white_full_convex.py:233–258` |
| M-side Bochner-PSD (SOC-relaxed)                                     | `lp_research_state/code/mside_bochner.py` |
| M-side Bochner via Schur lifting                                     | `lp_research_state/code/mside_bochner_schur.py` |
| M-side Bochner via Lasserre level-2                                  | `lp_research_state/code/mside_via_lasserre.py` |
| Fourier convention notes                                             | `lp_research_state/code/_fourier_convention_notes.md` |
| Together / Haugland equivalent formulation                           | `/private/tmp/together_repo/erdos-minimum-overlap/README.md:25–28` |
| Together's overlap functional computation                            | `lp_research_state/code/together_loader.py:431–458` |
| Empirical Ω vs M_T gap at `h*`                                      | `_fourier_convention_notes.md:255–263` |
| Why Together's `h*` is sub-optimal under Ω                           | `_fourier_convention_notes.md:283–296` |
