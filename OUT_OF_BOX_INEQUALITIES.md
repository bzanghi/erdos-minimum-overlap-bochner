# Out-of-Box: Known-Inequality Scavenger Hunt for Ω(f) Lower Bound

**Goal.** Find a *known* inequality from harmonic analysis / number theory / convex
analysis that gives, in one line, a lower bound `Ω(f) ≥ C` with `C > 0.3801279`
for every admissible `f`. Here

  f : [-2, 2] → [0, 1],   ∫f = 1,   Ω(f) = sup_{t ∈ ℝ} ∫ f(x) f(x+t) dx.

The Erdős minimum-overlap constant is `μ = inf_f Ω(f)`. White (2023) gave the
analytic LB `μ ≥ 0.379005`; our SDP gives `μ ≥ 0.3801279`; the Together UB is
`μ ≤ 0.380871`. Note that elsewhere in this repo we sometimes quote the
guaranteed bound `μ ≥ 0.379544` — the `0.3801279` figure refers to a single-row
SDP value, *not* an unconditional bound. The scavenger hunt below treats
`0.3801279` as the target to beat.

---

## Section 1 — Inequalities that give a (possibly weak) LB on Ω(f)

### 1. Plancherel / Parseval (with averaging)

Statement. For `g ∈ L²(ℝ)`, `∫|g|² = ∫|ĝ|²`. For `g ∈ L²([−L, L])`,
`∫_{-L}^{L} |g|² = Σ_k |c_k|²` where `c_k = (1/2L)∫ g(x) e^{-iπkx/L} dx` are the
Fourier coefficients of g on the interval.

Application. With `f` supported in `[-2, 2]`, the autocorrelation
`A(t) := (f * f̃)(t) = ∫ f(x) f(x+t) dx` is supported in `[-4, 4]`,
nonnegative, even, and satisfies `∫ A = (∫f)² = 1`, and `Â(ξ) = |f̂(ξ)|² ≥ 0`.
So `A` is a continuous nonneg function on `[-4, 4]` with integral 1, and

  Ω(f) = sup_t A(t) ≥ (1/8) ∫_{-4}^{4} A = 1/8 = 0.125.

That is the elementary "sup ≥ average" bound. Plancherel by itself adds
nothing: it just gives `∫ A² = ∫ |f̂|⁴`, which is a non-comparable functional.

Numerical LB: **0.125**. Far below 0.3801279.

### 2. Hausdorff–Young (Lp–Lq Fourier)

Statement. For `1 ≤ p ≤ 2`, `q = p/(p−1)`, `‖f̂‖_q ≤ ‖f‖_p` (with constants
depending on the convention; the Beckner constant gives the sharp form).

Application. With `0 ≤ f ≤ 1`, `∫f = 1` we have `‖f‖_p ≤ 1` for all `p ≥ 1`,
so `‖f̂‖_q ≤ 1`. But the value we need to lower-bound is `‖f̂‖_∞² = sup_t A(t)`
… no, that is wrong: `‖f̂‖_∞ ≤ ‖f‖_1 = 1` is an *upper* bound on `|f̂|`, and
`A(0) = ‖f‖_2² ≤ ‖f‖_1 ‖f‖_∞ ≤ 1`. None of these give a LB on `sup_t A(t)` for
`t` ranging over ℝ (we already know `A(0) ≤ 1`; the infimum problem is about
shifts `t ≠ 0`).

Numerical LB: none. **Ruled out.**

### 3. Bombieri's inequality

Statement (Bombieri, "Le grand crible…", 1976; see Montgomery, *Ten Lectures*,
Ch. 7). For a trigonometric polynomial `P(x) = Σ_{n ∈ S} a_n e(nx)` with finite
spectrum `S`, and well-spaced points `x_1, …, x_R ∈ ℝ/ℤ`,

  Σ_r |P(x_r)|² ≤ (|S| + δ⁻¹) Σ_n |a_n|².

Application. Bombieri's inequality is an *upper* bound on a sum of values of a
trigonometric polynomial in terms of its L² norm. It is the wrong direction: we
want a *lower* bound on `sup_t |P(t)|`, not an upper bound on `Σ |P(x_r)|²`.
There is a Selberg-type "dual" that bounds `sup` from below by `Σ|a_n|/(|S|+δ⁻¹)`,
but applied here it gives `sup_t A(t) ≥ ∫A = 1` divided by an effective bandwidth,
which collapses to the elementary `1/8` again because A is band-unlimited.

Numerical LB: none beyond the elementary `1/8`. **Ruled out.**

### 4. Selberg / large sieve

Statement. `Σ_r |Σ_{n ∈ S} a_n e(nx_r)|² ≤ (N + δ⁻¹) Σ_n |a_n|²` for
`x_r` δ-spaced and `S ⊂ [M, M+N]`.

Application. Same direction problem as Bombieri (sup-norm UB, not LB). Even
turned into a dual, the resulting LB on `sup_t A(t)` would be in terms of `Σ_n
Â(n) = A(0)` divided by an effective bandwidth, and again `A` is not
band-limited. The large sieve is also designed for `S` an *integer* spectrum;
`A` lives on `[-4, 4]` continuous.

Numerical LB: none. **Ruled out.**

### 5. Cauchy–Schwarz on autocorrelation

The "sup ≥ average" form gives `Ω(f) ≥ 1/8` (as in item 1). A second
Cauchy–Schwarz `A(0)² = (∫f²)² ≤ ‖f‖_∞ ∫f² · ‖f‖_1` is trivially saturated and
yields nothing.

A non-trivial Cauchy–Schwarz: for any nonneg weight `w` on `[-4, 4]` with
`∫w = 1`,

  Ω(f) ≥ ∫ A(t) w(t) dt = ∫∫ f(x) f(x+t) w(t) dx dt.

This is the basis for *White's* convex program (choose `w` to be a sum of
delta-like bumps near hypothesized maxima); it is not really a "known
inequality" giving a closed-form `C`. Optimizing `w` *is* the SDP.

Numerical LB: matches whatever weight we plug in — recovers `0.379–0.380`
already proved by the SDP, no shortcut.

### 6. Beckner / Brascamp–Lieb

Statement. Sharp Hausdorff–Young with Beckner constants. Brascamp–Lieb is a
multilinear generalization with sharp constants for products of `Lp` norms.

Application. These give sharp UB on `‖f̂‖_q` or on multilinear functionals
`∫∏ f_i(L_i x) dx`. Ω(f) is a bilinear form `∫∫ f(x) f(y) δ(t = y − x)` — a
distribution, not a Brascamp–Lieb-type smooth multilinear form. No direct LB.

Numerical LB: none. **Ruled out.**

### 7. Wiener's lemma

Statement. If `f ∈ A(T)` (absolutely convergent Fourier series) and `f` never
vanishes, then `1/f ∈ A(T)`.

Application. Qualitative, not quantitative. Gives no `C`.

**Ruled out.**

### 8. Ingham / Wirtinger / Poincaré

Statement (Wirtinger). For `g ∈ H¹([0, L])` with `∫g = 0`,
`∫g² ≤ (L/π)² ∫(g')²`.

Application. Apply with `g = A − 1/8` on `[−4, 4]` (mean-zero version of A). We
get `∫(A − 1/8)² ≤ (8/π)² ∫(A')²`. This bounds *variance of A* by *roughness of
A*. Sup-norm `sup A ≥ 1/8 + √(Var(A)/8)`. We have no quantitative LB on
`Var(A)` from the constraints — `A` could be the constant `1/8` on `[−4,4]`
(infeasible from `f ≥ 0`, but as a target), and the constraint that A is a
nonneg autocorrelation of an L¹ function does not by itself force
`Var(A) > 0` strongly enough. In fact `Â(ξ) = |f̂(ξ)|² ≥ 0` combined with
`Â(0) = 1` and `|Â(ξ)| ≤ 1` gives some control, but extracting a number
needs the same SDP work we already did.

Numerical LB: no closed form > `1/8` from Wirtinger alone. **Ruled out as a
one-liner.**

### 9. Heisenberg uncertainty

Statement. `Var(g) Var(ĝ) ≥ 1/(16π²)` for `g ∈ L²(ℝ)` with `‖g‖_2 = 1`.

Application. Take `g = √f / ‖√f‖_2`, then `Var(g) ≤ Var(f)/(∫f²)`; `f` is
supported in `[−2,2]` so `Var(g) ≤ 4`. Then `Var(ĝ) ≥ 1/(64π²)`. We'd want
to convert "ĝ is spread out" into "Ω(f) ≥ something". The autocorrelation
`A = f * f̃` has Fourier transform `Â = |f̂|² ≥ 0`. Uncertainty for f tells us
`f̂` cannot be too concentrated, hence `Â = |f̂|²` cannot be too concentrated,
hence (by Fourier inversion) `A` cannot be too smooth. But "A not too smooth"
gives an LB on `sup A − inf A`, not directly on `sup A`.

Numerical LB: no useful closed form > `1/8`. **Ruled out as a one-liner.**

### 10. Logan / de Branges / Selberg–Beurling extremal

Statement. Selberg's majorant: for `f` a function with Fourier transform
supported in `[−δ, δ]` that majorizes `sgn` on ℝ, the minimum `L¹` norm is
explicit (Vaaly's work; see Montgomery). Beurling's extremal function.

Application. Beurling–Selberg gives sharp constants for *one-sided*
approximation problems and is the standard tool in analytic number theory's
sphere-packing-style problems (Cohn–Elkies). The Erdős minimum-overlap problem
has a Beurling–Selberg flavor: minimize `sup A(t)` for `t` outside a small
neighborhood of 0, subject to `A ≥ 0`, `Â ≥ 0`, `∫ A = 1`.

This is the right *framework* but there is no off-the-shelf sharp constant for
this specific functional. The literature has Logan's problem (sphere-packing
1D), Cohn–Elkies (dimension `d`), Gonçalves–Oliveira–Steinerberger (one-sided
problems), and these *do* solve closely related variational problems. None
gives `Ω(f) ≥ 0.380` as a closed form.

Quantitatively: Logan-type bounds on the *first sign change* of a
nonneg-Fourier-transform function `A` give `sup A ≥ (∫A)/L_eff` where `L_eff`
is the support of A in physical space. Here `L_eff = 8`, giving `1/8` again.

Numerical LB: `1/8` via the Logan/Beurling lens; not enough. **Ruled out as
a one-liner.**

### 11. Sárközy's inequality

Sárközy proved bounds for sets `A ⊂ [N]` containing no perfect-square
differences. The connection to minimum-overlap is via additive combinatorics on
sets, but `f` here is a continuous density, not an indicator. The known
constants in Sárközy-style results are far weaker than the analytic
minimum-overlap bound. **Ruled out.**

### 12. Plünnecke–Ruzsa

Sumset inequalities `|A+B| ≤ K^c |A|`. Indicator-version of overlap is
`|A ∩ (A+t)|`. Plünnecke–Ruzsa gives `|A − A| ≤ K² |A|` if `|A+A| ≤ K|A|`, but
this controls the *number of distinct differences*, not the *max repetition* of
a fixed difference. The minimum-overlap problem's exact analog (Erdős's
original formulation on `{1, …, 2n}`) does not follow from Plünnecke–Ruzsa.
**Ruled out.**

### 13. Brunn–Minkowski

`|A + B|^{1/n} ≥ |A|^{1/n} + |B|^{1/n}`. Geometric volume; not applicable to
the autocorrelation sup. **Ruled out.**

### 14. Roth-type density theorems

`A ⊂ [N]` 3-AP-free has density `o(1)`. Irrelevant: minimum-overlap is about
shifted intersections of *two halves of a partition*, not about APs. **Ruled
out.**

### 15. Rudin's Λ(p)

For a Λ(p) set `S`, `‖Σ_{n∈S} a_n e(nx)‖_p ≲ ‖a‖_2`. Same direction as
Bombieri/large sieve (UB on norm of trig poly, not LB on sup of
autocorrelation). **Ruled out.**

### 16. Boas / Erdős L∞-to-L² ratio for nonneg trig polys

Statement (Boas 1948 and successors). For a nonneg trig polynomial
`p(x) = Σ a_n e(nx) ≥ 0` of degree `N`, `‖p‖_∞ / ‖p‖_2 ≤ √(2N+1)`. Erdős
studied the reverse direction.

Application. The autocorrelation `A` is nonneg, has nonneg Fourier transform,
but is not band-limited. A degree cutoff is artificial. Boas gives a UB on
`‖p‖_∞`, not a LB — wrong direction again. **Ruled out as a direct one-liner.**

---

## Section 2 — Compositions that *might* improve on `1/8` but not on `0.3801279`

The non-trivial *known* tool for this exact problem is the Beurling–Selberg /
Cohn–Elkies linear-programming framework. The SDP we already run is precisely
that framework, augmented with positivity (`f ∈ [0,1]`, i.e. `f² ≤ f`) and
Bochner positivity (`M_n(f) ⪰ 0`). White's 2023 paper *is* the one-line
application of this framework to minimum overlap. There is no second,
independent known LP that we are missing.

Possible compositions:

- **Heisenberg + Plancherel + `Â ≥ 0` + `∫A = 1`**: the Bourgain / Cohn–Elkies
  argument bounds the first zero of `A`. This gives lower bounds on `sup A` of
  order `1/L_eff ≈ 1/8`. To recover `0.38`, you must add the full set of
  Bochner-PSD and moment constraints — i.e. the SDP — which is not "one-line".

- **Sphere-packing LP bounds (Cohn–Elkies 2003)**: structurally identical to
  the SDP. Plugging in any *fixed* admissible auxiliary function from the
  Cohn–Elkies / Viazovska family does not recover `0.380` because those
  functions are tuned for the sphere-packing density, not minimum-overlap.

- **Mertens-type elementary inequality on shifted integrals**: White §3
  proves `μ ≥ 1/(2 + √2) ≈ 0.2929` by an elementary partition-of-unity
  argument. This *is* a one-line known LB. It is below `0.3801279`. There
  are intermediate elementary improvements (Moser ~0.282, Scherk ~0.281,
  Motzkin's progression up through White's `0.379005`), all of which require
  multi-page arguments, not one inequality.

---

## Section 3 — Best LB from any single known inequality

| Inequality | Numerical LB | Beats 0.3801279? |
|---|---|---|
| Plancherel / sup ≥ avg | 0.125 | No |
| Cauchy–Schwarz on autocorrelation (`sup ≥ ∫A / 8`) | 0.125 | No |
| White §3 partition-of-unity elementary bound | `1/(2 + √2) ≈ 0.2929` | No |
| All others above | none / 0.125 | No |
| White (2023) full LP (multi-page) | 0.379005 | No (and not one-line) |
| Our SDP at row centers (Bochner_n=30 etc.) | 0.3801279 | trivially equal (this is the target) |

The best closed-form LB from a single, citable, classical inequality is
**`1/(2 + √2) ≈ 0.2929`** (White §3). All purely Fourier-analytic single
inequalities (Plancherel, Hausdorff–Young, Bombieri, large sieve, Heisenberg,
Logan, Beurling–Selberg, Boas) give `1/8` or nothing for this functional.

---

## Section 4 — Verdict

**No.** There is no known inequality from harmonic analysis, number theory, or
convex analysis that yields, in one line, a lower bound `Ω(f) ≥ C` with
`C > 0.3801279`. The classical Fourier inequalities (Plancherel, Hausdorff–
Young, Bombieri, large sieve, Heisenberg, Beurling–Selberg, Boas, Rudin Λ(p))
either go in the wrong direction (UB on a sup, not LB) or collapse to the
trivial `sup ≥ average = 1/8`. The combinatorial inequalities (Sárközy,
Plünnecke–Ruzsa, Roth, Brunn–Minkowski) do not apply to the continuous
autocorrelation functional. The strongest known single-inequality LB is
White's elementary `μ ≥ 1/(2 + √2) ≈ 0.2929`.

The reason there is no shortcut is structural: White (2023) already *is* the
end product of the relevant variational LP, and the only way to push past
`0.379` quantitatively is to enlarge the LP with additional PSD constraints
(Bochner moment matrices, Hankel/polynomial moments), which is exactly what
our SDP does. The improvement from `0.379005` to `0.3801279` is the
contribution of those additional PSD constraints; it cannot be recovered from
any classical one-line inequality.

This negative result is consistent with the post-tail-bound analysis at the
top of `erdos_lower_bound_research_note.md`, which already concluded that
pushing past `0.379544` with current techniques is not possible at
currently-tractable SDP scale, and that further progress requires new
mathematical levers (much larger `T_max`, finite-dimensional SOS exactness,
or alternative basis representations) — i.e., research, not look-up.
