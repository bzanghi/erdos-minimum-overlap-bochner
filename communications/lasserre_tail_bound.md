# Rigorizing the Lasserre Level-2 Localizing Constraint via a Fejér–Riesz / Parseval Tail Bound

**Author:** Ben Zanghi
**Date:** 2026-05-10
**Status:** Mathematical write-up. No SDP solves were performed.

## 1. Setup and notation

Following White (Acta Arith. 2023, §3 and §5), the Erdős minimum-overlap convex
program optimizes over a function `f : [-1, 1] → [0, 1]` with `∫_{-1}^{1} f(x) dx = 1`.
Write the (period-2) Fourier expansion

  `f(x) = 1/2 + Σ_{k≥1} ( c_k cos(πkx) + d_k sin(πkx) )`,

so that with the convention
  `f̂(0) = 1/2,  f̂(k) = (c_k - i d_k)/2  (k ≥ 1),  f̂(-k) = conj f̂(k)`,
we have `|f̂(k)|² = (c_k² + d_k²)/4` for all `k ≥ 1`.

Three program facts that are used below:

- **(A) Range bounds.** `0 ≤ f ≤ 1`, so `‖f‖_∞ ≤ 1`, and `‖f‖_1 = 1`. Hence
  `‖f‖_2² = ∫ f² ≤ ‖f‖_∞ · ‖f‖_1 = 1`.
  By Plancherel, `Σ_{k∈ℤ} |f̂(k)|² = (1/2) ‖f‖_2² ≤ 1/2` (factor `1/2` from the `[-1,1]` measure normalization, equivalently `1/(2L)` with `L = 1`).
- **(B) Box bounds (White (5.10)).** `|c_k|, |d_k| ≤ 2/π`.
- **(C) Parseval / White (3.7), (5.11).** `Σ_{k≥1} (c_k² + d_k²) ≤ 1/2`,
  equivalently `Σ_{k≥1} |f̂(k)|² ≤ 1/8`, and (using `|f̂(0)|² = 1/4`)
  `Σ_{k∈ℤ} |f̂(k)|² ≤ 1/4 + 2 · 1/8 = 1/2`. (Consistent with (A).)

Define the truncation tail energy

  `E(T) := Σ_{k > T} (c_k² + d_k²) / 4 = Σ_{|k|>T} |f̂(k)|² · (1/2)`,

and the Parseval-bounded version

  `E(T) ≤ 1/8 − Σ_{1 ≤ k ≤ T} (c_k² + d_k²) / 4`.        (★)

The right side of (★) is **linear in the SDP variables** `(c_k², d_k²)` (which
are themselves `M_top[k,k]` and `M_top[T_max+k, T_max+k]` after the Lasserre
lift), so it is exactly the kind of expression cvxpy can use as a tail-bound
multiplier inside a PSD constraint.

## 2. The heuristic constraint as written in `lasserre.py`

The file `/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code/lasserre.py`
implements the level-2 Lasserre localizing constraint for `f - f² ≥ 0`:

  `Loc_truncated[j, k] := (f - f²)̂(j - k)`,  for `j, k = 0, …, T_loc`,

required to be Hermitian PSD. The Fourier coefficient `(f - f²)̂(m)` is
expanded as `f̂(m) − Σ_{n} f̂(n) f̂(m − n)`, but in `_f2_hat_re_im` the
convolution sum is **truncated**:

```python
n_lo = max(-T_max, m - T_max)
n_hi = min( T_max, m + T_max)
for n in range(n_lo, n_hi + 1):
    ...   # contributes f̂(n) f̂(m-n) only when |n|, |m-n| ≤ T_max
```

So the implemented matrix entry is

  `Loc_truncated[j, k] = f̂(m) − Σ_{|n|≤T_max,  |m-n|≤T_max} f̂(n) f̂(m-n)`,
  with `m = j − k`.

Equivalently, the omitted tail at lag `m` is

  `T_m := Σ_{|n|>T_max  OR  |m-n|>T_max} f̂(n) f̂(m-n)`,

so that `(f - f²)̂(m) = Loc_truncated_entry(m) − T_m`.
Because the tail `T_m` is **not bounded** in the implementation, the resulting
PSD requirement on `Loc_truncated` is heuristic: it tests `f² ≤ f` only against
trig polynomials of degree `≤ T_loc` whose Fourier-side bilinear form involves
only the in-band coefficients. The omitted high-frequency cross terms are silently dropped.

The Lasserre paragraph in the docstring is explicit about this: "for an SDP
UPPER bound on Ω we can include only the truncated sum (giving a TIGHTER
constraint, hence a LOOSER bound on Ω)" — i.e., the dropped tail `T_m` is
*not* a one-signed quantity, so the truncation is not conservatively safe in
the lower-bound direction.

## 3. The rigorous tail bound

We want a closed-form upper bound on `|T_m|` that is **affine in the SDP variables**, so that subtracting it as a diagonal correction gives a valid PSD constraint.

### 3.1. Diagonal lag (`m = 0`)

For `m = 0`, the tail is `T_0 = Σ_{|n|>T_max} f̂(n) f̂(−n) = Σ_{|n|>T_max} |f̂(n)|²`,
which is **real and nonnegative**. By Plancherel/Parseval,

  `T_0 = 2 Σ_{k>T_max} |f̂(k)|² = 2 · (1/2) Σ_{k>T_max} (c_k²+d_k²)/2`
       `= (1/2) Σ_{k>T_max} (c_k²+d_k²)`.

Wait — let me redo this carefully. `|f̂(k)|² = (c_k² + d_k²)/4` for `k ≥ 1`, and
`|f̂(−k)|² = |f̂(k)|²`, so

  `T_0 = Σ_{|n|>T_max} |f̂(n)|² = 2 Σ_{k>T_max} (c_k² + d_k²)/4 = (1/2) Σ_{k>T_max} (c_k²+d_k²)`.

By Parseval (C), `Σ_{k≥1}(c_k²+d_k²) ≤ 1/2`, so

  `T_0 ≤ (1/2) [ 1/2 − Σ_{1≤k≤T_max} (c_k² + d_k²) ] = 1/4 − (1/2) Σ_{1≤k≤T_max}(c_k²+d_k²)`.   (D₀)

(D₀) is **affine in the lifted bilinears** `c_k² = M_top[k,k]`,
`d_k² = M_top[T_max+k, T_max+k]`, which by the Schur lift are linear in
`M_top`. So `T_0` admits a tight, linear-in-program-variables upper bound.

This bound is **exact** (i.e. tight up to the slackness in Parseval (C)).
It is the only piece of the proof that does not require Cauchy–Schwarz.

### 3.2. Off-diagonal lag (`m ≠ 0`)

For `m ≠ 0`, the tail
  `T_m = Σ_{|n|>T_max  or  |m-n|>T_max} f̂(n) f̂(m-n)`
is no longer a positive sum, so we control `|T_m|` by Cauchy–Schwarz.
Split the bad index set
  `B_m := { n ∈ ℤ : |n|>T_max or |m-n|>T_max }`
into the union `B_m^A ∪ B_m^B`, where

  `B_m^A = { n : |n| > T_max }`,  `B_m^B = { n : |m − n| > T_max }` (some overlap allowed: the union bound is OK).

Then

  `|T_m| ≤ Σ_{n ∈ B_m^A} |f̂(n)| · |f̂(m-n)| + Σ_{n ∈ B_m^B} |f̂(n)| · |f̂(m-n)|`.

For the first sum, by Cauchy–Schwarz on `ℓ²(ℤ)`,

  `Σ_{|n|>T_max} |f̂(n)| · |f̂(m-n)|`
  `≤ ( Σ_{|n|>T_max} |f̂(n)|² )^{1/2} · ( Σ_{|n|>T_max} |f̂(m-n)|² )^{1/2}`
  `≤ ( Σ_{|n|>T_max} |f̂(n)|² )^{1/2} · ( Σ_{j ∈ ℤ} |f̂(j)|² )^{1/2}`
  `= sqrt(T_0) · sqrt(1/2)`        (using Plancherel for the second factor).

(Note: I dropped the second factor's restriction to `m-n ∈ B_m^A` and replaced it by all of `ℤ`; this only makes the bound larger, hence still valid.)

The second sum is symmetric: substituting `n' = m - n` gives the same expression. So

  `|T_m| ≤ 2 · sqrt(T_0) · sqrt(1/2) = sqrt(2) · sqrt(T_0) = sqrt(2 T_0)`.    (D_m)

Using (D₀),

  `|T_m| ≤ sqrt( 2 · [ 1/4 − (1/2) Σ_{1≤k≤T_max}(c_k²+d_k²) ] )`
        `= sqrt( 1/2 − Σ_{1≤k≤T_max}(c_k²+d_k²) )`.        (E)

This bound is uniform in `m`. Crucially, the right-hand side of (E) is
the square root of an affine expression in lifted bilinears
`c_k² + d_k² = M_top[k,k] + M_top[T_max+k, T_max+k]`. cvxpy can express this
via a slack scalar:

```python
S = 1/2 - cp.sum([M_top[k,k] + M_top[T_max+k, T_max+k] for k in range(1, T_max+1)])
# S is linear in M_top, and S ≥ 0 because of Parseval (C).
C = cp.Variable(nonneg=True)            # the tail-bound scalar
cons.append(cp.square(C) <= S)          # SOC: C ≤ sqrt(S)
```

Then `C ≥ |T_m|` for every `m ≠ 0`, and (using (D₀))

  `C₀ := 1/4 − (1/2) Σ_{k≤T_max}(c_k² + d_k²) ≥ T_0`     (linear, no SOC needed).

### 3.3. The rigorized localizing constraint

Define the **diagonal tail correction matrix**

  `C_{T_max} = diag( C₀, C, C, …, C )`     (size `(T_loc+1) × (T_loc+1)`).

Note `C₀` is in slot (0,0) and `C` (the SOC slack) sits on every other diagonal slot. The choice of the diagonal pattern is justified next.

**Claim (rigorous PSD constraint).** The constraint

  `Loc_rigorous := Loc_truncated  −  C_{T_max}  ⪰  0`        (R)

is a valid Lasserre level-2 localizing constraint (i.e., it is implied by `f² ≤ f` together with (A)–(C)).

**Proof.** Let `Δ := Loc_truncated − Loc_true`, so that `Δ[j,k] = T_{j-k}`.
By (D₀) and (E), each row `j` of `Δ` has `ℓ¹` norm bounded by

  `Σ_k |Δ[j,k]| = Σ_k |T_{j-k}| ≤ |T_0| + (T_loc) · max_{m≠0} |T_m| ≤ C₀ + T_loc · C`.

But we don't actually need a Gershgorin argument: the claim follows from a
sharper bound. We show `Δ ⪯ C_{T_max}` (in the Loewner order), i.e.,
`C_{T_max} − Δ ⪰ 0`, which combined with `Loc_true ⪰ 0` gives
`Loc_rigorous = Loc_truncated − C_{T_max} = Loc_true + (Δ − C_{T_max}) ⪰ Loc_true − (C_{T_max} − Δ) +  …`.
This is *not quite right* — let me be careful.

We have `Loc_truncated = Loc_true + Δ` where `Loc_true ⪰ 0` (Bochner, since
`f − f² ≥ 0`). So

  `Loc_rigorous = Loc_truncated − C_{T_max} = Loc_true + (Δ − C_{T_max})`.

For `Loc_rigorous ⪰ 0` we need `Δ − C_{T_max} ⪰ −Loc_true`, i.e., it suffices
that **`C_{T_max} ⪰ Δ` in the Loewner order**. (Then
`Loc_rigorous = Loc_true − (C_{T_max} − Δ) ⪰ 0` is **not** automatic — the
sum of a PSD and a non-PSD is not PSD in general.)

The honest statement of the rescue is therefore that **(R) is the constraint
we wish to hold**, and to *prove* it follows from `f − f² ≥ 0` we must
separately bound `C_{T_max} − Δ ⪰ 0` (which is what `C_{T_max}` is for).

For diagonal `Δ` (which `Δ` is **not** — `T_m` for `m ≠ 0` populates the off-diagonals), Gershgorin does the job:

  `Δ ⪯ ‖Δ‖_op · I  ≤  diag-of-row-sums I`.

But row sums depend on `T_loc · C` which is large. A sharper, rigorous version
uses the **Schur-product / Toeplitz structure of `Δ`**: `Δ` is a `(T_loc+1) ×
(T_loc+1)` Hermitian Toeplitz matrix with symbol `g(θ) = Σ_m T_m e^{imθ}`. By
the Bochner–Toeplitz theorem,

  `‖Δ‖_op = sup_θ |g(θ)| ≤ Σ_m |T_m|`.

The sum `Σ_{m=−T_loc}^{T_loc} |T_m|` is bounded by `|T_0| + 2 T_loc · max_m |T_m|`,
which — using (D₀) and (E) — gives an upper bound

  `‖Δ‖_op ≤ C₀ + 2 T_loc · sqrt( 1/2 − Σ_{k≤T_max}(c_k²+d_k²) )`.       (F)

So setting **`C_{T_max} = ‖Δ‖_op · I`** with `‖Δ‖_op` upper-bounded by (F) yields a rigorous PSD constraint via

  `Loc_rigorous = Loc_truncated − ( C₀ + 2 T_loc · C ) · I  ⪰ 0`,
  with `C ≥ sqrt( 1/2 − Σ_{k≤T_max}(c_k²+d_k²) )` enforced by SOC.

This is the cleanest statement and the one I recommend.

## 4. Quantitative forecast at White's `T_max = 30`

At White's optimum (`Ω* ≈ 0.379544`), the in-band Parseval slack
`1/2 − Σ_{1≤k≤30}(c_k²+d_k²)` is small but nonzero. From the empirically observed
coefficient decay `|c_k|, |d_k| ≲ 2/(π k²)` (heuristic only — not proven; comes from
inspecting solver output of the row-6 White table), the missing tail energy
`Σ_{k>30}(c_k²+d_k²)` is roughly

  `Σ_{k>30} 2 · (2/(π k²))² ≈ (8/π²) · Σ_{k>30} k^{-4} ≈ (8/π²) · (1/(3·30³)) ≈ 1.0 × 10⁻⁵`.

Then:

- `C₀ ≤ (1/2) · 1.0 × 10⁻⁵ ≈ 5 × 10⁻⁶` (negligible vs Loc diagonal ~ Ω/2 ≈ 0.19).
- `C ≤ sqrt(1.0 × 10⁻⁵) ≈ 3.2 × 10⁻³` (this is the off-diagonal magnitude bound).
- `‖Δ‖_op ≤ C₀ + 2 · 30 · C ≈ 5 × 10⁻⁶ + 60 · 3.2 × 10⁻³ ≈ 0.19`.

**This last number is catastrophically large** — it is the same order as the diagonal of `Loc` itself (`(f-f²)̂(0) ≈ Ω/2 ≈ 0.19`). The Toeplitz-operator-norm bound (F) is too loose to leave any rigorous gain.

The reason is that `Σ_{m=-T_loc}^{T_loc} |T_m|` includes `2 T_loc + 1 ≈ 61` terms, each as large as `3 × 10⁻³`. Cauchy–Schwarz at each lag is uniform, and the resulting `2 T_loc · C` factor dominates `C₀`.

A sharper rescue would replace (F) with an `ℓ²`-on-symbol estimate:
  `‖Δ‖_op ≤ ‖g‖_∞ ≤ ‖g‖_{H^{1/2+ε}}` and use that `Σ_m |T_m|² ≤ ?` is more tractable. But by Plancherel,
  `Σ_m |T_m|² = ∫ |g|² ≤ ∫ |g| · ‖g‖_∞`,
which is circular without an independent `‖g‖_∞` bound.

A more honest route is to bound the operator norm via **the Hilbert-Schmidt norm of `Δ`**, which would replace the `2 T_loc · C` term with `sqrt(2 T_loc + 1) · C`:

  `‖Δ‖_op ≤ ‖Δ‖_HS = sqrt( Σ_{j,k} |Δ[j,k]|² ) = sqrt( (2T_loc+1) · ‖T_·‖_2² )`,

still giving a `sqrt(60)` ≈ 8 factor and a final bound `8 · 3 × 10⁻³ ≈ 0.025`,
**which is also still ~100× larger than the +2.8 × 10⁻⁴ heuristic gain.**

### Forecast

The rescued contribution is **essentially zero** (or negative — it would
*degrade* the rigorous bound). At `T_max = 30`, the Toeplitz-Cauchy–Schwarz
machinery is too weak to certify the +2.8 × 10⁻⁴. To rigorize the Lasserre
gain one of the following must be true:

1. The `ℓ¹`-tail `Σ_{k>T_max} |c_k|` decays much faster than the `ℓ²`-tail
   would suggest (which is the case for analytic `f`, but White's program does not
   prove analyticity).
2. The constraint is reformulated to use the `ℓ²` tail directly rather than via Loc operator norm — e.g., by treating the tail as a separate slack absorbed into a different PSD block. I have not found such a formulation that closes.
3. `T_max` is increased substantially. To bring `‖Δ‖_op` below `2 × 10⁻⁴` would
   require `T_max` such that `60 · sqrt(missing tail)` is small — i.e.,
   `missing tail` ≪ `10⁻¹¹`, which under `k^{-4}` decay needs `T_max > 1000`. This is computationally feasible for the moment matrix (`(2·1000+1)² ≈ 4 × 10⁶` entries) but borderline for a level-2 lift PSD solve.

## 5. cvxpy implementation pseudocode

Given the forecast above, here is what an *attempted* rigorous version would look like. I include it for completeness; I do not claim it improves the bound.

```python
def add_lasserre2_constraint_rigorous(cons, c, d, T_max, T_loc):
    # ... (build M_top, border, link rows as before) ...

    # Parseval slack S in [0, 1/2], linear in M_top diagonal entries.
    S = cp.Constant(0.5)
    for k in range(1, T_max + 1):
        # M_top[k, k] = c_k^2 (post-lift); M_top[T_max+k, T_max+k] = d_k^2.
        S = S - M_top[k, k] - M_top[T_max + k, T_max + k]
    cons.append(S >= 0)   # implied by Parseval but tightens numerics

    # SOC slack: C_off >= sqrt(S), so C_off bounds |T_m| for every m != 0.
    C_off = cp.Variable(nonneg=True)
    cons.append(cp.square(C_off) <= S)

    # Diagonal tail at m=0:  T_0 <= 0.5 * S / (1 - 0) ... actually = (1/2) * (S - small).
    # Use C_0 := 0.5 * S as a safe affine upper bound on T_0.
    C_0 = 0.5 * S

    # Toeplitz operator-norm bound (Gershgorin/symbol-sup).
    # ||Delta||_op <= C_0 + 2 * T_loc * C_off.
    op_bound = C_0 + 2 * T_loc * C_off

    # Build the truncated Loc as before:
    # ... Loc_truncated = bmat([[Re_M, -Im_M], [Im_M, Re_M]]) ...

    # The rigorous PSD constraint:
    n_real = 2 * (T_loc + 1)
    cons.append(Loc_truncated - op_bound * np.eye(n_real) >> 0)
```

This is convex (the LMI `Loc_truncated − op_bound · I ⪰ 0` is jointly convex in
`(c, d, M_top, op_bound)` since `op_bound` enters linearly). The SOC constraint
`C_off² ≤ S` is also convex.

## 6. Hazards and adversarial-review honesty

- **The Toeplitz bound (F) is too loose.** As computed in §4, `2 T_loc · C ≈ 0.19`
  at `T_max = T_loc = 30`. This swamps the heuristic gain by 3 orders of magnitude.
  The HS-norm refinement still falls short by ~100×.
- **Sign convention.** I checked: `T_m := (Loc_truncated − Loc_true)[j, k]`
  with `m = j − k`. The sign in (R) (subtracting `C_{T_max}`) is correct
  *if* we want to lower-bound `Loc_true ⪰ 0`. Note `Δ` can be either sign,
  but `op_bound` is a uniform two-sided bound on `‖Δ‖_op`, so subtracting it
  from `Loc_truncated` certifies `Loc_true ⪰ Loc_truncated − op_bound · I − Δ`...
  actually, the certified statement is
    `Loc_truncated − op_bound · I ⪯ Loc_truncated − Δ = Loc_true`,
  i.e., `Loc_true ⪰ Loc_truncated − op_bound · I`,
  so requiring `Loc_truncated − op_bound · I ⪰ 0` is **a sufficient but not necessary** condition for `Loc_true ⪰ 0`. Fine for an LP/SDP relaxation: tighter than nothing, looser than the heuristic.
- **Numerical conditioning.** At large `T_max` the Schur-bordered matrix is
  `(2T_max + 2) × (2T_max + 2)`. CLARABEL handles ~6000-dim PSD blocks in ~1 hr;
  going to `T_max = 1000` (block dim ~2002) is feasible but slow.
- **The tail-bound matrix shape.** I argued for `op_bound · I`. A shape-aware
  alternative is `diag(C_0, C_off, C_off, …)`, exploiting that the diagonal
  lag is *much* smaller (`5 × 10⁻⁶`) than the off-diagonal bound. This would
  shave a factor of 2 at best, not the orders of magnitude needed.
- **Even-`f` reduction.** In the `assume_even=True` mode, `d_k = 0`, halving the
  `ℓ²` energy. But the `ℓ¹` issue (Toeplitz row sum) is unchanged.
- **Adversarial review.** Honest verdict: this write-up does **not** rescue the
  +2.8 × 10⁻⁴. A reviewer with experience in Toeplitz/Bochner machinery would
  spot the `2 T_loc` factor in (F) immediately and ask for a sharper symbol
  estimate. My honest forecast is that no `O(1)` symbol estimate exists: the
  fundamental issue is that we have only an `ℓ²` tail control, and `‖·‖_op ≥
  ‖·‖_{HS}/sqrt(n) ≥ ‖f̂‖_2`, so the operator norm of a length-`(T_loc+1)`
  Toeplitz with `ℓ²`-bounded symbol is at best `‖symbol‖_∞`, which we cannot
  improve without additional smoothness.

## 7. Summary

The rigorous tail bound goes through algebraically: (D₀), (E), (F), and the
cvxpy encoding in §5 are all correct and convex. But the *quantitative*
rescue fails by 2–3 orders of magnitude at White's parameter setting.
**The +2.8 × 10⁻⁴ heuristic gain cannot be rescued by Cauchy–Schwarz / Plancherel / Parseval alone at `T_max = 30`.**

To certify the gain rigorously, the path forward is **not** via tail bounds on
the existing truncation, but via either (a) a direct proof that the Lasserre
moment-matrix structure is exact for *low-degree* moments (which, since the SOS cone
at degree 4 in the Fourier coefficients is finite-dimensional, is provable but
requires a separate combinatorial argument), or (b) increasing `T_max` to the
point where the tail is negligible at the SDP precision (`T_max > 1000`).
Neither is in scope for this note.
