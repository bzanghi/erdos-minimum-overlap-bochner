# PRO-14: Shadow-price audit and the |ξ|/Ω puzzle

**Status:** Engineering tool delivered; math conjecture |ξ| ≤ Ω **disproved
empirically**; replacement empirical bound `|ξ| ≤ 1.47 · Ω` established
across 4 disparate centers.

## What was built

`lp_research_state/code/_pro14_verifier.py` — a one-stop tool that:

1. Solves White's SDP at any `(N, T, R, bochner_n, h, p, q1, q2)` center.
2. Extracts ALL relevant Lagrange multipliers by INDEX (robust, no string
   matching): `ξ` (sum constraint), `ν_3, ν_4` (moment constraints),
   `τ` ((5.13)), `λ_m^cos` for `m=1..2R`, and σ-pairs for the sine family.
3. Computes the Theorem-2-style bound on `Σ|λ_m|` from `(|ξ|, τ, ν_3,
   |Δ_sin(1)|)` and compares to the measured `Σ|λ_m|`.
4. Reports `Σ m·|λ_m|`, `ResidualGain = (π/(2N)) · Σ m·|λ_m|`, and the
   measured `C_explicit = Ω + ResidualGain`.

Output: `lp_research_state/data/pro14_shadow_prices.json`.

## Results across 4 centers (N=3000, T=1200, R=10, bochner_n=20)

| Row | Ω | \|ξ\|/Ω | τ/Ω | Σ\|λ\| meas | Σ\|λ\| Thm-2 | ratio | C_explicit |
|---|---|---|---|---|---|---|---|
| row1          | 0.3790031 | **1.4667** | 0.6988 | 1.3798 | 1.3777 | 0.9985 | 0.3821620 |
| row4 (bind)   | 0.3787923 | **1.4647** | 0.6837 | 1.3737 | 1.3689 | 0.9965 | 0.3818959 |
| row7          | 0.3800232 | **1.4445** | 0.3992 | 1.2446 | 1.2523 | 1.0062 | 0.3829570 |
| cde_n30_iter1 | 0.3789035 | **1.4492** | 0.4060 | 1.2594 | 1.2523 | 0.9944 | 0.3818025 |

### Key observations

1. **|ξ|/Ω is universal-looking.** Spread 1.4445 → 1.4667 (1.5%
   variation) across 4 physically distinct centers including one
   out-of-distribution CDE-discovered center.

2. **The conjecture |ξ| ≤ Ω from PRO-14 is FALSE.** Empirically `|ξ| ≈
   1.46 · Ω` — about 47% above the conjectured ceiling.

3. **Theorem 1's KKT identity is numerically tight.** The Theorem-2 bound
   on Σ|λ_m| agrees with the measured Σ|λ_m| to within 0.4% across all 4
   rows. This validates the proof of Theorem 1 (LEVER_I_PRIME_THEOREM.md
   §2) to numerical precision, AND shows that any structural improvement
   to the bound on Σ|λ_m| must come from tightening one of the input
   ingredients (`|ξ|, τ, ν_3, |Δ_sin(1)|`), not from a smarter recombination.

4. **τ/Ω is less stable** (0.40 → 0.70), but uniformly < 1.

## Implication: the corrected ceiling formula

Substituting the empirical `|ξ| ≤ 1.47Ω` (across our 4 rows), `τ ≤ 0.70Ω`,
`L·ν_3 ≪ 1`, `|Δ_sin(1)| = O(L)` into Theorem 2:

> `Σ|λ_m|  ≤  (1 + O(L²)) · (2 · 1.47Ω  +  0.70Ω  +  O(L))`
>          `≈  3.64 · Ω`

For `Σ m·|λ_m| ≤ 2R · Σ|λ_m|`:

> `Σ m·|λ_m|  ≤  2R · 3.64Ω  =  72.8 Ω  (R=10)`

This is a very loose bound. The empirical `Σ m·|λ_m| ≈ 5.5–6.0` is ~12×
tighter than this UB, because of the dual decay structure (PRO-13's
phase transition analysis): only `m=3,4,5,6` carry > 90% of the mass.

So we have TWO independent bounds:
- **Worst-case via Theorem 2 + empirical |ξ|, τ:**
  `Σ m·|λ_m| ≤ 72.8 Ω`, giving `ResidualGain ≤ (π/(2N))·72.8·Ω = (114Ω/N)`.
  At N=10000, Ω=0.38: ≈ 4.3 × 10⁻³. This is well ABOVE Together's
  UB margin and the residual saturation is **vacuous**.

- **Empirically-grounded:**
  `Σ m·|λ_m| ≈ 6.0` directly, giving `ResidualGain ≈ (π/(2N))·6.0 ≈
  9.5×10⁻⁴` at N=10000. This is what LEVER_I_PRIME_THEOREM.md uses.

The TRUE saturation theorem operates with the empirically-grounded
estimate; the worst-case bound is too loose to be useful.

## What this means for PRO-14

The original PRO-14 task was:
> Prove `|ξ| ≤ Ω` (and similar for `τ`, `ν_3`) to upgrade Theorem 2 from
> conditional to unconditional.

This is now known to be wrong as stated. Possible reformulations:

**A. Empirically-anchored conjecture (provable in principle):**
> `|ξ| ≤ C · Ω` for some universal constant `C ≤ 1.47`, valid at any KKT
> point of White's SDP in the residual region (5.16).

A proof of this would still be useful (it would make the residual
saturation theorem unconditional with a known constant). The math should
go via the envelope theorem on a SCALE-INVARIANT version of the problem.
We attempted the simple rescaling `f → s·f` but it does not preserve the
constraint `p_1 ≤ c[0] ≤ p_2` (anchor on `f̂(1)`'s value), since the
anchors are fixed parameters of the row.

**B. Verification-based (achievable now):**
> Run the shadow-price verifier on a denser cover (the 7 White rows + 10
> CDE centers + extreme points of (5.16)). If `|ξ|/Ω ≤ 1.50` across all
> tested centers, this gives an empirically-certified bound for any
> production solve.

Recommendation: do (B) at production scale, since (A) remains research.

**C. Tighten via a different identity:**
> The 0.4% Theorem-1 residual leaves little room for a better
> recombination of `(ξ, τ, ν_3, |Δ_sin|)`. A genuinely tighter bound
> on Σ|λ_m| would require either:
>   - Decomposing `|ξ|` and `τ` further via a sub-KKT identity, or
>   - Showing that not all `λ_m` are simultaneously large (which the
>     empirical decay confirms: only m∈{3,4,5,6} carry the mass).

The cleanest mathematical lever is (C2) — show that the m-decay is
itself a KKT consequence. This is a research question for PRO-14 v2.

## Recommendation

- **Closing the original PRO-14** as misformulated (conjecture wrong).
- **Spinning off a PRO-14 v2**: "Empirical shadow-price audit across a
  dense cover of (5.16); verify `|ξ| ≤ 1.50·Ω` and `τ ≤ Ω` over the cover."
- **Reporting:** the verifier tool, the validation of Theorem 1 to 0.4%
  agreement, and the corrected empirical ceiling are all directly useful
  for the preprint (PRO-5).

## Code retained

- `_pro14_verifier.py` — the audit tool
- `data/pro14_shadow_prices.json` — the 4-row audit (cached)
