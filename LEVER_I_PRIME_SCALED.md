# Lever I' Step D: Cell-Envelope Saturation Theorem at Scaled N

**Status:** N=15000 complete (4/4 rows). N=20000 in progress.

**Goal:** Test the corrected cell-envelope cosine + sine saturation bound at
larger `N` to see if (a) the multiplier sums `Σ m·λ` (cos) and `Σ m·σ` (sin)
remain near the N=3000 measured values, (b) C_explicit drops below
Together's UB at any tractable `N`.

**Predecessor:** [LEVER_I_PRIME_THEOREM.md](LEVER_I_PRIME_THEOREM.md), [SESSION_FINAL.md](SESSION_FINAL.md).

**Data:** [lp_research_state/data/lambda_m_scaled.json](lp_research_state/data/lambda_m_scaled.json) (incremental).
**Code:** [lp_research_state/code/_lever_i_prime_lambda_m_scaled.py](lp_research_state/code/_lever_i_prime_lambda_m_scaled.py).

---

## 1. Formula recap

Per Theorem 3 of LEVER_I_PRIME_THEOREM.md, the corrected cell-envelope
residual at scale `N` is

```
ResidualGain ≤ (π/(2N)) · Σ_m m·(λ_m^cos + σ_m^1 + σ_m^2)
             + (π²Ω/(3N³)) · Σ_m m³·(λ_m^cos + σ_m^1 + σ_m^2)
```

with `Ω ≈ 0.38`, `λ_m^cos` the cosine cell-envelope duals (line 182 of
white_full_convex.py), and `σ_m^{1,2}` the sine cell-envelope duals (lines
187, 190). `C_explicit(N) := 0.3801279 + ResidualGain(N)`; theorem
non-vacuous iff `C_explicit(N) < 0.380871`.

---

## 2. Measurements

### 2.1. N=3000 (baseline, from prior session)

| Row | `Σ m·λ` (cos) | `Σ m³·λ` | `Σ m·σ` (sin) | `Σ m³·σ` |
|---|---|---|---|---|
| row1 | 6.033 | 170.31 | 0.689 | 36.63 |
| row4 | 5.927 | 156.93 | 0.034 | 0.98 |
| row7 | 5.603 | 204.30 | **2.145** | **83.06** |
| cde_n30_iter1 | 5.537 | 187.55 | 0.000 | 0.00 |

### 2.2. N=15000 (this section)

| Row | `Ω` (SDP value) | `Σ m·λ` | `Σ m³·λ` | `Σ m·σ` | `Σ m³·σ` | solve time |
|---|---|---|---|---|---|---|
| row1 | 0.380110 | 6.222 | 189.22 | 0.993 | 55.71 | 68s |
| row4 | 0.379778 | 6.170 | 181.70 | 0.159 | 9.60 | 66s |
| **row7** | **0.381329** | 5.813 | 227.41 | **3.661** | **298.59** | 55s |
| cde_n30_iter1 | 0.379794 | 5.717 | 191.02 | 0.000 | 0.00 | 60s |

### 2.3. C_explicit at N=15000

| Row | cos residual | sin residual | combined | C_explicit | Margin to UB 0.380871 |
|---|---|---|---|---|---|
| row1 | 6.516 × 10⁻⁴ | 1.039 × 10⁻⁴ | 7.555 × 10⁻⁴ | 0.380883 | `-1.2 × 10⁻⁵` (just vacuous) |
| row4 | 6.461 × 10⁻⁴ | 1.668 × 10⁻⁵ | 6.628 × 10⁻⁴ | 0.380791 | `+8.0 × 10⁻⁵` **NON-VACUOUS** |
| **row7** | 6.088 × 10⁻⁴ | 3.834 × 10⁻⁴ | **9.922 × 10⁻⁴** | **0.381120** | `-2.5 × 10⁻⁴` (vacuous) |
| cde_n30_iter1 | 5.987 × 10⁻⁴ | 3.4 × 10⁻¹⁰ | 5.987 × 10⁻⁴ | 0.380727 | `+1.4 × 10⁻⁴` **NON-VACUOUS** |

**Sup-over-rows at N=15000:** row7 with `C_explicit = 0.381120`, theorem
**vacuous** by `2.5 × 10⁻⁴`.

---

## 3. Key findings at N=15000

### 3.1. Cosine multipliers remain row-stable, slightly grow with N

| Row | Σ m·λ (N=3000) | Σ m·λ (N=15000) | Δ |
|---|---|---|---|
| row1 | 6.033 | 6.222 | `+3.1%` |
| row4 | 5.927 | 6.170 | `+4.1%` |
| row7 | 5.603 | 5.813 | `+3.8%` |
| cde_n30_iter1 | 5.537 | 5.717 | `+3.3%` |

**Σ m·λ grows mildly (~3-4%) with `N`**, not constant as predicted. This
prediction error is small: it adds `~25 µm/N` to the cosine residual at any
`N`, a `0.4%` underestimate. The empirical sup `Σ m·λ ≤ 6.22` (row1 at
N=15000) is essentially the same as the N=3000 sup `6.03`.

### 3.2. Sine multipliers are NOT row-stable AND grow significantly with N

| Row | Σ m·σ (N=3000) | Σ m·σ (N=15000) | Δ |
|---|---|---|---|
| row1 | 0.689 | 0.993 | `+44%` |
| row4 | 0.034 | 0.159 | `+367%` |
| row7 | **2.145** | **3.661** | **+71%** |
| cde_n30_iter1 | 0.000 | 0.000 | flat |

The sine multipliers grow with `N` (except for cde where they stay zero,
suggesting cde-like centers don't even activate the sine constraint). **The
break-even N prediction from Step C using N=3000 multipliers
(`16,378`) was optimistic** because Σ m·σ grew faster than the prediction
assumed.

The row7 Σ m·σ growth from 2.14 → 3.66 (70%) is the largest contributor to
the persistent vacuous-ness at N=15000.

### 3.3. Per-row picture: 2/4 rows non-vacuous at N=15000

- **row4, cde_n30_iter1:** NON-VACUOUS (margins +8.0e-5, +1.4e-4)
- **row1:** essentially at break-even (margin -1.2e-5)
- **row7:** firmly vacuous (margin -2.5e-4)

The theorem holds row-by-row for row4 and cde. The row-sup theorem (the
statement we want for the full saturation argument) remains **vacuous** at
N=15000 because row7's sine multiplier is too large.

### 3.4. Updated break-even estimates using N=15000 multipliers

Using row7's N=15000 measurement (Σ m·(λ+σ) = 9.47):

> `N_break_row7 ≥ π · 9.47 / (2 · 7.43 × 10⁻⁴) ≈ 20,026`.

So N ≈ 20,000 should make row7 (and thereby the sup-row theorem)
non-vacuous, IF the multipliers don't keep growing as N increases further.

---

## 4. N=20000 results (in progress)

(To be populated when solves complete.)

| Row | Ω | Σ m·λ | Σ m·σ | combined | C_explicit | Margin |
|---|---|---|---|---|---|---|
| row7 | TBD | TBD | TBD | TBD | TBD | TBD |
| row1 | TBD | TBD | TBD | TBD | TBD | TBD |
| row4 | TBD | TBD | TBD | TBD | TBD | TBD |
| cde_n30_iter1 | TBD | TBD | TBD | TBD | TBD | TBD |

---

## 5. Memory and timing observations

| `N` | bochner_n | T | Solve time per row | Peak RAM | Status |
|---|---|---|---|---|---|
| 3,000 | 20 | 1200 | ~20s | ~0.8 GB | done (prior session) |
| 15,000 | 20 | 1200 | 55-68s | ~1.6 GB | **DONE (4/4)** |
| 20,000 | 20 | 1200 | est ~90-100s | est ~2.1 GB | in progress |

Memory scaling is sublinear (1.6 GB at N=15000 vs 0.8 GB at N=3000, factor
~2× for 5× cell count). Solve time per row is roughly linear in N.

A 4-row N=15000 sweep took ~31 minutes total (with overheads). A 4-row
N=20000 sweep should take ~40 minutes.

---

## 6. Honest summary so far

The first-pass prediction from Step C — break-even at N=16,378 with C_explicit
= 0.380871 at that N — assumed multiplier-sum scale-invariance. The actual
N=15000 data shows the multiplier sums **grow** with N, especially for the
sine family. The corrected break-even estimate is `~20,000`.

The theorem becomes non-vacuous row-by-row at N=15000 for 2 of 4 rows tested
(row4 and cde_n30_iter1), but remains vacuous at the sup (driven by row7).
The N=20000 run will test whether the multiplier growth saturates, in which
case the sup-row theorem becomes non-vacuous, or whether it keeps growing,
in which case break-even continues to drift upward.

The cleanest publishable statement now is:

> **(scaled saturation, observed)** At N=15000, the cell-envelope cosine +
> sine residual is row-dependent and ranges from `5.99 × 10⁻⁴` (cde) to
> `9.92 × 10⁻⁴` (row7). For 2 of the 4 representative rows, the framework
> saturation theorem at N=15000 is **non-vacuous** (C_explicit < 0.380871);
> for the row7 sup, it remains vacuous by `2.5 × 10⁻⁴`.

(Final statement pending N=20000 results.)
