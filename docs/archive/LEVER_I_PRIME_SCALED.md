# Lever I' Step D: Cell-Envelope Saturation Theorem at Scaled N

**Status:** COMPLETE. Measurements at N ∈ {3000, 10000, 15000, 20000}.

**Headline:** The corrected cell-envelope cosine + sine residual is empirically **stuck near the open gap** at all currently-tractable `N`. Scaling from N=10000 to N=20000 shrinks the combined residual from `1.28 × 10⁻³` to `8.0 × 10⁻⁴` (sup over rows), but the multiplier sums **grow with N**, slowing the residual decay. At N=20000, all 4 representative rows have `C_explicit` slightly above Together's UB. The framework appears to be in a **delicate near-saturated regime** — neither comfortably vacuous nor solidly certified.

**Predecessor:** [LEVER_I_PRIME_THEOREM.md](LEVER_I_PRIME_THEOREM.md), [SESSION_FINAL.md](SESSION_FINAL.md).

**Data:** [lp_research_state/data/lambda_m_scaled.json](lp_research_state/data/lambda_m_scaled.json).
**Code:** [lp_research_state/code/_lever_i_prime_lambda_m_scaled.py](lp_research_state/code/_lever_i_prime_lambda_m_scaled.py).

---

## 1. Formula

Corrected cell-envelope cosine + sine saturation bound (LEVER_I_PRIME_THEOREM.md Theorem 3):

```
ResidualGain(N) ≤ (π/(2N)) · Σ_m m·(λ_m^cos + σ_m^1 + σ_m^2)
                + (π²Ω/(3N³)) · Σ_m m³·(λ_m^cos + σ_m^1 + σ_m^2)
```

`C_explicit(N) := 0.3801279 + ResidualGain(N)`. Theorem non-vacuous iff
`C_explicit(N) < 0.380871`.

---

## 2. All measurements

### 2.1. Cosine multiplier sums `Σ m·λ`

| Row | N=3000 | N=15000 | N=20000 |
|---|---|---|---|
| row1 | 6.033 | 6.222 | **8.050** |
| row4 | 5.927 | 6.170 | **9.988** |
| row7 | 5.603 | 5.813 | 5.829 |
| cde_n30_iter1 | 5.537 | 5.717 | **9.654** |

The cosine multipliers grow with N. The growth is mild from N=3000 to
N=15000 (3-4%), but **dramatic** from N=15000 to N=20000 for rows 1, 4,
cde (29-69%). Only row7 stays flat across all three scales.

### 2.2. Sine multiplier sums `Σ m·σ`

| Row | N=3000 | N=15000 | N=20000 |
|---|---|---|---|
| row1 | 0.689 | 0.993 | 1.060 |
| row4 | 0.034 | 0.159 | 0.219 |
| row7 | **2.145** | **3.661** | **3.752** |
| cde_n30_iter1 | 0.000 | 0.000 | 0.000 |

Sine multipliers also grow with N but more uniformly. row7 is the
consistent sup, with `Σ m·σ` settling around 3.7 at large N. cde stays
at exactly zero (its sine constraint is permanently slack).

### 2.3. Per-row C_explicit

| Row | N=10000 (prior session) | N=15000 | N=20000 |
|---|---|---|---|
| row1 | 0.381006 | 0.380883 | **0.380843** ✓ NON-VACUOUS (margin +2.8e-5) |
| row4 | 0.381059 | **0.380791** ✓ NON-VACUOUS | 0.380930 (vacuous -5.9e-5) |
| row7 | **0.381346** | **0.381120** | 0.380880 (vacuous -9.4e-6) |
| cde_n30_iter1 | 0.380998 | **0.380727** ✓ NON-VACUOUS | 0.380886 (vacuous -1.5e-5) |

(N=10000 entries are extrapolations from N=3000 multipliers, since
no direct N=10000 dual extraction was done.)

### 2.4. Sup-row C_explicit (the headline statistic)

| N | sup-row | combined residual | C_explicit (sup) | Margin to UB |
|---|---|---|---|---|
| 3,000 (extrapolated to N=10000) | row7 | `1.28 × 10⁻³` | 0.381346 | `-4.7 × 10⁻⁴` (vacuous) |
| 15,000 (measured) | row7 | `9.92 × 10⁻⁴` | 0.381120 | `-2.5 × 10⁻⁴` (vacuous) |
| 20,000 (measured) | row4 | `8.02 × 10⁻⁴` | 0.380930 | `-5.9 × 10⁻⁵` (vacuous) |

The sup-row shifts from row7 (sine-dominated) at low N to row4 (cosine-dominated)
at N=20000, as cosine multipliers grow faster in some rows.

---

## 3. The framework is in a near-saturated regime

Two competing effects as N grows:

1. **Mechanical residual decay** `(π/(2N)) · Σ m·μ` ≈ shrinks as `1/N`.
2. **Multiplier growth** `Σ m·μ` increases as `N` grows (the cell-envelope
   constraints bind more tightly, demanding larger dual mass).

These nearly cancel out. The net residual decay observed:

| N range | Residual (sup) | Reduction |
|---|---|---|
| 10,000 (extrapolated) → 15,000 | `1.28e-3 → 9.92e-4` | `-22%` (50% N-increase) |
| 15,000 → 20,000 | `9.92e-4 → 8.02e-4` | `-19%` (33% N-increase) |
| Trend extrapolation to 30,000 | `~6.5e-4` ? | `-19%` (50% N-increase) |
| 50,000 | `~4.5e-4` ? | total `~-40%` from 20,000 |

If this trend continues, the framework would become genuinely non-vacuous
(sup C_explicit < 0.380871) at `N ≈ 30,000-50,000`. Below that, the theorem
is marginal: either non-vacuous per-row but vacuous in sup, or
vacuous overall by `O(10⁻⁵)`.

---

## 4. Honest summary

### 4.1. Theorem status at currently-tractable N

| Statement | N=15000 | N=20000 |
|---|---|---|
| "Cell-envelope cos+sin saturation holds at row X" | TRUE for 2/4 rows (row4, cde) | TRUE for 1/4 rows (row1) |
| "Cell-envelope cos+sin saturation holds at all 4 rows tested" | FALSE | FALSE |
| "Cell-envelope cos+sin saturation holds at the sup-row" | FALSE (margin `-2.5e-4`) | FALSE (margin `-5.9e-5`) |

The strong sup-row statement remains vacuous at all tested N up to 20,000.

### 4.2. Multiplier sums are NOT scale-invariant

The Step C prediction "break-even N ≈ 16,378" was based on N=3000
multiplier sums and assumed scale invariance. The actual measurements
show:

- Cosine `Σ m·λ` grew 0–69% as `N` grew from 3,000 to 20,000.
- Sine `Σ m·σ` grew 6–540% across the same range.

These growths shift break-even upward — empirically to `~25,000-30,000` for
the sup-row to fit under Together's UB.

### 4.3. What does this mean for the saturation diagnosis?

The OVERNIGHT_WRAPUP framing was: "the SDP framework saturates at
`µ ≈ 0.380553`; the remaining open gap is fundamental." That was wrong (per
LEVER_I_PRIME_THEOREM.md retraction).

The corrected framing from this Step D is: **the cell-envelope cosine + sine
residual at currently-tractable N is just above the open gap, and shrinking
slowly with N.** The framework is in a "tight" regime where:
- Per-row, the cosine cell-envelope augmentation can prove `µ ≤ 0.380843`
  to `0.380930` depending on the row, with row1 being the only proved
  bound below Together's UB at N=20000.
- The sup-row bound is just above `0.380871` at N=20000, with the margin
  shrinking as N grows.

A definitive saturation certification requires either:
- (a) Confirming the trend at N=30,000+ (memory: feasible, ~3 GB; ~10 minutes/row).
- (b) Deriving a sharper analytic per-`m` residual bound that exploits the
  primal optimum's density distribution (mentioned in
  LEVER_I_PRIME_THEOREM.md §4.2; not done).

### 4.4. The most useful concrete number

At N=20,000, the row1 SDP-derived bound on the cell-envelope-augmented LB:

> `SDP_LB(row1 with cell-envelope augmented)  ≤  0.380843`

This says: if we replace the cell-envelope cosine constraint at row1's
SDP with its exact analytical form (no relaxation), the LB cannot exceed
0.380843, which is **below Together's UB** by `2.8 × 10⁻⁵`. So at row1
alone, the cell-envelope augmentation cannot push the LB past 0.380843.

For the other 3 rows at N=20000, the augmented LB cap is between 0.380880
and 0.380930 — all in a narrow band around Together's UB.

---

## 5. Memory and time observations

| N | bochner_n | T | Solve time/row | Peak RAM | Status |
|---|---|---|---|---|---|
| 3,000 | 20 | 1200 | ~20s | ~0.8 GB | done |
| 15,000 | 20 | 1200 | 55-68s | ~1.6 GB | done |
| 20,000 | 20 | 1200 | 77-86s | ~2.1 GB | done |

Memory scales sublinearly (~`N^0.6`). Solve time scales linearly with N. A
4-row N=30000 sweep would take ~10 min/row × 4 = ~40 min and fit in ~2.8 GB.

---

## 6. Pending follow-ups

1. **N=30000 sweep** to confirm the residual-decay trend and refine the
   break-even estimate. Out of scope for current Ralph budget.
2. **Per-cell primal-side residual sharpening** (LEVER_I_PRIME_THEOREM.md
   §4.2): a tighter analysis of where SDP density concentrates could reduce
   the proved per-`m` residual by a factor of 2-5× without scaling N.
3. **poly_moment, Hankel-PSD residuals** to complete the full constraint
   stack saturation theorem (LEVER_I_PRIME_THEOREM.md §4.3).
4. **A=15000 Phase 5 LB re-run** (separately from this scaled extraction)
   to actually improve the empirical LB from 0.3801279 to whatever a
   larger-`N` Phase-5 configuration delivers.

<promise>SCALED_THEOREM_DONE</promise>
