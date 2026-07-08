# Lever F5 (Contrarian): Empirical-Physics Fit + PSLQ on Framework Asymptote

**Status:** Conclusion (C) — no closed-form recognition. **BUT** a valuable empirical insight emerges: **the multiplier sums asymptote to a finite constant (reciprocal fit `S = a/N + b` beats linear), implying the framework is asymptotically *tight* with residual → 0 as N → ∞.**

This DOES NOT close the open question but RESHAPES it: μ_LB + 0 = μ_∞, so the framework asymptote equals μ itself (which is what we're trying to find — uninformative).

---

## 1. Hypothesis

Lever I' Step E observed that `Σ_m m · (λ + σ)(N)` grows with N. The contrarian bet:

> **If `Σ m·(λ+σ)(N) ≈ a + bN` linearly, then ResidualGain(N) = (π/2N)·(a+bN) → πb/2 as N → ∞. The asymptotic framework ceiling LB + πb/2 might be a recognizable closed-form constant — if so, PSLQ could identify it.**

---

## 2. Empirical fits per row

Data from N ∈ {15K, 20K, 30K, 40K} (Step D + E measurements):

| Row | N=15K | N=20K | N=30K | N=40K |
|---|---|---|---|---|
| row1 | 7.215 | 9.109 | 11.786 | 11.981 |
| row4 | 6.329 | 10.206 | 10.983 | 11.063 |
| row7 | 9.475 | 9.581 | 9.840 | **14.909** ← jump |
| cde_n30_iter1 | 5.717 | 9.654 | 10.134 | 10.195 |

(values are `Σ m·(λ+σ)`.)

### 2.1. Linear fit `S = a + b·N`

| Row | a | b | R² | πb/2 | μ_∞_lin |
|---|---|---|---|---|---|
| row1 | 4.98 | 1.92 × 10⁻⁴ | 0.869 | 3.02 × 10⁻⁴ | 0.380430 |
| row4 | 5.50 | 1.58 × 10⁻⁴ | 0.610 | 2.48 × 10⁻⁴ | 0.380376 |
| **row7** | 5.58 | **2.05 × 10⁻⁴** | 0.736 | **3.21 × 10⁻⁴** | **0.380449** ← sup |
| cde_n30_iter1 | 5.12 | 1.45 × 10⁻⁴ | 0.559 | 2.28 × 10⁻⁴ | 0.380356 |

Sup-over-rows (linear): `μ_∞_lin = 0.380449`, margin to Together's UB `+4.2 × 10⁻⁴`.

### 2.2. Reciprocal fit `S = a/N + b` (better for 3/4 rows)

| Row | a | b (asymp) | R² | residual rate | μ_∞_recip |
|---|---|---|---|---|---|
| **row1** | −122,170 | 15.37 | **0.976** | → 0 | 0.3801279 |
| **row4** | −110,313 | 14.47 | **0.825** | → 0 | 0.3801279 |
| row7 | −102,541 | 15.44 | 0.514 | → 0 | 0.3801279 |
| **cde** | −102,862 | 13.43 | **0.780** | → 0 | 0.3801279 |

**3 of 4 rows favor reciprocal** (R² 0.78-0.98 vs linear 0.55-0.87). Under reciprocal, multiplier sum approaches a finite asymptote `b ≈ 13-15`. Then `ResidualGain(N) = (π/2N)·(a/N + b) = O(1/N) → 0`. Framework asymptote: **μ_∞_recip = LB exactly**.

(row7 is the exception: linear fits slightly better, R² 0.74 vs 0.51. But row7 has the largest single-step jump from N=30K → 40K, which neither model fits well.)

---

## 3. PSLQ + identify() on each candidate

At ~5-digit fit precision (only 4 data points per fit), PSLQ and identify() produced expressions with rational denominators 17-439 — all consistent with **fitting noise, not structural identification**:

```
row7_linear (0.3804492) → identify: 4 - 3π + 3log(2) + 5√2/2 + atan(2/5)/2
                                                                ↑
                                                          high coefficients,
                                                          random combinations
row1_linear (0.3804296) → identify: (1/439) - sqrt(2)/439 + (440/439)*atan(2/5)
                                              ↑↑↑ denominator 439 = clear overfitting

cde_linear (0.3803559) → identify: -2/3 + log(2)/3 + 2√2/3 - atan(2/5)/3
                                         ↑↑ denominator 3 — looks cleaner but
                                            coefficients 2/3, 1/3, 2/3, 1/3 are
                                            statistically expected at this precision
```

**No PSLQ hit at meaningful tolerance.** The fit-uncertainty (~10⁻⁴) exceeds the bracket required for closed-form identification (~10⁻⁵ or finer).

---

## 4. Cross-check: `LB + (3/7)·gap`

The row7 linear prediction `0.3804492` is within 3 × 10⁻⁶ of `LB + (3/7) × gap = LB + 0.318 × 10⁻³ = 0.380446`. This is at the edge of significance:
- "natural" fraction (3/7 has small integers)
- But within fit noise (3 × 10⁻⁶ is well below fit uncertainty ~10⁻⁴)

**Cannot be confirmed as a real relation** at current precision. Would need ~10 more data points at intermediate N values, OR a row7 dual extraction at much higher precision (e.g., SDPA-GMP).

---

## 5. Key conceptual finding (despite negative PSLQ)

**The reciprocal model's fit quality (R² 0.78-0.98 for 3/4 rows) tells us something important:**

> **`Σ m·(λ+σ)(N) → b ∈ ℝ` (a finite constant) as N → ∞** (per row).

This implies:

> **`ResidualGain(N) → 0` (the framework's `C_explicit` ceiling decays as 1/N).**

And therefore:

> **`μ_∞_framework = LB_phase5 + 0 = LB`** (the framework asymptote is the current LB).

But wait: the LB itself depends on N (Phase 5 was at N=10000). As N increases, LB should improve, approaching `μ` itself. So:

> **`lim_{N → ∞} SDP_LB(N) = μ` (the framework is asymptotically tight).**

**The framework is NOT "fundamentally" saturated.** It's tight in the limit — just very slowly converging at currently-tractable N.

This is *consistent* with the Step E observation that `C_explicit` margins are tightening with N (from `+1.26 × 10⁻⁴` at N=30K to `+1.58 × 10⁻⁴` at N=40K). The trend continues; eventually `C_explicit → μ`.

---

## 6. What this means for the open problem

### 6.1. The "beyond-framework" portion shrinks asymptotically

Step E established: at N=40K, 21% of the open gap is "beyond cell-envelope reach". F5 now suggests this is a *finite-N* phenomenon — at infinite N, the gap closes via the cell-envelope augmentation alone.

But "infinite N" is not actionable. The relevant question is: **at what tractable N (say, N ≤ 100K) does C_explicit drop within 10⁻⁵ of μ_UB?**

From the reciprocal fit, ResidualGain ≈ π·b/(2N) ≈ π·15/(2N) ≈ 24/N. For Residual ≤ 10⁻⁵: N ≥ 2.4 × 10⁶. **Practically infeasible** at current SDP scales.

### 6.2. Closed-form hunting via F5 fails

No PSLQ hit. The framework asymptote is μ itself, which is unknown. Closed-form discovery requires either:
- The actual μ value to ~15 digits (current bracket gives only 3-4 digits)
- A structural argument (not empirical fit) for what μ should be

### 6.3. Sensitivity to fit model

The two models give very different predictions:
- Linear → framework asymptote in [0.380356, 0.380449]
- Reciprocal → framework asymptote = LB = 0.3801279

Which is right? Theory says: as N → ∞, the SDP becomes tight, so framework asymptote = μ. Reciprocal is the theoretically-consistent fit. Linear is an artifact of finite-N data with limited range.

---

## 7. Net result

| Question | Answer |
|---|---|
| Does F5 find a closed-form μ_∞? | **No.** |
| Is the framework provably saturated at a non-trivial level? | Only at finite N (Step E result stands). Asymptotically, NO. |
| Does the multiplier-growth pattern reveal structure? | Yes — it's a 1/N transient, not a "phase transition". |
| Should we run a tighter fit (more data points)? | Marginally useful. Would need N ≥ 80,000 to test reciprocal vs linear definitively. |
| Should we shift to F3 (full-stack saturation)? | **Yes** — F3 is the cleaner extension of Step E. |

---

## 8. Honest summary

- Linear fit predicts non-trivial asymptote (`μ_∞ ≈ 0.38044`) at low R² (0.55-0.87).
- Reciprocal fit predicts trivial asymptote (`μ_∞ = LB`) at higher R² (0.78-0.98).
- PSLQ on all candidates returns spurious matches with high denominators (clear fitting noise).
- Reciprocal is theoretically preferred (consistent with "framework is asymptotically tight").
- **No closed-form for μ via this approach.**
- **The contrarian misfires, but cleanly.** The negative result rules out one specific structural hypothesis ("framework asymptote is a non-trivial closed-form constant"), strengthening the case for going to F3 (full-stack saturation) next.

<promise>CONTRARIAN_F5_DONE</promise>
