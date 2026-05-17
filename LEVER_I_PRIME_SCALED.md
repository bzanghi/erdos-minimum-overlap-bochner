# Lever I' Step D: Cell-Envelope Saturation Theorem at Scaled N

**Status:** IN PROGRESS — first solves at N=15000 running.

**Goal:** Test the corrected cell-envelope cosine + sine saturation bound at
larger `N` to see if (a) the multiplier sums `Σ m·λ` (cos) and `Σ m·σ` (sin)
remain near the N=3000 measured values, (b) C_explicit drops below
Together's UB at any tractable `N`.

**Predecessor:** [LEVER_I_PRIME_THEOREM.md](LEVER_I_PRIME_THEOREM.md), [SESSION_FINAL.md](SESSION_FINAL.md).

**Data:** [lp_research_state/data/lambda_m_scaled.json](lp_research_state/data/lambda_m_scaled.json) (incremental).
**Code:** [lp_research_state/code/_lever_i_prime_lambda_m_scaled.py](lp_research_state/code/_lever_i_prime_lambda_m_scaled.py).

---

## Formula recap

Per Theorem 3 of LEVER_I_PRIME_THEOREM.md, the corrected cell-envelope
residual at scale `N` is

```
ResidualGain ≤ (π/(2N)) · Σ_m m·λ_m^cos
             + (π²Ω/(3N³)) · Σ_m m³·λ_m^cos
             + (π/(2N)) · Σ_m m·(σ_m^1 + σ_m^2)
             + (π²Ω/(3N³)) · Σ_m m³·(σ_m^1 + σ_m^2)
```

with `Ω ≈ 0.38` and `(λ, σ)` the cell-envelope cosine/sine dual multipliers.

`C_explicit(N) := SDP_LB(Phase5) + ResidualGain(N)`. With `SDP_LB(Phase5) =
0.3801279`, the theorem is **non-vacuous at scale `N`** iff `C_explicit(N) <
0.380871`.

---

## Measurements at N=3000 (baseline; from prior session)

| Row | `Σ m·λ` (cos) | `Σ m³·λ` | `Σ m·σ` (sin) | `Σ m³·σ` |
|---|---|---|---|---|
| row1 | 6.033 | 170.31 | 0.689 | 36.63 |
| row4 | 5.927 | 156.93 | 0.034 | 0.98 |
| row7 | 5.603 | 204.30 | **2.145** | **83.06** |
| cde_n30_iter1 | 5.537 | 187.55 | 0.000 | 0.00 |
| **sup** | **6.033** | **204.30** | **2.145** | **83.06** |

**Open question:** Do these sums change at larger `N`? They COULD change
because the cell-envelope relaxation gets tighter as `L = 2/N` shrinks, which
shifts the SDP optimum and rebalances the multipliers.

---

## Measurements at N=15000 (in progress)

(Will be populated as solves complete.)

| Row | `Σ m·λ` | `Σ m³·λ` | `Σ m·σ` | `Σ m³·σ` | Status |
|---|---|---|---|---|---|
| row7 | TBD | TBD | TBD | TBD | running |
| row4 | TBD | TBD | TBD | TBD | pending |
| row1 | TBD | TBD | TBD | TBD | pending |
| cde_n30_iter1 | TBD | TBD | TBD | TBD | pending |

---

## Predicted C_explicit (using N=3000 sup-row values)

If the multiplier sums are scale-invariant (Σ m·λ_max stays at 6.03, Σ m·σ_max
stays at 2.145), then:

| `N` | Cosine residual | Sine residual | Combined | C_explicit | Non-vacuous? |
|---|---|---|---|---|---|
| 10,000 | `9.48e-4` | `3.37e-4` | `1.28e-3` | 0.381346 | NO (`+4.7e-4`) |
| 15,000 | `6.32e-4` | `2.25e-4` | `8.57e-4` | 0.380985 | NO (`+1.1e-4`) |
| 16,378 | `5.78e-4` | `2.06e-4` | `7.84e-4` | 0.380871 | break-even |
| 20,000 | `4.74e-4` | `1.69e-4` | `6.42e-4` | 0.380770 | YES (`-1.0e-4`) |

**Hypothesis to test:** at the measured `N=15000`, will the combined residual
match the prediction `8.6e-4` (vacuous) or will the multiplier sums drop?

If `Σ m·λ` drops to (say) 4.5 at N=15000 with Σ m·σ to 1.5, the residual
would be `(π/30000)·6 = 6.3e-4` and the theorem becomes non-vacuous at
N=15000 directly.

---

## Memory and timing observations

| `N` | bochner_n | T | Solve time (s) | Peak RAM (GB) | Status |
|---|---|---|---|---|---|
| 3,000 | 20 | 1200 | ~20 | ~0.8 | done (prior session) |
| 15,000 | 20 | 1200 | ≥ 100 | ≥ 1.6 | running |

---

## (To be completed when solves return)

- Per-row measurements at N=15000
- Per-row measurements at N=20000 (if memory allows)
- Comparison: do Σ m·λ, Σ m·σ scale with N, or stay constant?
- Final C_explicit at each tested N
- Yes/no answer: is the cos+sin saturation theorem non-vacuous at any
  tractable N?
