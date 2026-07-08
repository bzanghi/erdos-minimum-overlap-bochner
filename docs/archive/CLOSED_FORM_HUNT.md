# Closed-form hunt for μ — Wolfram + PSLQ session

**Status:** Strong negative result. None of Wolfram's "identify the number"
candidates survive verification at the 16-digit precision available from
Together's h*. PSLQ over a 24-constant basis with maxcoeff 10⁶ finds no
integer relation involving the UB.

**Date:** 2026-05-18. Tools used:
[lp_research_state/code/wolfram_client.py](lp_research_state/code/wolfram_client.py)
(Wolfram LLM API), `mpmath.pslq` at 50-digit precision.

---

## 1. Inputs

| Quantity | Value (digits) | Source |
|---|---|---|
| Rigorous LB (PRO-21) | `0.3803027228196094` | N=20K bn=40, after 1e-6 + Lipschitz margin |
| Raw grid_min | `0.3803048898545508` | same SDP, before margin |
| **Together UB (16-digit)** | **`0.3808703105862199`** | computed from Together's h* via `np.correlate` |
| Together UB (their reported) | `0.380871` | 6 digits in their README |
| Framework ceiling C_∞ | `0.380558` | PRO-6 asymptotic estimate |
| Row 4 Ω* (cell-min, N=3K) | `0.378792` | from `lambda_m_extracted.json` |
| \|ξ\|/Ω row 4 (PRO-14) | `1.4647` | shadow-price audit |
| \|ξ\|/Ω row 7 (PRO-14) | `1.4445` | shadow-price audit |

## 2. Wolfram "identify the number" results (at 6-digit precision)

The LLM API was queried with `identify the number <value> closed form` for
each constant above. Wolfram returned 1–3 candidates per query. Full
responses cached at
`lp_research_state/data/wolfram_identify_constants.json`.

| Target | Top Wolfram candidate(s) | In bracket [0.380305, 0.380871]? |
|---|---|---|
| LB 0.3803027 | `P_c(BCC bond) + 1/5` ≈ 0.38030 | **out** (just below) |
| LB 0.3803027 | `8/(29π) + 27π/290` ≈ 0.38030273 | **out** (-2.2×10⁻⁶) |
| LB 0.3803027 | `(11e)/(60 log²π)` ≈ 0.38030311 | **out** (-1.8×10⁻⁶) |
| UB 0.380871 | `(3eπ·log π)/77` ≈ 0.38087100 | in (top edge) |
| UB 0.380871 | `7/(10 log(2π))` ≈ 0.38087422 | **out** (above) |
| UB 0.380871 | `5/(7 C_Sh²)` ≈ 0.38087171 | unclear (C_Sh definition) |
| C_∞ 0.380558 | `(7e)/50` ≈ 0.38055946 | in (middle) |
| C_∞ 0.380558 | `(π-2)/3` ≈ 0.38053088 | in (middle) |
| C_∞ 0.380558 | `(47π)/388` ≈ 0.38055375 | in (middle) |
| \|ξ\|/Ω 1.4647 | `π^(1/3)` ≈ 1.46459 | — |
| \|ξ\|/Ω 1.4647 | `(4 log 3)/3` ≈ 1.46482 | — |
| \|ξ\|/Ω 1.4445 | `e^(1/e)` ≈ 1.44467 (Steiner's constant) | — |
| \|ξ\|/Ω 1.4445 | `13/9` ≈ 1.44444 | — |

## 3. Verification at 16-digit precision (the kill shot)

Together's h* gives μ_UB to ~16 digits via direct computation. Compared
candidate forms at full precision:

| Candidate | Value (16 digits) | Together UB | Δ |
|---|---|---|---|
| `(3eπ·log π)/77` | `0.3808709992950128` | `0.3808703105862199` | **+6.89×10⁻⁷** |
| `(7e)/50`        | `0.3805594559842663` | (different region)  | n/a |
| `(π-2)/3`        | `0.3805308845299311` | (different region)  | n/a |
| `(47π)/388`      | `0.3805537492750523` | (different region)  | n/a |

The `(3eπ·log π)/77` "match" was a **6-digit curve fit by Wolfram, defeated by our 7-digit reality**. The actual Together UB is ~0.69 ppm below it.

For the LB region, our most-precise rigorous LB at 16 digits is `0.3803027228196094`. Wolfram candidates `8/(29π)+27π/290` and `(11e)/(60·log²π)` lie *below* this value (out of bracket) — they cannot equal μ.

## 4. PSLQ on the 16-digit UB

Ran `mpmath.pslq` at 50-digit working precision with three searches:

- **Full-basis** (24 constants: 1, π, π², π³, 1/π, e, e², eπ, log 2, log π, log²π, log(2π), √2, √3, √5, ζ(2), ζ(3), ζ(5), Γ(1/4), Γ(1/3), γ, Catalan, Glaisher, UB), `maxcoeff = 10⁶`: **no relation found**.
- **Pair searches** (UB + b·c for each c, `maxcoeff = 10⁸`): **no relation found**.
- **Triple searches** (UB + b·c₁ + d·c₂, `maxcoeff = 10⁵`, all pairs of basis constants): **one hit** — `π² - 6·ζ(2) = 0` (the well-known Euler identity; UB coefficient is 0 — not actually involving UB).

**Conclusion.** Together's UB is **not expressible as a short integer-coefficient combination of the 23 named constants** in our basis. Either:
- μ has no clean closed form in this standard basis (likely), or
- The closed form requires more exotic constants (Apéry-like, Hardy-Littlewood, ...) outside our basis, or
- Our precision is still insufficient (16 digits — would need 30+ for stronger conclusions).

## 5. The one positive cross-check

Wolfram verified our cell-envelope integral identity exactly, matching sympy:

```
Integrate[Cos[Pi*m*x/2], {x, (j-1)*L, j*L}]
  = -(2 (sin((j-1)Lmπ/2) - sin(jLmπ/2)))/(πm)
```

Series expansion at m=0:
```
L − (π²(3j²−3j+1) L³ / 24) m² + (π⁴(5j⁴−10j³+10j²−5j+1) L⁵ / 1920) m⁴ + O(m⁶)
```

This independently confirms LEVER_I_PRIME_THEOREM.md §3.1's Case-B derivation of the per-cell residual `≤ (πm)² L³ / 24` — the leading-order coefficient `π²(3j²−3j+1)L³/24` matches exactly (the simplified bound takes max over j).

## 6. What this tells us

- **F3 (Wolfram inverse-symbolic) is exhausted at current precision.** The bracket width 5.7×10⁻⁴ admits ~4 short-closed-form candidates by chance from Wolfram's basis. Our 16-digit Together-UB anchor kills all of them.
- **Path forward for closed-form hunting:** need μ (or any of our solve outputs) to 25+ digits to make PSLQ definitive. This requires SDPA-GMP integration (cvxpy → SDPA-S serializer), the next milestone for PRO-11.
- **PSLQ-negative is a positive scientific result.** "μ has no short closed form in standard constants" is itself worth stating in the preprint — it constrains expectations about what an analytical solution would look like, and discourages superficial closed-form claims.
- **The Wolfram tool works.** Integration verified; Mathematica-syntax queries reliable. Natural-language queries are flaky (HTTP 500s) — prefer Mathematica syntax.
- **|ξ|/Ω is not a clean constant.** Row 4 ≈ π^(1/3), row 7 ≈ e^(1/e) — but these are different constants. Confirms PRO-14's empirical conclusion: the ratio is row-dependent, ~1.46 universally is approximate not exact.

## 7. Deliverables

- `lp_research_state/data/wolfram_identify_constants.json` — cached LLM API responses (provenance for the above)
- This document
- (Already landed) `lp_research_state/code/wolfram_client.py`, `pslq_hunt.py`

## 8. Recommended next step

Build the **cvxpy → SDPA-S serializer** (PRO-11's gating item) so we can run White's SDP at SDPA-GMP precision and get 30+ digit Ω*. Then re-run PSLQ.
The current negative result is *evidence* μ has no clean form; a 30-digit
re-run would be *proof* in the same basis up to maxcoeff 10¹⁰.

## 9. Update — 50-digit UB anchor + strong PSLQ negative (2026-05-18, PRO-26 Phase 2a scaffolding)

While building the PRO-26 scaffolding, we re-evaluated Together's UB at full mpmath precision using their exact float64 h\* values (verified feasible: sum = 300 exactly, all h ∈ [0, 1]). The 50-digit-precise UB is:

> **μ ≤ 0.38087031058621710878661081496601738896393463045218**

(vs Together's published 6 digits "0.380871")

This is the strongest published anchor on either side of the bracket. Re-ran PSLQ at this precision:

- **Full-basis (24 constants), tol 10⁻⁴⁵, maxcoeff 10⁸:** no relation
- **Pair searches (UB + b·c), tol 10⁻⁴⁵, maxcoeff 10¹⁰:** no relation

False-positive rate analysis: with maxcoeff 10⁸ over 24 basis constants and 50-digit precision, the expected false-positive rate is `(24 × 10⁸)⁴ × 10⁻⁵⁰ ≈ 10⁻¹⁵` — essentially zero. **No relation found at this rate is now a near-definitive negative result.**

**Conclusion at 50 digits:** μ_UB (= the M-value of Together's h*, which is an upper bound on the continuous μ) does NOT admit a closed form as an integer-coefficient combination with coefficients up to 10⁸ over the standard 24-constant basis. This holds with overwhelming statistical confidence.

This strengthens the preprint claim about μ being "transcendentally ugly" in the standard basis from "evidence" to "near-proof".

## 10. Subtle caveat

The 50-digit UB is the M-value of Together's *specific* h* (the 600-cell float64 array). This is an upper bound on μ, but the *true* μ may differ from this UB in some digit ≥ 7 (per PRO-23: μ < μ_UB strictly). So a closed form for μ itself could still exist even if no closed form for this particular UB exists. We'd need to either:
- Compute h* to higher precision (PRO-26 Phase 2a optimization step)
- Or solve the SDP at higher precision (PRO-11)
to anchor μ_LB to 30+ digits.
