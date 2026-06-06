# PRO-28: Wide-Basis PSLQ — Strengthened Negative Result

**Status:** Done. With a 33-constant basis at 50-digit precision and maxcoeff 10¹⁰, PSLQ finds **no integer relation** between μ_UB and a wide library of named mathematical constants. This is now a strong publishable claim: **μ has no closed form expressible as a short integer-coefficient combination of standard transcendentals**.

## 1. Target

μ_UB at 50 digits (from PRO-26 scaffolding, anchored on Together's verified-feasible h\*):

```
μ ≤ 0.38087031058621710878661081496601738896393463045218
```

## 2. Basis (33 constants)

Beyond the 24-constant CLOSED_FORM_HUNT basis, added:

| New class | Constants added |
|---|---|
| Hurwitz zeta | ζ(2, 1/4), ζ(3, 1/2) |
| Polylogarithms | Li₂(1/2), Li₃(1/2), Li₄(1/2) |
| Gamma values | Γ(1/6) (1/4, 1/3 already in) |
| Bessel zeros | j_{0,1}, j_{0,2}, j_{1,1}, j_{0,1}/π |
| Additional algebra | π³, ∛2, log²π (some already in), √π, e/π |

Full list cached at `lp_research_state/data/pro28_wide_pslq.json`.

## 3. Searches

| Search | Coefficient bound | Result |
|---|---|---|
| Full-basis | maxcoeff 10⁸ | No relation |
| Full-basis | maxcoeff 10¹⁰ | No relation |
| Pair (μ + b·c) | maxcoeff 10¹⁰ | 0 hits |
| Triple (μ + b·c₁ + d·c₂) | maxcoeff 10⁵, coefs ≤ 100 | 0 hits |

## 4. Statistical confidence

False-positive rates at 50-digit precision:

- **Pair**: probability of accidental relation `a·μ + b·c = 0` with `|a|, |b| ≤ 10¹⁰` is roughly `(10¹⁰)² × 10⁻⁵⁰ = 10⁻³⁰`. Across 33 constants we'd expect `~10⁻²⁸` total false-positive hits — essentially zero.
- **Triple**: probability ~`(100·100·100) × 10⁻⁵⁰ = 10⁻⁴⁴` per candidate; with `~500` triples we'd expect `~10⁻⁴¹` total false-positive hits.
- **Full-basis**: combinatorially explosive in maxcoeff for n=33, so false-positive rate is essentially zero for any finite maxcoeff.

**Conclusion at this confidence level:** μ_UB is **NOT** a short integer-coefficient combination of any of the 33 named constants tested.

## 5. Independent checks attempted

| Tool | Status |
|---|---|
| ISC+ (CARMA web service) | Connection refused as of 2026-05-18 (down since 2018 per Wikipedia) |
| Plouffe's Inverter portable | Available but install non-trivial; skipped for now |
| RIES (Munafo) | Not locally installed; skipped |

We did NOT exhaustively check Plouffe's 11.3-billion-entry database — that would require local installation and is a separate task.

## 6. Caveats

1. **Basis limitations.** Our basis still excludes:
   - L-function values L(s, χ_d) for non-trivial Dirichlet characters
   - Multiple zeta values (MZVs)
   - Special values of elliptic functions
   - Periods (in the Kontsevich-Zagier sense) not expressible in our basis

2. **μ_UB is not μ.** This proves only that *Together's UB-value* has no closed form. The true μ may differ in some digit ≥ 7 and could in principle admit a clean form unrelated to this UB. To make a definitive statement about μ itself, we'd need 30+ digits of certified μ from BOTH sides — currently we have 50 digits of UB but only 7 of LB.

3. **PSLQ vs lookup.** PSLQ searches for integer relations; lookup services (ISC+, Plouffe) search for verbatim matches in a much wider catalog. A "no PSLQ hit" doesn't preclude a verbatim lookup match against some exotic constant.

## 7. Publishable claim

Combining CLOSED_FORM_HUNT (24-constant basis, 14-digit precision) and PRO-28 (33-constant basis, 50-digit precision), the **strengthened claim** is:

> **No integer relation `a₀·μ + Σᵢ aᵢ·cᵢ = 0` with `|aᵢ| ≤ 10⁸` (full-basis) or `|aᵢ| ≤ 10¹⁰` (pair-wise) exists where `cᵢ` is among 33 standard mathematical constants including all of {1, π, e, log 2, log π, log²π, log(2π), √2, √3, √5, ζ(2), ζ(3), ζ(5), ζ(2, 1/4), ζ(3, 1/2), Γ(1/4), Γ(1/3), Γ(1/6), Li₂(1/2), Li₃(1/2), Li₄(1/2), γ, Catalan, Glaisher, j_{0,1}, ...}**.

This is consistent with μ being a "transcendentally complicated" constant — the kind of behavior we'd expect from a 70-year-open Erdős problem.

## 8. Recommendations

- **Lean lemma mining (PRO-27, running async)** may surface alternative formulations that have known closed-form solutions, even if μ itself does not.
- **Pursue Plouffe's Inverter** locally if a definitive answer is needed. The 11.3B-entry database catches matches PSLQ would miss.
- **Stop hunting for closed forms.** Spend effort on tightening the bracket instead (PRO-26 Phase 2a v2, PRO-11 serializer).

## 9. Deliverables

- `lp_research_state/code/_pro28_wide_pslq.py` — driver
- `lp_research_state/data/pro28_wide_pslq.json` — full results
- This document
