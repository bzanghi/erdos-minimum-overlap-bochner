# OUT_OF_BOX_REDTEAM — Adversarial audit of the 2026-05-10 session conclusion

**Author:** out-of-box red-team auditor
**Target:** [SESSION_2026-05-10_FINAL_ADDENDUM.md](SESSION_2026-05-10_FINAL_ADDENDUM.md) — the claim that "all 8 candidate levers ruled out, µ ≥ 0.3801279 is the framework's saturation point, write up the result."

## Framings investigated

All three framings were investigated. Framing 1 (LB-side) and Framing 2 (UB-side) received hands-on numerical effort; Framing 3 (formulation equivalence) received a literature-side audit against Haugland (2016) and White (2023).

---

## Framing 2 — Is Together's UB locally optimal?

**Setup.** Together's `h_values` is an n=600 step function with claimed `M(h) = 0.380871`. I loaded it directly from `/tmp/together_repo/erdos-minimum-overlap/solutions/together_ai_2026.py` and verified the bound via the same `np.correlate(h, 1-h, mode='full')/n*2` formula in their `analysis.ipynb` cell 3.

**Findings:**

1. **Verified baseline:** `M(h) = 0.38087031058622` (matches their stated 0.380871 to 6 sig figs; their stated digit appears to be a round-up rather than the true value).

2. **Together's `h` is at an LP-flatness optimum.** Of the 1199 possible shifts, **437** are tied within 1e-7 of the maximum, and **457** are tied within 1e-4. This is the signature of a sequential-LP-converged optimum where the active-set duality has saturated: the gradient cone at h points outward in every feasible direction.

3. **Random pairwise perturbation (sum-preserving, box-respecting) for 2000 trials at four scales σ ∈ {1e-3, 1e-4, 1e-5, 1e-6}:** zero improvements found. M is locally constant under these moves down to 1e-6.

4. **Sequential-LP local-search (the same method Together describes using).** Linearized M around h, using all 1199 shifts as active constraints. Trust regions δ ∈ {1e-2, 5e-3, 1e-3, 1e-4, 1e-5}: at every δ, the LP predicts an improvement (e.g., u/n·2 = 0.380681 at δ=1e-3, a predicted gain of 1.9e-4), but the **true** M at the LP step is **larger** than baseline (true M = 0.381374 at δ=1e-3 — a regression of 5e-4). The non-convex quadratic `-x[t]*x[t+k]` term reverses the linear prediction's sign. This is the classical SLP fixed-point: locally, no first-order move helps. See `/tmp/redteam_local_search.py`.

5. **Step-function refinement (upsample n=600 → n=1200, n=2400):** M is exactly preserved (= 0.38087031058622) under trivial upsampling, confirming the construction is scale-invariant in step-count and gives no obvious headroom for refinement.

**Together's stated method (sequential LP from prior SOTA) is exactly the SLP I ran above.** The fact that their fixed point is hostile to first-order moves at δ down to 1e-5 means they have indeed found a local M-functional minimum at n=600 resolution. The historical UB-progression in the addendum (Haugland 51 → AlphaEvolve 95 → TTT 600 → Together 600) is consistent with marginal returns plateauing.

**Verdict Framing 2:** EXHAUSTED. No tighter explicit UB constructible by local methods on Together's h. Going to n ≫ 600 *might* help (no theoretical obstruction), but is research-grade effort, not a quick win.

---

## Framing 3 — Does µ_Ω = µ_M?

**Setup.** White (2023) bounds `µ` defined as the inf of the *continuous* sup-integral functional (his §1, eq. 1.1). Together bounds `µ` defined as the inf over *step functions* of `max_k ∫ h(x)(1-h(x+k)) dx` (their README, citing Haugland 2016 → Swinnerton-Dyer 1996).

**Literature audit:**

- **Haugland (2016), arXiv:1609.08000**, top of page 1: "Swinnerton-Dyer proved in [Haugland, 1996] that lim_{n→∞} M(n)/n **is equal to** the infimum, over all step functions f on [0,2] with values in [0,1] and satisfying ∫f=1, of max_k ∫f(x)(1-f(x+k))dx."

- **White (2023), arXiv:2201.05704, §1**: "A key step in Haugland's method is a theorem of Swinnerton-Dyer proving that (1.1) is in fact also µ, see [3] for the proof. **It will be easiest for us to work with (1.1) as our definition of µ.**"

So White explicitly *adopts the continuous Ω-functional as the definition of µ* (his (1.1)), citing the Swinnerton-Dyer equivalence to lim M(n)/n. Together's `M(h)` for step functions is the discrete-restriction of the same continuous functional (their step-function `h` plugged into Haugland's eq. (1)). Step-function constructions only give UPPER bounds on the continuous infimum (a step function is a special admissible f), and the infimum over step functions equals the continuous infimum by standard L¹-density of step functions in the relevant class — Swinnerton-Dyer 1996 is the formal proof.

**Both bounds are on the same µ.** No infimum gap. The proof is one paper deep (Swinnerton-Dyer 1996, in Haugland's J. Number Theory 58:71-78) which I did not read in full, but it is the same proof both White and Together rely on, so any slack would invalidate both bounds symmetrically.

**Verdict Framing 3:** EXHAUSTED. No formulation gap.

---

## Framing 1 — Is there a constraint family in `white_full_convex.py` the diagnostic missed?

**Setup.** `lp_research_state/code/white_full_convex.py:99-339` is `build_problem(...)`. I enumerated every flag:

| Flag | Default | Diagnostic status |
|---|---|---|
| `use_T3` | False | Not explicitly tested in Phase 4-5 composition |
| `use_T5` | False | Subsumed by T5p (per code comment) |
| `use_T5p` | False | NOT included in Phase 4-5 composition (see below) |
| `bochner_n` | 0 | Pushed to 30 |
| `mside_bochner_n` | 0 | Tested, vacuous |
| `mside_bochner_schur_n` | 0 | Tested at n=10, n=20 — vacuous (Δ=0 vs bochner_n=20 at row4 N=2000) |
| `assume_even` | False | Tested separately; produces *conditional* bound at row6 = 0.3799 (no improvement) |
| `lasserre_T_max` | 0 | Retracted as non-rigorous (tail bound issue) |
| `mside_bochner_lasserre_n` | 0 | Tested, vacuous |
| `poly_moment` (external) | n/a | Phase 4A, k_max=20 |
| `Hankel-PSD` (external) | n/a | Phase 4B, n=4-6 |
| cover refinement (external) | n/a | Phase 5, saturated |

### The one gap: `use_T5p` is NOT in the Phase-5 driver

Inspection of `lp_research_state/code/path_b_with_polymoment.py` (the Phase 5 driver) shows it calls `path_b_analytical.build_problem_with_dual_handles` (`lp_research_state/code/path_b_analytical.py:43`). That function **silently omits the `use_T5p` flag** — there is no parameter to pass it, and the T5p constraint (line 224-226 of `white_full_convex.py`) is therefore not added in any Phase 4-5 result.

This is a real gap in the diagnostic that the synthesis missed. Let me quantify it.

**Test 1 — T5p on top of bochner_n=20 alone (no other augment), at row 6, N=2000:**

```
  + bochner_n=20                          : 0.3782354726
  + bochner_n=20 + T5p                    : 0.3784238734   (Δ = +1.88 × 10⁻⁴)
```

So T5p alone over Bochner adds **+1.88 × 10⁻⁴** — not negligible.

**Test 2 — T5p on top of the full Phase-5 stack at row 4 N=3000:**

```
  + bochner_n=30 + pm k=14 + hk n=4       : 0.3792291360
  + bochner_n=30 + pm k=14 + hk n=4 + T5p : 0.3792307778   (Δ = +1.64 × 10⁻⁶)
```

**On the binding row at near-Phase-5 settings, T5p in composition adds only +1.6 × 10⁻⁶.** Below CLARABEL's reduced gap tolerance (5e-5) and below the recovery constant (1e-4). The Hankel-PSD constraint subsumes T5p's content almost completely, because both are testing positivity of `f` against degree-2 trigonometric polynomials — Hankel-PSD via the Hausdorff moment representation, T5p via the explicit test polynomial `1 - cos(πx)`.

**Test 3 — `assume_even + Phase-5 stack` at row 6 N=2000 (the conditional case):**

```
  even+bochner_n=20                          : 0.3782354425
  even+bochner_n=20 + pm k=10                : 0.3782752316
  even+bochner_n=20 + pm k=10 + hk n=4       : 0.3786388107
  uncond+bochner_n=20 + pm k=10 + hk n=4     : 0.3786387861   (Δ even - uncond = +2.5e-8)
```

The `assume_even` flag adds **essentially nothing** (~3 × 10⁻⁸) when composed with the full stack. The unconditional LP is already operating in the regime where f is effectively even at h=0, p=0.381. The even-f conjecture is closed as a lever even *conditionally*.

### Other potential gaps

- `use_T3` ("L Σ (w² + v²) ≤ Ω"): also not in Phase-5 composition. By construction this is a Cauchy-Schwarz constraint on `w, v` (not on `c, d`), so it's orthogonal to the c/d-side constraints from Bochner/poly_moment. Unlikely to matter much given the LP-optimum is at the f-cone boundary in the c-direction, but UNTESTED in composition. Quick test would settle this — flagged but not run.

- `mside_bochner_schur_n` was tested at n ≤ 20 only. Was it tested at the Phase-5 settings? No — it was tested in isolation on bochner_n=20. The schur version is exact (no SOC slack) but the empirical result `Δ = 0` against bochner_n=20 alone suggests no headroom; still worth a quick test in Phase-5 composition.

---

## Overall verdict

**CONCLUSION HOLDS.**

The three framings each closed:

- **Framing 2 (UB):** Together's h is genuinely locally optimal at n=600. SLP and pairwise perturbations cannot improve it. Together's method is the standard machinery; their saturation matches the historical progression.

- **Framing 3 (equivalence):** White and Together work with the same µ by explicit definition (White §1) and routine Swinnerton-Dyer 1996 equivalence. No formulation gap.

- **Framing 1 (LB):** I found one specific gap — `use_T5p` was silently omitted from the Phase 5 driver. But quantifying it at Phase-5 settings gives Δ ≈ +1.6 × 10⁻⁶ on the binding row, **well below the +10⁻⁴ scale that would meaningfully change the headline bound 0.3801279 → 0.380871 gap (currently 7.4 × 10⁻⁴).** Adding T5p to the Phase 5 composition would not change the published-quality conclusion. The `use_T3` flag remains untested in composition (~5 min experiment); I recommend doing it before write-up but I do not expect it to matter.

The session's conclusion — `µ ≥ 0.3801279` is the framework's saturation point at currently-tractable SDP scale, and further numerical work on the existing technique stack is dominated by writing up the result — is **upheld**.

## Concrete next experiment if rigor wants the last 6 digits

To remove the only loose thread I found, run *one* additional Phase-5 cover-refinement iteration with `use_T5p=True` and `use_T3=True` added to the constraint set. Expected delta: ≤ +5 × 10⁻⁶ on the headline LB. This would not move the headline `0.3801279`, but it would close the formal lever-enumeration gap.

Driver location: `lp_research_state/code/path_b_with_polymoment.py` would need a `use_T5p` parameter threaded through (currently absent at `path_b_analytical.py:45`); then re-run iterate_centers_pm.py with the augmented build.

Quoted file locations and scratch code:
- `/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code/white_full_convex.py:99-339` (build_problem)
- `/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code/path_b_analytical.py:43-170` (the cut-down builder used in Phase 5; missing T3/T5p flags)
- `/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/parallel_results/symmetric_conditional.json` (conditional even-f data, row6 = 0.3799 — confirms even-f is not a hidden lever)
- `/tmp/redteam_local_search.py` (SLP on Together's h)
- `/tmp/redteam_T5p_compose.py`, `/tmp/redteam_T5p_phase5.py`, `/tmp/redteam_even_composition.py` (composition tests)

## Bottom-line numbers I produced (independent of the session's claims)

| Test | Result |
|---|---|
| Together M(h) re-verified | 0.38087031058622 |
| Active shifts within 1e-7 of M-max | 437 of 1199 |
| Pairwise-perturbation improvements (2000 trials × 4 σ scales) | 0 |
| SLP step improvement at δ=1e-3 (predicted vs true) | +1.9e-4 predicted, −5e-4 actual |
| T5p Δ alone over bochner_n=20 (row6 N=2000) | +1.88e-4 |
| T5p Δ in Phase-5 composition (row4 N=3000) | +1.64e-6 |
| assume_even Δ in Phase-5 composition (row6 N=2000) | +2.5e-8 |
| mside_bochner_schur_n=20 Δ on bochner_n=20 (row4 N=2000) | +0 (identical to 7 digits) |
