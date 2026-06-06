# Together-as-Primal SDP Diagnostic

**Date:** 2026-05-10
**Predecessor:** [CDE_PHASE1_RESULT.md](CDE_PHASE1_RESULT.md) — µ ≥ 0.3801279 (saturated CDE stack)
**Together (2026):** [github.com/togethercomputer/erdos-minimum-overlap](https://github.com/togethercomputer/erdos-minimum-overlap) — µ ≤ 0.380871
**Open gap:** [0.3801279, 0.380871], width 7.43 × 10⁻⁴

## TL;DR

The diagnostic plugged Together's piecewise-constant minimizer h* into White's SDP and read off slacks per constraint family. Three structural findings:

1. **Together's h* maps cleanly into White's variable space via even reflection** (`f_even(x) = h(|x|)/2`). The asymmetric direct embedding is *infeasible* in White's SDP — rejected by the Bochner-PSD constraint on (1 − f), λ_min M_n(1−f_direct) = −0.882.
2. **The CDE Phase 1-5 augmentation stack (Bochner + poly_moment + Hankel) is satisfied at f_even with substantial slack** — none of the augmentations is binding on Together's certificate. The binding constraint family must be the *cell-kernel cos/sin autocorrelation envelope* in [`white_full_convex.py:176-190`](lp_research_state/code/white_full_convex.py).
3. **The gap between Phase-5's SDP-optimal f̃ and Together's f_even is 99.9% low-frequency (smooth structural deviation), NOT high-frequency (Gibbs ringing).** This is an empirical reversal of the prior hypothesis from `findings.md` and `probe.py` — the CDE stack already damped the Gibbs problem (recovered f̃ has min −0.97 vs the historical probe's −3.78).

**Recommendation: pivot to lever D — structural restriction.** The next attack is not "more PSD constraints" or "alternative basis"; it is "force f̃ to lie near a piecewise-constant family with O(1) breakpoints, justified by a discrete-to-continuous limit theorem on M(n)." Lever A (Lukács SOS) is unlikely to help. Lever C (combinatorial M(n)) remains independently valuable as an upper-bound attack.

## Q1: Is Together's f* feasible in our SDP encoding?

**f_even** (even reflection, `f(x) = h(|x|)/2` on [−2, 2]): yes, all constraints satisfied with margin.

| Constraint family | Quantity | Value at f_even | Notes |
|---|---|---|---|
| Bochner-PSD n=30 (f)       | λ_min M_n(f)     | 6.27 × 10⁻⁵ | PSD ✓ |
| Bochner-PSD n=30 (1−f)     | λ_min M_n(1−f)   | 5.66 × 10⁻² | PSD ✓ (large margin) |
| Poly-moment k=14 (min over k) | min slack    | 1.24 × 10⁻² | 34× the k=14 tail bound |
| Hankel-PSD n=6 (H)         | λ_min H_n        | 3.04 × 10⁻⁵ | PSD ✓ |
| Hankel-PSD n=6 (A)         | λ_min A_n        | 1.43 × 10⁻⁴ | PSD ✓ |
| Row 4 box on f̂(1)          | c[0] vs [0.3875, 0.3875] | c[0] = −2.31 × 10⁻⁴, d[0] = 4.07 × 10⁻¹⁷ | **NOT IN BOX** — requires box override |

The row-box mismatch is informative: Together's f does not "live in" row 4 of White's residual-region cover. White's 7 rows (and the CDE-discovered centers) were placed at *the SDP's optimizer's binding points*, not at where the true Erdős optimizer lives. The two are structurally different.

**f_direct** (asymmetric, `f(x) = h(x)` on [0, 2], 0 elsewhere): infeasible. Bochner-PSD on (1 − f_direct) has λ_min = −0.882, decisively rejecting it as a valid 0 ≤ f ≤ 1 function. The direct embedding concentrates all mass on [0, 2] with f = 1 in many cells; this is incompatible with the SDP's symmetric/L¹-1 framework. (The pre-Bochner trigonometric envelope already flags it: |d̂(1)| = 0.807 > 2/π ≈ 0.637; Σ(c² + d²) = 1.05 > 0.5.)

## Q2: What is Ω(f*) in our encoding?

For f_even (with row 4's c[0], d[0] bounds overridden to match f_even's f̂(1) = (−2.3 × 10⁻⁴, 0)):

- **SDP-encoded Ω(f_even, pinned) = 0.459311**
- f_even's TRUE autocorrelation = 0.387337
- White Phase 5 SDP optimum (row 4, no pin) = 0.380128
- Together's M(h*) = 0.380871 (M_T functional, not Ω_W)

The SDP's certified upper bound on autocorrelation at f_even is **0.072 above f_even's actual autocorrelation**. The SDP's cell-kernel envelope is not tight on this primal point; the Bochner / poly_moment / Hankel augmentations do not sharpen it in the directions f_even probes.

For f_direct: infeasible (24 s solve), per Q1.

## Q3: Which constraint has the most slack at f_even?

Poly-moment (k=2): slack = 0.292, **5662× the tail bound**. The poly-moment constraint family is the LEAST binding by a wide margin — almost completely inactive on Together's certificate.

Full poly-moment slack table at f_even:

| k | slack | tail bound | slack / tail |
|---|---|---|---|
| 2  | 2.92 × 10⁻¹ | 5.16 × 10⁻⁵ | 5662 |
| 4  | 1.29 × 10⁻¹ | 1.03 × 10⁻⁴ | 1254 |
| 6  | 6.96 × 10⁻² | 1.55 × 10⁻⁴ |  450 |
| 8  | 4.18 × 10⁻² | 2.06 × 10⁻⁴ |  202 |
| 10 | 2.68 × 10⁻² | 2.58 × 10⁻⁴ |  104 |
| 12 | 1.79 × 10⁻² | 3.10 × 10⁻⁴ |   58 |
| 14 | 1.24 × 10⁻² | 3.61 × 10⁻⁴ |   34 |

Even the tightest poly-moment row (k=14) has 34× the tail bound and 1.4× the L² trust threshold (0.0089) — useful for the SDP's own optimum, irrelevant for Together's primal.

## Q4: Which constraint is binding (the active lever)?

**None of the explicit augmentation families.** All slacks at f_even are above the L² trust threshold (≈ 0.009 from Parseval-exact T=4000 tail):

- Bochner-PSD: smallest margin is λ_min M_n(f) = 6.3 × 10⁻⁵ (small but positive — PSD).
- Poly-moment: smallest slack is 1.24 × 10⁻² at k=14.
- Hankel-PSD: smallest margin is λ_min H_n = 3.0 × 10⁻⁵ (PSD with margin).

The binding constraints must therefore be the *cell-kernel cos/sin autocorrelation envelope* ([`white_full_convex.py:176-190`](lp_research_state/code/white_full_convex.py)) — the original White constraints, NOT the CDE additions. Those encode: for each m ∈ {1, 2, ..., 2R}, `Ω ≥ Σ_j (w_j + v_j) · envelope_m(j)` where envelope_m is built from cos(πmx/2) cell-bound integrals.

This is the constraint family that produces Ω = 0.459 when (c, d) are pinned at f_even — it is the *core* White machinery, untouched by the CDE.

## Q5: Gap function structure

`g(x) := f̃(x) − f_even(x)`, with f̃ recovered from Phase 5's SDP-optimal (c, d) on a 10 000-point grid over [−2, 2].

- max|g| = **1.218**
- L²(g) = **0.589**
- **low-band fraction = 0.99902** (first 250/5001 rfft modes)
- **high-band fraction = 0.00098**
- argmin at x = **1.986**, g = **−0.971** (right boundary, where f̃ dips below f_even = 0)
- argmax at x = **0.145**, g = **+1.218** (interior, where f̃ spikes above f_even's plateau ≈ 0.5)

Visual: see [`lp_research_state/data/together_gap_plot.png`](lp_research_state/data/together_gap_plot.png).

**The gap is smooth low-frequency structural deviation, NOT Gibbs oscillation.** This refutes the prior hypothesis from `findings.md` and `probe.py` that the Fourier basis is structurally inadequate. The CDE stack already substantially damped the Gibbs problem (f̃ min = −0.97 vs historical −3.78 in pre-CDE probes). What remains is a coherent, low-frequency disagreement between f̃ and f_even — most of g's energy lives in the first 5 % of the FFT band.

## Q6: The call

**Lever D — structural restriction theorem.**

Rationale:

- The augmentation stack (Bochner + poly_moment + Hankel) does not bind on Together's f. Adding more PSD families along the same lines will continue to over-engineer that part of the SDP without affecting the binding constraint.
- The remaining gap is structurally low-frequency. A constraint family that forces f̃ near a piecewise-constant low-breakpoint family (Together-like) would attack the binding constraint directly.
- The natural mathematical object: a *limit theorem* of the form "the optimal f for the continuous Erdős problem is in the closure of piecewise-constant densities with O(1) breakpoints on a uniform grid." Such a theorem would follow from a discrete-to-continuous limit argument on M(n), and once established it would let us *restrict the SDP* to that family, collapsing the gap.
- Lever A (Lukács SOS / alternative basis) is unlikely to help — the gap is no longer Gibbs-dominated.
- Lever C (combinatorial M(n) push past Haugland's n=43) remains independently valuable as a way to tighten the upper bound and as evidence supporting the structural theorem (the finite-n optimizers should converge to Together's piecewise-constant shape).

**Concrete next steps (next session):**

1. Survey the literature: discrete-to-continuous limit theorems for set-overlap problems on Z. Specifically: is there an analogue of compactness for {densities of partitions} → {weak limits} that picks out piecewise-constant minimizers?
2. Empirically: compute M(n) optimizers for n = 20, 30, 40 (where exact values are known). Are the optimizers always piecewise-constant on a uniform grid? With how many breakpoints?
3. If the structural restriction can be proven: implement it as an additional constraint family in `white_full_convex.py` and re-run Phase 5.

## Open conceptual question

White's Ω and Together's M are different functionals that share µ as a limit. The SDP-encoded Ω(f_even) = 0.459 has no direct relationship to Together's M(f_even) = 0.380871 — they're measuring different things on the same input. A bolder reframing of the lower bound would be: *encode Together's M-functional directly as an SDP, in place of (or in parallel to) White's Ω-functional.* This is potentially a fifth lever (E) worth investigating in a follow-up — see if M-side Bochner (`mside_bochner_n` in `white_full_convex.py`) can be sharpened along the same CDE pattern that worked for f-side.

## Reproducing

```
cd /Users/benzanghi/Documents/Claude/Projects/Erdos
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import aggregate_results; aggregate_results()"
```

## Raw data

- `lp_research_state/data/together_f_star.json` — parsed h* (both Together and White conventions)
- `lp_research_state/data/together_f_star_fourier_even.npz` — projected f_even Fourier coefs (T=4000)
- `lp_research_state/data/together_f_star_fourier_direct.npz` — projected f_direct Fourier coefs
- `lp_research_state/data/row4_phase5_primal.npz` — Phase 5 SDP-optimal (c, d) at row 4
- `lp_research_state/data/row4_f_tilde.npz` — LP-optimal f̃ on 10 000-pt grid
- `lp_research_state/data/together_gap_function.npz` — f̃, f_even, gap on the grid
- `lp_research_state/data/together_diagnostic_results.json` — every numeric result above
- `lp_research_state/data/together_gap_plot.png` — three-panel visual

## Code

- `lp_research_state/code/together_loader.py` — h* loader + two White embeddings
- `lp_research_state/code/together_diagnostic.py` — Fourier projection, constraint diagnostics, SDP solve, gap analysis, `aggregate_results()`
- `lp_research_state/code/_together_projection_independent.py` — quadrature cross-check (≥10 digits)
- `lp_research_state/code/_fourier_convention_notes.md` — White's convention reference
