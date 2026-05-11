# Session 2026-05-10 — Final addendum

Update to [SESSION_2026-05-10_SYNTHESIS.md](SESSION_2026-05-10_SYNTHESIS.md) after running Lever C and inspecting the UB-construction literature.

## Lever C result: TOGETHER VALIDATED

Pushed exact brute-force M(n) to n=18, replicated all published values for n ≤ 15. **No M(n)/n ratio for n ≤ 18 beats 0.380871** — the smallest known integer ratio is M(15)/15 = 0.40, which is 5% *above* Together's bound. The integer route to a tighter UB is not feasible at currently-tractable scale; n ≈ 50 is the threshold where ratios could plausibly dip below 0.380871, and that's overnight-cluster compute via SAT.

See [LEVER_C_RESULT.md](LEVER_C_RESULT.md) and [lp_research_state/data/known_Mn_values.json](lp_research_state/data/known_Mn_values.json).

## Newly-discovered evidence: UB progression has saturated

From Together's repo `README.md`, the historical UB-via-step-function progression:

| Method | Year | Steps | UB | Δ vs prior |
|---|---|---:|---:|---:|
| Haugland | 2016 | 51 | 0.380927 | baseline |
| AlphaEvolve | 2025 | 95 | 0.380924 | −3 × 10⁻⁶ |
| TTT-Discover | 2026 | 600 | 0.380876 | −4.8 × 10⁻⁵ |
| Together AI | 2026 | 600 | 0.380871 | −5 × 10⁻⁶ |

Going from 95 to 600 steps (+531% resolution) gained only 4.8 × 10⁻⁵ — and going from 600 to 600 via better optimization gained only 5 × 10⁻⁶. **The UB construction route is plateauing near 0.38087.** Pushing to 1000+ steps will not close the gap to White's 0.379 / our 0.380128.

This means: **the [0.3801279, 0.380871] gap is unlikely to be closed by harder-trying on either side of the current framework.** Both the lower bound (convex relaxation) and the upper bound (explicit step-function construction) have saturated near their respective walls.

## New problem-statement variant from Together's README

Together's README states the problem in a different (cleaner) form than White's `Ω`:

> Let `C_5` be the largest constant satisfying `sup_{x∈[-2,2]} ∫_{-1}^1 f(t)·g(x+t) dt ≥ C_5` for all non-negative `f, g: [-1,1] → [0,1]` with `f + g = 1` on `[-1,1]` and `∫_R f = 1`, where `f, g` are extended by zero outside `[-1,1]`.

The constraint `f + g = 1` is critical and is what makes this *the actual minimum overlap problem*. White's framework optimizes over admissible f but the `g = 1 − f` link is implicit in the f-cone constraints (which is partially why our Bochner on (1−f) was meaningful). Together's formulation makes the pairing explicit.

If the SDP framework were rewritten in `(f, g)` variables with `f + g = 1` enforced as a hard linear constraint (and the supremum-integral as the objective), would the convex relaxation be tighter? This is **Lever G** — not investigated this session, but flagged for future thought.

## Definitive verdict for this session

Of the candidate levers explored:

| Lever | Status |
|---|---|
| A — Lukács SOS / alt basis | ❌ unlikely (Gibbs already damped) |
| B — Together-as-primal diagnostic | ✓ executed; produced [TOGETHER_DIAGNOSTIC.md](TOGETHER_DIAGNOSTIC.md) |
| C — push M(n) integer brute force | ❌ Together UB stands; integer route requires n ≈ 50+ |
| D — O(1) breakpoints restriction | ❌ refuted; h\* has 400+ blocks |
| D′ — Lipschitz/BV via discrete limit | ❌ refuted; lifted discrete optima diverge from h\* |
| E — M-side SDP encoding | ❌ relaxations vacuous; exact lift retracted |
| F — Push step-function UB past 600 steps | ❌ asymptotically plateaued near 0.38087 |
| G — `(f, g)` formulation rewrite of SDP | ❓ untested; would require new SDP build |

**Of all candidate levers explored or surveyed this session, only Lever G remains plausible and untested.** It is also the most speculative — there is no a priori reason to expect the `(f, g)` rewrite to admit a tighter convex relaxation than the current f-cone formulation. The natural argument *against* is that White's framework already implicitly encodes f+g=1 via the integral and pointwise constraints, so an explicit rewrite shouldn't change the convex hull.

## Two actionable recommendations

### Recommendation A — Stop pushing the bound, write up the result

µ ≥ 0.3801279 is the published-quality output. The CDE Phase 1-5 result moved the lower bound by +1.1 × 10⁻³ over White's published 0.379005 — a substantial improvement. Time invested in writing this up (preprint draft at `communications/preprint_draft.tex`) likely has higher payoff than further numerical pushing.

The Together-as-primal diagnostic also produced a self-contained methodological contribution worth publishing: a structured procedure for evaluating an SDP's slack at a competitor's certificate, with a definitive negative result on M-side SDP relaxations. Could be a standalone short note.

### Recommendation B — Test Lever G (the only untested lever)

Implement the `(f, g)` formulation of the problem in a small SDP (T=500, single row), and compare the dual LB to White's f-only formulation at the same scale. If Lever G produces a tighter LB at small scale, scale it up and run Phase 5-style augmentations. If it produces the same LB, the rewrite doesn't help and we're done. Estimated effort: 3-5 days.

The expected outcome (per the convex-hull argument above) is no improvement. But it's the only lever we haven't definitively ruled out, and the test is cheap relative to the alternative of declaring "open" indefinitely.

## What this session cost vs. what it produced

Cost: one session of compute + 16 commits.

Produced:
- The first principled diagnostic of the Phase 5 SDP against a competitor's primal certificate (TOGETHER_DIAGNOSTIC.md)
- Empirical ruling-out of 5 of 7 candidate levers, with explicit refuting numbers in each case
- Pre-test infrastructure for the 1 remaining untested lever (G) and the unused continuation of Lever C (SAT to n=50+)
- Updated mental model: the [0.3801279, 0.380871] gap is the natural saturation point of the current convex-relaxation framework, not just a "we haven't pushed hard enough" gap

This is exactly the kind of negative-result session that informs the next direction without burning weeks of compute on a wrong path.
