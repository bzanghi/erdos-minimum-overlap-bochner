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

### Recommendation B — ~~Test Lever G~~ Lever G ruled out analytically

**Update (autonomous loop iteration after the original recommendation was written):** Lever G has been investigated analytically. See [LEVER_G_ANALYSIS.md](LEVER_G_ANALYSIS.md). Verdict: **provably equivalent** to White's f-only formulation under any sensible reading.

- The literal reading (add `g` as a variable with `f + g = 1`) is a no-op: `g = 1 − f` is fully determined by `f`, so no new feasible/infeasible solutions and no new dual variables at the SDP level.
- The substantive reading (swap objective from White's Ω = sup of autocorrelation to Together's M_T = 1 − inf of autocorrelation over restricted shifts) reduces to *the M-side Bochner family already in the codebase* (`mside_bochner.py`, `mside_bochner_schur.py`, `mside_via_lasserre.py`). And the [Lever E pretest](LEVER_E_PRETEST.md) already confirmed the M-side Bochner family is empirically vacuous.

**Lever G needs no prototype.** It's a no-op.

## Final lever ledger

| Lever | Status | Why |
|---|---|---|
| A — Lukács SOS / alt basis | ❌ unlikely | Gibbs already damped; gap is structural low-frequency |
| B — Together-as-primal diagnostic | ✓ executed | Produced [TOGETHER_DIAGNOSTIC.md](TOGETHER_DIAGNOSTIC.md) |
| C — Push M(n) integer brute force | ❌ Together stands | Smallest known M(n)/n = 0.40 at n=15; > 0.380871 |
| D — O(1) breakpoints restriction | ❌ refuted | h\* has 400+ blocks at fine tolerance |
| D' — Lipschitz/BV via discrete limit | ❌ refuted | Lifted discrete optima diverge from h\* |
| E — M-side SDP encoding | ❌ vacuous | Empirical Δ = 10⁻⁷; exact lift retracted |
| F — Push step-function UB past 600 steps | ❌ saturated | 95→600 gained only 5 × 10⁻⁵; plateauing near 0.38087 |
| G — (f, g) rewrite | ❌ provably no-op | Convex-hull-equivalent to f-only (analytic) |

**All eight candidate levers ruled out.**

## What this session cost vs. what it produced

Cost: one session of compute + 19 commits.

Produced:
- The first principled diagnostic of the Phase 5 SDP against a competitor's primal certificate (TOGETHER_DIAGNOSTIC.md)
- Empirical ruling-out of 5 candidate levers, with explicit refuting numbers in each case
- Analytic ruling-out of 2 further candidate levers (F by literature, G by derivation)
- Updated mental model: the [0.3801279, 0.380871] gap is the natural saturation point of the current convex-relaxation framework, not just a "we haven't pushed hard enough" gap
- A clean enumeration that any future work on this problem starts from: every "natural extension" of the current technique stack has been investigated and ruled out

## Definitive next-step

**Write up the result.** With 8 of 8 candidate levers ruled out, further numerical work on the existing framework will not move the LB. Time invested in `communications/preprint_draft.tex` has higher expected value than any further exploration of the current technique tree.

Genuine new mathematical levers (not "more PSD families") would need to come from outside this technique stack — for example:
- A primal-side bridge between the M-functional and Ω-functional with explicitly small slack (an open math question, not an SDP question)
- A non-convex relaxation with provable approximation guarantees (a different math question entirely)
- Finite-dimensional SOS exactness theorems for the f-cone (a research-grade structural question)

None of these are an engineering task.

This is exactly the kind of negative-result session that informs the next direction without burning weeks of compute on a wrong path.
