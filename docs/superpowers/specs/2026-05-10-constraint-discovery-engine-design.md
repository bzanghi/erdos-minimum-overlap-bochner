# Constraint Discovery Engine — design

**Date:** 2026-05-10
**Status:** Phase-1 spec (probe + manual candidates). Phase-2 (autonomous overnight loop) seeded later.
**Context:** White's framework reduces "improve µ" to "find a new valid convex constraint." Current rigorous bound `µ ≥ 0.379544` (Bochner-PSD + ellipse extension). Lasserre-2 withdrawn (tail kills gain). Goal: build tooling that systematically discovers new valid constraints rather than relying on artisanal mathematical insight.

## Target outcome

User selected **D (maximum ambition)**: aim for closed-form µ or gap < 10⁻⁶. Realistic delivery is incremental rigorous improvement; the tool is structured so that even if no single discovery closes the gap, the search infrastructure compounds over time and across analogous problems.

## Mechanism (the meta-thesis)

White's program has the property that **any new valid convex constraint can only tighten the bound**. So "improve µ" is operationally: discover convex inequalities `c(f, M, c, d, w, v, …) ≤ 0` that

1. Are mathematically valid (rigorously provable),
2. Are convex in the LP's existing variables (or composable via an SDP / SOC lifting),
3. Are not already implied by existing constraints (cut a face of the current feasible region),
4. Survive the ellipse-extension argument (uniform across the residual region).

The Lasserre attempt failed property 1 once truncation was honest. The Bochner constraint satisfied all four. The Engine is the systematic search for more.

## Three-layer architecture

```
┌──────────────────────────────────────────────────────────┐
│  Outer loop (creative, Claude-driven)                    │
│  - propose candidate inequalities by analogy / theory    │
│  - sketch validity proofs                                │
│  - parallel agents brainstorm across families            │
└──────────────────────────────────────────────────────────┘
              ↓ candidate constraint spec ↓
┌──────────────────────────────────────────────────────────┐
│  Middle loop (registry, engineering)                     │
│  - typed catalog of constraints with provenance          │
│  - composition / dedup against existing constraints      │
│  - cutting-power measurement per row × per scale         │
└──────────────────────────────────────────────────────────┘
              ↓ measured candidates ↓
┌──────────────────────────────────────────────────────────┐
│  Inner loop (mechanical, automated)                      │
│  - probe: reconstruct f̃ from LP optimum                  │
│  - test physical f-properties; flag violations           │
│  - violation → suggest a constraint family               │
│  - run SDP with the candidate; record Δµ                 │
└──────────────────────────────────────────────────────────┘
```

The inner loop is the "cutting-plane oracle." Given (c*, d*, w*, v*) from the current LP, synthesize f̃(x) = ½ + Σ c*_k cos(πkx) + d*_k sin(πkx) and test against a battery of properties known to hold for any feasible f. Each *violation* is a hint that a constraint expressing that property would cut.

## Phase 1 (this session)

**Goal:** prove the probe works and produce ≥ 1 measured candidate constraint with non-zero cutting power on row 4.

Concrete steps:
1. Solve the Bochner-augmented LP at row 4 with modest scale (N=2000, T=800, R=10, bochner_n=15).
2. Reconstruct f̃ from (c*, d*).
3. Compute a violation panel:
   - `min_x f̃(x)`, `max_x f̃(x)` (pointwise positivity / boundedness)
   - `∫ f̃²` vs `∫ f̃` (relates to `f² ≤ f`)
   - `min eig` of higher-level moment matrices than `bochner_n` (does Bochner cut deeper?)
   - sign / monotonicity of `M̃ = f̃ * (1 − f̃)`
   - `f̃(x) · (1 − f̃(x))` integrated against test polynomials
4. For the most-violated property, propose a convex constraint targeting it.
5. Encode in cvxpy, re-solve, measure ΔΩ*.
6. If ΔΩ* > 0 with a clean derivation, queue for ellipse-extension testing at scale.

## Phase 2 (overnight, parallel agents)

Once Phase 1 yields ≥ 1 working candidate:
1. Push the candidate into `experiments_queue.json` with cron_runner integration.
2. Sweep `(N, T, n)` to characterize ceiling.
3. Run `path_b_*.py` with the augmented dual to test ellipse coverage.
4. Spawn parallel Claude agents to brainstorm sibling candidates (varied basis, varied test measure).
5. Document survivors in `findings.md`.

## Decisions

- **Scope:** focus on this problem; the framework is general but won't be packaged for reuse until it produces ≥ 1 µ improvement.
- **Rigor floor:** every constraint added must have a written validity argument. No heuristic constraints in the LP, ever — that was the lesson of the Lasserre retraction.
- **Verification policy:** any candidate that produces an improvement must be re-encoded by an independent script before being claimed; matches the project's existing cross-verification convention.

## Out of scope (Phase 1)

- Lean/Coq formalization (Approach D from the brainstorm)
- Combinatorial brute force on `M(n)` (Approach A)
- Alternative bases beyond Fourier (Approach B's harder variants — wavelet/RKHS)

These are queued for later phases if Phase 1 surfaces signal.
