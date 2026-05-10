# Together-as-Primal SDP Diagnostic — Design

**Date:** 2026-05-10
**Status:** design (pre-implementation)
**Author:** session 2026-05-10
**Successor to:** `2026-05-10-constraint-discovery-engine-design.md` (the CDE stack, now saturated at µ ≥ 0.3801279)

## Background

The Constraint Discovery Engine stack {Bochner-PSD + poly_moment + Hankel-PSD + iterative cover refinement} on the Fourier basis at SDP scale (N=10000, T=4000, bochner_n=30) has empirically saturated. Phases 4–5 added only +6 × 10⁻⁵ over Phase 3. Both findings.md's joint (N→∞, n→∞) extrapolation analysis and the probe diagnostic (which found `f̃ < −3.78` in the recovered density) point to the *Fourier basis itself* as the structural blocker, via Gibbs oscillations the SDP cannot see through with PSD constraints alone.

The current gap is `[0.3801279, 0.380871]`, width 7.43 × 10⁻⁴. Together's upper bound 0.380871 comes from an explicit 600-step piecewise-constant minimizer `f*` (the Together repo). No prior session has fed that `f*` into our SDP. We have only compared *values*.

## Goal

**Treat Together's `f*` as a primal candidate in our SDP, read off which constraints are tight or slack at `f*`, and produce a structural memo that selects the next attack vector.** The deliverable is *direction*, not a new µ value.

The next attack vector will be one of:
- **A — Lukács SOS for `f² ≤ f`**: encode the quadratic feasibility constraint via Fejér-Riesz / Toeplitz lift with a rigorous truncation tail bound (the rigorous descendant of the withdrawn Lasserre attempt).
- **C — Combinatorial M(n) at scale**: push exact M(n) past Haugland's n=43 using modern SAT/IP/branch-and-bound; attacks the gap from above.
- **D — Structural restriction theorem**: a theorem of the form "the optimal f is piecewise constant on a uniform grid" or similar structural fact, provably restricting the SDP search space.

The diagnostic chooses among A, C, D based on evidence rather than guess.

## Protocol

### Step 1 — Fetch and verify Together's `f*`

- Locate `f*` in `github.com/togethercomputer/erdos-minimum-overlap`. Expect either a CSV/JSON file of breakpoints+values or a script that regenerates it.
- Parse into a canonical representation: numpy arrays `(breakpoints, values)` for a step function on `[0, 1]` (or whatever interval Together uses).
- **Verify their value 0.380871 directly** from `f*` using the overlap functional as Together formulates it (independent of our SDP). Document any normalization differences vs. White's formulation in writing — both signs and scales.
- If the repo does not expose `f*` in raw form, regenerate by running their solver locally. If that fails (e.g., dependency issues, GPU requirement), reconstruct from their paper's description.

**Output:** `together_f_star.json` containing breakpoints, values, source provenance, and Together's claimed value re-computed by us.

### Step 2 — Project `f*` into our Fourier basis

- Read White's Fourier convention from `lp_research_state/code/white_full_convex.py`. The expected form is `(c_k, d_k)` for k = 1..T as cosine/sine coefficients on the relevant interval.
- Compute `(c_k, d_k)` from `f*` analytically — a step function's Fourier coefficients have closed-form per breakpoint. No numerical quadrature unless cross-checked.
- Compute and **record a rigorous bound** on the truncation tail `‖f* − f*_T‖_{L²}` and `‖f* − f*_T‖_{L¹}` for T = 4000. This is the irreducible representation error of `f*` in our basis.
- **Cross-verification:** independent re-implementation in `_together_projection_independent.py`. Both implementations must agree on every `(c_k, d_k)` to ≥10 significant digits.

**Output:** `together_f_star_fourier.npz` (the truncated coefficients) plus written truncation-error bound.

### Step 3 — Evaluate constraints at `f*`

For row 4 (the binding White center) and the binding CDE center from Phase 5:

- Substitute the projected `(c_k, d_k)` for the SDP's decision variables in `build_problem(...)` at that center. Compute, with the constraints from current Phase 5 parameters:
  - **Bochner-PSD:** `λ_min(M_n(f*))` and `λ_min(M_n(1 − f*))` at `n = 30`
  - **Poly-moment:** for each `k ∈ {2, 4, …, 14}`, the slack `m_{2k}(f*) − (−tail_k)`
  - **Hankel-PSD:** `λ_min` of the Hankel block at `n = 6`
  - **Objective:** `Ω(f*)` in our SDP encoding
- Tabulate: per constraint, is `f*` strictly feasible, tight, or infeasible? If infeasible, by how much (and from which truncation: ours via T=4000, or theirs via 600-step quantization)?

The critical comparison: does `Ω(f*)` in our SDP encoding equal Together's claimed 0.380871, or is it lower?
- **If equal (within ~10⁻⁵):** our encoding correctly captures `f*`'s value but the optimizer of our SDP is a *different* function (the Gibbs-oscillating one). Gap is structural: a constraint that restricts to "f*-like" functions would close it.
- **If our Ω(f*) is lower than 0.380871:** our objective is encoded weaker than the true overlap functional, OR our truncation discards information `f*` exploits. Either case is a finding.
- **If `f*` is infeasible in our SDP:** our relaxation has a constraint that rules out the true optimizer — a bug-grade finding. Highly informative; would require constraint audit.

**Output:** `together_diagnostic_results.json` with per-center, per-constraint slack tables.

### Step 4 — Dual comparison

- Solve our Phase 5 SDP at row 4 normally. Extract the dual `ψ` (the test-function multiplier on the primal constraints) via the `dual_extractor.py` machinery.
- Recover our LP-optimal `f̃` from the primal Fourier coefficients (this is the Gibbs-oscillating density from the probe).
- Compute the *gap function* `g(x) = f̃(x) − f*(x)` on a dense grid (`x ∈ [0, 1]`, M ≈ 10⁴ points).
- Plot and inspect `g(x)`:
  - If `g` is dominated by high-frequency oscillation: **Gibbs-only** failure → A (Lukács SOS) or alternative-basis is the lever.
  - If `g` has low-frequency, spatially-localized structure: there is a *structural* difference between our optimizer and Together's → D (restriction theorem) is the lever.
  - If `g` is mostly localized where `f̃` is near 0 or 1: the relaxation of `f ∈ [0, 1]` is the weak point → revisit Bochner on `f` and `1 − f` jointly.
- Compute `⟨ψ, constraint_slacks(f*)⟩` — does the dual certify Together's value, or is it blind to the constraint family that makes `f*` feasible?

**Output:** `together_gap_function.npz` (the function `g` sampled densely) plus an annotated plot.

### Step 5 — Structural diagnosis memo

Write `TOGETHER_DIAGNOSTIC.md` answering, with numbers and plots:

- Q1: Is `f*` feasible in our SDP encoding? With what slack on each constraint?
- Q2: What is `Ω(f*)` in our encoding vs. Together's 0.380871?
- Q3: Which constraint has the *most* slack at `f*` (over-engineered, candidate to remove)?
- Q4: Which constraint has the *least* slack (the binding one — that's the family to strengthen)?
- Q5: Is `g(x) = f̃(x) − f*(x)` Gibbs-dominated, structurally-localized, or boundary-localized?
- Q6: **The call.** Next lever: A, C, or D, with a one-paragraph justification.

## Code

- `lp_research_state/code/together_loader.py` — Step 1 driver
- `lp_research_state/code/together_diagnostic.py` — Steps 2–4 driver
- `lp_research_state/code/_together_projection_independent.py` — independent re-implementation of Step 2 (project convention: `_` prefix = single-purpose script, not a library)

## Verification discipline

Project rule: every new computation gets independent cross-verification to ≥6 significant digits (≥10 for things that *should* be bit-equal).
- Together's stated value vs. our direct re-computation from `f*`: ≥6 digits.
- Fourier projection: ≥10 digits between the two implementations.
- `Ω(f*)` computed via direct overlap functional vs. via SDP-substituted: ≥6 digits.

## Stopping criteria

- If Step 1 takes > 60 minutes of wall time: switch to reconstructing `f*` from Together's paper description; if that also stalls, write up Step 1 as a "we cannot verify Together's certificate from public materials" finding and pivot the session to lever A directly.
- If Step 3 shows `f*` is *infeasible* in our SDP (any constraint rejected with slack > truncation error): stop and write up immediately. This is a major finding by itself — either our relaxation is wrongly tight or Together's certificate is approximate.
- If Step 3 shows `f*` is feasible and `Ω(f*) ≈ 0.380871` to within ~10⁻⁵: the diagnosis is "gap is structural", memo selects D, and the *next* session is the structural-theorem hunt.
- If Step 3 shows `f*` is feasible but `Ω(f*)` is materially lower than 0.380871: there is an SDP-encoding issue; investigate before drawing conclusions about levers.

## Risks

- **Repo content unknown.** Together's repo may not expose `f*` raw. Mitigation: reconstruct from paper; failing that, write up the obstruction.
- **Normalization differences.** Together and White may use different problem normalizations. Mitigation: Step 1 verifies their value end-to-end before Step 2 runs.
- **Truncation artifacts.** Projecting a 600-step function onto T = 4000 Fourier modes introduces error. Mitigation: rigorous tail bound recorded; re-run at T = 8000 if T = 4000 results are within the tail-bound noise.
- **Quantization artifacts in `f*` itself.** If their certificate is approximate (their solver's epsilon), our diagnostic inherits that error. Mitigation: memo reports the gap between their *stated* value and our *recomputed* value alongside everything else.
- **Inconclusive diagnostic.** It is possible the answer is "all three constraints are equally slack, no clean signal". Mitigation: that is itself a finding — would suggest the gap is genuinely a combination, and the right next move is A (most likely to be tractable) by default, with explicit caveat.

## What this is NOT

- Not a push for a new µ value. If µ improves it's a bonus, not the goal.
- Not a rewrite of `white_full_convex.py`. The diagnostic only *evaluates* the existing program with `f*` substituted.
- Not a commitment to lever A, C, or D. The memo's recommendation may be reversed in the next session if execution surfaces new information.

## Success criteria

A committed `TOGETHER_DIAGNOSTIC.md` that answers Q1–Q6 with numbers and a defensible next-lever recommendation. Reproducible via committed code + data. Cross-verified to project standard.
