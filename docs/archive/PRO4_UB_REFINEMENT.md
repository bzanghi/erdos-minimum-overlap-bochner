# PRO-4: Together UB refinement — NEGATIVE RESULT

**Status:** Closed. Together's h* is a tight local optimum of the discrete
piecewise-constant minimax problem at n = 600 (and 1200) cells. Local
refinement does not push the UB below 0.380871.

## Motivation

PRO-23 found that Together's h* has a residual ε ≈ 7.6 × 10⁻³ when tested
against the continuous KKT functional equation `Σγ · [h(x+t) + h(x-t)] = κ`.
We interpreted this as evidence that "Together's UB 0.380871 is not tight"
and `μ < 0.380871` strictly. PRO-4 was to verify this empirically by
attempting local refinement of h*.

## Methods tried

### 1. LP-minimax steepest descent (`_pro4_refine_together.py`)

At each iteration, solve
```
min t   s.t.   g_j · d ≤ t  for j ∈ S (active set)
              Σ d = 0
              d_i ≥ 0 if h_i = 0,  d_i ≤ 0 if h_i = 1
              ||d||_∞ ≤ 1
```
where g_j = ∇M(jL). Line-search for α along d.

**Result:** LP returns t* ≈ −0.19 (substantial descent direction in linearized
sense), but line-search caps α at ~10⁻¹⁰. Cumulative improvement after
20 iterations: 3 × 10⁻¹³.

**Why:** any nontrivial step in d makes a near-active shift (within ~10⁻⁹
of max) overtake the active set's reduction. Active set is *too dense* to
admit a meaningful step.

### 2. Smoothed log-sum-exp descent (`_pro4_refine_v2.py`)

Surrogate `f_τ(h) = (1/τ) log Σ_j exp(τ·M(jL))`. Projected gradient with
τ ramped 10⁴ → 10⁷.

**Result:** ||grad|| = 0.024 consistently (nonzero), but every line-search
step rejected from α = 10⁻¹² downward. Zero net improvement.

### 3. Upsample to 1200 cells + LP descent (`_pro4_upsample_refine.py`)

Upsampling preserves M at the original (even) shifts; odd shifts have
M = 0.380870 also. Total feasible-set dimension doubles.

**Result:** Same behavior as 600-cell case. Cumulative improvement ~10⁻¹¹
over 26 iterations, then stalls.

## Reconciliation with PRO-23

The KKT residual 7.6 × 10⁻³ found in PRO-23 was at the *continuous-problem*
level (KKT functional equation for h : [0, 2] → [0, 1] in L²). What it
measures is the discrepancy between the discrete h* and the unknown
continuous h\*\* that satisfies the analytic KKT equation.

The discrete h* IS optimal (or nearly so) within the 600-cell
piecewise-constant ansatz: any descent direction is defeated by the
density of near-active shifts. PRO-23's residual does NOT directly
quantify the gap |μ − M(h\*)|.

The byproduct claim `μ < 0.380871 strictly` is still likely true (because
the continuous infimum should be ≤ the discrete minimum), but quantifying
the gap rigorously is unsolved.

## What would unblock further UB refinement

Beyond local descent, one of:
- **Finer discretization** (n ≥ 10⁴) with global optimization (SA or
  evolutionary methods). Together's published h* used ~600 cells, but
  their internal pipeline involved annealing — replicating that at higher
  resolution is plausibly the right path. Cost: 10–100 GPU-hr.
- **Smoothing the ansatz** to allow `h` to be piecewise-linear or
  piecewise-polynomial, expanding the feasible class beyond piecewise
  constant. Would change the discrete optimum.
- **Analytical solution** of the KKT functional equation directly (the
  unfinished Step 4 of PRO-23). Blocked on having a tight enough h* to
  start from — chicken-and-egg.

## Status and recommendation

PRO-4 closed as negative result. The "Together UB has slack" claim from
PRO-23 stands as a *qualitative* statement but does not translate to a
tractable LB-side push at our discretization.

The published UB 0.380871 is robust against ~10⁻¹⁰-level local
perturbations. Anyone trying to claim `μ < 0.380871` strictly will need
to either reproduce + extend Together's global-optimization pipeline or
solve the continuous problem analytically.

## Code retained

- `_pro4_verify_together.py` — verifies M(h*) = 0.3808703106
- `_pro4_refine_together.py` — LP-minimax descent prototype
- `_pro4_refine_v2.py` — smoothed log-sum-exp descent prototype
- `_pro4_upsample_refine.py` — 1200-cell upsample + LP descent

All three refinement scripts are *demonstration of stuck-at-local-min*
artifacts; not to be used as production tools.
