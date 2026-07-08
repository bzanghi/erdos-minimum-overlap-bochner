# PRO-34: Grid-refinement descent test — Together's h* is stationary at n = 1200 too

**Status:** Complete. **Question:** PRO-33 established h* is a numerically exact KKT
point of the n=600 problem, so any remaining slack in `µ ≤ 0.380871` must come from
grid refinement. Does cell-doubling to n = 1200 open a first-order descent direction?

**Answer: No.**

## Method

- Cell-double h* (600 → 1200 cells, value-preserving). Sanity: `M_max` at n=1200 is
  bit-identical to n=600 (`0.380870310586220`) — doubling preserves the functional. ✓
- Solve the minimax linearization LP at n=1200: variables (δ, u), constraints
  `M_j + ∇M_j·δ ≤ u` on all near-active lags, `Σδ = 0`, box `0 ≤ h+δ ≤ 1`, trust
  radius `|δ| ≤ 10⁻⁴`.
- **Lag restriction (validity-preserving):** a step of radius r can change any `M_j`
  by at most `(2/n)·Σ_k|∇M_j,k|·r ≤ 4r`, so only lags with `M_j ≥ M_max − (4r + 10⁻⁴)`
  can become binding — 946 of 2399 lags. Post-check: omitted-lag ceiling `0.380619 < u*`. ✓
- **Solver note:** `method='highs'` (simplex) effectively hangs on this LP
  (>45 CPU-min, killed) due to the ~900-way degenerate near-tie structure;
  `method='highs-ipm'` solves it in **7 s**. Use IPM for all minimax LPs on this
  problem family.

## Result (`_pro34_refine_n1200.py` superseded by `scratchpad/n1200_diag` method)

| Quantity | n = 600 | n = 1200 |
|---|---|---|
| LP-certified max first-order gain (r=10⁻⁴) | +1.94 × 10⁻¹⁰ | **+1.94 × 10⁻¹⁰** |
| True M at LP step (exact eval) | no gain (≤ solver noise) | no gain (−3.4 × 10⁻⁹) |

The certified maximal first-order improvement at n=1200 equals the n=600 value to
3 significant digits — i.e. the 600 new degrees of freedom contribute **nothing**.
h* is first-order stationary under 2× grid refinement.

## Second-order analysis (n = 1200)

Tangent space: null space of the 860 active-lag gradients (tol 10⁻⁹) + mass
constraint, with box-active cells frozen — only **54-dimensional** (812 free cells,
rank 758). First-order cone directions are already excluded by the LP, so descent
would require a tangent direction δ with `Q_j(δ) = −(2/n)Σ_i δ_i δ_{i+j} < 0` for
every active lag j.

1. **Random probe:** 300 random tangent directions, exact quadratic M evaluated at
   t ∈ ±[10⁻³, 3×10⁻²]: **no descent**.
2. **PSD certificate attempt:** does some convex combination `Σγ_j Hess_j` become PSD
   on the tangent space? Supergradient ascent of `λ_min(Σγ B_j)` over the simplex
   peaks at **−5.5 × 10⁻⁶ < 0** — no certificate. By minimax duality this means a PSD
   *mixture* X of tangent directions exists with `tr(B_j X) < 0` for all j — but a
   mixture is not a direction.
3. **Rank-1 hunt:** smoothed-minimax (log-sum-exp) descent of `max_j Q_j(δ)` over the
   54-sphere, 30 restarts × 1600 iterations: best achievable common curvature is
   **+3.1 × 10⁻⁵ > 0** — every direction found has at least one active lag curving up.
   The gap between the dual value (−5.5×10⁻⁶) and the rank-1 optimum (+3.1×10⁻⁵) is
   the standard SDP-relaxation lifting gap, not evidence of a real descent direction.

**Conclusion: h\* is (numerically) a strict second-order local minimum at n = 1200.**

## Interpretation

Together's h* is not merely a converged local optimum of their n=600 search — it is
first-order stationary *as a continuum candidate* (tested through the n=1200 embedding).
Combined with PRO-33 this **strengthens** the status of `0.380871` as a serious
candidate for µ itself, and sharpens the project's structural dilemma:

- The LB framework is provably capped at `C_explicit = 0.380713` (Lever I′ saturation
  theorem, N=40000).
- If h*'s basin is globally optimal and µ ≈ 0.380871, then **no amount of work inside
  the current LB framework can close the gap** — the missing 1.6 × 10⁻⁴ between the
  framework ceiling and µ requires a qualitatively new LB architecture.
- If µ < 0.380713 (inside framework reach), then the UB extremizer lives in a
  *different basin* than h*, and UB progress requires global moves (basin hopping,
  combinatorial restructuring of the block pattern), not local refinement.

Either way, the "refine Together's grid" lever (the natural reading of PRO-4) is
**dead at first order in the h* basin** — worth knowing before anyone spends compute
on n = 2400 SLP runs.

## Follow-on directions (ranked)

1. **KKT-continuation to high precision (new best candidate for a µ value).** h*
   satisfies the corrected KKT system to 1.3×10⁻⁸ with LP-extracted multipliers
   (γ on 391 lags, λ* = −3.798×10⁻⁴). Newton on the corrected system
   (h, γ, λ) at n = 1200–4800 could polish the candidate to ~10⁻¹⁴, giving a
   high-precision µ candidate for PSLQ closed-form identification — the
   MU_HIGH_PRECISION PSLQ attempt failed for lack of precision in the *bracket*,
   but a KKT-polished stationary value is a much sharper target.
2. **Basin-diversity search on the UB side.** All known constructions (Haugland 51 →
   AlphaEvolve → TTT → Together 600) may share one basin. Structured restarts (e.g.
   random block patterns with the observed 28%/62%/10% lower/interior/upper cell
   fractions, or symmetrized/antisymmetrized seeds) at n=600 with the IPM-SLP would
   test whether 0.380871 is basin-unique.
3. **Beyond-framework LB architecture** — the only path to the LB side's missing
   1.6×10⁻⁴ if µ ≈ 0.380871. The saturation theorem now doubles as a *design spec*:
   any new architecture must not factor through cell-envelope + Bochner-PSD duals.
