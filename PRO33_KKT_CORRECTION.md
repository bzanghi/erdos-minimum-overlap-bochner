# PRO-33: Correction — Together's h* IS a KKT point; PRO-23's functional equation dropped edge terms

**Status:** Complete (this session). **Severity: preprint-affecting.** The claim
"Together's h* is not a tight KKT point" (PRO-23, [LEVER_FUNCTIONAL_EQUATION.md](LEVER_FUNCTIONAL_EQUATION.md),
propagated to `communications/preprint_draft.tex` §"KKT functional equation") is **refuted**.
With the correctly derived stationarity system, h* satisfies the KKT conditions to
**1.26 × 10⁻⁸** (interior) and **3 × 10⁻¹¹** (boundary sign conditions) — fully consistent
with Together's stated 10⁻⁹ optimization precision.

## 1. The bug in PRO-23's functional equation

PRO-23 tested the equation

```
Σ_{t∈S} γ_t · [h*(x+t) + h*(x−t)] = κ   on interior cells
```

This form implicitly assumes `x+t` and `x−t` always lie inside the support `[0, 2]`.
They do not: the active set S extends to |t| ≤ 256 cells (of 600), so for roughly
**85% of cells** at least one active shift falls off the support edge. The correct
per-cell gradient of `M_j(h) = (2/n) Σ_i h_i (1−h_{i+j})` (with h zero-extended,
`g = 1−h` supported on the array only) is

```
∂M_j/∂h_k = (2/n) · [ g_{k+j}·1{k+j ∈ [0,n)}  −  h_{k−j}·1{k−j ∈ [0,n)} ]
```

The term `g_{k+j} = 1 − h_{k+j}` contributes a **k-dependent indicator** `1{k+j in range}`
that PRO-23's symmetrized form treats as the constant 1. The omitted edge terms are O(1)
per cell — exactly the size of PRO-23's observed 7.6 × 10⁻³–5.6 × 10⁻² "residual".
Tested directly: the PRO-23 form evaluated with the *true* optimal multipliers has
interior spread **0.157** — the equation itself is wrong, not h*.

## 2. The correct test (LP-dual extraction)

Instead of guessing (γ, κ) via a threshold-classified feasibility LP, extract the
exact multipliers from the minimax linearization LP:

```
min u   s.t.  M_j(h*) + ∇M_j·δ ≤ u  (all 1199 signed lags),
              Σδ = 0,  −h* ≤ δ ≤ 1−h*,  |δ| ≤ r
```

(HiGHS, r = 10⁻⁴.) Results at n = 600 (`code: _pro33_kkt_correct.py`):

| Quantity | Value |
|---|---|
| LP-predicted first-order gain | **1.94 × 10⁻¹⁰** (≈ solver floor → h* is first-order optimal) |
| γ support | 391 lags, ≈-symmetric (asym 7.9 × 10⁻³), peak lags ±109..±124 |
| λ* (mass-constraint multiplier) | −3.798 × 10⁻⁴ |
| Interior stationarity residual (383 cells) | max **1.26 × 10⁻⁸** |
| Lower-cell sign condition (164 cells, need grad ≥ λ*) | worst violation **2.8 × 10⁻¹¹** |
| Upper-cell sign condition (53 cells, need grad ≤ λ*) | worst violation **2.9 × 10⁻¹¹** |

**Conclusion: h* is a numerically exact KKT point of the n=600 problem.**

## 3. Red-team SLP claim also corrected (both directions)

[OUT_OF_BOX_REDTEAM.md](OUT_OF_BOX_REDTEAM.md) reported the SLP "predicts an improvement
(e.g., 1.9 × 10⁻⁴ at δ=10⁻³) but true M regresses". Since M is *exactly quadratic* in h,
the per-lag linearization error is bounded by `(2/n)‖δ‖₂² ≤ 2r²·(‖δ‖₂²/nr²)` — at
r = 10⁻⁵ that is ≤ 2 × 10⁻¹¹, so a genuine 10⁻⁶-scale predicted LP gain could never
regress. The correct LP predicts gain ≈ 2 × 10⁻¹⁰, not 1.9 × 10⁻⁴. The red-team's
**conclusion** (h* locally optimal at n=600) was right; its **evidence** (large predicted
gain reversed by curvature) was an artifact of the same missing-edge-term/derivation bug class.

## 4. Consequences

1. **Preprint correction required.** The paragraph "KKT functional equation: Together's
   h* is not tight" and the display `µ < 0.380871 strictly` (justified via KKT residual)
   are unsupported and must be replaced by the opposite statement: h* is KKT-stationary
   at n = 600; strictness of `µ < 0.380871` (if claimed) needs a different argument
   (e.g., a strictly better finer-grid construction — see PRO-34).
2. **PRO-23 Step 4 (analytical solve) is un-blocked** — the corrected functional
   equation *is* satisfied by h*, and the exact discrete multipliers γ are now available
   (`data: _pro33_gamma.json` if regenerated). An analytical/continuum solve of the
   corrected equation `Σ_j γ_j[(1−h(x+j))·1{x+j∈supp} − h(x−j)·1{x−j∈supp}] = λ` on
   {0 < h < 1} is a live direction again.
3. **Recurring bug-class note** (third instance): truncated/idealized identity minus
   its boundary/tail terms = overclaim. Lasserre tail (retracted), positive-lag sup
   (fixed), now KKT edge terms (corrected). Any derived identity must be validated
   against a brute-force gradient before use as evidence.

## 5. Companion experiment: grid refinement (PRO-34)

h* being KKT-tight at n=600 relocates the UB slack question to the *grid*: is the
cell-doubled h* still stationary at n = 1200, or do the new degrees of freedom admit
descent? See [PRO34_UB_REFINEMENT.md](PRO34_UB_REFINEMENT.md).
