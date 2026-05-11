# Lever decision memo: discrete M(n) optimizers vs Together's h*

**Date:** 2026-05-10
**Question:** Does the lifted continuous density `f_n` from discrete M(n) optimizers
converge in shape to Together's continuous minimizer `h*` as n grows?

## TL;DR — Verdict: DIVERGE (with one caveat)

- L¹ distance between lifted-and-aligned `f_n` and Together's `h_fold` plateaus
  near **0.25** across n = 4..18 — no decreasing trend.
- Pearson correlation hovers around **+0.5–0.65**, never approaching 1.
- KS distance plateaus near **0.07** — stable, not shrinking.
- `M(n)/n` at n = 18 is **0.444**, still far above µ ≈ 0.379 — n is small enough
  that the discrete optimizer is not yet in any asymptotic regime.
- `h*` is genuinely **fractional**: 58% of cells in `(0.05, 0.95)`, 233 distinct
  values; while discrete `f_n = 1_{A*}` is necessarily a `{0, 1}` indicator on a
  `2n`-cell grid. Even at `n = 18`, the discrete problem cannot **represent**
  `h*`'s fractional richness.

Verdict: **DIVERGE**.
**Caveat:** the divergence is genuine but may be partly artifactual due to
discrete vs continuous representational mismatch. Going to n = 30+ would
tell us whether `M(n)/n` keeps refining its shape toward h* or remains stuck.

## Next-step pivot: **E (M-side SDP)** as the primary lever, with optional supporting work

Because the discrete sequence cannot represent h*'s fractional structure even
in principle, the "would `f_n → h*`?" question is undecidable from purely
discrete optimization at tractable n. The right way to test whether h* is
asymptotically optimal is the **M-side SDP** — encode `h: [0, 2] → [0, 1]`
with fractional values directly and ask the dual whether any improvement on
Together's `0.380871` is possible, mirroring the Ω-side Bochner-PSD machinery.

Concrete plan:

1. **M-side SDP encoding (~ 2 weeks).** Write `mside_sdp.py` analogous to
   `white_full_convex.py` but minimizing Together's `M(h) = sup_k ∫ h(x)(1 −
   h(x+k)) dx` instead of White's `Ω(f) = sup_t ∫ f(x) f(x+t) dx`. Reuse
   Bochner moment-matrix infrastructure from `bochner.py`. Sanity check:
   at `h = h*`, the SDP value should equal `0.380871` to working precision.
2. **Cross-check Together's h*.** Run a single solve at moderate N (say
   N = 1200) to get a rigorous dual lower bound on M-side; compare to
   the 0.379544 from Ω-side. If the M-side bound matches or exceeds, then
   the two formulations are tight at the µ-limit and h* IS near-optimal —
   in that case our DIVERGE verdict here would re-route to a "the discrete
   sequence at n ≤ 18 is too small to see the shape" reading and we'd
   resurrect Pivot D' (Lipschitz/Sobolev) on the M side.
3. **If M-side gives a meaningfully different lower bound from Ω-side,**
   then the two formulations have different optimal `h`, requiring separate
   augmentation strategies. In that case the Ω-side ceiling of 0.379544 is
   not the right number to push — we should report `min(M-side LB, Ω-side LB)`
   as the unified bound and pursue augmentation on the looser side.

Effort estimate: **2 weeks for prototype**, 1 month if Bochner-on-M-side needs
new identities (the convolution `h * (1 − h)` is different from the
auto-correlation `f * f`; the Fourier moment encoding will require fresh derivation).

**Optional supporting computation (low priority): push M(n) further** — extend
to n ≈ 22 via better B&B (took 11 minutes at n = 18 single-threaded; n = 22
likely an overnight run). Useful as supporting evidence but does NOT change the
strategic direction since the representational mismatch is the dominant issue.

## Evidence

### Discrete M(n) sequence (n = 2..18, brute-force exact)

| n | M(n) | M(n)/n | blocks | time (s) |
|---|------|--------|--------|----------|
| 2 | 1 | 0.5000 | 3 | < 0.01 |
| 3 | 2 | 0.6667 | 4 | < 0.01 |
| 4 | 2 | 0.5000 | 5 | < 0.01 |
| 5 | 3 | 0.6000 | 4 | < 0.01 |
| 6 | 3 | 0.5000 | 7 | < 0.01 |
| 7 | 3 | 0.4286 | 5 | < 0.01 |
| 8 | 4 | 0.5000 | 7 | 0.02 |
| 9 | 4 | 0.4444 | 7 | 0.03 |
| 10 | 5 | 0.5000 | 10 | 0.21 |
| 11 | 5 | 0.4545 | 9 | 0.28 |
| 12 | 5 | 0.4167 | 7 | 0.40 |
| 13 | 6 | 0.4615 | 9 | 3.4 |
| 14 | 6 | 0.4286 | 7 | 4.6 |
| 15 | 6 | 0.4000 | 13 | 8.0 |
| 16 | 7 | 0.4375 | 9 | 54 |
| 17 | 7 | 0.4118 | 9 | 66 |
| 18 | 8 | 0.4444 | 11 | 664 |

M(n)/n trends downward but very slowly; at n = 18 we're at 0.444 vs target µ ≈ 0.379.

### Convergence metrics (smoothed f_n vs h_fold, best alignment over shift × complement)

| n | L¹ (raw) | L¹ (aligned) | Pearson | KS |
|---|----------|--------------|---------|----|
| 2 | 0.361 | 0.281 | +0.495 | 0.062 |
| 3 | 0.172 | 0.172 | +0.780 | 0.037 |
| 4 | 0.285 | 0.237 | +0.624 | 0.039 |
| 5 | 0.215 | 0.215 | +0.650 | 0.082 |
| 6 | 0.236 | 0.212 | +0.670 | 0.041 |
| 7 | 0.368 | 0.245 | +0.600 | 0.068 |
| 8 | 0.306 | 0.241 | +0.588 | 0.054 |
| 9 | 0.336 | 0.259 | +0.519 | 0.072 |
| 10 | 0.277 | 0.265 | +0.533 | 0.071 |
| 11 | 0.326 | 0.240 | +0.527 | 0.080 |
| 12 | 0.395 | 0.257 | +0.569 | 0.072 |
| 13 | 0.330 | 0.260 | +0.567 | 0.063 |
| 14 | 0.406 | 0.253 | +0.584 | 0.077 |
| 15 | 0.392 | 0.248 | +0.546 | 0.071 |
| 16 | 0.388 | 0.276 | +0.541 | 0.068 |
| 17 | 0.399 | 0.268 | +0.561 | 0.078 |
| 18 | 0.337 | 0.267 | +0.513 | 0.056 |

Take-home: no decreasing trend on L¹(aligned) from n ≥ 4. Pearson is bounded
away from 1, plateau around 0.55. KS plateau ~0.07.

### Plot

See `lp_research_state/data/fn_vs_h_convergence.png`. Top: overlay of f_n
(smoothed, n = 4, 8, 12, 14, 16) with h_fold (red) and h restricted to [0,1]
(dashed black). Bottom: convergence metrics vs n.

The visual confirms: lifted f_n is essentially `{0, 1}`-step (with smoothing
edges) at any scale below 1/n; h_fold is gently fractional. No alignment
can make them close because they live in different function classes.

### Why the lifted f_n cannot match h*

`h*` has 58% of its cells strictly fractional, with 233 distinct values out
of 600. The discrete optimizer's lifted indicator can only take values 0
or 1 on the 2n-grid (or, after smoothing at width 1/n, the same with
softened edges occupying ~1 cell each side). At n = 18, the "smoothed
indicator" still has a strongly bimodal value distribution near {0, 1};
fractional values constitute only the boundary cells. To represent the
30%+ of mid-range fractional structure in h*, we'd need either continuous
weights (which is what the M-side SDP would give) or much finer
discretization than n ≤ 18 supports.

## Anti-pattern check

- Did we soften the verdict? No — DIVERGE is the right call.
- Did we account for symmetries? Yes — cyclic shift + complement alignment was
  performed; without it correlations were negative (-0.5 to -0.1) due to
  arbitrary orientation choice.
- Did we check `h*`'s structure? Yes — confirmed it is genuinely fractional;
  the discrete vs continuous mismatch is the dominant obstruction.

## Artifacts

- `lp_research_state/code/_brute_force_Mn_extended.py` — B&B + reflection
  symmetry, validated against existing n ≤ 12 results.
- `lp_research_state/code/_lifted_density_compare.py` — comparison machinery.
- `lp_research_state/data/Mn_optimizers_large.json` — n = 2..18.
- `lp_research_state/data/lifted_densities.npz` — sampled densities.
- `lp_research_state/data/lifted_comparison_table.json` — metrics table.
- `lp_research_state/data/fn_vs_h_convergence.png` — diagnostic plot.
