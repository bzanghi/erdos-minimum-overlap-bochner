# Constraint Discovery Engine — Phase 1, 2, 3 Results

**Date:** 2026-05-10
**Status:** rigorous, reproducible, under identical rigor convention as the published
`path_b_rigorous.json` (margin = 1e-6, Lipschitz grid error ≈ 2.15 × 10⁻⁶).

## Headline

| | µ ≥ | Δ vs White | Δ vs prior |
|---|---|---|---|
| White (2023) | 0.379005 | — | — |
| Bochner-PSD + ellipse, 7 White centers (prior published) | 0.379544 | +5.39 × 10⁻⁴ | — |
| **Phase 1**: cover refinement (n=20 + 4 CDE centers) | 0.379620 | +6.15 × 10⁻⁴ | +7.57 × 10⁻⁵ |
| **Phase 2**: cover refinement at bochner_n=30 (7 White n=20 + 5 CDE n=30) | 0.379879 | +8.74 × 10⁻⁴ | +3.35 × 10⁻⁴ |
| **Phase 3**: poly-moment + n=30 + 12-center cover | **0.380067** | **+1.062 × 10⁻³** | **+5.23 × 10⁻⁴** |

### Open gap (post-Phase-3)
- Lower: `µ ≥ 0.380067` (this work)
- Upper: `µ ≤ 0.380871` (Together Computer, March 2026)
- **Open width: ~8.0 × 10⁻⁴**, down from 1.3 × 10⁻³ at the start of this session — about 40% of the gap closed in one session.

## Phase 1 — cover refinement (committed earlier)

White's 7 centers were chosen in 2023 for *his* unaugmented bound. Once Bochner-PSD is
added, the dual landscape shifts. The CDE reads off the binding point of the current
envelope cover, solves a new SDP center there at full scale, extracts duals via path-B
machinery, and adds the new ellipse. Iterating converges in 4 steps. **+7.57 × 10⁻⁵.**

## Phase 2 — cover refinement at bochner_n=30

Same loop, new centers solved at `bochner_n=30` instead of `n=20` (existing 7 White
centers kept at `n=20`). Saturates after 5 iterations.

| iter | new center (h_c, p_c) | V_c (primal) | grid_min after | Δ |
|---|---|---|---|---|
| 0 (baseline) | — | — | 0.3795476 | — |
| 1 | (0.00000, 0.39417) | 0.3799355 | 0.3796814 | +1.34 × 10⁻⁴ |
| 2 | (0.00343, 0.38417) | 0.3799422 | 0.3798356 | +1.54 × 10⁻⁴ |
| 3 | (0.00004, 0.39015) | 0.3798862 | 0.3798676 | +0.32 × 10⁻⁴ |
| 4 | (0.01222, 0.39075) | 0.3801721 | 0.3798771 | +0.10 × 10⁻⁴ |
| 5 | (0.00814, 0.38955) | 0.3799726 | 0.3798813 | +0.04 × 10⁻⁴ |

Under uniform margin: **µ ≥ 0.3798792** (rigorous). Binding at (h=0.00513, p=0.38862),
witness cde_n30_iter5.

## Phase 3 — polynomial moments (in progress)

**A brand-new constraint family.** The probe showed the LP-optimal `f̃` (reconstructed
from Fourier coefficients) had `m_10 = -0.0072` — impossible for any `f ≥ 0` since
`m_{2k} = ∫x^{2k}f ≥ 0` always holds. Adding `m_{2k} ≥ -tail_bound_k` as a linear
constraint on `(c, d)` rigorously cuts.

Validity: `m_k = (1/2)α_0^(k) + Σ_{j=1..T} (c_j α_j^(k) + d_j β_j^(k)) + tail_k`. The
tail `|tail_k| ≤ (2/π)Σ_{j>T}(|α_j^(k)| + |β_j^(k)|)`. By integration-by-parts recurrence,
`|α_j^(k)|, |β_j^(k)| ≤ 2k/(πj)²` leading term for even k, giving `tail_k ≤ O(k/T)`.
At T=4000 and k≤14, tail bounds are 5e-5 to 3.6e-4 — small.

**Single-row test (row 4, full scale, bochner_n=20)**:
- Baseline: Ω* = 0.3796534
- + poly_moment (k=2,4,...,14): **Ω* = 0.3800915**
- **ΔΩ* = +4.38 × 10⁻⁴**

**Full-cover Phase 3** (12 centers re-solved at bochner_n=30 + poly_moment k_max=14):
- Per-center V_c range: 0.3801–0.3817 (all centers ≥ 0.38)
- grid_min = 0.3800695
- Lipschitz eps_grid = 2.16 × 10⁻⁶
- **µ ≥ 0.3800673** rigorously
- Binding at (h=0.00453, p=0.39215), witness `cde_n30_iter3` (a CDE-discovered center)

Per-center V_c contributions (Phase 3, sorted):

| center | V_c | h_c | p_c |
|---|---|---|---|
| row7 | 0.3816760 | 0.030 | 0.375 |
| row3 | 0.3807412 | 0.020 | 0.375 |
| row5 | 0.3804245 | 0.000 | 0.400 |
| row1 | 0.3804241 | 0.015 | 0.381 |
| row2 | 0.3804120 | 0.015 | 0.385 |
| cde_n30_iter4 | 0.3803548 | 0.012 | 0.391 |
| row6 | 0.3801852 | 0.000 | 0.381 |
| cde_n30_iter5 | 0.3801556 | 0.008 | 0.390 |
| cde_n30_iter1 | 0.3801413 | 0.000 | 0.394 |
| cde_n30_iter2 | 0.3801333 | 0.003 | 0.384 |
| row4 | 0.3801027 | 0.004 | 0.388 |
| cde_n30_iter3 | 0.3800792 | 0.000 | 0.390 |

## Code

- [lp_research_state/code/probe.py](lp_research_state/code/probe.py) — cutting-plane diagnostic probe
- [lp_research_state/code/iterate_centers.py](lp_research_state/code/iterate_centers.py) — iterative cover refinement
- [lp_research_state/code/cde_evaluate.py](lp_research_state/code/cde_evaluate.py) — Phase 1 rigorous LB
- [lp_research_state/code/cde_phase2_eval.py](lp_research_state/code/cde_phase2_eval.py) — Phase 2 rigorous LB
- [lp_research_state/code/hankel_probe.py](lp_research_state/code/hankel_probe.py) — Hankel-PSD violation probe
- [lp_research_state/code/poly_moment.py](lp_research_state/code/poly_moment.py) — polynomial-moment constraint encoder
- [lp_research_state/code/test_poly_moment.py](lp_research_state/code/test_poly_moment.py) — single-row poly-moment test
- [lp_research_state/code/path_b_with_polymoment.py](lp_research_state/code/path_b_with_polymoment.py) — Phase 3 full driver

## Reproducing

```bash
# Phase 1
.venv/bin/python lp_research_state/code/iterate_centers.py --max_iters 4 \
    --N 10000 --T 4000 --R 10 --bochner_n 20 --min_delta 5e-7

# Phase 2
.venv/bin/python lp_research_state/code/iterate_centers.py --max_iters 5 \
    --N 10000 --T 4000 --R 10 --bochner_n 30 --min_delta 5e-7 \
    --out lp_research_state/parallel_results/cde_iter_n30.json

# Phase 3
.venv/bin/python lp_research_state/code/path_b_with_polymoment.py \
    --N 10000 --T 4000 --R 10 --bochner_n 30 --pm_k_max 14

# Evaluate
.venv/bin/python lp_research_state/code/cde_phase2_eval.py
```

## Methodological note

White's framework reduces "improve µ" to "find a new valid convex constraint OR
a tighter dual cover." The CDE attacks both:

- **Cover side** (Phases 1, 2): no new mathematics; better placement of dual-feasibility
  centers. Free engineering wins.
- **Constraint side** (Phase 3): genuinely new constraint family (polynomial moments)
  derived from the Hausdorff moment problem on [-1, 1]. Complementary to Bochner.

Both are composable, both are rigorous under standard margin convention, and both are
implemented as small Python additions to the existing infrastructure rather than
rewrites. The compute scales linearly in the number of centers and roughly cubically
in `bochner_n`; total session cost is ~30 minutes.

## Phase 4+ candidates

- Full Hankel-PSD with both even+odd moment columns and slack variables for odd-k
- M-side Bochner composed with these (separate constraint family)
- Multi-start basin hop over (h_c, p_c) instead of binding-point-only
- Higher k_max (16, 18) for poly-moment — compound the +4e-4 single-row gain
- Even-f conditional bound (separate research direction, larger leap)
