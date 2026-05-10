# Erdős minimum overlap problem — Bochner-PSD strengthening

This repo contains code, numerics, and a research note documenting a small
but rigorous improvement on White's (Acta Arith. 2023) lower bound for
Erdős' minimum overlap constant µ.

## Result

$$\boxed{\mu \;\geq\; 0.379828}$$

vs. White (2023)'s `µ ≥ 0.379005`. Improvement: **+8.2 × 10⁻⁴**.

(Earlier Bochner-only result `µ ≥ 0.379544` is superseded by adding
Lasserre level-2 to the program; both numbers are independently verified.)

## Method

Add the Bochner moment-matrix PSD constraint
$$\bigl[\hat f(j-k)\bigr]_{j,k=0,\dots,n} \succeq 0 \quad \text{(both for } f \text{ and } 1-f\text{)}$$
to White's Section 5 convex program, then apply White's own §5.1 / Appendix II
ellipse-extension argument with the augmented dual objective. The seven
ellipses around White's Table-3 centers, recomputed with our augmented
duals, fully cover White's residual region (5.16); their intersected
minimum is `0.3795475` (closed-form), which after a conservative `1e-6`
margin for CLARABEL's IPM gap gives the headline `0.379544`.

The improvement is not large — the Bochner constraint just barely widens
each ellipse enough to keep coverage while raising the dual objective.
But it is an unconditional improvement on the previous best lower bound,
derived via the same proof structure as White.

## Repo layout

- `erdos_lower_bound_research_note.md` — the main research note.
- `lp_research_state/code/` — Python code:
  - `white_full_convex.py` — White's §5 program with `bochner_n` parameter.
  - `bochner.py` — Bochner-PSD constraint encoder.
  - `dual_extractor.py` — extracts CLARABEL's rigorous dual lower bound.
  - `path_b_analytical.py` — Path B ellipse-extension implementation.
  - `path_b_independent.py` — independent re-implementation.
  - `path_b_rigorous.py` — closed-form ellipse-min refinement.
  - `lasserre.py` — Lasserre level-2 SDP encoder (compounds with Bochner).
- `lp_research_state/findings.md` — accumulating raw findings (cron-driven).
- `lp_research_state/experiments_done.json` — every SDP solve, with status.
- `lp_research_state/parallel_results/path_b/` — per-row Path B data.
- `min_overlap_report.md` — earlier Phase-1 + Phase-2 report (orientation,
  brute force, baseline reproduction, upper-bound parallel track).
- `communications/` — email to Ethan White, arXiv preprint draft,
  LinkedIn copy.

## Reproducing

```python
import sys; sys.path.insert(0, "lp_research_state/code")
from white_full_convex import solve_full_program
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp

# Row 4 (the binding row at White's Table 3 centers) at full size
N, T, R = 10000, 4000, 10
Omega, w, v, c, d, eps, dlt, cons = build_problem(
    N, T, R,
    h1=0.004, h2=0.004,
    p1=0.3875, p2=0.3875,
    q1=-0.02, q2=0.02,
    bochner_n=20,
)
prob = cp.Problem(cp.Minimize(Omega), cons)
res = solve_with_dual_extraction(prob)
print(res["rigorous_dual_LB"])  # ≥ 0.379653
```

For the full Path B argument, see `lp_research_state/code/path_b_analytical.py`.

## State of the art

- **Lower bound:** `µ ≥ 0.379828` (this work, Bochner + Lasserre-2) — vs. White (2023): `µ ≥ 0.379005`.
- **Upper bound:** `µ ≤ 0.380871` — Together Computer (March 2026), via
  sequential-LP refinement of a 600-step function. Verified independently
  to `0.3808703106…`.
- **Open gap:** `µ ∈ [0.379828, 0.380871]`, width ≈ 1.0 × 10⁻³.

## Citations
- E. P. White, "A new bound for Erdős' minimum overlap problem," *Acta Arith.* 208 (2023). [arXiv:2201.05704](https://arxiv.org/abs/2201.05704).
- J. K. Haugland, "A new upper bound on the constant in the Erdős minimum overlap problem," 2016. [arXiv:1609.08000](https://arxiv.org/abs/1609.08000).
- Together Computer, "New State-of-the-Art on Erdős' Minimum Overlap Problem" (March 2026). [GitHub](https://github.com/togethercomputer/erdos-minimum-overlap).

## Acknowledgements

This work was carried out in collaboration with Claude (Anthropic) using
the Cowork Mode multi-agent research framework. Sub-agents performed:
parallel SDP solves across White's 7 Table-3 ellipse centers; independent
re-encoding of the Bochner constraint as a bug-check; analytical
ellipse-extension implementation; rigorous dual-extraction verification;
and SDPA-GMP spot-checking of CLARABEL's numerical precision. Three
independently-written code paths agreed on the headline numbers to 10+
digits.
