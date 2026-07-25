# Erdős minimum overlap problem — improved lower bound on µ

Code, numerics, and research notes for a rigorous improvement on White's
(Acta Arith. 2023) lower bound for the Erdős minimum-overlap constant µ, built
by augmenting White's §5 convex program with Bochner moment-matrix PSD
constraints and polynomial-moment cuts, then certifying coverage over White's
full parameter space.

## Result

$$\boxed{\mu \;\geq\; 0.3802946}\quad\text{over White's \emph{entire} }(E(M),\,c_1,\,d_1)\text{ parameter space}$$

vs. White (2023)'s `µ ≥ 0.379005` and the current best *published* bound `µ ≥ 0.37912` (Kim & Pilanci, ICML 2026) — an **unconditional, full-space improvement of +1.17 × 10⁻³**, with **no White-published number used in the bound** (every region is certified by our own augmented dual cover).

- **Core residual region (5.16):** `µ ≥ 0.3802946016`, the **binding** region. All 12 anchors carry a Jansson-Chaykin-Keil interval-arithmetic certificate (`_jansson_reanchor.py`), with the duals read from the same solve that produced the certificate. This replaced the earlier `primal − 1e-5` convention, which was a haircut on a solver-reported value rather than a theorem — and the certified anchors came out *above* it, so the fix also raised the bound (0.3802838 → 0.3802946).
- **Full-space promotion (verified, re-run 2026-07-25):** the augmented dual cover (121 dual-feasible centers) clears the core value over all 18 of White's Table-2 "outside" regions, so the full-space minimum equals the core. Tightest outside region R6 at 0.3803090. See [`lp_research_state/FULLSPACE_VERIFICATION.md`](lp_research_state/FULLSPACE_VERIFICATION.md).
- **Upper bound, certified here:** `µ ≤ 0.380859056651254094841565810818` (exact rational, all `2n−1` signed lags in integer arithmetic). Note the widely-cited `0.380856` **is not a bound** — see [`MINIMUM_OVERLAP_STATE_2026-07-25b.md`](MINIMUM_OVERLAP_STATE_2026-07-25b.md) §2.

**Honest caveats (these travel with the bound):**
- It is **load-bearing on the polynomial-moment cuts** (`pm_k_max=20`), which are rigorous as of the 2026-05-22 tail-bound fix (see [`lp_research_state/findings.md`](lp_research_state/findings.md)), and on a set of fresh "promotion" centers in regions R16/R17 (with the 12 core anchors alone, those corners fall to 0.3802561, −2.8 × 10⁻⁵ below target).
- **Margins are thin** (binding outside region R16 clears by +1.2 × 10⁻⁴). A margin-hardening re-solve at N ≥ 24000 + Farkas certificates for the (non-load-bearing) infeasibility exclusions is in progress (PRO-44).
- A prior Bochner-only headline was `µ ≥ 0.379544`; an earlier Lasserre-level-2 extension was **retracted** (truncated moment expansion without a tail bound). Both lessons are recorded in the research note.

**Author validation.** E. P. White (the author of the program we augment) confirmed (2026-05-31) that the Bochner-PSD constraint is "a valid constraint to add," and supplied two corrections to his published program: constraints 5.6/5.7 should have a `4` (not `8`) in the RHS numerator — **applied** (`mside_sin_coeff=4.0`; impact verified neutral, PRO-43) — and 5.8/5.9 should use `2m−1`, which our code already did.

## Method

1. **Base program** — White's §5 Fourier-analytic convex program ([`lp_research_state/code/white_full_convex.py`](lp_research_state/code/white_full_convex.py), `build_problem(...)`).
2. **Bochner-PSD augmentation** — add `M_n(f) ⪰ 0` and `M_n(1−f) ⪰ 0` (`bochner_n`); the rigorous core improvement.
3. **Polynomial-moment cuts** — `m_{2k} ≥ −tail_bound_k` from the Hausdorff moment theorem, with an analytic tail remainder ([`poly_moment.py`](lp_research_state/code/poly_moment.py)).
4. **Dual cover + ellipse extension** — each center's dual objective is a globally-valid lower bound; the cover is `max_c Φ_c(h,p,q)` ([`path_b_analytical.py`](lp_research_state/code/path_b_analytical.py)).
5. **Full-space certification** — rigorous box-min via grid + Lipschitz `eps_grid`, with **adaptive subdivision** to control `eps_grid` on White's wide outside regions ([`_fullspace_eval.py`](lp_research_state/code/_fullspace_eval.py)).

Verification convention: independent re-implementations agreeing to 10+ digits, and `rigorous_dual_LB = value − last_gap` (dual extraction), not unit tests.

## Reproducing

See **[`REPRODUCE.md`](REPRODUCE.md)** for step-by-step recipes (core headline, full-space verification, the 8→4 correction check). Quick smoke (binding row 4):

```python
import sys; sys.path.insert(0, "lp_research_state/code")
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp
Omega, w, v, c, d, eps, dlt, cons = build_problem(
    10000, 4000, 10, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02, bochner_n=20)
res = solve_with_dual_extraction(cp.Problem(cp.Minimize(Omega), cons))
print(res["rigorous_dual_LB"])  # ≥ 0.379653 (Bochner-only, single row)
```

## State of the art

- **Lower bound (this work):** `µ ≥ 0.380284`, full-space, verified. vs White (2023) `0.379005`.
- **Upper bound:** `µ ≤ 0.380871` — Together Computer (March 2026); verified independently to `0.3808703106…`.
- **Open gap:** `µ ∈ [0.380284, 0.380871]`, width ≈ **5.87 × 10⁻⁴**.

## Key documents

- [`lp_research_state/findings.md`](lp_research_state/findings.md) — rolling research ledger (leading line = latest result).
- [`lp_research_state/FULLSPACE_VERIFICATION.md`](lp_research_state/FULLSPACE_VERIFICATION.md) — the full-space promotion verification of record.
- [`lp_research_state/WHITE_EMAIL_CORRECTION.md`](lp_research_state/WHITE_EMAIL_CORRECTION.md) — the 8→4 correction analysis.
- [`erdos_lower_bound_research_note.md`](erdos_lower_bound_research_note.md) — the main research note.
- [`communications/`](communications/) — correspondence with E. P. White, preprint draft.

## Citations
- E. P. White, "A new bound for Erdős' minimum overlap problem," *Acta Arith.* 208 (2023). [arXiv:2201.05704](https://arxiv.org/abs/2201.05704).
- J. K. Haugland, "A new upper bound on the constant in the Erdős minimum overlap problem," 2016. [arXiv:1609.08000](https://arxiv.org/abs/1609.08000).
- Together Computer, "New State-of-the-Art on Erdős' Minimum Overlap Problem" (March 2026). [GitHub](https://github.com/togethercomputer/erdos-minimum-overlap).

## Acknowledgements

Carried out in collaboration with Claude (Anthropic). Independent sub-agent
re-implementations cross-checked every headline number to 10+ digits; results
are reported with their load-bearing dependencies stated explicitly.
