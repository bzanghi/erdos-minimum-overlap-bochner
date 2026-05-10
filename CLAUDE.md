# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A **mathematical research project** on the Erdős minimum overlap constant `µ`, not a software product. The core artifact is the research note ([erdos_lower_bound_research_note.md](erdos_lower_bound_research_note.md)) and a reproducible numerical proof of `µ ≥ 0.379544` — a `+5.4 × 10⁻⁴` improvement over White (Acta Arith. 2023). The Python code exists to *generate* and *cross-verify* that numeric.

There is no build system, no test runner, no lint config. Verification is done by **independent re-implementations agreeing to 10+ digits**, not by unit tests.

## Reproducing the headline result

```bash
cd lp_research_state/code
python3 -c "
import sys; sys.path.insert(0, '.')
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp
# Row 4 — the binding row at White's Table-3 centers
Omega, w, v, c, d, eps, dlt, cons = build_problem(
    10000, 4000, 10, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02,
    bochner_n=20,
)
prob = cp.Problem(cp.Minimize(Omega), cons)
res = solve_with_dual_extraction(prob)
print('rigorous LB:', res['rigorous_dual_LB'])  # expect ≥ 0.379653
"
```

The full ellipse-extension argument (which converts the 7 single-point row results into an unconditional bound on `µ`) is in `lp_research_state/code/path_b_rigorous.py` and `path_b_independent.py`.

## Architecture: the SDP and what gets bolted onto it

The whole project orbits one file: **[lp_research_state/code/white_full_convex.py](lp_research_state/code/white_full_convex.py)**. Its `build_problem(...)` constructs White's §5 Fourier-analytic convex program in cvxpy, with augmentations toggled by parameters:

- `bochner_n=20` — adds the Bochner moment-matrix PSD constraint `M_n(f) ⪰ 0` and `M_n(1−f) ⪰ 0`. This is the rigorous improvement.
- `mside_bochner_n=k` — adds an SOC-relaxed M-side Bochner constraint (also rigorous).
- `use_T5p=True` — enables tightening 5' (`f² ≤ f` tested against `1 − cos(πx)`).
- `use_T5`, `use_T3` — older tightenings, mostly subsumed.

Helper modules:
- **`dual_extractor.py`** — runs the SDP with `verbose=True`, parses CLARABEL's iteration log, and returns `rigorous_dual_LB = reported_value − last_gap`. **This is the project's central epistemic trick.** CLARABEL's status `optimal_inaccurate` is a labeling artifact (true gaps are ~10⁻⁷); the dual objective at any iteration with small dual residual is a rigorous LB. Always prefer `rigorous_dual_LB` over `value` when reporting bounds — it recovers ~+1 × 10⁻⁴ over the safety convention.
- **`bochner.py`** vs **`bochner_independent.py`** — the encoded constraint and an independent re-implementation by a separate agent who didn't see the first. They agree bit-for-bit; if you modify one, expect to verify against the other.
- **`path_b_analytical.py`**, **`path_b_rigorous.py`**, **`path_b_independent.py`** — three independently-written implementations of White's Section 5.1 / Appendix II ellipse-extension argument. They agree to 10+ digits; same cross-check policy applies.
- **`lasserre.py`** — Lasserre level-2 augmentation. **Documented as non-rigorous and withdrawn.** Truncates `(f²)̂(m)` without a tail bound; `communications/lasserre_tail_bound.md` derives the natural Fejér-Riesz tail bound and shows it kills the gain quantitatively at currently-tractable `T_max`. Don't quote any Lasserre-augmented value as a rigorous bound.

## The 7 "rows"

White (2023) covers the residual parameter region `(5.16)` with 7 ellipses; we run the SDP at each ellipse *center*. Codenames `row1..row7` refer to those centers (defined in `lp_research_state/cron_runner.py:POINTS_BY_ROW`). **Row 4** (`(h, p, q) = (0.004, 0.3875, ±0.02)`) is the *binding* row — it produces the MIN over rows. Row 5 is memory-heavier (use `bochner_n=15` instead of 20 to fit in 4 GB).

## Critical caveat — read before claiming any bound

**A single-point SDP solve at row centers is NOT a rigorous bound on `µ`.** White's `µ ≥ 0.379005` covers the full `(h, p, q)` ranges via the §5.1 ellipse-extension argument. To convert per-row dual objectives into an unconditional lower bound on `µ`, you must run the `path_b_*` machinery, which checks that the 7 ellipses (recomputed with augmented duals) still cover the residual region. The retracted Lasserre claim aside, this is the most common source of overclaiming in this repo's history. **`findings.md` distinguishes "MIN over 7 row centers" from "rigorous bound on `µ`" — preserve that distinction in any new analysis.**

## The cron-driven experiment loop

`lp_research_state/cron_runner.py` is a single-experiment driver. The workflow:

1. `lp_research_state/experiments_queue.json` lists pending experiments (kind, params, priority).
2. Run `python3 lp_research_state/cron_runner.py` to execute the highest-priority pending one. Use `--dry-run` to preview without solving.
3. Results append to `experiments_done.json`. The driver dedupes via `is_done(...)` matching on `(kind, N, T, R, row, bochner_n, ...)`.
4. **`lp_research_state/findings.md`** is the human-curated rolling log. The leading line is updated each run with the latest finding; older findings are demoted into the body. Treat it as a research lab notebook, not documentation.

Supported experiment kinds: `lp_run`, `lp_run_bochner`, `lp_run_bochner_dual`, `lp_run_bochner_sweep`, `lp_run_mside_bochner`. Other kinds (`infra`, `alpha_sweep`) are flagged for manual handling.

The driver is path-tolerant: it auto-resolves `LP_STATE_DIR` / `LP_CODE_DIR` across several mount points (legacy session paths, the canonical repo path). When working from a new mount, set those env vars rather than editing the resolution list.

## Files prefixed `_` are throwaways

`code/_*.py` files (`_lasserre_scan.py`, `_run_one_rigorous.py`, etc.) are one-shot scratchpads that wired up specific cron sweeps. They're kept for provenance but should not be imported as libraries. New scratchpads should follow the same `_` convention.

## When proposing improvements

Pushing the rigorous bound past `0.379544` with this program's current technique set (Bochner-PSD + ellipse extension) is **not possible** at currently-tractable SDP scale — see the post-tail-bound update at the top of `erdos_lower_bound_research_note.md`. New levers (much larger `T_max > 1000`, finite-dimensional SOS exactness, alternative basis representations) are research questions, not engineering tasks. Don't propose more Lasserre level / Bochner level scans as a path forward without first reading that section.

## When making changes to the SDP encoding

If you modify `bochner.py`, `white_full_convex.py`'s constraint section, or `path_b_*.py`:
1. Run the corresponding `_independent` file's path on at least one row.
2. Compare to 10+ significant digits. Any disagreement is a bug, not a rounding artifact.
3. Spot-check with SDPA-GMP at small `N` if the change touches numerical precision.

## External resources

- White (2023) Acta Arith.: [arXiv:2201.05704](https://arxiv.org/abs/2201.05704) — the program being augmented.
- Together Computer (March 2026) upper bound: [GitHub](https://github.com/togethercomputer/erdos-minimum-overlap) — gives `µ ≤ 0.380871`. Open gap is `[0.379544, 0.380871]`, width ~1.3 × 10⁻³.
- Preprint draft: [communications/preprint_draft.tex](communications/preprint_draft.tex). Email draft to E. P. White: [communications/email_to_ethan_white.md](communications/email_to_ethan_white.md).
