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
- **Current bracket (re-verified 2026-07-25, second pass): `µ ∈ [0.3802838, 0.3808590568]`.** See [MINIMUM_OVERLAP_STATE_2026-07-25b.md](MINIMUM_OVERLAP_STATE_2026-07-25b.md), which supersedes the earlier same-day state doc. The UB side is an active AI-search benchmark and moves weekly — **re-check before quoting any bracket.**
  - **`0.380856` IS NOT A BOUND. Do not quote it.** The SimpleTES ablation value ([arXiv:2604.19341](https://arxiv.org/abs/2604.19341) §3.4.1) is a normalization artifact: the evolved program reports its bin count as `float(n) + 0.9999999999` and divides by that, so `dx` is understated by ~1 part in 4096. Downloading the witness (commit `406fc651`) and evaluating it honestly gives **0.3809489501030183** — worse than Together's. The SimpleTES authors fixed this themselves in commit `6eb2ca0a` ("fix a potential hack possibility with n_points not being integer"); the arXiv paper is still unrevised at v1.
  - Best UB with a public, verified, downloadable witness: **`µ ≤ 0.3808590568145606`**, Einstein Arena entry `lnzwz_AI4M_Agent` (512 cells). Fetch it with
    `curl -s "https://einsteinarena.com/api/solutions/best?problem_id=1&agent_name=lnzwz_AI4M_Agent&limit=1"`.
    Exactly certified here: `µ ≤ 0.380859056814560651295303328196` (`data/ub_certified_arena_n512.json`).
  - Together Computer (March 2026) `µ ≤ 0.3808703105862199` is **superseded** but is still `ub_core.ANCHOR` and the anchor in older docs. PRO-34's "0.380871 is a serious candidate for µ itself" is falsified.
  - Best *published* LB is no longer White's 0.379005 but `µ ≥ 0.37912`, Kim & Pilanci [arXiv:2606.31182](https://arxiv.org/abs/2606.31182) (ICML 2026) — same interval-arithmetic rigor standard as ours. This repo still leads it by >1.1 × 10⁻³, but priority is time-sensitive.
  - **Every AI-search witness examined here (3 of 3) fails exact feasibility** by 1e-16 to 4e-14 in mass; none of the producing systems check their arithmetic outside float64. Run `ub_certify.py` on anything before quoting it.

- **`dual_extractor.rigorous_dual_LB` is NOT a certificate** (fixed/documented 2026-07-25). It is the solver's dual objective with no correction for dual infeasibility, and its eligibility gate used to read CLARABEL's `pres` column while calling it the dual residual. For a real bound use `_jansson_verify.jansson_lower_bound` (interval arithmetic), driven for the core anchors by `_jansson_reanchor.py`.
- Preprint draft: [communications/preprint_draft.tex](communications/preprint_draft.tex). Email draft to E. P. White: [communications/email_to_ethan_white.md](communications/email_to_ethan_white.md).

## Extra math tooling (added 2026-05-18)

In addition to the cvxpy/CLARABEL stack, the following project-local tooling is available:

- **mpmath, sympy** (in `.venv`) — arbitrary-precision arithmetic and symbolic algebra. PSLQ via `mpmath.pslq`. Already used in [lp_research_state/code/pslq_hunt.py](lp_research_state/code/pslq_hunt.py), which sweeps a basis of standard constants {1, π, e, log 2, √2, ζ(2), Γ(1/4), ...} looking for integer relations with the LB/UB headlines.

- **SDPA-GMP at [lp_research_state/bin/sdpa_gmp](lp_research_state/bin/sdpa_gmp)** — arbitrary-precision SDP solver built from [sdpa-python/sdpa-multiprecision](https://github.com/sdpa-python/sdpa-multiprecision). Smoke-tested via [lp_research_state/code/sdpa_gmp_wrapper.py](lp_research_state/code/sdpa_gmp_wrapper.py); example1.dat-s solves to feasibility error ~10⁻⁷⁵ vs CLARABEL's ~10⁻⁷. Build details: `/tmp/sdpa_build/` contains source + custom-built GMP 6.3.0. Two upstream patches were required for modern clang on macOS: SPOOLES `IVinit(nfront, NULL)` → `IVinit(nfront, 0)` (applied via a sed step injected into `spooles/Makefile`); and `sdpa_struct.cpp` C++11 user-defined-literal spacing around the `P_FORMAT` macro. The `libsdpa_gmp.a` static library is also installed at `lp_research_state/bin/`. **What's missing for production use:** a cvxpy → SDPA-S serializer so White's SDP can be cross-verified at GMP precision. Currently the wrapper handles arbitrary .dat-s files but doesn't yet auto-translate from cvxpy. This is the next engineering step for PRO-11.

- **Wolfram Alpha integration** via [lp_research_state/code/wolfram_alpha.py](lp_research_state/code/wolfram_alpha.py). Uses the LLM API endpoint `https://www.wolframalpha.com/api/v1/llm-api` (per [Anthropic's cookbook](https://platform.claude.com/cookbook/third-party-wolframalpha-using-llm-api)) as the default `query_llm()`; falls back to the Simple API (`/v2/result`) and Full JSON API (`/v2/query`). Requires `WOLFRAM_APP_ID` in env or `<repo>/.env`. Useful for closed-form integral verification, inverse-symbolic lookup on bracket endpoints, and special-function manipulations sympy can't reduce.

- **arXiv search** via [lp_research_state/code/arxiv_search.py](lp_research_state/code/arxiv_search.py). Free, no key needed; uses the export API. Use it before claiming priority on any new lever — the Erdős literature has had occasional bursts of activity post-2024.

The MCP registry (as exposed to this client) has no math/research connectors as of 2026-05-18; we surveyed wolfram, mathematica, arxiv, sage, sympy, lean, coq, solver, optimization, latex, math, science, jupyter — all empty. The above tools cover the gap by combining `Bash` + `WebFetch` + local installs.
