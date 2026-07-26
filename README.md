# Erdős minimum overlap problem — improved lower bound on µ

Code, numerics, and research notes for a rigorous improvement on White's
(Acta Arith. 2023) lower bound for the Erdős minimum-overlap constant µ, built
by augmenting White's §5 convex program with Bochner moment-matrix PSD
constraints and polynomial-moment cuts, then certifying coverage over White's
full parameter space.

## Result

$$\boxed{\mu \;\geq\; 0.3803954}\quad\text{over White's \emph{entire} }(E(M),\,c_1,\,d_1)\text{ parameter space}$$

vs. White (2023)'s `µ ≥ 0.379005` and the current best *published* bound `µ ≥ 0.37912` (Kim & Pilanci, ICML 2026) — an **unconditional, full-space improvement of +1.28 × 10⁻³** over the published record (`+1.39 × 10⁻³` over White), with **no White-published number used in the bound** (every region is certified by our own augmented dual cover).

- **Core residual region (5.16):** `µ ≥ 0.3803954`. All 12 anchors carry a Jansson-Chaykin-Keil interval-arithmetic certificate (`_jansson_reanchor.py`) at `N=48000`, with the duals read from the same solve that produced the certificate. This replaced the earlier `primal − 1e-5` convention, which was a haircut on a solver-reported value rather than a theorem — and the certified anchors came out *above* it, so the fix raised the bound rather than costing anything (0.3802838 → 0.3802946 → 0.3803954).
- **Full-space promotion (re-run 2026-07-26):** the augmented dual cover clears `0.3803954` over all 18 of White's Table-2 "outside" regions, so the **binding region is the core** — R16 0.3803961, R9 0.3803979, R17 0.3804045, R6 0.3804601, R7 0.3805539. Margins on the three tightest are thin (+7.5e-7, +2.6e-6, +9.2e-6). Reproduce in one command:

  ```bash
  cd lp_research_state/code && LP_DUALEXT=../parallel_results/dualext_reanchored_N48000.json LP_TARGET=0.3803954 ../../.venv/bin/python _fs_recompute.py
  ```

  See [`MINIMUM_OVERLAP_STATE_2026-07-25b.md`](MINIMUM_OVERLAP_STATE_2026-07-25b.md) §7–§8.
- **Upper bound, certified here:** `µ ≤ 0.380859056614806899090596051448` (exact rational, all `2n−1` signed lags in integer arithmetic). Note the widely-cited `0.380856` **is not a bound** — see [`MINIMUM_OVERLAP_STATE_2026-07-25b.md`](MINIMUM_OVERLAP_STATE_2026-07-25b.md) §2.

**Honest caveats (these travel with the bound):**
- It is **load-bearing on the polynomial-moment cuts** (`pm_k_max=20`), which are rigorous as of the 2026-05-22 tail-bound fix (see [`lp_research_state/findings.md`](lp_research_state/findings.md)), and on a set of fresh "promotion" centers in regions R16/R17 (with the 12 core anchors alone, those corners fall to 0.3802561, −2.8 × 10⁻⁵ below target).
- **Margins are thin**: the core binds, and the three tightest outside regions clear it by only +7.5 × 10⁻⁷ (R16), +2.6 × 10⁻⁶ (R9) and +9.2 × 10⁻⁶ (R17). Any further core gain runs into them immediately, and must re-run R16 and R9 at the raised target in the same change. Farkas certificates for the (non-load-bearing) infeasibility exclusions remain outstanding.
- **Region floors are target-limited, not infima.** The adaptive evaluators stop as soon as a sub-box clears their target, so each floor means "at least this". The core floor is also grid-dependent: `0.3803953504` at `n_grid=801`, but `0.3803953255` at 401 — always quote the resolution with the number.
- A prior Bochner-only headline was `µ ≥ 0.379544`; an earlier Lasserre-level-2 extension was **retracted** (truncated moment expansion without a tail bound). Both lessons are recorded in the research note.

**Author validation.** E. P. White (the author of the program we augment) confirmed (2026-05-31) that the Bochner-PSD constraint is "a valid constraint to add," and supplied two corrections to his published program: constraints 5.6/5.7 should have a `4` (not `8`) in the RHS numerator — **applied** (`mside_sin_coeff=4.0`; impact verified neutral, PRO-43) — and 5.8/5.9 should use `2m−1`, which our code already did.

## Method

1. **Base program** — White's §5 Fourier-analytic convex program ([`lp_research_state/code/white_full_convex.py`](lp_research_state/code/white_full_convex.py), `build_problem(...)`).
2. **Bochner-PSD augmentation** — add `M_n(f) ⪰ 0` and `M_n(1−f) ⪰ 0` (`bochner_n`); the rigorous core improvement.
3. **Polynomial-moment cuts** — `m_{2k} ≥ −tail_bound_k` from the Hausdorff moment theorem, with an analytic tail remainder ([`poly_moment.py`](lp_research_state/code/poly_moment.py)).
4. **Dual cover + ellipse extension** — each center's dual objective is a globally-valid lower bound; the cover is `max_c Φ_c(h,p,q)` ([`path_b_analytical.py`](lp_research_state/code/path_b_analytical.py)).
5. **Full-space certification** — rigorous box-min via grid + Lipschitz `eps_grid`, with **adaptive subdivision** to control `eps_grid` on White's wide outside regions ([`_fullspace_eval.py`](lp_research_state/code/_fullspace_eval.py)).

Verification convention: independent re-implementations agreeing to 10+ digits, plus a
Jansson-Chaykin-Keil a-posteriori certificate in directed-rounding interval arithmetic
([`_jansson_verify.py`](lp_research_state/code/_jansson_verify.py)) for every anchor in the
binding region. **`dual_extractor.rigorous_dual_LB` is NOT a certificate** — it is the
solver's dual objective with no correction for dual infeasibility, and its eligibility gate
formerly read CLARABEL's `pres` column while calling it the dual residual (both documented and
fixed 2026-07-25). Use it to steer a search, never to state a bound.

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

- **Lower bound (this work):** `µ ≥ 0.3803954`, full-space, every core anchor certified in
  interval arithmetic. Best *published* LB is `0.37912` (Kim & Pilanci, ICML 2026), which
  displaced White's `0.379005`.
- **Upper bound:** `µ ≤ 0.380859056614806899090596051448`, certified here in exact rational
  arithmetic from the Einstein Arena `lnzwz_AI4M_Agent` witness (n=512, publicly
  downloadable) after a structured basin-hop polish.
- **Do not quote `0.380856`.** It is a normalization artifact, not a bound; the honest value
  of that construction is `0.3809490`. See
  [`MINIMUM_OVERLAP_STATE_2026-07-25b.md`](MINIMUM_OVERLAP_STATE_2026-07-25b.md) §2.
- **Open gap:** `µ ∈ [0.3803954, 0.3808591]`, width ≈ **5.64 × 10⁻⁴**.
- The UB side is an active AI-search benchmark that moves weekly — re-check before quoting.

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
