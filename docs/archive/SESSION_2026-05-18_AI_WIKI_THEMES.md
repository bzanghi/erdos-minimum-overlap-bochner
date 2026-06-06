# SESSION 2026-05-18 — AI-wiki themes applied to min-overlap

**Headline:** Four engineering-shaped levers tested, all null. Together's `µ ≤ 0.380871` is plausibly tight on UB side (equioscillation signature); the CDE Phase 5 "alternative basis to suppress Gibbs" hypothesis is empirically refuted at the dual level; M-side Bochner — SOC and Schur — gives ≤ 10⁻⁷ benefit and is empirically dead. **No bound moves.** Remaining levers are mathematician-shaped, not engineering-shaped.

---

## What triggered the session

User pointed at the Erdős-problems wiki page on AI contributions
([teorth/erdosproblems/wiki/AI-contributions-to-Erdős-problems](https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erd%C5%91s-problems)) and asked which themes might apply to the minimum overlap problem.

Themes extracted from the wiki:

1. Lean formalization dominates fully-AI solutions — doesn't transfer to a real-constant problem like µ.
2. **AlphaEvolve = construction search**, repeatedly closes UB gaps. Problem [36] on the wiki page is min-overlap itself; AlphaEvolve already improved Haugland (2016), and Together then improved AlphaEvolve.
3. Numerical discovery → symbolic certification split.
4. Counterexample-vs-proof asymmetry favors AI on counterexamples.
5. Multi-tool ensemble for hard problems.
6. Literature synthesis remains weak.

The actionable proposal was theme (2): **stop attacking only from below; run construction search on the UB side.**

---

## PRO-24: UB-side equioscillation finding

**Question:** Can the local refiner (`refiner.py`, `basin_hop.py`, `run_search.py` in `lp_research_state/code/`) push past Together's `µ ≤ 0.380871`?

**Method.**
1. SLP polishing (`basin_hop.slp_step`, *true* minimax linearization, not smoothed surrogate) from Together's h\* at native N=600 across trust radii ∈ {0.01, 0.003}.
2. Upscale to N=1200 by `np.repeat(h_tog, 2)`, re-run SLP.
3. Active-set analysis of Together's h\*.

**Result.**
- SLP at N=600: no descent, Δ ≤ 10⁻¹⁰.
- SLP at N=1200: no descent, Δ ≤ 10⁻¹⁰.
- **457 of 2(N−1)=1198 J^±_j modes lie within 10⁻¹⁰ of M(h\*)=0.3808703106**, corresponding to ≈230 distinct shifts each tied at both `+` and `−` signs.
- h\* is highly asymmetric (`|h - h_rev|_∞ = 0.5316`), so the `+` and `−` activations are genuinely independent — not a symmetry artifact.

**Function-class follow-up (same probe step):**
- Gaussian convolution of h\* with σ ∈ {0.5, 1, 2, 3, 5, 10, 20} cells then clip+rescale: every σ raises M, monotonically from +7.2 × 10⁻⁵ at σ=0.5 cells to +4.5 × 10⁻³ at σ=20.
- Block-average downsample to N ∈ {300, 200, 150, 100, 60, 40, 30, 20, 10}: also monotone uphill, +3.2 × 10⁻⁴ at N=300 growing to +1.4 × 10⁻² at N=10.

**Interpretation.** Equioscillation across ≈230 distinct active shifts is the KKT-saturation signature of a true minimax optimum, not a search artifact. Together's h\* sits at a local minimum in both smoothness and refinement directions. The AI-wiki "AlphaEvolve closes UB construction gaps" theme is fully spent on min-overlap at the piecewise-constant function class.

**Verdict.** The whole 5.7 × 10⁻⁴ open gap [0.3803, 0.380871] is structurally LB-side.

---

## PRO-24b: alternative-basis hypothesis refuted at the dual level

**Question.** CDE Phase 5 (the LB-side saturation analysis at `µ ≥ 0.3801279`) proposed *alternative basis (wavelet / Chebyshev to suppress Gibbs)* as a candidate lever. Worth a rebuild?

**Method.** Solve at row 4 (h=0.004, p=0.3875, q ∈ [−0.02, 0.02]), N=2000, T=800, R=10, `bochner_n=20`. Extract the 42×42 PSD dual matrix Z on the `f ≥ 0` Bochner-Toeplitz constraint. Compute per-lag mass `mass[ℓ] = Σ_{|j−k|=ℓ} |Z_{jk}|` to see where in frequency the binding multipliers concentrate.

**Result.**

| lag range | mass fraction |
|---|---|
| 0–5 | **78.3%** |
| 6–15 | 20.1% |
| 16–20 | **1.6%** |

The `1 − f ≥ 0` Bochner constraint dual is essentially zero (`||A||_F = 3.5 × 10⁻¹⁰`) — that side of the family is inactive.

**Interpretation.** Gibbs lives at high Fourier frequencies; the binding dual mass lives at low Fourier lags. A basis change (Chebyshev / wavelet) re-expresses the high-frequency tail but cannot tighten constraints whose multipliers live at lags 0–5. **The CDE Phase 5 "alternative basis" lever is refuted by the dual structure.**

**Redirection.** The binding bottleneck is the *constraint family at low frequencies*, not the basis. Most plausible un-mined option in this lineage: **M-side Bochner via exact Schur lifting** (`mside_bochner_schur_n` in `white_full_convex.py:281-292`).

---

## PRO-24c: M-side Bochner Schur probe — clean null

**Question.** Does M-side Bochner — in the SOC formulation (`mside_bochner_n`) or the exact-Schur formulation (`mside_bochner_schur_n`) — provide headroom on top of f-side Bochner?

**Method.** Same row 4 / N=2000 / T=800 / R=10 scale as PRO-24b. 18 runs total:
- 2 baselines (no Bochner; f-side `bochner_n=20`)
- 8 M-side-alone (SOC × {3,5,10,20} and Schur × {3,5,10,20})
- 8 combined f-side n=20 + M-side at the same 8 settings.

**Result.**

| config | Ω\* | Δ |
|---|---|---|
| baseline (no Bochner) | 0.3762765238 | — |
| f-side `bochner_n=20` | 0.3782041527 | +1.93 × 10⁻³ |
| SOC M-side alone, n_M ∈ {3,5,10,20} | 0.3762765..3 | ≤ 3 × 10⁻⁸ vs baseline |
| Schur M-side alone, n_M ∈ {3,5,10,20} | 0.3762765..6 | ≤ 1 × 10⁻⁷ vs baseline |
| **f-side n=20 + SOC M-side n_M ∈ {3,5,10,20}** | 0.37820414..6 | **≤ 5 × 10⁻⁸ vs f-side-only** |
| **f-side n=20 + Schur M-side n_M ∈ {3,5,10,20}** | 0.37820413..5 | **≤ 5 × 10⁻⁸ vs f-side-only** |

CLARABEL noise floor at this scale is ~5 × 10⁻⁵; all M-side movement is 3 orders of magnitude below that.

**Mechanism.** The SOC slack `U_m ≥ |f̂(m)|²` (and the equivalent Schur 3×3 form, which defines the *same* convex set per the docstring of `mside_bochner_schur.py`) is one-sided. The LP optimum picks `U_m` exactly to make T_relax PSD trivially, so the M-Toeplitz PSD adds essentially no information beyond f-side Bochner. The 2026-05-10 expectation of "≤ +1 × 10⁻⁴" was an order of magnitude too generous; empirical reality is ≤ 10⁻⁷.

**The only variant that could salvage this lever** is `mside_bochner_lasserre_n` — exact-bilinear M-side using Lasserre-lifted moments, replacing the inequality with equality. That requires `lasserre_T_max > 0`, gated on the Fejér-Riesz tail-bound problem that CLAUDE.md flags as the Lasserre family's known failure mode.

---

## State of the bounds (unchanged)

| Quantity | Value | Source |
|---|---|---|
| **LB** (rigorous, Phase 5, post-margin) | `µ ≥ 0.3803027` | PRO-21, `bochner_n=40` at N=20K |
| **Framework ceiling** (asymptotic) | `µ ≥ 0.380558` | PRO-6 complementarity proof |
| **UB** (Together March-2026 certificate) | `µ ≤ 0.380871` | now plausibly tight (PRO-24) |
| **Open gap** | 5.84 × 10⁻⁴ | unchanged |
| Framework-attainable portion of gap | [0.3803027, 0.380558] = 2.55 × 10⁻⁴ (44%) | |
| Beyond-framework portion of gap | [0.380558, 0.380871] = 3.13 × 10⁻⁴ (54%) | |

---

## Levers remaining after this session

All engineering-shaped levers within the present framework now have empirical null results. What's left is mathematician-shaped:

| Lever | What it is | Why we didn't pull it this session |
|---|---|---|
| **Even-f conditional bound** | Prove or refute that the µ-optimizer can be taken even. Halves SDP dimension if true; unlocks much larger T. | Real theorem, not a code probe. Deserves its own brainstorm session. |
| **Lasserre level-2 with corrected tail bound** | Open arithmetic question whether the Fejér-Riesz tail bound, applied correctly, leaves any net gain at tractable T_max. | `communications/lasserre_tail_bound.md` derives the bound and shows it kills the gain at currently-tractable scale. Re-examination needs careful arithmetic. |
| **Combinatorial M(n) at n ≥ 50 via SAT** | Haugland (2016) reached n ≤ 43; SAT + symmetry tricks. Would beat Together's UB iff M(n)/n drops below 0.380871 at some accessible n. | Weeks of compute; per LEVER_C, all known M(n)/n for n ≤ 18 are above 0.40. |
| **Stop and write up the +5.4 × 10⁻⁴ result** | The preprint draft already exists at `communications/preprint_draft.tex`. | The honest call if no further idea materializes. |

---

## Files touched this session

- [`lp_research_state/findings.md`](lp_research_state/findings.md) — PRO-24, PRO-24b, PRO-24c folded into the rolling lab notebook entry.
- This file (`SESSION_2026-05-18_AI_WIKI_THEMES.md`) — standalone session writeup.

**No code changes.** All probes used existing tooling (`together_loader`, `fast_eval`, `basin_hop.slp_step`, `white_full_convex.build_problem`, `mside_bochner.py`, `mside_bochner_schur.py`).

---

## One-line takeaway

The Erdős-problems wiki's "AlphaEvolve closes UB gaps" theme has already been played on min-overlap (`alphaevolve_2025.py` → Together), and the lever is spent. The five-minute probe established this; the follow-on probes established that the *internal* CDE Phase 5 candidate levers — alternative basis, M-side Bochner — are also empirically dead. The remaining open gap is mathematician-territory.
