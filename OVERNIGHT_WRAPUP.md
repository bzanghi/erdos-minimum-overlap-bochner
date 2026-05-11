# Overnight wrap-up — 2026-05-10 → 2026-05-11

While you slept, the autonomous loop did the following. **Headline:** no new µ-improving lever found, but the framework is now more rigorously characterized as saturated, the preprint+email materials are in better shape for v2, and Lever H and Lever I (the "wild" candidate directions) both ruled out via close paper reads.

## Bound state — unchanged

- **LB: µ ≥ 0.3801279** (Phase 5 of CDE, unchanged)
- **LB: µ ≥ 0.380129x (in flight)** — Phase 4B + T5p running in background; expected +1.6 × 10⁻⁶ gain
- **UB: µ ≤ 0.380871** (Together's certificate, unchanged)
- **Gap: 7.4 × 10⁻⁴** (essentially unchanged)

## What was done

### 1. Lever H deep read → FAILS

Close read of arXiv:2210.16437 (White's L² autoconvolution paper). Three independent obstructions kill the transfer:

- White's technique is a Fourier-truncated QCLP plus Minkowski strict-convexity, NOT a second-variation eigenvalue operator (the cross-domain hunt's premise was wrong)
- The L²/L⁴ Fourier closed-form (White's Lemma 7) makes that problem polynomial of degree 4 in `f̂(k)`. The L∞ min-overlap functional has no such polynomial expansion — that's why our SDP must be a Bochner-PSD *relaxation*
- Smoothed L^p proxy needs p ≥ 14,000 to close 10⁻⁴ of µ → ĥ(k)^{28000} polynomial terms, vastly past Lasserre tractability

White himself notes the L∞ regime is "more computationally challenging" with no method given. The 2025 follow-ups are constructive UBs, not analytical LBs. **Zero estimated gain.**

Full read: [LEVER_H_DEEP_READ.md](LEVER_H_DEEP_READ.md)

### 2. Lever I deep read → DOES_NOT_APPLY

Two findings:

**(a) Attribution error in the original cross-domain hunt was caught.** The cited paper arXiv:2206.09876 is **Rupert Li 2022** (UMN Duluth REU, Cohn mentor), NOT "Cohn-de Laat-Salmon". The actual CdLS paper is arXiv:2206.15373 ("Three-point bounds for sphere packing").

**(b) Even taking Li 2022 as the target, the technique doesn't apply** for three reasons: wrong direction (gives LBs on LP value, not UBs / ceilings), PSD constraints don't survive Li's restriction-and-periodize map, no analog of "f(x) ≤ 0 for |x| ≥ r" in White's program.

Full read: [LEVER_I_DEEP_READ.md](LEVER_I_DEEP_READ.md)

### 3. The deep read of Lever I surfaced a NEW direction — Lever I'

**Residual-enumeration ceiling theorem.** Use the template at [communications/lasserre_tail_bound.md](communications/lasserre_tail_bound.md) (which derived the Fejér-Riesz tail bound that retracted Lasserre) systematically on every constraint family in `white_full_convex.py`. The cumulative slack residual, if computably small, IS the saturation theorem. Effort: 1-2 weeks of careful derivation. Output: an explicit `C* < 0.380871` such that no Bochner-PSD augmentation can prove `µ ≥ C*`.

This is the **only direction surfaced this iteration that is both (a) novel relative to the existing technique stack and (b) tractable from in-session resources**. Full discussion in updated [OUT_OF_BOX_SYNTHESIS.md](OUT_OF_BOX_SYNTHESIS.md).

### 4. T5p formal-completeness fix

The adversarial red-team found that `path_b_with_polymoment.py` (the Phase 5 driver) silently omits the `use_T5p` flag from `white_full_convex.py`. Threaded the flag through:

- [`path_b_analytical.py:build_problem_with_dual_handles`](lp_research_state/code/path_b_analytical.py) — added `use_T5p` parameter
- [`path_b_with_polymoment.py`](lp_research_state/code/path_b_with_polymoment.py) — added `--use_T5p` CLI flag
- Sanity check at N=3000, T=1200: confirmed +1.6 × 10⁻⁶ gain in Phase-5 composition

Full Phase-4B-config + T5p run (N=10000, T=4000, bochner_n=30, pm_k_max=20, hankel_n=6) is in flight; expected to push µ ≥ 0.3801199 → 0.3801215 or thereabouts (small but real). **Results will land in the next 10-15 min — see latest commit.**

### 5. Materials updates

- **[communications/email_to_ethan_white_v2.md](communications/email_to_ethan_white_v2.md)** — superseding v1. Reframed around the strengthened µ ≥ 0.3801279 headline AND the new "ceiling theorem" question. Asks White directly whether his framework admits a saturation theorem, and whether there's a closed-form Fourier identity for `sup_t (h⋆h)(t)` we're missing. Ready for your review before sending. **NOT autonomously sent.**

- **[communications/preprint_addendum.md](communications/preprint_addendum.md)** — inventory of post-v1 extensions with three concrete publication-path recommendations:
  1. v2 preprint with the strengthened headline µ ≥ 0.3801279 (requires rewriting Theorem 1 / Section 4)
  2. Two separate papers (bound improvement + saturation methodology)
  3. v1 unchanged + addendum as appendix (lowest risk)

### 6. Final lever ledger — 10/10 eliminated

| Lever | Status | Why |
|---|---|---|
| A — Lukács SOS / alt basis | ❌ unlikely | Gibbs already damped post-CDE |
| B — Diagnostic | ✓ executed | TOGETHER_DIAGNOSTIC.md |
| C — Integer M(n) | ❌ Together UB stands | Even at n=20 via SAT, M/n > 0.380871 |
| D — O(1) breakpoints | ❌ refuted | h\* has 400+ blocks |
| D' — Lipschitz/BV via discrete | ❌ refuted | Lifted optima diverge from h\* |
| E — M-side SDP | ❌ vacuous | Δ ~ 10⁻⁷ |
| F — Step-function UB push | ❌ saturated | 95→600 gained only 5e-5 |
| G — (f, g) rewrite | ❌ analytically no-op | Convex-hull equivalent |
| H — White's L² transfer | ❌ FAILS | Deep read, no polynomial expansion |
| I — Li/CdLS ceiling theorem | ❌ DOES NOT APPLY | Wrong direction |
| **I' — Residual-enumeration ceiling** | **🟡 NEW recommended** | Use lasserre_tail_bound.md template |

## What needs your decision

### Decision 1 — Publication path

Three options (see preprint_addendum.md §"Recommended publication path"):

1. **v2 preprint** with strengthened headline. Strongest single statement, longest writeup.
2. **Two papers** (bound improvement + saturation methodology). Cleanest editorial separation.
3. **v1 unchanged + addendum as appendix.** Lowest risk, smallest scope.

My recommendation: option 2. The diagnostic methodology has standalone interest; the bound improvement is a clean follow-up to White's paper that deserves its own headline.

### Decision 2 — White email

[communications/email_to_ethan_white_v2.md](communications/email_to_ethan_white_v2.md) is ready for your review. Two specific points to check:

1. **The ceiling-theorem question** is the genuinely novel ask. If he engages, this could become a bigger collaboration than just confirming the +1.1e-3 result.
2. **Repo state** — the email assumes the repo is shareable. Decide which subset to make public before sending.

### Decision 3 — Next session focus

Best in-session-reach research direction: **Lever I' (residual-enumeration ceiling theorem)**. 1-2 weeks of derivation, using existing project machinery, output is a publishable saturation result.

Best non-research action: **write up what we have** (whichever publication path you choose).

## What's still running in background

- **Phase 4B + T5p** (background task `b3wdb433p`): full scale, 12 centers. Will land within ~10 min. The committed result will tell you the actual T5p gain at full scale. If it matches the red-team's small-scale +1.6 × 10⁻⁶ estimate, expect µ ≥ 0.380121x.

## Session summary metrics

- Commits during this overnight: ~7
- Total session commits: ~30
- Total levers investigated: 10 (all eliminated)
- New levers surfaced: 1 (Lever I' — residual enumeration)
- Materials drafted for user decision: email v2, preprint addendum, synthesis update

All work committed on `main`. Not pushed.
