# Out-of-box wild investigation — synthesis

Result of 4 parallel "go wild" investigations launched after the standard diagnostic chain had ruled out all 8 candidate levers (A through G). Each agent took an unconventional angle on the Erdős minimum-overlap problem and reported back independently.

## Bottom line — UPDATED after deep paper reads

Initial impression: two genuinely new actionable directions (Levers H, I). After deep paper reads of both: **both fail.** Final outcome of the wild iteration: **ten of ten candidate levers eliminated**, three valuable confirmations of the diagnostic, and a refined understanding of WHY the framework is saturated.

The deep reads are committed at [LEVER_H_DEEP_READ.md](LEVER_H_DEEP_READ.md) and [LEVER_I_DEEP_READ.md](LEVER_I_DEEP_READ.md). Their verdicts:

| # | Investigation | Outcome |
|---|---|---|
| 1 | Cross-domain analogy hunt | Initial lead, both deep reads FAILED (see below) |
| 2 | Inequality scavenger hunt | Confirmation: 16 classical inequalities scanned, none beats 0.3801279 |
| 3 | Adversarial red-team audit | Confirmation + minor finding: T5p missing from Phase-5 driver (+1.64 × 10⁻⁶) |
| 4 | SAT-based M(n) extension | New data: extended exact M(n) table from Haugland's n=15 to **n=20** |
| 5 | Lever H deep read (arXiv:2210.16437) | **FAILS**: L² autoconvolution technique doesn't transfer (no polynomial-in-Fourier-coefs expansion of sup) |
| 6 | Lever I deep read (arXiv:2206.09876) | **DOES_NOT_APPLY**: paper is Li 2022 (not CdLS — attribution caught and corrected); technique gives LBs on LP value, wrong direction for ceiling theorem |

## Lever H — RULED OUT after deep read

The initial cross-domain hunt flagged this as a "strong lead." The deep read of arXiv:2210.16437 ([LEVER_H_DEEP_READ.md](LEVER_H_DEEP_READ.md)) overturned that judgment with three independent obstructions:

1. The "second-variation eigenvalue operator" premise was wrong. White's actual technique is a Fourier-truncated quadratically-constrained linear program plus a Minkowski-equality-case strict-convexity uniqueness argument.
2. The L² → L⁴-in-Fourier closed-form (White's Lemma 7) makes the L² problem polynomial of degree 4 in `f̂(k)`. The L∞ min-overlap functional has *no such polynomial expansion* — that's the structural reason our SDP must be a *relaxation* via Bochner-PSD rather than a direct polynomial encoding.
3. The smoothed proxy `F_p(h) = ‖h⋆h‖_p` only converges at rate `(log M)/p`. Closing 10⁻⁴ of µ requires `p ≥ 14,000`, producing `ĥ(k)^{28000}` polynomial terms — vastly past Lasserre tractability.

White himself notes the L∞ regime is "more computationally challenging" in arXiv:2210.16437 §2 and gives no method for it. His 2025 follow-ups (arXiv:2508.02803, arXiv:2506.16750) are constructive UB constructions, not analytical LB transfers.

**Verdict: zero estimated gain.** The technique class is fundamentally different.

## Lever I — RULED OUT, with attribution error caught

The cross-domain hunt cited "Cohn-de Laat-Salmon arXiv:2206.09876". The deep read ([LEVER_I_DEEP_READ.md](LEVER_I_DEEP_READ.md)) caught that **arXiv:2206.09876 is Rupert Li 2022 (sole author, UMN Duluth REU)**, NOT Cohn-de Laat-Salmon. The actual CdLS paper is arXiv:2206.15373 ("Three-point bounds for sphere packing").

Evaluating Li 2022 (the paper actually cited): three blockers prevent applying it as a saturation theorem:
1. **Wrong direction.** Li's Corollary 3.3 says continuous Cohn-Elkies LP ≥ discrete LP — gives LBs on LP value, not the UB / ceiling we'd need.
2. **PSD constraints don't survive Li's restriction-and-periodize map.** Toeplitz PSD constraints at Bochner level n=30 don't commute with sampling.
3. **No analog of "f(x) ≤ 0 for |x| ≥ r"** in White's program. White uses cell-bound LP constraints + Fourier moment constraints, not a free function with pointwise sign constraints. Li's m ≥ 2r hypothesis has no Erdős-side analog.

**Verdict: technique does not apply.** No saturation ceiling theorem available from this source.

## Lever I' — RECOMMENDED alternative (in-session reach)

The Lever I deep read suggested a self-contained alternative: enumerate every augmentation in `white_full_convex.py` and derive Lasserre-tail-bound-style residuals for each. If all residuals are computably small, the cumulative bound IS the saturation theorem.

Template exists at [communications/lasserre_tail_bound.md](communications/lasserre_tail_bound.md). Effort estimate: 1-2 weeks of careful derivation. Output: an explicit `C* < 0.380871` such that no Bochner-PSD augmentation can prove `µ ≥ C*`.

This is the **only direction surfaced this session that is both (a) novel relative to the existing technique stack and (b) tractable from in-session resources**. The original "transfer external machinery" framings (H, I) both failed; the residual-enumeration approach uses the project's own machinery.

## ORIGINAL "strong lead" — Lever H: transfer White's L² autoconvolution machinery (HISTORICAL — RULED OUT, see above)

[OUT_OF_BOX_CROSS_DOMAIN.md](OUT_OF_BOX_CROSS_DOMAIN.md) identified that **White himself has an "almost-tight" L² autoconvolution result** ([arXiv:2210.16437](https://arxiv.org/abs/2210.16437)) using the same Fourier toolkit on a structurally analogous extremal problem. He achieves uniqueness of minimizer to 0.0014% precision — a level of analytical control completely absent from the minimum-overlap state of the art.

If his L² uniqueness / eigenvalue technique transfers (perhaps via a smoothed `L² + ε·L∞` proxy) to the min-overlap setting, this could:
- Pin down the minimizer's structure analytically (rather than empirically via SDP)
- Yield a sharp LB via the same machinery that worked for autoconvolution

**Concrete first action:** Read arXiv:2210.16437 §§3-5 carefully. White has a section on uniqueness via a positive-definiteness argument on the second-variation operator. If the analogous operator for min-overlap admits the same treatment, the framework genuinely transfers.

**Concrete second action:** The project already has [communications/email_to_ethan_white.md](communications/email_to_ethan_white.md). Adapt and send — ask directly whether the L² technique applies. White is the author of both papers; he is the right person to ask. Low effort, potentially decisive answer.

This is **Lever H** — the only lever surfaced this session that hasn't been investigated or ruled out.

## The publishable-negative lead — Lever I: prove the saturation theorem

The same cross-domain hunt also flagged **Cohn-de Laat-Salmon** ([arXiv:2206.09876](https://arxiv.org/abs/2206.09876), Adv. Math. 2024). Their discrete-reduction dual-bound machinery is designed precisely to prove duality gaps for SDP hierarchies on extremal problems.

If their framework can be adapted to White's Bochner-PSD hierarchy, it would convert this session's empirical finding ("the CDE Phase 5 stack is saturated") into a **theorem** that 0.379544 (or some related explicit value) is a proved ceiling of the Bochner-PSD-plus-ellipse technique class.

**Why publishable:** A theorem of the form "no PSD-based augmentation of White's program can prove µ ≥ C* for C* > 0.380X" is independently valuable. It tells the community that further work must be qualitatively different, and it pins down the exact ceiling.

This is **Lever I** — a converse to the saturation diagnostic.

## What we confirmed (no change to the bound)

### Lever 2 — inequalities don't help

[OUT_OF_BOX_INEQUALITIES.md](OUT_OF_BOX_INEQUALITIES.md) confirmed: 16 classical inequalities (Plancherel, Hausdorff-Young, Bombieri, Selberg, Beckner, Wiener, Ingham/Wirtinger, Heisenberg, Beurling-Selberg/Logan, Sárközy, Plünnecke-Ruzsa, Brunn-Minkowski, Roth, Rudin Λ(p), Boas, Cauchy-Schwarz) scanned. The strongest LB from a single classical inequality is **1/(2 + √2) ≈ 0.2929** (White §3 partition-of-unity). All other natural Fourier-analytic one-liners collapse to the trivial 1/8 = 0.125.

The structural reason this matters: White's framework is the *terminal point* of the LP-duality side of classical inequality reasoning. Pushing further isn't about finding a sharper inequality — it's about finding a different mathematical framework.

### Lever 3 — red-team audit holds

[OUT_OF_BOX_REDTEAM.md](OUT_OF_BOX_REDTEAM.md): all three adversarial framings failed to break the diagnostic.

- **LB side**: One minor gap (`use_T5p` flag silently omitted from Phase-5 driver `path_b_with_polymoment.py`). T5p adds +1.88 × 10⁻⁴ to bochner_n=20 alone but only **+1.64 × 10⁻⁶** in the full Phase-5 composition. Below the headline-bound threshold but worth closing for completeness.
- **UB side**: Together's h\* is genuinely locally optimal. 2000 random perturbations across 4 σ-scales found zero improvement. SLP linearizations *predict* improvement but true M *regresses* on each iteration. Their plateau is real.
- **Formulation side**: No gap. White §1 explicitly defines µ as the continuous Ω-functional; equivalence to lim M(n)/n via Swinnerton-Dyer 1996, in Haugland's J. Number Theory paper.

### Lever 4 — SAT extension

[OUT_OF_BOX_SAT_MN.md](OUT_OF_BOX_SAT_MN.md): SAT-based exact M(n) computation reached **n = 20** in 90 seconds (n = 21 exceeded the 10-min cutoff).

**New exact values not in Haugland's published table:**
- M(19) = 8 (ratio 0.4211)
- M(20) = 8 (ratio 0.4000, ties n=15 as smallest known integer ratio)

This is a small but genuine novelty — the published table stops at n=15. Worth flagging in the writeup as supporting data, even though it doesn't improve the UB on µ. The exponential scaling wall is real: n=53 (the threshold where ⌈0.380871n⌉/n drops below 0.4) is astronomical via SAT.

## FINAL lever ledger (10 levers eliminated + 1 new direction recommended)

| Lever | Status | Why / Where |
|---|---|---|
| A — Lukács SOS / alt basis | ❌ unlikely | Gibbs already damped post-CDE |
| B — Diagnostic | ✓ executed | [TOGETHER_DIAGNOSTIC.md](TOGETHER_DIAGNOSTIC.md) |
| C — Integer M(n) | ❌ Together UB stands | Smallest known M(n)/n = 0.40; even at SAT-reachable n=20, > 0.380871 |
| D — O(1) breakpoints | ❌ refuted | h\* has 400+ blocks |
| D' — Lipschitz/BV via discrete | ❌ refuted | Lifted discrete optima diverge from h\* |
| E — M-side SDP | ❌ vacuous | Relaxations empirically inactive |
| F — More step-function steps | ❌ saturated | 95→600 gained only 5e-5 |
| G — (f, g) rewrite | ❌ analytically no-op | Convex-hull equivalent |
| H — White's L² autoconvolution transfer | ❌ FAILS | Deep read found no polynomial expansion of sup-functional; L^p proxy needs p≥14000 |
| I — Li/CdLS saturation theorem | ❌ DOES NOT APPLY | Wrong direction; PSD constraints don't survive Li's map |
| **I' — Residual-enumeration ceiling theorem** | **🟡 NEW, recommended** | Use [lasserre_tail_bound.md](communications/lasserre_tail_bound.md) template on every augmentation; 1-2 weeks effort |

## Recommended next steps, prioritized — UPDATED

### Action 1 (DONE during overnight): Close the T5p gap

Threaded `use_T5p` through `path_b_analytical.build_problem_with_dual_handles` and `path_b_with_polymoment.py`. Sanity-check at N=3000 confirms +1.6 × 10⁻⁶ gain. Full Phase 5 + T5p run at N=10000 currently running; will commit results when complete.

### Action 2 (DONE during overnight): Email v2 draft to Ethan White

Reframed [communications/email_to_ethan_white_v2.md](communications/email_to_ethan_white_v2.md) (superseding the v1 draft) around:
- The strengthened headline µ ≥ 0.3801279
- The diagnostic finding that the cell-kernel envelope is binding, not the augmentations
- The "ceiling theorem" question — does White expect a saturation theorem in his framework?
- The closed-form sup-functional question — is there an identity for `sup_t (h⋆h)(t)` we're missing?

Ready for user review before sending. NOT autonomous-sent.

### Action 3 (DONE during overnight): Preprint addendum

Written at [communications/preprint_addendum.md](communications/preprint_addendum.md). Documents the post-v1 extensions (CDE Phase 1-5, Together-as-primal diagnostic, SAT M(n) extension, saturation lever ledger) with three publication-path recommendations:
- v2 preprint (strongest, requires rewrite of Theorem 1 / Section 4)
- Two separate papers (clean editorial separation between bound improvement and saturation methodology)
- v1 unchanged + addendum as appendix (lowest risk)

User-decision pending on which path to take.

### Action 4 (NEW recommended direction): Residual-enumeration ceiling theorem (Lever I')

The deep read of Lever I surfaced an in-session-reach approach to the saturation theorem question. Use the template at [communications/lasserre_tail_bound.md](communications/lasserre_tail_bound.md) (which already worked out the Fejér-Riesz tail bound that retracted Lasserre) and apply systematically to every constraint family in `white_full_convex.py`. The cumulative residual bound, if computable, IS the saturation theorem.

Effort: 1-2 weeks of careful derivation. Output: an explicit `C*` such that no Bochner-PSD augmentation can prove `µ ≥ C*`. Even partial results are publishable.

This is the **single recommended direction** for any future session that wants to push this problem further. Previous "transfer external machinery" approaches (H, I) both failed; the residual-enumeration uses only the project's own toolkit.

## What this session's "wild" iteration achieved

- Surfaced 2 new candidate levers (H, I) not previously considered, both with concrete arXiv leads and specific technical hooks
- Cross-verified the diagnostic's conclusion via 3 independent adversarial probes (all failed to break it)
- Extended Haugland's published M(n) table by 5 values via SAT
- Identified a trivial formal-completeness fix (T5p in Phase-5 composition)

This is exactly the kind of "out-of-box" search you can't substitute for during a normal diagnostic chain — it requires explicitly seeking weakness in your own framework and explicitly looking outside your immediate technique tree.

## Files committed this iteration

- [OUT_OF_BOX_CROSS_DOMAIN.md](OUT_OF_BOX_CROSS_DOMAIN.md) — cross-domain analogy hunt
- [OUT_OF_BOX_INEQUALITIES.md](OUT_OF_BOX_INEQUALITIES.md) — inequality scavenger hunt
- [OUT_OF_BOX_REDTEAM.md](OUT_OF_BOX_REDTEAM.md) — adversarial red-team audit
- [OUT_OF_BOX_SAT_MN.md](OUT_OF_BOX_SAT_MN.md) — SAT-based M(n) extension
- [lp_research_state/code/_sat_Mn.py](lp_research_state/code/_sat_Mn.py) — SAT solver code
- [lp_research_state/data/Mn_sat_results.json](lp_research_state/data/Mn_sat_results.json) — n=15..20 data
- [OUT_OF_BOX_SYNTHESIS.md](OUT_OF_BOX_SYNTHESIS.md) — this synthesis
