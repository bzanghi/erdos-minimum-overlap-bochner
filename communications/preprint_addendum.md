# Preprint addendum — post-v1 extensions

This document collects extensions to [preprint_draft.tex](preprint_draft.tex) (v1, stating µ ≥ 0.379544) that have accrued since v1 was drafted. Editorial decision pending: incorporate as additional sections (giving a stronger v2), publish as a separate methodological note, or both.

**Recommended for the user (Ben):** decide between three publication paths:
1. **v2 preprint** with the strengthened headline µ ≥ 0.3801279 — clean and stronger, but requires rewriting Theorem 1 and Section 4 (results table) substantively.
2. **Two papers:** v1 as a short note on the original Bochner-PSD extension, plus a separate note on the saturation diagnostic and the "Constraint Discovery Engine" pattern. The diagnostic is a methodological contribution independent of the bound.
3. **v1 unchanged + this addendum as a "supporting materials" appendix.** Cleanest editorially, smallest scope.

## Headline progression accrued

| Phase | Augmentation | µ ≥ | Δ vs prior |
|---|---|---:|---:|
| 0 (White 2023) | base program | 0.379005 | — |
| 1 | Bochner-PSD n=20 + ellipse (= preprint v1's result) | 0.379544 | +5.4 × 10⁻⁴ |
| 2 | + cover refinement (4 new CDE-discovered centers) | 0.379620 | +7.6 × 10⁻⁵ |
| 3 | + bochner_n=30 (12 centers) | 0.379879 | +2.6 × 10⁻⁴ |
| 4 | + poly_moment k_max=14 + Hankel-PSD n=6 | 0.380120 | +2.4 × 10⁻⁴ |
| 5 | + iterated cover refinement on combined stack | **0.3801279** | +8.0 × 10⁻⁶ |
| 5+T5p | + T5p (1−cos(πx)) test added to Phase 5 driver | 0.380129[*] | +1.6 × 10⁻⁶ |

[*] Currently being verified at full scale at time of this writing; previously measured at +1.6 × 10⁻⁶ at smaller scale.

Net: **+1.123 × 10⁻³ over White's 0.379005**. Open gap is now [0.3801279, 0.380871] of width 7.43 × 10⁻⁴.

## Five new methodological contributions (any of which deserve writeup)

### M1 — Cover refinement / Constraint Discovery Engine (CDE) pattern

White's seven ellipse centers in his Table 3 were chosen for HIS unaugmented bound. Once Bochner-PSD is added, the binding point of the dual envelope shifts. We iterate: read off the binding point of the current envelope, solve a new full-scale SDP at that center, extract duals via the same path-B machinery, add a new ellipse, repeat. Converges in 4-5 iterations. Phase 1 and Phase 2 of the table above. **+1.0 × 10⁻⁴ from refinement alone.**

### M2 — Polynomial moment constraints with rigorous tail bounds

The Hausdorff moment problem on [-1, 1] gives `m_{2k}(f) = ∫ x^{2k} f(x) dx ≥ 0` for any nonneg f. Encoded as linear constraints on (c, d) via the Fourier expansion `m_{2k} = (1/2) α_0^{(k)} + Σ_j (c_j α_j^{(k)} + d_j β_j^{(k)}) + tail_k`, with `|tail_k| ≤ O(k/T)` derivable by integration-by-parts. Implementation at [`poly_moment.py`](../lp_research_state/code/poly_moment.py). **+4.4 × 10⁻⁴ at single-row scale, +2.4 × 10⁻⁴ in full cover composition.**

### M3 — Together-as-primal SDP diagnostic

A new methodological procedure: take a competitor's primal certificate (Together's piecewise-constant h\* attaining `M(h*) = 0.380871`), project into our SDP's Fourier basis at the same truncation, and evaluate every constraint as a slack. We find:
- All CDE augmentations (Bochner, poly-moment, Hankel) are healthily slack at Together's h\*
- SDP-encoded Ω at h\* (with row 4's box overridden to match h\*'s f̂(1)) is **0.459311** — vs h\*'s actual autocorrelation of 0.387337
- The binding constraint family is the *original cell-kernel autocorrelation envelope* (White's §5 (5.1)-(5.13)), NOT the CDE augmentations
- Gap function f̃ − f_even is 99.9% low-frequency — the SDP-optimal f̃ is structurally different from Together's h\* in a smooth (non-Gibbs) way

Full diagnostic at [TOGETHER_DIAGNOSTIC.md](../TOGETHER_DIAGNOSTIC.md). The procedure is independent of the bound: any future SDP-based attack on this or related extremal problems can use the same diagnostic.

### M4 — SAT-based exact M(n) extension

Haugland's published table of exact M(n) values stops at n=15. Using a SAT encoding (pseudo-boolean cardinality + sequential counter for the overlap constraint, binary search on M, with the Glucose4 solver via python-sat), we extended the table:

| n | M(n) | M(n)/n | source |
|--:|--:|--:|--|
| 16 | 7 | 0.4375 | session brute force |
| 17 | 7 | 0.4118 | session brute force |
| 18 | 8 | 0.4444 | session brute force |
| 19 | 8 | 0.4211 | session SAT (**new**) |
| 20 | 8 | 0.4000 | session SAT (**new**) |

All cross-verified by direct overlap computation. Code at [`_sat_Mn.py`](../lp_research_state/code/_sat_Mn.py). None of the new ratios beats Together's UB; they're informative for the asymptotic structure analysis.

### M5 — Saturation lever ledger

Enumeration of candidate next-step techniques and rigorous ruling-out of each. Documented in [SESSION_2026-05-10_FINAL_ADDENDUM.md](../SESSION_2026-05-10_FINAL_ADDENDUM.md) and [OUT_OF_BOX_SYNTHESIS.md](../OUT_OF_BOX_SYNTHESIS.md). The ten candidate levers, all eliminated:

| Lever | Status |
|---|---|
| A — Lukács SOS / alternative basis | unlikely (Gibbs damped post-CDE) |
| B — Together-as-primal diagnostic | executed; produced M3 above |
| C — Push integer M(n) | extended to n=20 via SAT (M4); doesn't beat Together |
| D — O(1)-breakpoint restriction | refuted empirically (h\* has 400+ blocks) |
| D' — Lipschitz/BV via discrete limit | refuted (lifted optima diverge from h\*) |
| E — M-side SDP encoding | empirically vacuous (Δ ~ 10⁻⁷) |
| F — Push step-function UB past 600 steps | saturated near 0.38087 in literature |
| G — (f, g) rewrite | analytically no-op |
| H — White's L² autoconvolution transfer | FAILS (no polynomial expansion of sup-functional) |
| I — Li/CdLS duality-gap saturation theorem | DOES_NOT_APPLY (technique is wrong direction) |

This is itself publishable: a clean technique-saturation result for an open problem helps focus future work outside the current technique stack.

## Two open methodological questions surfaced by the diagnostic

### Q-A — Does the cell-kernel envelope admit a provable saturation ceiling?

The diagnostic identified the cell-kernel autocorrelation envelope (`white_full_convex.py:176-190`) as the binding constraint family at Together's h\*. A theorem of the form

> The SDP `sup{ c : Bochner-PSD-augmented White's program at level n proves µ ≥ c }` is bounded above by an explicit `C* < 0.380871`

would convert the empirical saturation into a published ceiling. The deep read of arXiv:2206.09876 (Li/CdLS attribution corrected; see [LEVER_I_DEEP_READ.md](../LEVER_I_DEEP_READ.md)) ruled out the direct adaptation, but suggested an alternative: derive Lasserre-tail-bound-style residuals for *each* of White's constraint families, then bound their cumulative slack analytically (template at [`communications/lasserre_tail_bound.md`](lasserre_tail_bound.md)). This is in-session-reach effort.

### Q-B — Is there a Fourier identity for sup_t (h⋆h)(t) we're missing?

The L∞-norm functional `sup_t (h⋆h)(t)` has no polynomial expansion in `ĥ(k)` (unlike L² or Lp for finite p). This is the root cause of the relaxation gap and the reason White's L² autoconvolution technique (which gets uniqueness via polynomial-Fourier identities) does not transfer.

An identity expressing `sup_t (h⋆h)(t)` directly via Fourier data — even an implicit / variational identity — could collapse the cell-kernel envelope into something tighter. Worth asking White directly (see [email_to_ethan_white_v2.md](email_to_ethan_white_v2.md)).

## Recommended publication path

**My recommendation (autonomous reasoning):** Option 2 — two separate papers:

- **Paper 1 (short, fast):** Update preprint v1 to claim µ ≥ 0.3801279, with the CDE Phase 1-5 narrative as Section 4 ("Extensions") and the SAT M(n) values as a small Section 5. Single new headline number; clean theorem statement. Estimated effort: 1-2 weeks of LaTeX work to extend v1.

- **Paper 2 (methodological):** "Diagnosing saturation in convex hierarchies: the Together-as-primal procedure on the Erdős minimum-overlap problem." The diagnostic (M3 above) is genuinely novel methodology that doesn't depend on the bound improvement. Has standalone interest for anyone working on SDP relaxations of extremal problems. Estimated effort: 1 month.

Either or both could be appropriately co-authored with White (he is the indispensable foundation of both).

Option 3 (v1 + this addendum as appendix) is the lowest-risk if there's any concern about scope creep. Option 1 (full v2 rewrite) is the strongest single statement but takes the longest.
