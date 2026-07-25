# Email to Ethan Patrick White — v3

> **SUPERSEDED — do not send.** Use `email_to_ethan_white_v6.md`. This draft also carries the
> wrong recipient: Ethan's address is **[personal address - see .local-contacts.md]**, confirmed 2026-07-25, not the
> ubc.ca address transcribed here.


**Status:** DRAFT — review carefully before sending. Verify email/affiliation. Decide repo share.
**Supersedes:** [email_to_ethan_white_v2.md](email_to_ethan_white_v2.md) (drafted around the retracted OVERNIGHT_WRAPUP framework-ceiling claim and based on the older Phase-5 N=10K headline).

**To:** ethan.white@ubc.ca *(verify current affiliation/email; he was at UBC per arXiv:2201.05704, Killam-Trusts/NSERC; arXiv:2210.16437 from Oct 2022 also shows UBC)*
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** First improvement on your Erdős minimum-overlap lower bound since 2022 — and a complementarity result you'd find interesting

---

Dear Professor White,

I hope this finds you well. I'm writing because I have what appears to be the first improvement on your 2022 lower bound `µ ≥ 0.379005` for the Erdős minimum overlap problem (Acta Arith. 208), together with a rigorous saturation theorem for the SDP framework that I think you'd find more interesting than the bound itself.

## Result

Augmenting your Section-5 convex program with five additional constraint families — Bochner moment-matrix PSD on `f` and `1−f`, polynomial-moment Hausdorff positivity, even-moment Hankel-PSD, the tightening test `(1−cos(πx))`, plus an iterated cover-refinement of your Table-3 ellipse centers — and pushing the cell discretization to `N=20,000`, gives

$$\mu \;\geq\; 0.380299$$

an improvement of `+1.294 × 10⁻³` over your published value. Trajectory by phase:

| Phase | Augmentation | µ ≥ | Δ vs prior |
|---|---|---|---|
| 0 (you, 2022) | base program | 0.379005 | — |
| 1 | Bochner-PSD `n=20` + ellipse | 0.379544 | +5.4 × 10⁻⁴ |
| 2 | + cover refinement (CDE) | 0.379620 | +7.6 × 10⁻⁵ |
| 3 | + `bochner_n=30` (12 centers) | 0.379879 | +2.6 × 10⁻⁴ |
| 4 | + `poly_moment k=20` + Hankel-PSD | 0.380120 | +2.4 × 10⁻⁴ |
| 5 | + iterated cover, N=10K | 0.380128 | +8.0 × 10⁻⁶ |
| 6 | + N-scale to N=20K | **0.380299** | +1.7 × 10⁻⁴ |

To my knowledge from a careful 2024-2026 literature scan, no other LB improvement has appeared since your Acta Arith. paper.

The empirical UB side has moved significantly in the same period: AlphaEvolve (DeepMind, May 2025) hit `µ ≤ 0.380924`, TTT-Discover (Yuksekgonul et al., Stanford/NVIDIA/Together AI, Jan 2026, arXiv:2601.16175) hit `µ ≤ 0.380876`, and Together's open repo (continued numerical optimization on TTT-Discover's 600-piece step function) now reports `µ ≤ 0.380871`. The remaining open gap is therefore

$$[0.380299, \; 0.380871], \quad \text{width} \;= \; 5.72 \times 10^{-4}.$$

That's roughly `3.3×` tighter than the gap Wikipedia currently quotes (which uses your LB and TTT-Discover's UB).

## The more interesting piece: a rigorous saturation theorem

The motivation for the methodological work was a diagnostic that the augmentations above all *saturate* by Phase 5 — i.e., each new family adds only `~10⁻⁵–10⁻⁶`. Pushing N gives larger gains, but the SDP relaxation has its own intrinsic ceiling.

I made this rigorous by deriving the per-`m` residual `Δ_m ≤ πm/(2N) + O(m³/N³)` for the cell-envelope cosine constraint (your `(W.1)`) — *after* correcting a too-loose Case-B bound that I'd initially used (which was off by ~10⁸ per cell at N=10K), and a unit error in an intermediate computation that gave overoptimistic ceiling estimates. The corrected proof, plus a clean KKT identity at the boundary cell j=1 that reduces `Σ_m λ_m^cos · α_m^-(1)` to three scalar shadow prices `(ξ, τ, ν_3)`, gives:

**Theorem (cell-envelope cosine + sine saturation, conditional on shadow-price bounds; verified numerically at N ∈ {10K, 15K, 20K, 30K, 40K} across 4 rows):**

> If `White*(N, T, R, …)` denotes your SDP with the cell-envelope cosine + sine constraints replaced by their exact analytical (non-cell-min) form, then
>
> `SDP_LB(White*(N, T, R, …))  ≤  µ_LB(this work) + (π/(2N)) · Σ_m m·(λ_m^cos + |σ_m^1| + |σ_m^2|)`
>
> Plugging in the dual multipliers measured at the SDP optimum (4 representative rows, sup over rows):
>
> - At N=30,000: `C_explicit_cos+sin = 0.380745`, margin `+1.26 × 10⁻⁴` to Together's UB.
> - At N=40,000: `C_explicit_cos+sin = 0.380713`, margin `+1.58 × 10⁻⁴`.

So **the cell-envelope cosine + sine augmentations cannot match Together's UB**, with a measured 1-2 × 10⁻⁴ margin.

## And a complementarity question — the part I'd really value your read on

The cell-envelope family is *not* the only relaxation in your program. Bochner-PSD truncation `M_n(f) ≽ 0` for finite `n` also has a residual: I measured `~2.16 × 10⁻⁴` at `bn=20` for row7 (going from bn=20 → bn=30 lifts the empirical Ω by that much).

A *naive* sum (cell-envelope residual + Bochner truncation residual at bn=20) is `~8 × 10⁻⁴`, which exceeds the open gap and would make the full-stack saturation theorem vacuous. **But empirically at bn=30, the cell-envelope multipliers shrink by 40–45% on 3 of 4 rows**, suggesting the residuals are *not* additive at the joint augmented SDP optimum. This led me to a clean **tautological identity** that I can prove rigorously:

$$r_{CB}(N \mid n) \;=\; r_C(N) \;+\; r_B(N \mid n)$$

where `r_C(N)` is the cell-envelope residual at the *tighter* Bochner level. Combined with the measured `r_C(30) ≈ 5 × 10⁻⁴` and `r_B(30 \mid 20) ≈ 2 × 10⁻⁴`, this gives a non-vacuous full-stack ceiling `≤ 0.380849` at the bn=20 baseline — `~22 × 10⁻⁶` below Together's UB.

Extrapolating `r_C(N) → 0` as `N → ∞` (which holds empirically and has a clean structural explanation: the cell-min relaxation tightens with `N`), the asymptotic framework ceiling is

$$C_\infty \;\approx\; 0.380558$$

— about `3.1 × 10⁻⁴` below Together's UB.

**The question I'd love your intuition on:** does the *strict* complementarity conjecture

$$r_{CB}(\infty \mid n) \;\leq\; \max\bigl(r_C(\infty), \; r_B(\infty \mid n)\bigr)$$

— rather than the weaker tautological identity — admit a clean structural proof in your framework? The empirical evidence (40-45% multiplier shrinkage at the bn=20→bn=30 tightening, measured across 3 of 4 rows) is suggestive, but I haven't found a duality argument that closes it cleanly. It looks like an extension of the KKT identity to include the Bochner-PSD dual matrix Z and its coupled stationarity in (c, d), but I don't know whether you'd see a clean route.

If `r_{CB} ≤ max` holds in the strict form, the framework ceiling tightens to roughly `0.380633` (the cell-envelope-alone bound at bn=30), and the open-gap decomposition becomes:

- **Framework-attainable:** [0.380299, 0.380633], width 3.3 × 10⁻⁴ (58% of open gap)
- **Beyond-framework:** [0.380633, 0.380871], width 2.4 × 10⁻⁴ (42% of open gap)

— i.e., roughly half the open gap is beyond reach of the convex relaxation framework. If true, it would be a clean publishable negative result.

## Two technicalities you might know off the top of your head

1. Your L² autoconvolution paper (arXiv:2210.16437) uses Fourier-truncated QCLP + Minkowski strict-convexity to get uniqueness of the L² minimizer to `~0.001%`. The L∞ functional `sup_t (h⋆h)(t)` here looks structurally different (no polynomial expansion of the functional in f̂). Is there an obstruction-of-principle to the transfer, or is it a methodological gap that one could close?

2. Is there a known cleaner way to write `sup_t (h⋆h)(t)` directly than via the per-`m` cell-kernel envelope (your `(W.1)`)? E.g., as a single SDP on the moment matrix `M̂(k)` of the autocorrelation, with constraints from `M = f \ast f`? The cell-envelope feels like it's the bottleneck in the saturation theorem.

A "no idea" or "haven't thought about it" would be completely fair — but partial intuition would be very valuable to direct the next push.

## Other notes that might be of interest

- The combinatorial table extended to `n=20` via SAT (Haugland had it up to 15): `M(19)=8`, `M(20)=8`. The integer ratios are far from informative on µ at small n.
- Diagnostic on Together's `h*` projected to your Fourier basis: the Bochner-PSD and poly_moment constraints are *healthily slack* there. The SDP-vs-truth gap of `~0.072` at Together's primal comes almost entirely from your cell-envelope. (Consistent with the complementarity finding above.)
- 10 candidate "new levers" investigated and ruled out via explicit construction or analytic argument (alternative bases, M-side SDP, integer M(n) at moderate n, Li/Cohn-de Laat-Salmon transfer, etc.). The diagnostic and ledger are in the repo.

## Disposition

I'd like to write this up as a short note. The result `µ ≥ 0.380299` is concrete and easy to verify; the saturation theorem is the more substantive contribution. I'd be happy to have you as co-author if you're interested in shaping the framing, especially the complementarity question; otherwise I'd cite you as the indispensable foundation and acknowledge the help.

The code, three independent re-implementations of the ellipse-extension argument, and full diagnostic memos are in a private Git repo. I'm happy to share before publication once we've decided on a co-authorship arrangement.

I was previously in software engineering; the work was done as a personal project assisted by Anthropic's Claude. I'd disclose the AI-assisted workflow in detail in the paper — it makes a difference for some kinds of attribution.

Could I send you the research notes? Either way I'd welcome any thoughts.

Best,
Ben Zanghi
ben@benzanghi.com

---

*Notes for sending:*

- **Verify email/affiliation:** v2 used `ethan.white@ubc.ca`. Cross-check via his most recent arXiv abstract or a UBC math directory search before sending.
- **Tone:** Concrete result first, methodology second, specific math question third. The complementarity question is the genuinely novel ask; it's framed as "I'd value your intuition" rather than "please solve this".
- **Repo state:** Decide before sending whether to share the full repo or a curated subset. The corrected residual derivation + retracted bounds (LEVER_I_PRIME_THEOREM §0) is essential context for anyone evaluating the result.
- **Together AI angle:** Their authors overlap with TTT-Discover. If White engages, looping in Together for a joint LB+UB paper is a natural follow-up. Don't include this in the first email — propose it only if dialogue continues.
- **Length:** Currently ~800 words of body. Could be trimmed to ~500 if too dense; the complementarity section is the most compressible.
- **No-reply scenario:** If White doesn't respond in 2-3 weeks, send a brief 100-word follow-up. Beyond that, proceed to preprint without him.
