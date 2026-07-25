# Email to Ethan Patrick White — v4

> **SUPERSEDED — do not send.** Use `email_to_ethan_white_v6.md`. This draft also carries the
> wrong recipient: Ethan's address is **[personal address - see .local-contacts.md]**, confirmed 2026-07-25, not the
> ubc.ca address transcribed here.


**Status:** DRAFT — review carefully before sending. Verify email/affiliation. Decide repo share.
**Supersedes:** [email_to_ethan_white_v3.md](email_to_ethan_white_v3.md) (used the pre-correction headline 0.380299 and called the asymptotic ceiling a theorem).
**Changes in v4:** headline corrected to `µ ≥ 0.3802973` (poly-moment tail-bound fix, 2026-05-22); production config now `bochner_n=40`; open gap `5.74 × 10⁻⁴`; asymptotic ceiling `C_∞` framed as a conjecture, not a theorem; added one line on the planned arbitrary-precision (SDPA-GMP) certification.

**To:** ethan.white@ubc.ca *(verify current affiliation/email; he was at UBC per arXiv:2201.05704, Killam-Trusts/NSERC; arXiv:2210.16437 from Oct 2022 also shows UBC)*
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** First improvement on your Erdős minimum-overlap lower bound since 2022 — and a complementarity result you'd find interesting

---

Dear Professor White,

I hope this finds you well. I'm writing because I have what appears to be the first improvement on your 2022 lower bound `µ ≥ 0.379005` for the Erdős minimum overlap problem (Acta Arith. 208), together with a framework-saturation result that I think you'd find more interesting than the bound itself.

## Result

Augmenting your Section-5 convex program with five additional constraint families — Bochner moment-matrix PSD on `f` and `1−f`, polynomial-moment Hausdorff positivity, even-moment Hankel-PSD, the tightening test against `(1−cos(πx))`, plus an iterated cover-refinement of your Table-3 ellipse centers — and pushing the cell discretization to `N=20,000` with Bochner truncation `n=40`, gives

$$\mu \;\geq\; 0.3802973$$

an improvement of `+1.292 × 10⁻³` over your published value. Trajectory by phase:

| Phase | Augmentation | µ ≥ | Δ vs prior |
|---|---|---|---|
| 0 (you, 2022) | base program | 0.379005 | — |
| 1 | Bochner-PSD `n=20` + ellipse | 0.379544 | +5.4 × 10⁻⁴ |
| 2 | + cover refinement (CDE) | 0.379620 | +7.6 × 10⁻⁵ |
| 3 | + `bochner_n=30` (12 centers) | 0.379879 | +2.6 × 10⁻⁴ |
| 4 | + `poly_moment k=20` + Hankel-PSD | 0.380120 | +2.4 × 10⁻⁴ |
| 5 | + iterated cover, N=10K | 0.380128 | +8.0 × 10⁻⁶ |
| 6 | + N=20K, `bochner_n=40` | **0.3802973** | +1.7 × 10⁻⁴ |

(One honest note on the numbers: an adversarial self-review caught that one of the five families — the polynomial-moment cuts — was using a Fourier-truncation tail bound that omitted the analytic remainder of an infinite sum, making it slightly too tight. The corrected, provably-rigorous bound lowered the headline by `5.4 × 10⁻⁶`, from a previously-computed `0.3803027` to the `0.3802973` above. Same failure mode that made me retract an earlier Lasserre-level-2 attempt; I mention it because you'll appreciate how easy this trap is.)

To my knowledge from a careful 2024–2026 literature scan, no other LB improvement has appeared since your Acta Arith. paper.

The empirical UB side has moved significantly in the same period: AlphaEvolve (DeepMind, May 2025) hit `µ ≤ 0.380924`, TTT-Discover (Yuksekgonul et al., Stanford/NVIDIA/Together AI, Jan 2026, arXiv:2601.16175) hit `µ ≤ 0.380876`, and Together's open repo (continued numerical optimization on TTT-Discover's 600-piece step function) now reports `µ ≤ 0.380871`. The remaining open gap is therefore

$$[0.3802973, \; 0.380871], \quad \text{width} \;= \; 5.74 \times 10^{-4}.$$

That's roughly `3.3×` tighter than the gap Wikipedia currently quotes (which uses your LB and TTT-Discover's UB).

## The more interesting piece: how far the framework can go

The motivation for the methodological work was a diagnostic that the augmentations above all *saturate* by Phase 5 — i.e., each new family adds only `~10⁻⁵–10⁻⁶`. Pushing N gives larger gains, but the SDP relaxation has its own intrinsic ceiling.

I made the fixed-baseline part rigorous by deriving the per-`m` residual `Δ_m ≤ πm/(2N) + O(m³/N³)` for the cell-envelope cosine constraint (your `(W.1)`) — *after* correcting a too-loose Case-B bound that I'd initially used (off by ~10⁸ per cell at N=10K) and a unit error in an intermediate computation. The corrected proof, plus a clean KKT identity at the boundary cell j=1 that reduces `Σ_m λ_m^cos · α_m^-(1)` to three scalar shadow prices `(ξ, τ, ν_3)`, gives:

**Theorem (cell-envelope cosine + sine ceiling; verified numerically at N ∈ {10K, 15K, 20K, 30K, 40K} across 4 rows):**

> If `White*(N, T, R, …)` denotes your SDP with the cell-envelope cosine + sine constraints replaced by their exact analytical (non-cell-min) form, then
>
> `SDP_LB(White*(N, T, R, …))  ≤  SDP_LB(White(...))  +  (π/(2N)) · Σ_m m·(λ_m^cos + |σ_m^1| + |σ_m^2|)`
>
> Plugging in the dual multipliers measured at the SDP optimum (sup over 4 representative rows):
>
> - At N=30,000: cell-envelope ceiling `≈ 0.380745`, margin `+1.26 × 10⁻⁴` to Together's UB.
> - At N=40,000: `≈ 0.380713`, margin `+1.58 × 10⁻⁴`.

So **the cell-envelope cosine + sine augmentations cannot match Together's UB**, with a measured 1–2 × 10⁻⁴ margin.

Combining the cell-envelope and Bochner-truncation residuals via an exact set identity `r_CB(N | n) = r_C(N) + r_B(N | n)` (the joint augmented region is just `K_N ∩ F_C`) gives a non-vacuous full-stack ceiling `≤ 0.380849` at the bn=20 baseline — `~2.2 × 10⁻⁵` below Together's UB. Extrapolating `r_C(N) → 0` as `N → ∞` (empirically a `1/N` law, though a Richardson fit gives exponent `0.94 ± 0.02`, hinting at log-corrections) suggests an asymptotic ceiling

$$C_\infty \;\approx\; 0.380558$$

— about `3.1 × 10⁻⁴` below Together's UB. I'm careful to call this last step a **conjecture**, since it rests on the empirical `N → ∞` extrapolation rather than a proof.

## And a complementarity question — the part I'd really value your read on

**The question:** does the *strict* complementarity conjecture

$$r_{CB}(\infty \mid n) \;\leq\; \max\bigl(r_C(\infty), \; r_B(\infty \mid n)\bigr)$$

— rather than the weaker exact identity — admit a clean structural proof in your framework? The empirical evidence (40–45% cell-envelope multiplier shrinkage at the bn=20→bn=30 tightening, measured across 3 of 4 rows) is suggestive, but I haven't found a duality argument that closes it. It looks like an extension of the KKT identity to include the Bochner-PSD dual matrix Z and its coupled stationarity in (c, d), but I don't know whether you'd see a clean route.

If `r_CB ≤ max` holds in the strict form, the framework ceiling tightens to roughly `0.380633`, and roughly half the open gap is provably beyond reach of the convex relaxation framework — a clean publishable negative result.

## Two technicalities you might know off the top of your head

1. Your L² autoconvolution paper (arXiv:2210.16437) uses Fourier-truncated QCLP + Minkowski strict-convexity to get uniqueness of the L² minimizer. The L∞ functional `sup_t (h⋆h)(t)` here looks structurally different (no polynomial expansion of the functional in f̂). Is there an obstruction-of-principle to the transfer, or is it a methodological gap one could close?

2. Is there a known cleaner way to write `sup_t (h⋆h)(t)` directly than via the per-`m` cell-kernel envelope (your `(W.1)`)? I tried a direct sup-t SDP and it produced an *invalid* lower bound, which suggests the cell-envelope is load-bearing for validity, not just convenience.

A "no idea" or "haven't thought about it" would be completely fair — but partial intuition would be very valuable to direct the next push.

## Other notes that might be of interest

- The combinatorial table extended to `n=20` via SAT (you had it up to 15): `M(19)=8`, `M(20)=8`. The integer ratios are far from informative on µ at small n.
- Diagnostic on Together's `h*` projected to your Fourier basis: the Bochner-PSD and poly_moment constraints are *healthily slack* there. The SDP-vs-truth gap of `~0.072` at Together's primal comes almost entirely from your cell-envelope (consistent with the complementarity finding above).
- ~10 candidate "new levers" investigated and ruled out via explicit construction or analytic argument (alternative bases, M-side SDP, integer M(n) at moderate n, transfer arguments, etc.). The ledger is in the repo.

## On rigor / next step

The bound is currently certified via CLARABEL (a primal-dual interior-point solver) with rigorous dual-objective extraction; the actual duality gaps are `10⁻⁶–10⁻⁸`, far inside my `2 × 10⁻⁶` margin. For a fully solver-independent certificate I'm planning to re-solve the binding row(s) at arbitrary precision (SDPA-GMP, already built and smoke-tested to `~10⁻⁷⁵`); the only missing piece is a cvxpy→SDPA-S serializer. I'd flag this honestly in any writeup.

## Disposition

I'd like to write this up as a short note. The result `µ ≥ 0.3802973` is concrete and easy to verify; the saturation result is the more substantive contribution. I'd be happy to have you as co-author if you're interested in shaping the framing — especially the complementarity question — otherwise I'd cite you as the indispensable foundation and acknowledge the help.

The code, three independent re-implementations of the ellipse-extension argument, and full diagnostic memos are in a private Git repo. I'm happy to share before publication once we've decided on a co-authorship arrangement.

I was previously in software engineering; the work was done as a personal project assisted by Anthropic's Claude. I'd disclose the AI-assisted workflow in detail in the paper — it makes a difference for some kinds of attribution.

Could I send you the research notes? Either way I'd welcome any thoughts.

Best,
Ben Zanghi
ben@benzanghi.com

---

*Notes for sending:*

- **Verify email/affiliation** before sending (UBC math directory or his latest arXiv abstract).
- **Repo state:** decide whether to share the full repo or a curated subset. The corrected residual derivation + the retracted/corrected bounds are essential context for anyone evaluating the result — don't hide them; they're a credibility asset, not a liability.
- **Together AI angle:** their authors overlap with TTT-Discover. If White engages, a joint LB+UB paper is a natural follow-up — propose only if dialogue continues.
- **No-reply scenario:** brief 100-word follow-up after 2–3 weeks; beyond that, proceed to preprint without him.
