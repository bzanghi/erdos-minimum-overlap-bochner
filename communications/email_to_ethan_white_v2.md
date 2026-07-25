# Email to Ethan Patrick White — v2

> **SUPERSEDED — do not send.** Use `email_to_ethan_white_v6.md`. This draft also carries the
> wrong recipient: Ethan's address is **[personal address - see .local-contacts.md]**, confirmed 2026-07-25, not the
> ubc.ca address transcribed here.


**Status:** DRAFT — review and adjust salutation/details before sending.
**Supersedes:** [email_to_ethan_white.md](email_to_ethan_white.md) (drafted when µ ≥ 0.379544 was the result).

**To:** ethan.white@ubc.ca *(verify current affiliation/email before sending — UBC was Killam-Trusts/NSERC affiliation per arXiv:2201.05704)*
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** Bochner-PSD strengthening of your minimum-overlap LB to µ ≥ 0.380128 — and a question about a possible ceiling

---

Dear Professor White,

I hope this finds you well. I'm writing because I have a rigorous improvement to your 2023 lower bound `µ ≥ 0.379005` for the Erdős minimum-overlap problem (Acta Arith. 208), and a more interesting follow-up question about whether the technique class admits a provable ceiling.

## Result

Augmenting your Section-5 convex program with a sequence of *Bochner moment-matrix PSD* constraints — the Hermitian Toeplitz constraints `[f̂(j-k)]_{j,k=0..n} ⪰ 0` for both `f` and `1−f`, together with polynomial-moment constraints `m_{2k}(f) ≥ −tail_k`, a Hankel-PSD on those moments, and an iterated cover-refinement of your Table-3 ellipse centers — gives

$$\mu \;\ge\; 0.3801279$$

an improvement of `+1.123 × 10⁻³` over your published value, with your same proof structure (§5.1 / Appendix II ellipse extension), augmented duals, and slightly enlarged ellipses. Together's recent UB is `µ ≤ 0.380871`, so the remaining open gap is now `[0.3801279, 0.380871]` of width `7.4 × 10⁻⁴`.

Five compositional steps:

| Phase | Augmentation | µ ≥ | Δ vs prior |
|---|---|---|---|
| 0 (you, 2023) | base program | 0.379005 | — |
| 1 | Bochner-PSD n=20 + ellipse | 0.379544 | +5.4 × 10⁻⁴ |
| 2 | + cover refinement (4 new centers) | 0.379620 | +7.6 × 10⁻⁵ |
| 3 | + bochner_n=30 (12 centers) | 0.379879 | +2.6 × 10⁻⁴ |
| 4 | + poly_moment k=14 + Hankel-PSD n=6 | 0.380120 | +2.4 × 10⁻⁴ |
| 5 | + iterated cover on the combined stack | **0.380128** | +8.0 × 10⁻⁶ |

The Phase 4–5 increments are tiny — strongly suggesting the technique class is saturating, which is the more interesting question below.

## What's NOT in this number

A withdrawn Lasserre-level-2 attempt that initially looked like `+8.2 × 10⁻⁴` failed adversarial review when I worked out the Fejér-Riesz tail bound rigorously (notes in repo at `communications/lasserre_tail_bound.md`): the natural rigorization kills the gain quantitatively at currently-tractable T_max. M-side Bochner via SOC or Schur relaxation is empirically vacuous in the codebase (slack absorbs constraint content). The headline 0.380128 is the verified ceiling of the Bochner-PSD + poly-moment + Hankel + cover-refinement family on your program.

## Diagnostic: where the remaining gap lives

I ran a focused diagnostic plugging Together's piecewise-constant minimizer `h*` (their 600-step certificate that attains `µ ≤ 0.380871`) into the augmented SDP as a *primal candidate*:

- At Together's `h*` projected into your Fourier basis, **all of the CDE augmentations are healthily slack**. Bochner-PSD `λ_min M_30(f) = 6 × 10⁻⁵`, poly-moment slacks are 30×–5700× the tail bound, Hankel-PSD is positive with margin.
- The SDP's encoded Ω evaluates to `0.459` when (c, d) are pinned at Together's `h*` — `+0.07` above `h*`'s own actual autocorrelation `0.3873`.
- That excess can only come from your *cell-kernel cos/sin autocorrelation envelope* (the bounds at `white_full_convex.py:176-190` in my repo, mirroring §5 of your paper). The Bochner-PSD + cover-refinement augmentations do not sharpen those envelopes in the directions Together's `h*` probes.

So the empirical ceiling of `0.380128` appears to come from your cell-kernel envelope, not from the relaxation of `f ∈ [0, 1]`.

## The question I'd love your read on

**Do you have any intuition about whether the cell-kernel autocorrelation envelope admits a provable upper bound on µ via your framework?**

Concretely: is there a closed-form Fourier identity, or a known duality-gap theorem (Cohn-de Laat-Salmon type, Beurling-Selberg type), that would pin down

$$C^* \;:=\; \sup\bigl\{\, c \;:\; \text{augmented Bochner-PSD + cell-envelope SDP can certify } \mu \ge c \,\bigr\}$$

— either constructively, or as a strict-duality-gap statement?

If `C* < 0.380871`, that would prove (a) Together's UB cannot be matched by any Bochner-PSD augmentation of your program, and (b) closing the remaining gap requires a qualitatively different technique. That's a clean publishable saturation theorem.

Two follow-up technicalities you might know off the top of your head:

- Your L² autoconvolution paper (arXiv:2210.16437) gives uniqueness of the L² minimizer to ~0.001% via Fourier-truncated QCLP + Minkowski strict-convexity. The L∞ functional `sup_t (h⋆h)(t)` here looks structurally different (no polynomial-in-Fourier-coefs expansion); is there an obstruction-of-principle to the transfer, or just a methodological gap?
- Is there a known cleaner way to write `sup_t (h⋆h)(t)` directly (rather than via the per-`m` cell-kernel envelope) that I'm missing?

I'd take a "no idea" or "haven't thought about it" as completely fair — but a partial intuition would be very valuable to direct the next push.

## Other notes that might be of interest

- The even-`f*` conditional bound is `µ_even ≥ 0.379904` — i.e., your §6 evenness assumption tightens the result barely at all beyond Bochner-PSD itself, suggesting the augmented optimum is close to even.
- Combinatorial M(n): I extended Haugland's exact-value table from `n ≤ 15` to `n=20` via SAT (M(19)=8, M(20)=8 — both new). The integer ratios at small n are far from informative on µ; the gap-closing route empirically goes through continuous step-function constructions, not exact M(n).
- The technique enumeration this session ruled out 10 distinct candidate "next levers" via explicit construction or analytic argument. The diagnostic, the lever ledger, and the saturation evidence are all in the repo.

## Disposition

I'd like to write this up as a short note, with you as co-author if you're interested in commenting or refining, or acknowledged as the indispensable foundation otherwise. The code, three independent re-implementations, and full diagnostic memos are in a Git repo I'd be happy to share before any publication.

The contribution here is honestly modest — "noticed Bochner is missing from your program, verified the extension carefully, ran a careful saturation diagnostic." The structural ideas are yours. But the saturation diagnosis (and the implicit conjecture that the Bochner+cell-envelope class has a ceiling strictly below Together's UB) seems like the more interesting follow-on question, and you'd know whether that's a plausible ceiling-theorem direction.

I'm a software engineer by day; the work was done as a personal project assisted by Anthropic's Claude. I'd be happy to disclose the AI-assisted workflow in detail if useful — it makes a difference for some kinds of attribution.

Could I send you the repo and the research note? Either way I'd welcome any thoughts.

Best,
Ben Zanghi
ben@benzanghi.com

---

*Notes for sending:*

- **Email/affiliation:** Verify current. arXiv:2201.05704 shows UBC + Killam Trusts/NSERC support 2022; recent talks (UBC March 2023) suggest he was still there. A quick search for "Ethan Patrick White" + current institution should give the right address. Failing that, the math department directory at UBC.
- **Tone:** Intentionally formal-but-warm. Lead with the result, ask the more interesting question second, offer co-authorship, acknowledge foundation.
- **Repo state:** Don't send before the Git repo is public OR you've decided which subset of it to share. He should be able to inspect the code immediately when he opens the email.
- **Saturation framing:** The "ceiling theorem" question is the genuinely novel ask; if he engages with that, it could become a much more interesting collaboration than just confirming the +1.1e-3 result. Lead the body with the bound, lead the *question* with the saturation.
