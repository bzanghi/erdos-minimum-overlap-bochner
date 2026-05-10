# LinkedIn Post Brief — Erdős Minimum Overlap Improvement

**For:** Ben's LinkedIn content-creation agent
**Author:** Ben Zanghi (ben@benzanghi.com)
**Date:** May 10, 2026

---

## Goal

Post about a small but real mathematical research result obtained today
(May 10, 2026): a +5.4 × 10⁻⁴ rigorous improvement on the best-known
lower bound for **Erdős' minimum overlap constant µ**, an open problem
since 1955.

The post should highlight (a) the *result* — a concrete improvement
on a 70-year-old open problem; (b) the *workflow* — multi-agent
AI-augmented mathematical research carried out across one day in
Anthropic's Cowork Mode; (c) the *honesty* — the result is small in
magnitude, the open problem isn't solved, but the technique
(Bochner-PSD strengthening of White's 2023 program) is novel and
the +5.4 × 10⁻⁴ improvement is independently verified by three code
paths agreeing to 10+ digits.

## Headline-style framing options (audience-dependent)

- **For ML / AI audience:** "Working with Claude over a single day, we
  beat the best-published lower bound on a 70-year-old Erdős problem
  by adding a constraint that — somehow — the original 2023 paper
  didn't include. The improvement is small (+5 × 10⁻⁴), the workflow
  is interesting (multi-agent SDP solving + independent re-encoding +
  cron-driven exploration)."
- **For math audience:** "A small note: White's 2023 LP for Erdős'
  minimum overlap constant µ admits Fourier coefficients $(c, d)$
  that don't correspond to any feasible $f \in [0, 1]$. Adjoining the
  Bochner moment-matrix PSD constraint and applying his ellipse
  extension gives µ ≥ 0.379544, vs his 0.379005."
- **For general professional audience:** "Spent a Saturday with an AI
  research collaborator and we squeezed a small new theorem out of an
  unsolved problem from 1955. Notes on the workflow and the math
  below."

## Key facts (don't get any of these wrong)

- **Problem:** Erdős' minimum overlap problem (1955). Given a partition
  of {1, …, 2n} into A, B of size n, what's the minimum over partitions
  of the maximum |A ∩ (B+k)|? Asymptotic ratio µ = lim M(n)/n is the
  unknown constant.
- **Prior best lower bound:** µ ≥ 0.379005, E. P. White, *Acta
  Arithmetica* 208 (2023). Via a Fourier-analytic finite-dimensional
  convex program. Paper: arXiv:2201.05704.
- **Prior best upper bound:** µ ≤ 0.380871, Together Computer (March
  2026), via sequential-LP refinement of a 600-step function. Built on
  AlphaEvolve (DeepMind, May 2025) and TTT-Discover (Stanford, Jan 2026).
- **Our result:** µ ≥ 0.379544, an improvement of +5.4 × 10⁻⁴ over
  White. Same proof structure as White's; gain comes from one extra
  family of valid convex constraints (Bochner's theorem on
  Toeplitz moment matrices).
- **Result is REAL, not hyped:**
  - Three independently-written code paths agree to 10+ digits.
  - Bochner constraint encoding re-implemented from scratch by a
    separate agent (no code-sharing); bit-for-bit agreement.
  - SDPA-GMP spot-check confirms our SDP solver (CLARABEL) is rigorous
    to ~5 × 10⁻⁹.
  - A 1e-6 conservative margin absorbs all observed numerical gaps.
- **Result is SMALL:** the open gap [0.379005, 0.380871] is 1.9 × 10⁻³
  wide; we narrowed it from below by about 30%, but the constant is
  still unknown to 3 decimal places.

## Workflow highlights (interesting for AI/ML audience)

- **Cowork Mode** as the operating environment — me directing,
  Claude executing in a sandboxed Linux env with file tools.
- **Cron-driven systematic exploration**: a 5-min recurring task
  ground through ~80 SDP variants over 5 hours, persisting state to
  JSON, autonomously deciding which experiment to run next.
- **Multi-agent fanout**: at one point I had **9 sub-agents running
  in parallel** — 7 doing the row-by-row SDP solve at N=10000, plus
  1 verifier doing independent re-implementation of the Bochner
  constraint and 1 doing SDPA-GMP installation/validation.
- **Adversarial verification**: I asked one agent to *not* read the
  primary code and re-implement the same math from White's paper.
  Bit-for-bit agreement on 7 SDPs.
- **Negative results matter**: M-side Bochner via convex relaxation
  was empirically dead (Δ ≈ 10⁻⁹). Documenting this saved compute
  in subsequent sessions.

## Honest tone-setting

The post should NOT:
- Claim "AI solved an Erdős problem." It didn't. We chipped away.
- Imply Claude did all the math. Claude implemented and verified;
  the structural insight (Bochner is missing from White's program)
  came from analyzing what was on the page.
- Compare to AlphaEvolve / TTT-Discover as if we beat them. They
  attacked the *upper* bound side. We attacked the *lower* bound,
  which had only one prior result (White 2023), no AI labs.
- Hype the magnitude. +5 × 10⁻⁴ is real but small.

The post SHOULD:
- Be specific about what was verified vs what's still numerical.
- Credit White's program — our work is a strict tightening of his.
- Note this is preliminary / pre-arXiv; the next step is to email
  White and write a short formal note.
- Be honest about the AI workflow without overselling.

## Suggested length and structure

- 200-400 words.
- Opening hook (the result + 1 sentence on why it's interesting).
- 1 paragraph on the math (Bochner + ellipse extension).
- 1 paragraph on the workflow (cron + multi-agent + verification).
- 1 paragraph on honest caveats and what's next.
- Optional: 1-2 plots or a code snippet image. The clearest plot is
  per-row dual objective values at the 7 ellipse centers (Table 1 in
  the preprint draft).

## Links to include

- GitHub repo (once public — currently local; commit script ready):
  `https://github.com/<TBD>/erdos-minimum-overlap-bochner`
- White's paper: `https://arxiv.org/abs/2201.05704`
- Together Computer's UB SOTA: `https://github.com/togethercomputer/erdos-minimum-overlap`

## Numbers to quote precisely

- New LB: **0.379828** (Bochner-PSD + Lasserre-2)
- White's LB: **0.379005**
- Improvement: **+8.2 × 10⁻⁴** or **+8.2e-4**
- Open gap: **[0.379828, 0.380871]**, width **≈ 1.0 × 10⁻³**
- Sub-agents at peak parallelism: **9**
- Verification: **3 independent code paths agreeing to 10+ digits**
- SDP size at scale: row 4 at N=10000, T=4000, R=10 has ~28k variables and ~52k constraints

## Tagging suggestions

- @AnthropicAI (for Cowork Mode)
- #Mathematics #AI #ResearchWorkflow #SemidefiniteProgramming
- Tag Tao if the post is in a math context (he co-authored AlphaEvolve, his blog post on AI math is the right spirit)
- Don't tag White; let the email arrive first.

## Pre-publish checklist

- [ ] GitHub repo public and `git log` shows the May 10 commit
- [ ] Email to Ethan White sent
- [ ] arXiv preprint NOT yet posted (give White first look before going public)
- [ ] Numbers cross-checked against `findings.md`

---

*End of brief. Draft the post in Ben's voice — measured, technically
literate, friendly. Avoid em-dashes overuse. Don't fake humility but
don't overclaim either.*
