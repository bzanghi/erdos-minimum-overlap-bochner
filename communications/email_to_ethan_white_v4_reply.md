# Email to Ethan White — v4 (REPLY to his 2026-07 response)

**Status:** DRAFT — Ben reviews before sending.
**Supersedes:** [email_to_ethan_white_v3.md](email_to_ethan_white_v3.md) (written as a cold intro; Ethan has since replied, so this is a threaded reply).
**Context:** Ethan confirmed the Bochner-PSD constraint "seems valid," suggested interval-average-of-f variables as an alternative, and reported two errata in the published program: (5.6)/(5.7) numerator 8 → 4, and (5.8)/(5.9) RHS m → 2m−1.

---

Dear Ethan,

Thank you for the reply, and especially for the two corrections — my code is a direct transcription of your Section 5, so both pointed at specific lines. The (5.8)/(5.9) indexing I'd already transcribed as 2m−1 (taken from the LHS), so no effect. The (5.6)/(5.7) factor of 8 I had inherited; I've adopted the 4 and re-run the full pipeline on the corrected program. At the binding center the change moves the optimum by only +1.8 × 10⁻⁶ — consistent with your "not material" finding — and in the fortunate direction: the corrected program is slightly tighter, so earlier reported values were valid in any case.

On average-values-of-f variables: I think your instinct is right — the PSD constraint packages f, 1−f ∈ [0, 1] spectrally without new variables. Interestingly, projecting the best known upper-bound construction into your variable space shows the f-side constraints are all slack; essentially the whole relaxation gap sits in the cell-kernel envelope on the M side.

One update that may interest you: the augmented program (PSD blocks + polynomial-moment cuts with tail bounds + your 5.1 extension, N = 20,000, corrected constants) now gives **µ ≥ 0.380284**, certified over your entire parameter space rather than only the (5.16) residual region. And on the other side, Together AI's 600-piece step function (µ ≤ 0.380871) turns out to be a numerically exact KKT point of the discretized min-max problem — and stays stationary, to first and second order, under refinement to 1200 cells. So 0.380871 looks like a serious candidate for µ itself. Its optimum is strangely flat (the active set of shifts has positive density) — does that degeneracy match your intuition for the true extremizer, or is it suspicious to you?

For honesty: I'm a software engineer, not a mathematician — this is a personal curiosity project done in close collaboration with an AI assistant (Anthropic's Claude), with everything load-bearing cross-checked mechanically (independent re-implementations, multiple solvers, dual certificates). I'm writing it up as a short note with that workflow disclosed; I'd be glad to send you the draft, and gladder still if you wanted to weigh in.

Thanks again — corrections like these are exactly what's hardest to catch from outside.

Best wishes,
Ben

---

*Notes for sending:*

- Reply in-thread to Ethan's email (keeps context; no need to re-introduce).
- The +1.8 × 10⁻⁶ figure is from row 4, N=3000, T=1200, bochner_n=20, CLARABEL (`_pro35_erratum_ab.py`). Post-merge with main (2026-06-01 session): `mside_sin_coeff=4.0` is the code default and the pipeline has been re-run on the corrected program; the quoted µ ≥ 0.380284 is main's unconditional full-space headline (core-region value is 0.3802973).
- The "healthily slack f-side / binding M-side envelope" claim is from TOGETHER_DIAGNOSTIC.md.
- The KKT/second-order claims are PRO33_KKT_CORRECTION.md and PRO34_UB_REFINEMENT.md; the self-correction is disclosed deliberately (credibility, and he'll appreciate the edge-term subtlety).
- Trimmed to ~350 words per Ben. Cut: the ceiling-theorem detail (now just implied), the co-authorship formalities (folded into "weigh in"), the self-correction disclosure (kept for the note itself), and v3's complementarity question (follow-up material).
- Per Ben: the engineer/AI/personal-project disclosure is deliberate and stays — one sentence version.
