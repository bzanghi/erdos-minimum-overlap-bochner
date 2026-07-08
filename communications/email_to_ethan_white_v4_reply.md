# Email to Ethan White — v4 (REPLY to his 2026-07 response)

**Status:** DRAFT — Ben reviews before sending.
**Supersedes:** [email_to_ethan_white_v3.md](email_to_ethan_white_v3.md) (written as a cold intro; Ethan has since replied, so this is a threaded reply).
**Context:** Ethan confirmed the Bochner-PSD constraint "seems valid," suggested interval-average-of-f variables as an alternative, and reported two errata in the published program: (5.6)/(5.7) numerator 8 → 4, and (5.8)/(5.9) RHS m → 2m−1.

---

Dear Ethan,

Thank you for the generous reply — and especially for the two corrections. Since our program is a direct transcription of your Section 5, both errata pointed straight at specific lines of our code, and chasing them down was immediately useful:

**On (5.8)/(5.9):** our implementation already used 2m−1 in the right-hand sides — when transcribing we took the index from the left-hand side δ_{2m−1}/ε_{2m−1}, so this erratum turns out not to affect us (the corrected form is the looser one, so it's the valid direction in any case).

**On (5.6)/(5.7):** here we had faithfully inherited the paper's 8. Re-solving at our binding center (your row 4, with the Bochner-PSD augmentation at n=20) with the corrected 4, the optimum moves by only +1.8 × 10⁻⁶ — consistent with your "not material" finding, and in the fortunate direction: the corrected program is (very slightly) *tighter*, so the values we've been reporting remain valid lower bounds, and re-running the full pipeline with the 4 should marginally improve them. That full re-run is queued; I'll use the corrected constant everywhere going forward.

**On average-values-of-f variables:** I think your instinct about why the PSD constraint helps is exactly right. The moment matrix encodes "f and 1−f are nonnegative" spectrally, with no new variables beyond the existing Fourier coefficients — it's the Bochner/Herglotz criterion truncated to an n×n principal minor. We did also test the interval-average route indirectly: a diagnostic that projects the best known upper-bound construction into your variable space shows the f-side constraints (including the PSD blocks) are healthily slack there, and essentially all of the relaxation gap sits in the cell-kernel envelope on the M side — which suggests average-f variables would add little beyond what the PSD blocks capture, and that the M-side envelope is where the framework's slack lives.

Two updates since I last wrote, in case they're of interest:

1. With the PSD blocks plus polynomial-moment (Hausdorff) cuts with analytic tail bounds, a Hankel-PSD family, an iterated refinement of your Table-3 ellipse cover, and N = 20,000, the bound through your Section 5.1 extension argument now stands at **µ ≥ 0.380302** (Mosek, genuine `optimal` status at all 12 centers, duality gaps ≤ 5 × 10⁻¹¹). Against the best current construction-side value (µ ≤ 0.380871, Together AI's 600-piece step function), the open window is ~5.7 × 10⁻⁴.

2. Possibly the more interesting development: we derived a computable ceiling for this technique class — using the measured dual multipliers of the cell-envelope family, no augmentation of the cosine+sine envelope can prove beyond ≈ 0.38071 at any tested scale. And on the construction side, Together's step function turns out to be a numerically exact KKT point of the discretized min-max problem (stationarity residual ~10⁻⁸), and it stays stationary — to first *and* second order — when the grid is refined from 600 to 1200 cells. (We initially concluded the opposite from a mis-derived stationarity equation that dropped the domain-edge terms; the corrected test reversed it.) Taken together: 0.380871 looks like a serious candidate for µ itself, and if that's right, closing the gap from below provably requires something outside this convex-relaxation framework. Do you have any instinct on whether the true extremizer should be expected at such a "flat" degenerate optimum — the active set of shifts has positive density — or is that degeneracy itself suspicious to you?

I'm writing this up as a short note (the bound plus the ceiling theorem); I'd be glad to send you the draft and the research notes, and happier still if you wanted to weigh in on the framing — with co-authorship open if any of it grows into something you'd want your name on, and a grateful acknowledgment otherwise.

Thanks again for taking the time, and for the errata — exactly the kind of detail that would have been very hard to catch from outside.

Best wishes,
Ben

---

*Notes for sending:*

- Reply in-thread to Ethan's email (keeps context; no need to re-introduce).
- The +1.8 × 10⁻⁶ figure is from row 4, N=3000, T=1200, bochner_n=20, CLARABEL (`_pro35_erratum_ab.py`); if the full-scale corrected re-run lands before sending, update the sentence with the headline-level number.
- The "healthily slack f-side / binding M-side envelope" claim is from TOGETHER_DIAGNOSTIC.md.
- The KKT/second-order claims are PRO33_KKT_CORRECTION.md and PRO34_UB_REFINEMENT.md; the self-correction is disclosed deliberately (credibility, and he'll appreciate the edge-term subtlety).
- Kept to ~600 words; the complementarity question from v3 is *dropped* here to keep the reply focused — it can go in a follow-up or the note itself.
