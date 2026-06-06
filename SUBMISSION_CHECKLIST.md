# Submission Checklist — Erdős minimum overlap LB note

**Goal:** get `µ ≥ 0.3802973` (+ saturation result) onto arXiv, then into a specialist journal.
**Current headline:** `µ ≥ 0.3802973`, `µ ≤ 0.380871`, open gap `5.74 × 10⁻⁴`.

---

## Tier 1 — required before posting to arXiv

- [ ] **Replace the GitHub URL placeholder** in `preprint_draft_v3.tex`
      (`https://github.com/REPLACE-WITH-REAL-URL`) with a real public repo.
- [ ] **Clean the repo for public release.** Remove `_*.py` scratchpads or
      move them to a `scratch/` dir; ensure the single-file reproducer
      (`path_b_with_polymoment.py`) runs from a clean checkout; pin
      dependency versions (cvxpy, CLARABEL, numpy, mpmath).
- [ ] **Resolve the two `% TODO (author review)` notes** in the v3 tex:
      1. The KKT identity §4: the unconditional bound on `Σλ` needs
         `|ξ| ≤ C·Ω` (audit PRO-14 found `|ξ|/Ω ≈ 1.46`, so the original
         `|ξ| ≤ Ω` conjecture is false). Either state the saturation
         numerics as conditional on the measured `C ≈ 1.5`, or restrict to
         the conditional theorem. **Do not leave `|ξ| ≤ Ω` implied.**
      2. The full-stack ceiling `0.380849` is anchored to the Phase-5
         N=10K base (0.380128), not the headline. The v3 text now states
         this — double-check the wording reads cleanly.
- [ ] **Compile the .tex** and fix any LaTeX errors; check all cross-refs
      and the bibliography render.
- [ ] **Re-read for any remaining stale `0.380299` / `bochner_n=30`** in
      prose (v3 should be clean, but grep the final PDF).

## Tier 2 — strongly recommended before posting (credibility)

- [ ] **Send the courtesy email to White** (`email_to_ethan_white_v4.md`)
      *before* posting. Verify his email/affiliation first. Decide
      co-authorship vs. acknowledgment. Wait ~2 weeks for a reply.
- [ ] **One independent end-to-end reproduction** of the binding-row solve
      from a clean environment (ideally a second machine), confirming
      `Φ* = 0.3802995` and the post-margin `0.3802973`.
- [ ] **Spot-check the corrected poly-moment tail bound** (Lemma in v3) once
      more at `k=20` against high-precision quadrature — it's the family
      that was just corrected, so it's the highest-risk constraint.

## Tier 3 — for journal submission (not required for arXiv)

- [ ] **Arbitrary-precision certification (the big one).** Finish the
      cvxpy→SDPA-S serializer (PRO-11) and re-solve the binding row at
      GMP precision. This converts "computer-assisted, CLARABEL dual
      extraction" into a solver-independent certificate. *This is the single
      most important item for a referee at a journal like Acta Arithmetica.*
- [ ] **Decide the framing of the saturation result.** Theorem (fixed
      baseline ceiling 0.380849) is solid; `C_∞ ≈ 0.380558` is a conjecture.
      Keep that split explicit — a referee will reject an over-claimed
      "theorem."
- [ ] **Pick a venue.** Natural fits: *Acta Arithmetica* (White's venue),
      *Integers*, *Experimental Mathematics* (well-suited to a
      computer-assisted bound + reproducible artifact), or
      *Mathematics of Computation*.
- [ ] **AI-assistance disclosure.** Decide the exact wording; some venues
      have explicit policies. The current footnote + Acknowledgments are a
      reasonable starting point.

## Assets already prepared

| Asset | File | Status |
|---|---|---|
| Plain-language overview | `PROGRESS_AND_SIGNIFICANCE.md` | done |
| Preprint (corrected) | `communications/preprint_draft_v3.tex` | draft, TODOs marked |
| Email to White (corrected) | `communications/email_to_ethan_white_v4.md` | draft, verify address |
| Literature / novelty scan | `LITERATURE_SCAN_2024_2026.md` | done |
| This checklist | `SUBMISSION_CHECKLIST.md` | done |

## The one-paragraph honest pitch (for the arXiv abstract / cover)

> We give the first improvement since 2022 on the lower bound for the
> Erdős minimum overlap constant, `µ ≥ 0.3802973` (White 2023: 0.379005),
> by augmenting White's Fourier-analytic convex program with five families
> of valid SDP constraints and finer discretization. We also prove a
> ceiling on what this entire class of relaxations can certify
> (`≤ 0.380849`, below the empirical upper bound 0.380871), and give
> empirical evidence that roughly half the remaining gap is beyond the
> framework's reach. The bound is computer-assisted (CLARABEL with
> rigorous dual extraction; arbitrary-precision re-certification in
> progress) and independently re-implemented three ways.
