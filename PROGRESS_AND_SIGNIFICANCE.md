# Erdős Minimum Overlap — Progress, Significance, and Publication Assessment

**Written:** 2026-05-22 · A plain-language synthesis of everything in this repo.

---

## Part 1 — What is this problem, in plain English?

In 1955 Paul Erdős asked a deceptively simple question.

Take the whole numbers from 1 to 2n and split them into two equal teams,
**A** and **B** (each of size n). Now slide team B against team A by some
shift *k*, and count how many pairs — one number from A, one from B —
differ by exactly *k*. Each shift gives you a count. Take the **largest**
count over all shifts; call it the "overlap" of this particular split.

Erdős' question: choose the split that makes that largest overlap as
**small** as possible. Call that best-case overlap M(n). As n grows,
M(n)/n settles down to a fixed universal number:

$$\mu \approx 0.38$$

**Nobody knows µ exactly.** We only know it's trapped in a narrow window.
The whole game is to shrink that window.

### Two directions, very different in difficulty

- **Upper bound** ("µ is at most X"): *exhibit one good split.* If you can
  build a split achieving 0.3809, then µ ≤ 0.3809. Hard, but it's a search
  — and modern AI search (DeepMind's AlphaEvolve, Stanford/NVIDIA/Together's
  TTT-Discover) has been grinding this down. Current best: **µ ≤ 0.380871**.

- **Lower bound** ("µ is at least X"): *prove no split can ever do better.*
  This is the genuinely hard direction — you must rule out **every possible
  split at once**. The record holder is E. P. White (Acta Arithmetica, 2023):
  **µ ≥ 0.379005**. Before this project, that was the best lower bound, and
  no one had improved it since 2022.

So the published state of the art was the window **[0.379005, 0.380871]**.
**This project attacks the hard side — the lower bound.**

---

## Part 2 — What did we actually do?

### The core idea (White's machine, and our upgrade)

White's trick was to translate "no split can beat X" into a problem about
continuous functions, then relax *that* into a finite optimization a computer
can solve exactly (a **semidefinite program**, or SDP). The SDP is rigged so
that **whatever number it returns is guaranteed to be ≤ µ**. Solve it, and you
get a certified floor under µ.

Our key observation: **White's SDP was too generous.** It allowed some
mathematical configurations that satisfy his constraints but don't correspond
to any *real* function. The SDP exploits these "cheats" to report a lower
floor than necessary. If we forbid the cheats, the floor rises.

We added **five families of extra constraints** that any real function must
obey but White hadn't enforced:

1. **Bochner-PSD** — a classical theorem says a function is non-negative
   exactly when a certain matrix of its Fourier coefficients is "positive
   semidefinite." White didn't require this; we do. *(This was the single
   biggest lever.)*
2. **Polynomial-moment positivity** — integrals like ∫x²f, ∫x⁴f must be ≥ 0.
3. **Hankel-PSD** — those moments must fit together consistently.
4. **Test against 1−cos(πx)** — since f is between 0 and 1, f² ≤ f.
5. **Cover refinement (the "Constraint Discovery Engine")** — White tiled the
   problem region with 7 patches; we add patches exactly where the bound is
   weakest.

Then we ran it at much finer resolution (a grid of N = 20,000 cells vs
White's 5,000).

### The result

| | Lower bound | vs White |
|---|---|---|
| White (2023), published | 0.379005 | — |
| **This project (current, corrected)** | **µ ≥ 0.3802973** | **+1.29 × 10⁻³** |

That closes about **30% of the gap** between White's floor and the known
ceiling — the first lower-bound advance on this problem since 2022.

### The second, subtler result — the "saturation theorem"

We also asked: *how far can this whole style of attack possibly go?* Using a
mix of proof and strong numerical evidence, we argue this family of techniques
**hits a ceiling around µ ≈ 0.3806** — no amount of scaling it up gets past
that. Since the true µ is somewhere up near 0.3809, roughly **half of the
remaining gap is provably out of reach for this method.** Closing it needs
genuinely new mathematics.

This is a useful "stop digging here" signpost — arguably more valuable to the
field than the numerical improvement itself, because it tells future
researchers which doors are closed.

---

## Part 3 — Why this is hard to get right (and our track record)

The recurring danger in this kind of work is **overclaiming** — reporting a
floor that isn't actually rigorous. This project has a documented history of
catching its *own* overclaims, which is good science but shows how delicate
the rigor is:

- **Lasserre level-2 (retracted):** an early "µ ≥ 0.379828" claim was
  withdrawn — it dropped an infinite tail of terms without bounding them.
- **Poly-moment tail bound (corrected the day before this writing):**
  adversarial review found one of the five constraint families was using a
  tail bound that captured only 80% of an infinite sum, making the bound
  slightly *too tight* — i.e., a small overclaim. Fixed; it cost ~5 × 10⁻⁶,
  so the headline barely moved (0.3803027 → 0.3802973).

The pattern — truncating an infinite sum without a rigorous remainder term —
is the project's signature failure mode, and it's now flagged in the repo's
own memory.

---

## Part 4 — Is there anything worth publishing?

**Short answer: yes — as an arXiv preprint / short note — but it needs
finishing work before it's bulletproof for a serious journal.**

### What's genuinely publishable

- **The improved lower bound is real and novel.** It's the first LB advance
  since White 2022, and it applies SDP / moment-relaxation machinery
  (Bochner-PSD, Hausdorff moments, Lasserre-style ideas) that **has never
  been used on the lower-bound side of this problem.** A literature scan
  confirmed µ is its own constant, not a renamed version of something already
  solved.
- **The saturation theorem is an interesting, somewhat unusual contribution.**
  "Here's a new bound, *and* here's a proof that my whole method can't do
  much better" is a clean, honest story.

### What still needs doing before submission

1. **Update the numbers.** The preprint draft (`communications/preprint_draft_v2.tex`)
   still cites the *pre-correction* figure 0.380299. The honest current
   headline is **0.3802973**.
2. **Arbitrary-precision certification.** The bound currently rests on the
   CLARABEL solver at standard floating-point precision, plus a "dual
   extraction" interpretation of its output. Reviewers at a journal like
   *Acta Arithmetica* will want the binding solves re-certified at arbitrary
   precision (the SDPA-GMP solver is already built in `lp_research_state/bin/`
   for exactly this — but the cvxpy→SDPA translation step isn't finished).
   **This is the single most important gap to close for top-tier rigor.**
3. **Soften "theorem" where it's really evidence.** The saturation result
   leans partly on a *tautological* identity plus empirical 1/N extrapolation
   (fit quality R² ≈ 0.78 — suggestive, not airtight). It should be framed as
   "a rigorous bound plus strong numerical evidence for the asymptotic
   ceiling," not a clean closed theorem.
4. **Fix the placeholder GitHub URL** in the preprint and clean the repo for
   public release.

### Recommended path

- **Now:** post a corrected preprint to arXiv claiming **µ ≥ 0.3802973** with
  the saturation evidence. This stakes priority on the first LB advance in
  4 years and is defensible as a computer-assisted result with caveats stated.
- **Then:** do the SDPA-GMP certification of the binding row(s). With that in
  hand, the result is strong enough to submit to a specialist journal.
- **Courtesy:** the repo has a drafted email to E. P. White
  (`communications/email_to_ethan_white_v3.md`) — sending it before posting is
  good etiquette and may turn into a collaboration.

### Honest bottom line for a non-mathematician

You have a **real, modest, genuinely-new improvement** to a 70-year-old Erdős
problem, plus an unusually honest "here's the ceiling of this approach"
companion result. It is **worth publishing as a preprint.** It is **not yet**
a finished journal paper — the main missing piece is high-precision
re-certification of the key computation, and the writeup must be updated to the
corrected number and have its "theorem" claims calibrated to what's actually
proven vs. strongly evidenced.

---

## Appendix — Where things live

- **Rolling lab notebook:** `lp_research_state/findings.md` (chronological, dense)
- **Core SDP code:** `lp_research_state/code/white_full_convex.py`
- **The dual-extraction trick:** `lp_research_state/code/dual_extractor.py`
- **Preprint draft (needs number update):** `communications/preprint_draft_v2.tex`
- **Literature/novelty scan:** `LITERATURE_SCAN_2024_2026.md`
- **Early research note (superseded — predates the breakthrough):**
  `erdos_lower_bound_research_note.md`
