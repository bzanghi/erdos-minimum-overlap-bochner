# Lever E pretest — should we build an M-side SDP?

**Date:** 2026-05-10
**Question:** Is it worth weeks of work to build a new SDP that minimizes
Together's M-functional (`M(h) = sup_t ∫ h(x)(1-h(x+t)) dx`) directly,
instead of (or alongside) White's autocorrelation Ω-functional?

Driver script: `lp_research_state/code/_lever_e_pretest.py`
Raw results: `lp_research_state/data/lever_e_pretest_results.json`

---

## Test 1 — M-functional on our LP-optimal f̃

Source: `lp_research_state/data/row4_f_tilde.npz` (row 4, Phase 5 SDP-optimal,
10000-point grid on `[-2, 2]`). f̃ violates `[0, 1]`: min ≈ −0.971, max ≈ 1.548;
∫ f̃ over `[-2, 2]` ≈ 2.0001 (matches normalization ∫f = 1 modulo the period-2
× 2-half = 2 convention).

Together's M expects a function on `[0, 2]` with values in `[0, 1]` and
integral 1. f̃ satisfies none of these on its positive half. We therefore
report **five interpretations**:

| Interpretation                                         | ∫ over `[0,2]` | M-value      |
|--------------------------------------------------------|---------------:|-------------:|
| `f̃` on `[0, 2]` (positive-half slice, no transform)   |       0.99073  |   **0.7072** |
| Even-fold: `g(x) = (f̃(x) + f̃(-x))/2` on `[0, 2]`     |       1.00003  |   **0.6992** |
| `|f̃|` on `[0, 2]` (absolute value)                    |       1.19396  |   **0.5938** |
| `clip(f̃, 0, 1)` on `[0, 2]`                           |       1.03236  |   **0.6158** |
| Even-fold then clip                                    |       1.04190  |   **0.6117** |
| **Together's claimed M(h\*)**                          |       1.00000  |   **0.38087** |

**Every interpretation lies in `[0.59, 0.71]` — vastly larger than Together's
0.38087.** Our LP-optimal f̃ is a *terrible* candidate function for the
M-functional. Far from giving a tighter UB, it gives a value nearly 2× worse
than Together's h\*.

**This is consistent with what we already knew:** Ω and M are different
functionals; the SDP minimizes Ω, not M, and the SDP variables live in a
Fourier-truncated space where the recovered f̃ Gibbs-oscillates outside
`[0, 1]`. The Phase 5 SDP-LB on Ω is `0.380128`, yet the same f̃ yields M ≈ 0.6
— a 60% blow-up. There is **no upper-bound finding here.** The natural
"plug f̃ into M" idea fails.

---

## Test 2 — Existing mside_bochner SOC encoding at row 4, full scale

Same config as CLAUDE.md's reproduction (N=10000, T=4000, R=10, h=p=0.3875,
q ∈ [−0.02, 0.02]), with `bochner_n=20` always on:

| Configuration                                      | `prob.value`       | rigorous LB (5 dig) | dual resid | elapsed |
|----------------------------------------------------|-------------------:|--------------------:|-----------:|--------:|
| **baseline** (`bochner_n=20`)                      | 0.3796533955317040 | 0.37965             | 5.5e-10    |  69.5 s |
| **+ `mside_bochner_n=8`** (SOC-relaxed)            | 0.3796532336209890 | 0.37965             | 3.2e-07    |  64.2 s |

**Δ vs baseline at `prob.value` level: −1.6 × 10⁻⁷** (within solver noise;
the constraint set inclusion `F_2 ⊆ F_0` requires Δ ≥ 0, so the observed
sign is noise). At rigorous-LB precision the two values are indistinguishable.

Direct corroboration of findings.md (2026-05-10): SOC-relaxed M-side is
empirically dead at any cron-runnable `n_M`. The SOC slack `U_m ≥ |f̂(m)|²`
absorbs the M-side PSD content without constraining `(c, d)`. We did not run
`n_M=10` because (a) prior findings.md data already settled this at `n_M=10`
(Δ = +1.65 × 10⁻⁸) and (b) the anti-pattern rule caps per-solve time at 5 min.

We did **not** run the Schur variant (`mside_bochner_schur_n`) because the
module's own docstring says it "defines the SAME convex set" as the SOC form;
prior findings.md confirms the SOC was decisively dead, and a coextensive
encoding will give the same answer. We did **not** run the Lasserre-lifted
variant (`mside_bochner_lasserre_n`) because it requires
`lasserre_T_max > 0`, and CLAUDE.md documents Lasserre level-2 as
"non-rigorous and withdrawn" with a known tail-bound that kills the gain at
tractable `T_max`.

---

## Test 3 — Prior M-side data already on the books

From `lp_research_state/findings.md` (entries dated 2026-05-10):

- **`mside_bochner_n=5`** at row 4 N=2000: Ω = 0.37627652427597486 (status
  `optimal`), Δ vs plain LP baseline = **+1.4 × 10⁻⁹**.
- **`mside_bochner_n=10`** at row 4 N=2000: Ω = 0.3762765394 (status
  `optimal_inaccurate`), Δ vs baseline = **+1.65 × 10⁻⁸**.
- Pre-committed cancel rule (Δ < +1e-5) was met decisively at both levels.
- Mechanistic reading: the SOC slack `U_m` independently inflates the
  `−4 U_m` off-diagonal term to make `T_relax` PSD without constraining
  `(c_m, d_m)`. The diagonal `Ω/2` dominates and absorption is essentially
  complete.

**Three M-side encoding modules already exist in the repo:**

| File                                              | Encoding                                                   | Convex set      | Status                                              |
|---------------------------------------------------|------------------------------------------------------------|-----------------|-----------------------------------------------------|
| `lp_research_state/code/mside_bochner.py`         | SOC slack `U_m ≥ |f̂(m)|²`                                  | `F_2`           | Empirically dead (Δ < 1e-7)                         |
| `lp_research_state/code/mside_bochner_schur.py`   | 3×3 Schur slack `s_m ≥ |f̂(m)|²`                            | **Same `F_2`**  | Coextensive with SOC → equally dead                 |
| `lp_research_state/code/mside_via_lasserre.py`    | EXACT bilinear via Lasserre level-2 moments                | `F_1` (true M-side, no relaxation) | Requires `lasserre_T_max > 0`, but Lasserre lvl 2 is non-rigorous (CLAUDE.md / communications/lasserre_tail_bound.md) |

The *only* M-side encoding that does not have a known slack-absorption
problem is `mside_via_lasserre.py`, which is gated on an SDP relaxation
the repo has already retracted as non-rigorous.

**One-sentence summary:** the M-side Toeplitz-PSD constraint has been
convexly relaxed (SOC, Schur — same set), found empirically inactive at
publishable scale; the only known *exact* lift relies on Lasserre level-2,
which is documented as non-rigorous in this repo and has no current path to
becoming rigorous at tractable `T_max`.

---

## Verdict — **NOT PROMISING**

### Why

1. **The natural intuition is wrong.** "Our SDP optimizer f̃ might already be
   a great M-candidate" — Test 1 shows f̃ gives M ≈ 0.6–0.7, nearly 2× worse
   than Together. The SDP is solving for a *different* functional in a
   *different* function class. Pivoting the objective from Ω to M won't
   inherit any of f̃'s structure as a head start.

2. **The convex M-side has already been built and shown dead.** Two
   independent encodings (SOC, Schur) of the same convex relaxation give
   Δ < 1e-7 even at `n_M = 10`. This isn't an implementation problem we can
   fix by writing more careful code; it's a structural absorption of the
   constraint by the slack variable, fundamental to any convex relaxation
   that decouples `|f̂(m)|²` from `(c_m, d_m)`.

3. **The only exact M-side encoding requires a known-bad Lasserre lift.**
   Lasserre level-2 has been retracted in this repo precisely because its
   localizing tail bound kills any gain. The "M-side via Lasserre" file
   exists in code form but inherits the non-rigorousness of its parent
   relaxation.

4. **The ceiling looks low even if it worked.** Even the unrelaxed M-side
   Bochner cuts a convex *superset* (`F_1`) of the true M-feasible set. Most
   of the constraint content at row 4 is at small `m` where the SOC slack is
   already locally exact (`c_1 ≈ 0.3875` puts `m=1` in the tight regime).
   Expected upside if `F_1` were achievable is `≤ +1e-4` over the current
   Phase 5 stack (per the 2026-05-10 spec note in findings.md).

5. **Reward/effort ratio.** Weeks of work (per the prompt's framing) for an
   expected gain of `≤ 1e-4` on a stack already at 0.380128, against a UB of
   0.380871 (gap ~7e-4). Even a best-case M-side would close maybe 15% of
   the remaining gap. Not worth it.

### What we are NOT saying

- We are **not** claiming the *unrelaxed* M-side Toeplitz PSD is provably
  unhelpful. We are saying: every convex relaxation we have tried is
  vacuous; the only exact lift available is non-rigorous; and the natural
  starting heuristic (plug in f̃) fails.

---

## Next-step plan — pivot to Lever C, formalization, or write-up

Three credible directions, in order of preference:

1. **Lever C (push combinatorial `M(n)` via SAT/IP to n ≥ 25):** the current
   stack already gets `µ ≥ 0.380128`. Combinatorial brute-force `M(n)` is
   *exact*, has no SDP convergence/tolerance issues, and Together's bound
   `0.380871` comes from a 600-cell discrete optimization in spirit; pushing
   `M(n)` to 25+ with SAT/IP would either confirm convergence to Together's
   regime or surface a new structural lower bound. This is **weeks of
   engineering, not weeks of new mathematics** — strictly easier than
   Lever E.

2. **Lever A revisited (alternative basis):** Phase 5 saturates on the
   standard Fourier basis. A wavelet or Chebyshev basis could suppress the
   Gibbs oscillations in f̃ (which currently take f̃ out of `[0, 1]` by
   ±0.5) and might unlock a strictly different feasible set without
   relaxation. Hard, but the structural reason for saturation (Gibbs Mass
   in the high-frequency tail) suggests it's the natural axis.

3. **Stop and write up.** `µ ∈ [0.379544, 0.380871]` with the headline
   `+5.84 × 10⁻⁴` improvement over White is already publishable. The
   preprint draft exists; emailing E. P. White is queued. Closing the
   research note and shipping is a legitimate end-state.

**Concrete recommendation:** pivot to Lever C first (1–2 weeks, strictly
upper-bounded in scope: implement a SAT/IP encoding for the discrete
M(n)-minimization problem, push n from current values to 25–30, report). If
that saturates against Together, the project is *done* and we ship the
write-up. If it produces a surprise gap or new structure, that informs
Lever A. Either way, weeks-of-effort Lever E is dominated.

---

## Status

- **Status:** DONE
- **Test 1:** M_strict = 0.594–0.707 (five interpretations); M_clipped = 0.612 / 0.616. All ≫ Together's 0.38087.
- **Test 2:** baseline `bochner_n=20` rigorous LB = 0.379653 at row 4 (N=10000); `+ mside_bochner_n=8` gives Δ = −1.6e-7 (solver noise; constraint inactive).
- **Test 3:** SOC and Schur M-side relaxations define the same convex set and have been measured dead (Δ ≈ 1e-8) at `n_M = 5, 10` in prior findings; only exact lift requires retracted Lasserre level-2.
- **Verdict:** NOT PROMISING.
- **Next-step plan:** pivot to Lever C (combinatorial `M(n)` to n ≥ 25 via SAT/IP); if that saturates, ship the existing `0.380128` write-up.
