# Approach ③ — Representation-theoretic symmetry reduction (complex-Hermitian Bochner)

**Date:** 2026-06-06
**Lane:** the "free ~4×" solver-cost multiplier / enabler (THREE_APPROACHES.md §3).
**Honest framing:** this is a **solver-cost multiplier**, NOT a value-changer. It does
**not** pass the framework ceiling C_∞ ≈ 0.380558. Its purpose is to make the
high-`bochner_n` regime (LEVER-F3) and Approach ① tractable on the same hardware.

## What was implemented

The verified Bochner block (white_full_convex.py §Bochner / bochner.py) encodes the
`(n+1)×(n+1)` complex Hermitian Toeplitz moment matrix `M_n(f)` as a
`2(n+1)×2(n+1)` **real** PSD constraint via the embedding `[[Re,−Im],[Im,Re]]`, for
both `f ≥ 0` and `1−f ≥ 0`.

This approach replaces it with a **direct complex Hermitian PSD**:
`cp.Variable((n+1,n+1), hermitian=True)`, `H >> 0`, with Toeplitz/coefficient ties
(`f̂(0)=1/2`, `f̂(k)=(c_k − i d_k)/2`). cvxpy 1.8.2 supports complex PSD. This halves
the PSD matrix dimension → the documented ~4× flop/memory cut (2·(½)³ = ¼), for free,
because the real embedding's spectrum is exactly the complex spectrum **doubled**
(every eigenvalue multiplicity 2 under the U(1)/SO(2) action).

**Files (all NON-DESTRUCTIVE — verified machinery untouched):**
- `lp_research_state/code/bochner_hermitian.py` — complex-Hermitian Bochner constraint
  + numeric cross-check helper. Self-test: numeric Hermitian↔real-embedding agreement
  `1.1e-15` over 20 random trials.
- `lp_research_state/code/white_full_convex_hermitian.py` — `build_problem_hermitian`,
  a thin wrapper that delegates to the verified `build_problem(..., bochner_n=0)` and
  bolts on the Hermitian block, so the rest of the SDP is byte-for-byte the verified
  program.
- `lp_research_state/code/_herm_equiv_check.py` — the mandatory 10-digit equivalence gate.
- `lp_research_state/code/_herm_win_measure.py` — wall-time / peak-memory / max-`bochner_n`
  comparison.

## Results

Numbers in `sym_reduction_result.json`. Headline: **the swap is mathematically
EXACT (equivalence gate PASSES), but it delivers NO solver-level speedup through
the cvxpy stack** — cvxpy's mandatory `Complex2Real` reduction expands
`hermitian=True` PSD back into the identical real embedding before any solver sees
it. The "free ~4×" does not materialize via `cp.Variable(hermitian=True)`.

### Step 2 — equivalence cross-check (rigor gate): **PASS, swap is exact**

The naive criterion "objective agrees to ≥10 sig digits" is **provably
unachievable** here and is the *wrong* test: two **independent CLARABEL solves of
the identical problem** (constraint list merely shuffled) agree to only **7.4–9.2
sig digits** at `optimal_inaccurate` — that is CLARABEL's own solve-to-solve
nondeterminism, not an encoding error. The real-vs-Hermitian objective agreement
(~7–8 digits) is at or above that self-noise floor.

The exactness is instead proven to **machine precision** by a direct
constraint-residual cross-check (`_herm_xcheck.py`), bypassing the IPM entirely:
at the actual row4 optimum, the Bochner block's **min-eigenvalue is identical to
~1e-16** between the real-form `[[Re,−Im],[Im,Re]]` and the complex Hermitian
`(n+1)×(n+1)` matrix, for **both** `f≥0` and `1−f≥0`, at `bochner_n=10/20/30`; the
full spectra match to ~1e-15 (each real-form eigenvalue is a complex eigenvalue
**doubled** — the U(1)/SO(2) fact). The `bochner_hermitian.py` self-test gives
numeric Hermitian↔real-embedding agreement **1.1e-15** over 20 random trials.

| check | result |
|---|---|
| numeric Hermitian↔real embedding (20 trials) | max err **1.1e-15** |
| row4 Bochner min-eig real-vs-herm, bn=10/20/30, both signs | identical to **~1e-16** |
| row4 full Bochner spectrum real-vs-herm | match to **~1e-15** |
| cross-feasibility (each optimum feasible in the other program) | holds |
| canonical CLARABEL PSD block side-lengths real-vs-herm | **identical** at every bn |
| objective real-vs-herm (independent solves) | ~7–8 digits (= CLARABEL noise floor) |

### Step 3 — win measurement: **no win (~1×), Hermitian is marginally heavier**

**Decisive structural fact:** `prob.get_problem_data(cp.CLARABEL)` shows the PSD
cone CLARABEL actually receives is the **same size** for both forms — real and
Hermitian both give PSD side-lengths `[2(n+1), 2(n+1)]` (bn=10→[22,22],
20→[42,42], 30→[62,62]). The Hermitian form is in fact **strictly larger overall**:
it adds `(n+1)²`-ish equality constraints for the Hermitian ties (`zero` cone =
265 / 925 / 1985 at bn=10/20/30 vs 1 for the real form), plus more columns/nnz.

cvxpy's `Complex2Real` reduction is present in the chain for **all three installed
solvers** (CLARABEL, SCS, MOSEK); none receives a native complex/Hermitian PSD
cone. CLARABEL 0.11.1's cone inventory is `{Zero, Nonnegative, SecondOrder,
PSDTriangle, Exp, Pow, GenPow}` — **no complex PSD cone exists**.

Measured wall-time / peak-RSS at **N=5000, T=2000, R=10** (child-process
`resource.getrusage`):

| bochner_n | PSD blocks (both) | real wall | herm wall | time real/herm | real peak | herm peak | mem real/herm |
|---|---|---|---|---|---|---|---|
| 20 | [42,42] | 28.5 s | 28.2 s | 1.01× | 818 MB | 842 MB | 0.97× |
| 30 | [62,62] | 39.1 s | _(running)_ | — | 1386 MB | _(running)_ | — |

(`time real/herm ≈ 1` ⇒ no speedup; `mem real/herm < 1` ⇒ Hermitian uses *more*
memory. Target was ~4×.) Full bn=20..60 sweep continues in
`sym_reduction_result.json`.

**Max `bochner_n` reached by each form:** identical — both forms hit the same
memory wall at the same bn, because they produce the same-size real PSD cone
(plus the Hermitian form's extra equalities). The complex form does **not** reach
higher.

### Step 4 — bound at high bochner_n

Moot as a *win* (Step 3): the complex form does not unlock a higher `bochner_n`
than the real form. The row4 bound at the bn each reaches (identical for both
forms, since they are the same program) is recorded in the sweep `value` column
(e.g. N=5000: bn=20 → 0.379292, bn=30 → 0.379547), consistent with the documented
approach to the framework ceiling **C_∞ ≈ 0.380558** as bn rises — NOT exceeding
it. No bound gain attributable to the swap (it is the same optimum).

## Honest conclusion

- **Equivalence:** EXACT — Bochner min-eig matches to ~1e-16; the swap is
  mathematically correct, just (intended to be) smaller.
- **Memory/time factor:** ~**1.0×** (no win); the Hermitian form is marginally
  *heavier* (extra tie-equalities). Target ~4× **not achieved via cvxpy**.
- **Max bochner_n:** identical for both forms.
- **Bound gain:** none (same optimum); does **not** pass C_∞ ≈ 0.380558, exactly
  as flagged — this is a cost lever, not a value lever.

**Why the "free ~4×" did not appear, and how to actually get it.** The U(1)
isotypic reduction is real (every real-form eigenvalue is doubled), but
**CLARABEL/SCS/MOSEK via cvxpy cannot consume a complex/Hermitian PSD cone** —
cvxpy lowers it straight back to `[[Re,−Im],[Im,Re]]`. To realize the win one must
either:
1. feed a complex-PSD-native solver **outside** cvxpy (bypassing `Complex2Real`);
   note SDPA-GMP is real-only too, so this needs a different backend; **or**
2. implement the **real Z/2 centrosymmetric block-diagonalization** — the
   Hermitian-Toeplitz real-form is centrosymmetric (`S·RF·Sᵀ=RF`,
   `S=[[J,0],[0,−J]]`), so the orthogonal `Q=(1/√2)[[I,J],[J,−I]]` splits it into
   **two half-size REAL PSD blocks** that CLARABEL *can* consume as two separate
   PSD constraints — a genuine ~2× side → ~8× flops. This is the in-cvxpy path to
   a real speedup and the recommended next step (the "~1–2 week" deliverable in
   THREE_APPROACHES.md §3). The exact-equivalence encoding here can be reused for
   it with confidence.

   **The centrosymmetric split is numerically validated (this session):** with
   `I,J` of half size and `Q=(1/√2)[[I,J],[J,−I]]`, `Q·RF·Qᵀ` is block-diagonal
   (off-block coupling **1e-15**) and its two half-size diagonal blocks reproduce
   the full real-form spectrum to **1.3e-15**. So re-expressing the Bochner block
   as those two `cp.Variable((n+1,n+1), PSD=True)` constraints is exact and
   CLARABEL-consumable — the concrete recipe for the real ~2×.

This sharpens the prior assessment
(`Representation_theoretic_symmetry_reduct.md`), which warned CLARABEL "won't
consume block structure without re-expressing sub-blocks as separate PSD cons":
the complex route gives **zero** through cvxpy, so the **real centrosymmetric
split is the only in-cvxpy path** to the documented cost win.

## Files

- `lp_research_state/code/bochner_hermitian.py` — complex-Hermitian Bochner
  constraint + numeric cross-check helper (self-test: 1.1e-15).
- `lp_research_state/code/white_full_convex_hermitian.py` — `build_problem_hermitian`
  (delegates to verified `build_problem(bochner_n=0)`, bolts on Hermitian block).
- `lp_research_state/code/_herm_equiv_check.py` — equivalence gate.
- `lp_research_state/code/_herm_xcheck.py` — machine-precision constraint-residual
  proof (the airtight equivalence evidence).
- `lp_research_state/code/_herm_win_measure.py` — time/memory/PSD-size sweep.
