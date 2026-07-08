# Deep-Dive Assessment — Representation-theoretic symmetry reduction of White's SDP

**Date:** 2026-06-06
**Status:** PROPOSAL, adversarially vetted. No results claimed.
**Verdict:** **PROMISING** (as an *enabler*), with one headline correction: it **cannot beat** the
`C_∞ ≈ 0.380558` two-point ceiling — it can at best *approach* it by making the F3-projected
`bochner_n ≈ 60–80` regime tractable. Roughly half of the realized speedup is obtainable by a
trivial complex-SDP swap that needs no representation theory at all.

Read alongside: `docs/NEW_APPROACHES/NA_optimization_hierarchies.md` (Approach 2 — the source
pitch), `docs/archive/LEVER_F3_FULL_SATURATION.md` (the `+4–6e-4 → ~0.38058` projection this
would realize), `docs/archive/LEVER_SUPT_DIRECT.md` / PRO-22 (load-bearing cell-envelope = the
ceiling's origin), `lp_research_state/code/white_full_convex.py` (the program),
`lp_research_state/code/bochner.py` (the dominant PSD block), and
`lp_research_state/code/symmetric_push.py` (the repo's only — *conditional* — symmetry use).

---

## 0. One-paragraph bottom line

The technique is real, textbook-standard, and correctly cited. The decisive question is **how large
the actual symmetry group is on the cost-dominant block**, and I answered it numerically: the
Bochner constraint — the only PSD block whose size grows with the level you want to push — carries a
`U(1)` (complex-Hermitian-embedding) symmetry worth a **free ~4× on PSD flops**, plus a `Z/2`
centro/persymmetry worth **another ~2×**. That is ~4–8× on the PSD work, *exactly* the honest worst
case the source memo flagged ("if the only symmetry is `Z/2` the split is ~2×"). Crucially: (a) the
`U(1)` half of that is already available by swapping cvxpy's real-embedding for a complex Hermitian
variable — **one line, zero representation theory**; (b) the `nonneg ≈ N` cell-envelope cone, which
does *not* shrink under this group, becomes the new floor; (c) even a perfect reduction realizes the
F3 projection (`µ ≳ 0.38058`) but **stops below `C_∞ ≈ 0.380558`** because that ceiling is a property
of the degree-2 / two-point moment cone, not of solver cost (PRO-22). So: a genuine, near-term
engineering win that **buys headroom for other levers**, mis-pitched as a way to *break* the
saturation.

---

## 1. FEASIBILITY with available tools — HIGH

Pure linear algebra + convex modeling; every dependency is already in the repo.

| Need | Available? | Where |
|---|---|---|
| Form group action, averaging projector `P_G`, commutant | yes | numpy in `.venv` (cvxpy 1.8.2) |
| Block-diagonalize a matrix *-algebra numerically | yes | implement de Klerk's routine (Math. Program. 2011) in ~150 lines, or `mpmath`/`sympy` for exact small cases |
| Re-solve reduced PSD blocks | yes | cvxpy/CLARABEL; **complex** variables already supported (`cp.Variable(..., complex=True)`, `H >> 0`) |
| 10-digit reduced-vs-unreduced cross-check | yes | repo's standard discipline; `dual_extractor.py` gives `rigorous_dual_LB` either way |
| GMP-precision spot check at small N | yes | `lp_research_state/bin/sdpa_gmp` |

**Cost estimate.** First decisive experiment: **~1 hour** (most of which I already executed below —
see §3). A *working* reduced solve that actually consumes block structure to gain wall-clock at
`bn ≈ 60–80`: **~1–2 weeks**, because CLARABEL does not accept a user-supplied block-diagonal
decomposition — you must either (i) feed each isotypic sub-block as a separate, smaller cvxpy PSD
constraint (straightforward for the `U(1)` and `Z/2` factors; this is the real deliverable), or
(ii) export to a solver that natively exploits low-rank/clustered structure
(SDPA, or `ClusteredLowRankSolver.jl`, which is purpose-built for exactly these Fourier SDPs — see
sphere-packing docs). The svec-layout / scaled-PSD bookkeeping against cvxpy's canonicalization is
finicky (same class of work as the rigor track's PSD-unpack), so budget independent re-derivation
per the `_independent` convention.

**Verdict:** feasible at low–moderate effort. No new theory, no new external dependency required for
the `U(1)+Z/2` realization.

---

## 2. PRIOR ART — technique is standard; this *application* is new

**Who invented / uses it (citations verified live):**
- **de Klerk, Pasechnik, Schrijver**, *Reduction of symmetric SDPs using the regular
  ∗-representation*, Math. Program. 109 (2007) 613–624. The commutant/centralizer reduction; matrix
  order drops to **the number of orbits of the group action** (`homepages.cwi.nl/~lex/files/symm.pdf`).
  This orbit-count fact is the single most important number for this proposal (see §4).
- **Bachoc, Vallentin**, *New upper bounds for kissing numbers from SDP*, JAMS (2008); **Vallentin**,
  *Symmetry in semidefinite programs*, arXiv:0706.4233 — the canonical recipe for
  Cohn–Elkies-type **Fourier** SDPs, the closest structural analogue to White's program.
- **Gatermann, Parrilo**, *Symmetry groups, SDPs, and SOS*, J. Pure Appl. Algebra (2004);
  **de Klerk**, *Numerical block diagonalization of matrix ∗-algebras*, Math. Program. (2011) — the
  algorithmic block-diagonalizer.
- **Invariant SDPs** survey: arXiv:1007.2905. **Complex→real reductions:** arXiv:2307.11599
  ("a more efficient reformulation of complex SDP as real SDP"), PICOS complex-SDP docs — these
  document the *free* `U(1)` factor (§3.1).

**Has any piece been tried here?**
- **Only `lp_research_state/code/symmetric_push.py`** — and it is a *different, weaker* thing: it
  *assumes* `f` even (`d_k = 0, dlt = 0, v_j = w_j`), which (i) is **conditional** on the open
  even-`f` conjecture and is therefore *not an unconditional bound on `µ`*, and (ii) is feasible only
  for rows 5 and 6 (`h = 0`); rows with `h > 0` contradict even-`f`. Its header says exactly this.
  That is a *manual variable drop on the optimizer*, not an *isotypic decomposition exploiting
  constraint symmetry*. The proposed reduction is **unconditional** (it uses the symmetry of the
  *constraint kernels*, true for all 7 rows regardless of what the optimizer does).
- A grep for `representation / regular.representation / block.diag / isotypic / irrep / commutant`
  across `lp_research_state/code/` returns nothing. Confirmed: the systematic reduction is unbuilt.
- **Application novelty:** the Erdős minimum-overlap / autoconvolution literature (White
  arXiv:2201.05704; Matolcsi–Vinuesa; the 2025 autoconvolution note arXiv:2508.02803) uses Fourier +
  convex optimization but **no symmetry block-diagonalization**. So priority on *this application* is
  clean, even though the *method* is decades old.

**Verdict:** mature, well-sourced technique; genuinely new on this problem; the repo's lone prior use
is conditional and orthogonal.

---

## 3. FIRST EXPERIMENT — RAN IT (cheap, analytic/numeric, no heavy solve)

The cheapest decisive signal is: *what symmetry actually survives on the cost-dominant PSD block, and
how big is the implied reduction?* I computed this directly (numpy, `<1 s`).

### 3.1. The Bochner real-form is automatically 2× redundant (free `U(1)`)

`bochner.py` encodes the complex Hermitian moment matrix `M_n(f)` (size `n+1`) as the real-symmetric
`real_form = [[Re, -Im], [Im, Re]]` of size `2(n+1)`. For a **generic** Hermitian Toeplitz `M`
(`n=6`), the `14×14` real form has **every eigenvalue of multiplicity exactly 2**:

```
eig(real_form) = [-2.149(×2), 0.046(×2), 0.114(×2), 0.742(×2), 0.848(×2), 1.381(×2), 2.519(×2)]
```

This is the standard fact (arXiv:2307.11599; PICOS complex SDP): the real embedding of an `m×m`
Hermitian PSD constraint *always* block-diagonalizes back into **two identical `m×m` complex blocks**
under the `SO(2)`/`U(1)` action `[[0,-I],[I,0]]`. Interior-point cost on a dense PSD block scales
≈ (side)³ in the Schur complement, so replacing one `2(n+1)` real block by one `(n+1)` **complex**
block is `2·(1/2)³ = 1/4` the flops on that constraint → **a free ~4× on the Bochner PSD work**.

> **Actionable, today, no representation theory:** swap
> `cons.append(real_form >> 0)` for a complex Hermitian `cp.Variable((n+1,n+1), hermitian=True)`
> (or `complex=True` with `M == M.conj().T`) in `bochner.py` / `white_full_convex.py`'s Bochner
> section. Cross-check `rigorous_dual_LB` to 10 digits, then re-attempt the binding **row 4** at the
> first `bochner_n` that previously OOM'd. **This one change captures the larger of the two symmetry
> factors and is the recommended *actual* first experiment.**

### 3.2. The Hermitian Toeplitz block carries a real `Z/2` centro/persymmetry

For a Hermitian Toeplitz `M`, index-reversal `J = antidiag(1)` gives `J M J = Mᵀ = conj(M)`. I
verified on the real form that `S · real_form · Sᵀ = real_form` for `S = [[J,0],[0,-J]]` (reversal +
sign flip on the Im block) — i.e. a genuine `Z/2` symmetry **of the constraint itself** (not an
optimizer assumption). Centrosymmetric matrices split into a `⌈m/2⌉` + `⌊m/2⌋` pair of blocks → a
further **~2×**. This is the de Klerk–Pasechnik–Schrijver "order = #orbits" statement: `Z/2` on
`n+1` indices has `≈ (n+1)/2` orbits.

### 3.3. Where the cost actually lives (cone census, `N=300, T=200`)

`get_problem_data` (CLARABEL canonicalization), three Bochner levels:

| `bochner_n` | PSD block sides | PSD svec total | nonneg cone | soc cones |
|---|---|---|---|---|
| 12 | `[26, 26]` | 702 | **2529** | 40×`[3]` + 2×`[202]` |
| 30 | `[62, 62]` | 3906 | **2529** | 40×`[3]` + 2×`[202]` |
| 60 | `[122, 122]` | 15006 | **2529** | 40×`[3]` + 2×`[202]` |

Two readings, both load-bearing for the verdict:
- **The premise is correct:** the PSD blocks are exactly the `2(n+1)` real forms, they scale `O(n²)`
  in svec (≈`O(n³)` in IPM flops), and they *do* dominate as you push the level. So reducing them is
  the right target.
- **The honest ceiling on the win:** the `nonneg = 2529` cone (cell-envelope, ≈`N`-sized) and the two
  `soc=[202]` blocks are **invariant under the Bochner group and do not shrink**. The reflection
  `x ↦ 2−x` *does* act on `(w, v)` (it is the `v ↔ w` swap), giving the cell cone its own `Z/2` — but
  that is again only ~2×, and at large `N` with modest `bochner_n` the cell cone is the wall, not the
  PSD block. The reduction helps precisely in the regime the proposal targets (large `bochner_n`),
  and helps little elsewhere.

**Finding from the first experiment:** the realizable symmetry is `U(1) × Z/2 ≈ 4–8×` **on the PSD
blocks only**; ~4× of it is a free complex-SDP swap. This is concrete, decisive, and bounds the
upside before any heavy solve.

---

## 4. WHY IT COULD (PARTIALLY) BEAT THE SATURATION — and the sharp correction

**The honest positive case.** `LEVER_F3_FULL_SATURATION.md` §5.3 projects that pushing
`bochner_n ≥ 40–50` plus a tighter cell-envelope yields `+4–6×10⁻⁴`, *plausibly* `µ ≥ 0.38058`, and
flags it as "research, not engineering" because the scale is intractable in 4 GB. Symmetry reduction
**does not change the relaxation** — it changes the cost of solving it. A 4–8× cut on the PSD blocks
turns `bochner_n = 60–80` at the binding **row 4** from OOM into a run, directly realizing the F3
projection. Because it multiplies the reach of *every already-validated rigorous lever* (Bochner-PSD,
ellipse extension), it is the highest-leverage *enabler* on the board, and its first experiment is the
cheapest decisive one. That is the real, defensible value.

**The correction the source pitch understates.** The workflow goal names a target: "push the LOWER
bound past the conjectured `~0.380558` SDP-framework ceiling." **This approach cannot do that.**
PRO-22 (`LEVER_SUPT_DIRECT.md`) and `NA_optimization_hierarchies.md` establish that `C_∞ ≈ 0.380558`
is a property of the **degree-2 / two-point moment cone** (the cell-envelope Parseval link
constrains only `|f̂(m)|²` of a single density), not a solver-cost artifact. Bochner-PSD at
`bochner_n → ∞` is *still inside that cone*. So symmetry reduction can take `µ_LB` *up to ~the
ceiling faster*, but **the limit it approaches is `C_∞`, not beyond it** — the source memo's own §5.2
and its risk #2 say exactly this ("even a perfect 10× still only buys `µ ≈ 0.38058` … it cannot pass
the `C_∞` 2-point ceiling, only approach it"). Transcending the ceiling requires a *different cone*
(higher-arity moments — Approach 1 in the same memo), which symmetry reduction merely *enables at
scale*, not replaces.

**Net:** "beat the saturation" is **false as literally stated**; "realize the F3 `+4–6e-4` and *reach*
the ceiling cheaply, then enable the cone-changing lever that *can* transcend it" is **true and
valuable**. Verdict tracks the true statement.

---

## 5. RISKS / why it might fail (adversarial)

1. **Small group ⇒ small win (the dominant risk, now quantified).** The realized symmetry is
   `U(1) × Z/2 ≈ 4–8×` on the PSD blocks (§3), **not** the order-of-magnitude the "large-`T_max`
   barrier" rhetorically needs. There is no hidden dihedral/translation group of large order: the
   Toeplitz structure gives persymmetry (`Z/2`), the complex embedding gives `U(1)`, and that is the
   inventory. A 4–8× is a real but *modest* enabler.

2. **Half the win is free and needs none of this machinery.** The `U(1)` factor (~4× of the ~4–8×) is
   a one-line complex-SDP swap (§3.1). If the *only* delivered gain is that swap, the
   "representation-theoretic / regular ∗-representation" framing is overkill — the honest deliverable
   is "use complex Hermitian Bochner blocks + exploit centrosymmetry," not a commutant computation.
   (This is a *scoping* risk, not a correctness one: the win is real, the billing is inflated.)

3. **The cell-envelope cone is the floor and is barely touched.** The `nonneg ≈ N` and `soc=[202]`
   cones (§3.3) carry only the `v ↔ w` reflection `Z/2`. At large `N` / modest `bochner_n` they
   dominate and the PSD reduction buys little wall-clock there. The reduction pays off *only* in the
   high-`bochner_n` regime — fortunately the regime the proposal wants, but the gain is not uniform.

4. **CLARABEL won't consume block structure for free.** It is a general IPM; to *realize* the speedup
   you must either re-express each isotypic sub-block as a separate smaller cvxpy PSD constraint
   (doable for `U(1)+Z/2`; this is the actual work) or export to a structure-aware solver
   (SDPA-GMP exists in-repo for precision but not for block-exploitation; `ClusteredLowRankSolver.jl`
   would be a new external dependency). Without that plumbing, forming the commutant changes nothing
   the solver sees.

5. **Bookkeeping/rigor surface.** svec scaling, the `[[Re,-Im],[Im,Re]]` ↔ complex map, and matching
   cvxpy canonicalization are exactly the kind of off-by-`√2`/sign traps that silently corrupt a
   `rigorous_dual_LB`. Mitigate with the repo's `_independent` re-derivation + 10-digit cross-check +
   one SDPA-GMP spot check at small `N`. (Standard, but non-zero.)

6. **Ceiling, restated as a risk:** even with a perfect reduction and `bochner_n = 200`, the bound
   asymptotes to `C_∞ ≈ 0.380558 < µ_UB = 0.380871`. The gap does **not** close from this lever alone.
   A successful run that lands at, say, `µ_LB ≈ 0.38056` would be a real `+~3e-4` headline improvement
   over the current `0.380284` **and still leave a `~3e-4` gap** — a strong outcome, not a closure.

---

## 6. Recommendation

**Do the §3.1 experiment first (complex Hermitian Bochner swap), because I have already shown it
captures the larger symmetry factor for one line of code.** Then add the centrosymmetry split. If the
combined wall-clock cut lets **row 4** reach `bochner_n ≈ 60–80` within memory, read the new
`rigorous_dual_LB`; the realistic prize is the F3-projected `µ_LB ≈ 0.3805–0.38058`. Treat this as the
**enabler for `NA_optimization_hierarchies.md` Approach 1** (the 3-point measure hierarchy that can
actually transcend `C_∞`), per that memo's own "Combination worth flagging: 2 ⊕ 1." Do **not** market
it as breaking the saturation — market it as reaching the ceiling cheaply and unlocking the lever that
breaks it.

---

## Appendix — exact reproduction of the first experiment

```python
# Bochner real-form is 2x redundant (free U(1)); Hermitian Toeplitz has a Z/2 centro-symmetry.
import numpy as np
n = 6; rng = np.random.default_rng(0)
c = rng.standard_normal(n); d = rng.standard_normal(n)
Re = np.eye(n+1)*0.5; Im = np.zeros((n+1, n+1))
for j in range(n+1):
    for k in range(n+1):
        ell = j-k
        if ell:
            a = abs(ell); Re[j, k] = c[a-1]/2
            Im[j, k] = (-d[a-1]/2) if ell > 0 else (d[a-1]/2)
RF = np.block([[Re, -Im], [Im, Re]])
print(np.round(np.linalg.eigvalsh(RF), 3))          # every eigenvalue doubled -> free 2x (4x flops)
J = np.fliplr(np.eye(n+1))
S = np.block([[J, 0*J], [0*J, -J]])
print(np.allclose(S @ RF @ S.T, RF))                # True -> constraint-level Z/2 (another ~2x)
```

```python
# Where the cost lives (run from lp_research_state/code with ../../.venv/bin/python):
import cvxpy as cp
from white_full_convex import build_problem
for bn in (12, 30, 60):
    O, w, v, c, d, e, dl, cons = build_problem(300, 200, 10, 0.004, 0.004,
                                               0.3875, 0.3875, -0.02, 0.02, bochner_n=bn)
    data, _, _ = cp.Problem(cp.Minimize(O), cons).get_problem_data(solver=cp.CLARABEL)
    print(bn, 'PSD sides', list(data['dims'].psd), 'nonneg', data['dims'].nonneg)
# -> PSD sides [2(bn+1)]*2 grow O(bn^2); nonneg stays 2529 (cell-envelope floor, untouched by Bochner group)
```
