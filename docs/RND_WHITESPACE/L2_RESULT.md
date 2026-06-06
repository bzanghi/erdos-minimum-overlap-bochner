# L2 — Jansson rigorous a-posteriori bound: RESULT

**Date:** 2026-06-06
**Author:** Claude (machine-assisted), L2 thrust
**Script:** [`lp_research_state/code/_jansson_verify.py`](../../lp_research_state/code/_jansson_verify.py)
 (+ driver [`_jansson_ladder.py`](../../lp_research_state/code/_jansson_ladder.py))
**Machine results:** [`L2_RESULT.json`](L2_RESULT.json)

---

## Verdict: **A — VERIFIED** (at the tested binding centers, N ≤ 3000)

The Jansson–Chaykin–Keil (2007) Thm 3.2 rigorous a-posteriori lower bound,
implemented end-to-end in directed-rounding interval arithmetic (`mpmath.iv`) on
the *exact* CLARABEL canonical conic data, **certifies the augmented-SDP optimum
at the binding centers and the rigorous bound `p_lo` stays above White's
0.379005** — the rigorous correction is negligible. The "critical risk" the
charter flagged (a slightly-negative interval λ_min on the Bochner PSD block
eating the margin) was **resolved**: it was an artifact of *unpivoted* interval
Cholesky, not a real infeasibility. With symmetric-diagonal-pivoted interval
LDLᵀ, the verified λ_min of the Bochner dual block is **positive** (the block is
rigorously inside the PSD cone), so its penalty contribution is exactly zero.

> Scope (honest): this certifies **`SDP_opt(center) ≥ p_lo`** at the two tested
> centers — NOT `µ ≥ p_lo`. The latter additionally needs (a) the constraint
> data to be a valid relaxation of the overlap problem (it is, by White's
> derivation), and (b) the ellipse/cover lift from the 7 centers to the full
> (h,p,q) region (`path_b_*`). Those are out of scope here; flagged below.

---

## What was built

`_jansson_verify.py` — a reusable, unit-tested verifier. For an augmented
`build_problem(...)` instance at a center it:

1. Pulls the **exact solver-form data** the way the solver consumed it:
   `prob.get_problem_data(cp.CLARABEL)` → `(A, b, c, dims)`, then
   `chain.solve_via_data(...)` → raw `(x, z, s, obj_val, obj_val_dual, r_prim,
   r_dual)`. Self-checks `c@x == obj_val` (= Ω*) and `−b@z == obj_val_dual` to
   machine precision (verifies the canonical objective IS the modeled Ω).
2. Forms the **dual defect** `D = c + Aᵀz` (CLARABEL convention: dual-feasible
   ⇔ `D = 0`; `z ∈ K*`).
3. **Rigorous per-cone lower bound** `d_j` on the cone-λ_min of `z` (PSD blocks:
   verified symmetric-pivoted interval LDLᵀ; nonneg: min coord; SOC: `t−‖x‖` with
   `‖x‖` interval-upper-bounded). Unit-tested to enclose-from-below numpy
   `eigvalsh` on 60 random matrices (worst gap 1e-13).
4. **Finite primal box bounds** `x̄`: model columns mapped to their exact
   canonical columns via cvxpy's authoritative `var_offsets` (NOT a guessed
   id-sort — the solver interleaves `eps/dlt` between `v` and `c`; a built-in
   assertion `|x_i| ≤ x̄_i` on mapped columns guards the mapping). Auxiliary lift
   columns get a **provable structural cap** (every lift is a squared model
   expression: `sum_squares(c,d) ≤ 0.5`, `a_m² ≤` an interval-computed constant),
   optionally tightened by interval bound-propagation.
5. **Assembles `p_lo = −bᵀz + Σⱼ min(0,dⱼ)·s̄ⱼ − Σᵢ |Dᵢ|·x̄ᵢ`** entirely in
   interval arithmetic and returns the lower endpoint.

### The inequality (proved in the script docstring)
For any primal-feasible `x` (so `s = b − Ax ∈ K`):
`cᵀx = −bᵀz + zᵀs + Dᵀx`, then `zᵀs ≥ Σ min(0,λ_min^{K*}(z_j))·trace(s_j)`
and `Dᵀx ≥ −Σ|Dᵢ|x̄ᵢ`. Taking `inf` over feasible `x` gives
`optimum ≥ p_lo`. All terms evaluated with directed rounding.

---

## The decisive technical finding: pivoting fixes the Bochner λ_min

The f≥0 Bochner moment matrix is **near-singular at the optimum** (complementary
slackness — the optimal `f` makes `M_n(f)` rank-deficient; condition number ~1e9
at N=3000, bn=16). Plain (unpivoted) interval Cholesky of the dual block fails at
the small pivot — accumulated interval width (~2.5e-9) swamps it — and spuriously
reports λ_min ≈ −3.3e-6, which (× trace 17) would inflate the penalty to ~−5.6e-5
and *appear* to eat much of the margin. **Increasing arithmetic precision (mpmath
dps 30→120) does NOT help** — it's an instability of the unpivoted factorization,
not a precision shortfall.

**Symmetric diagonal pivoting** (move the largest remaining diagonal to the pivot;
a symmetric permutation is congruent, preserving inertia) makes the interval LDLᵀ
**tight**: it certifies λ_min ≥ 0.999 × numpy's λ_min = **+3.7e-12 > 0**. So the
Bochner dual block is rigorously PSD and contributes **zero** penalty. numpy
`eigvalsh` independently confirms λ_min ≈ +4e-12 (the block genuinely IS PSD).

**Robustness — the `zᵀs` term is exactly 0, so the result is parameter-free.**
At every N, ALL cone blocks of the dual `z` are verified to have **positive**
cone-λ_min (e.g. N=2000: 1 nonneg + 42 SOC + 2 PSD blocks, min `d_lower` = 3e-12
> 0). So `z ∈ K*` rigorously, `min(0,dⱼ) = 0` for every block, and the entire
`zᵀs` penalty vanishes **independent of the slack bound `s̄ⱼ` / `slack_infl`**.
The only surviving correction is the dual-stationarity term `Dᵀx` (~1e-7), which
is intrinsic to the defect `‖c+Aᵀz‖` and the O(1) primal box — it does not depend
on any tunable. `D` itself is enclosed in interval arithmetic (sparse `Aᵀz`
matvec in `mpmath.iv`), so even the matvec rounding is bounded.

> Per-center caveat: that `p_lo > 0.379005` holds at these **2 of 7** centers at
> N=3000 is NOT yet `µ ≥ 0.379005` — that needs all 7 centers AND the cover step.
> What is established is the *mechanism*: the rigorous interval certificate
> reproduces the per-center SDP value to ~2e-7, so it will carry whatever bound
> the production-N solves produce, with a provably negligible correction.

---

## Per-N results (slack_infl = 1; pm_k_max = 14; all PSD blocks verified PSD)

**Binding center row4 = (h=0.004, p=0.3875, q∈[−0.02,0.02]):**

| N | bn | prob.value (Ω*) | p_lo (rigorous) | penalty | p_lo − 0.379005 (White) | p_lo − 0.380284 |
|---|---|---|---|---|---|---|
| 300  | 6  | 0.36954989 | 0.36954989 | −8.1e-10 | −9.46e-3 | −1.07e-2 |
| 1000 | 10 | 0.37662382 | 0.37662351 | −2.6e-7  | −2.38e-3 | −3.66e-3 |
| 2000 | 12 | 0.37848719 | 0.37848690 | −2.4e-7  | −5.18e-4 | −1.80e-3 |
| 3000 | 16 | 0.37916388 | **0.37916366** | −1.9e-7 | **+1.59e-4** | −1.12e-3 |

**Headline binding center cde_n30_iter3 = (h≈4.5e-5, p=0.39015, q∈[−0.02,0.02]):**

| N | bn | prob.value (Ω*) | p_lo (rigorous) | penalty | p_lo − 0.379005 (White) | p_lo − 0.380284 |
|---|---|---|---|---|---|---|
| 300  | 6  | 0.36975245 | 0.36975244 | −2.5e-9 | −9.25e-3 | −1.05e-2 |
| 1000 | 10 | 0.37672719 | 0.37672705 | −1.2e-7 | −2.28e-3 | −3.56e-3 |
| 2000 | 12 | 0.37852748 | 0.37852722 | −2.1e-7 | −4.78e-4 | −1.76e-3 |
| 3000 | 16 | 0.37917988 | **0.37917936** | −4.5e-7 | **+1.74e-4** | −1.10e-3 |

**Reading the table.**
- `p_lo` tracks Ω* to within the penalty (~2–4.5e-7): the verified interval
  certificate essentially **equals** the SDP value — the rigorous correction
  costs nothing.
- The penalty at every N is dominated by the tiny `Dᵀx` term (~1e-7…1e-9); the
  Bochner `zᵀs` term is **exactly 0** because the verified (pivoted) λ_min of
  both Bochner PSD dual blocks is **positive** at every N (column "PSDlmin>0?" =
  True for all 8 runs in the JSON).
- Both centers clear White only at **N=3000** (+1.6e-4 / +1.7e-4). At N≤2000 the
  *discretization* (coarse N underestimates the true SDP value) keeps Ω* — and
  hence `p_lo` — below White. This is expected: the published bound uses
  N=10000–24000. **The margin lost at small N is discretization, NOT the Jansson
  penalty** (which never exceeds 5e-7). So at production N the verified `p_lo`
  will sit ~5e-7 below the (larger) Ω* and comfortably above White.
- Margin to the 0.380284 headline is still negative at N≤3000 for the same
  discretization reason; resolving it needs production-N solves (separate step).

---

## What this unlocks / next step

- The Jansson mechanism **works** for this program: the rigor upgrade
  ("numerically-certified log-parse" → "interval-arithmetic certificate") costs
  essentially nothing in bound strength, *provided the PSD λ_min is computed with
  a pivoted verified factorization*. The `dual_extractor.py` log-parse can now be
  replaced by (or cross-checked against) this certificate.
- **Next step (to make it a bound on µ, not just on the center SDP):** feed each
  center's verified `p_lo` through the existing rigorous ellipse/cover machinery
  (`path_b_with_polymoment.py` / `path_b_rigorous.py`), replacing the
  `value − margin` per-center input with the interval-certified `p_lo`. That
  lifts the per-center certificates to a verified full-region `µ ≥ …`. Production
  N (10000–24000) solves are a separate step (this charter capped N ≤ 3000); the
  verifier is N-agnostic and runs on whatever solve is produced.
- This is the verified front-end for the L1 (rational-certificate) finish.
