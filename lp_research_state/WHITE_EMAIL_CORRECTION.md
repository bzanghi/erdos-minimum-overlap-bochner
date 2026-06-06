# White email correction (2026-05-31): constraints 5.6/5.7 coefficient 8 → 4

**Status: RESOLVED — correctness/provenance fix applied; headline NOT changed.**
**Verdict: NEUTRAL at the binding center** (the constraint is effectively slack where it
governs the published bound). Where it *is* binding (row7 center, off the binding row) the
old coefficient 8 was *conservative*, not an overclaim. No center shows a material negative
delta, so the published core-region headline stands.

---

## 1. White's email (the trigger)

E. P. White (author of the program we augment, White 2023, Acta Arith.,
[arXiv:2201.05704](https://arxiv.org/abs/2201.05704)) wrote, after reviewing our
Bochner-PSD strengthening:

1. **Validation.** He confirmed the Bochner moment-matrix PSD constraint
   `M_n(f) ⪰ 0` / `M_n(1−f) ⪰ 0` is a legitimate, rigorous strengthening of his §5
   convex program.
2. **Correction A — constraints 5.6/5.7.** The RHS of constraints 5.6/5.7
   "have an 8 in the numerator on the RHS, it should be a 4."
3. **Correction B — constraints 5.8/5.9.** These "should use `2m−1`" (not `m`) on the RHS.
4. **Materiality.** "My calculations show that this does not make a material difference to
   the bound."

---

## 2. Our code status (before vs. after)

### Correction A (5.6/5.7 — the imaginary/sine cell-consistency constraint)

The canonical site is **`lp_research_state/code/white_full_convex.py:188`** (now line ~192
after the edit), inside `build_problem(...)`:

```python
# BEFORE (buggy):
rhs = -(8.0 / (m * np.pi)) * sin_pi_half_m * bm
# AFTER (corrected, parametrized):
rhs = -(mside_sin_coeff / (m * np.pi)) * sin_pi_half_m * bm   # mside_sin_coeff default 4.0
```

The **real-part sibling** of this constraint (line 184/185,
`rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am`) already correctly used **4**. The 5.6/5.7
constraint is exactly the imaginary half of the same complex identity (Lemma 2:
`M̂(m) = a_m f̂(m) − 4|f̂(m)|²`, with `a_m = (4/(mπ)) sin(mπ/2)` and
`f̂(m) = (c_m − i d_m)/2`), so by cos/sin symmetry the multiplier must be **4 on both**. The
hardcoded **8** doubled the linear-in-data term and is not derivable from `M(x) ≥ 0`. See §3.

**All hardcoded-8 occurrences fixed (7 active sites + 1 transitive inheritor):**

| # | File | Function | Fix |
|---|------|----------|-----|
| 1 | `white_full_convex.py:188` | `build_problem` | **Parametrized** `mside_sin_coeff` (default 4.0); pass `8.0` to reproduce old behavior |
| 2 | `white_full_convex_exact.py:161` | `build_problem_exact` | literal `8.0 → 4.0` + citation comment |
| 3 | `path_b_independent.py:144` | `solve_one_center` | literal `8.0 → 4.0` + citation comment |
| 4 | `path_b_analytical.py:105` | `build_problem_with_dual_handles` | literal `8.0 → 4.0` + citation comment |
| 5 | `path_b_lasserre.py:97` | `build_problem_with_dual_handles_BL` | literal `8.0 → 4.0` (withdrawn-Lasserre branch) |
| 6 | `symmetric_push.py:117` | `build_problem_even_with_handles` | literal `8.0 → 4.0` + citation comment |
| 7 | `_run_lasserre3_test.py:80` | `build_program` | literal `8.0 → 4.0` (throwaway scratchpad) |
| — | `path_b_rigorous.py` | (imports #4) | **transitive** — inherits fix from `path_b_analytical.py`; no own copy |

Verified clean (no buggy-8, no 5.6/5.7 sine constraint at all): `bochner.py`,
`bochner_independent.py`, `poly_moment.py`, `white_full_convex_supt.py`,
`ai_discovery/dsl.py`. A repo-wide grep for `8.0 / (m * np.pi)` now returns **zero** matches.
`py_compile` passes on all edited files; `import white_full_convex` + a smoke
`build_problem(...)` both succeed (project venv `/Users/benzanghi/Documents/Claude/Projects/Erdos/.venv`).

### Correction B (5.8/5.9 — the eps/dlt tail-bound constraints) — ALREADY CORRECT

Our code **already uses `2m−1`**, matching White's note. Verified in
`white_full_convex.py:196-200`: the tail-bound loop iterates `m` over `1..R`, sets
`m_odd = 2*m-1`, and passes `m_odd` (not `m`) into both `tail_bound_eps(m_odd, T)` and
`tail_bound_delta(m_odd, T)`; the helper bodies (lines 91-96) consume only `m_odd`. The LHS
slot `eps[m-1]`/`dlt[m-1]` holds the `(2(m-1)+1) = (2m-1)`-th odd Fourier coefficient, so LHS
and RHS indexing match. **No change needed.** (The plain-`m` denominators at lines 184/192 are
correct: those loops run over the *true* frequency `1..2R`, so `m` is the actual frequency, not
a half-index.)

---

## 3. Independent derivation (confirming the coefficient is 4)

From Lemma 2 as encoded at `white_full_convex.py:261`, `M̂(m) = a_m f̂(m) − 4|f̂(m)|²` with
`a_m = (4/(mπ)) sin(mπ/2)` and `f̂(m) = (c_m − i d_m)/2`. Splitting into real/imaginary parts:

- `Re M̂(m) = (2 sin(mπ/2)/(mπ)) c_m − (c_m² + d_m²)`
- `Im M̂(m) = −(2 sin(mπ/2)/(mπ)) d_m`

The linear normalization is **identical in magnitude** — `2 sin(mπ/2)/(mπ)` on `c_m` resp.
`d_m` — by cos/sin symmetry. With the code's scaling (`am = c/2`, `bm = d/2` on the even branch),
`(4 s/(mπ))·(c/2) = (2 s/(mπ)) c_m` reproduces `Re M̂` (line 184/185, coeff 4 — correct), and the
sin block needs the *same* coeff 4 to reproduce `|Im M̂|`. The literal `(8 s/(mπ))·(d/2) =
(4 s/(mπ)) d_m` is **exactly double** and is **not** derivable from `M(x) ≥ 0`. (sympy-verified.)

**Provenance.** `git log --follow white_full_convex.py` returns a *single* commit
(`5344383`, 2026-05-10, "Bochner-PSD strengthening … µ ≥ 0.379544"); `git blame` attributes the
whole rhs/cons block to that boundary commit. **No commit ever changed 4 → 8** and no comment
justifies an 8. It is an inadvertent transcription artifact baked into the first encoding of
§5 and copy-pasted into every independent re-derivation. `rationale_for_8_found = none`. The
correction confirms White's email exactly.

---

## 4. Rigor direction (why this is not an overclaim)

The 5.6/5.7 constraint is the two-sided band
`(L/2)(b₋@w − b₊@v) ≤ rhs ≤ (L/2)(b₊@w − b₋@v)` with `rhs = −(coeff/(mπ)) sin(mπ/2)·bm`.

- Coefficient **4 < 8** ⇒ `|rhs|` **smaller** ⇒ the two-sided band is **narrower** ⇒ the
  relaxation is **tighter but still valid** ⇒ the SDP minimum can only **rise or hold**.
- Therefore the old coefficient 8 was a **valid-but-looser (conservative)** constraint, **not**
  an overclaim. Switching to the correct 4 can only **improve or leave unchanged** the bound.

This is the opposite failure mode from the retired Lasserre / poly-moment tail-bound traps
(where an under-sized bound made a cut *too tight* and *inflated* the min). Here the corrected
value is the *tighter* one, so there is no risk of a retroactive reduction. The only question
the numerics decide is *how much* the bound improves — and at the binding center, the answer is
"immaterially," exactly matching White's remark.

---

## 5. Per-center 8-vs-4 measurement table

All solves use `dual_extractor.solve_with_dual_extraction` (CLARABEL), with identical solver
config per center for the two coefficients (so the *difference* is the only variable). These are
single-point row-center solves at a **light config (N≈2000, T≈800, bochner_n=20)** — the
*difference* is meaningful; the absolute values are NOT certified µ bounds nor production numbers.

`delta := lb_coeff4 − lb_coeff8` (rigorous_dual_LB); `Δ(value) := prob.value_4 − prob.value_8`
(higher resolution).

| Center | (h, p) | binding? | lb (coeff 8) | lb (coeff 4) | delta (LB) | Δ(value) | reading |
|--------|--------|----------|--------------|--------------|------------|----------|---------|
| **cde_n30_iter3** (governs headline) | (0.00392, 0.39225) | binding region | 0.37829 | 0.37829 | 0.0 † | **−2.74e-8** | **NEUTRAL / slack** |
| **row4** (binding row) | (0.004, 0.3875) | binding region | 0.3782 | 0.3782 | 0.0 † | **−3.63e-8** | **NEUTRAL / slack** |
| row7 | (0.03, 0.375) | NOT binding | 0.37924 | 0.38084 | **+1.60e-3** | +1.5969e-3 | constraint BINDING; 8 was conservative |

† `rigorous_dual_LB` is parsed from CLARABEL's iteration log, which prints `dcost` to ~5 sig
figs, so both coeffs round to the same value and `delta_LB = 0` carries **no information below
~1e-5**. The informative signal at the binding centers is full-precision `prob.value`.

**Interpretation of the binding centers.** At both binding-region centers the magnitude is
`|Δ| ≈ 2.7e-8 … 3.6e-8`, which is (a) far below the `1e-6` neutrality threshold and (b) below
CLARABEL's own last-gap floor (~9e-8; both solves `optimal_inaccurate`, dual residual ~1e-10).
The tiny **negative** sign is convergence-point noise at the residual scale, **not** a real
overclaim. (If it were a genuine overclaim the magnitude would have to exceed the solver's own
gap floor, and the sign would persist as the solve is tightened — neither holds.)

**Interpretation of row7.** `delta = +1.60e-3` is strongly positive, far above `+1e-6`, so the
sine constraint *is* binding at this off-binding-row center. Coeff 4 gives a **higher** bound,
confirming the rigor direction in §4 (8 was conservative). row7 is **not** the binding row for
the headline (row4 is), so this large local gain does not move the certified core bound —
consistent with White's "not material" remark.

---

## 6. Cross-check result (independent re-implementation agreement)

`path_b_independent.py` rebuilds the 5.6/5.7 constraint itself, inline in `solve_one_center()`
(line 145), and its constant has *already* been independently corrected to `4.0` with the
2026-05-31 citation. Running the canonical `build_problem` against the independent
`solve_one_center` in one process, **corrected coeff = 4**:

- Config A (N=600, T=240, R=8, bn=6): main = indep = `0.372255512283854`, abs_diff = `0.0`
  (bit-for-bit).
- Config B (N=1000, T=400, R=10, bn=8): main = indep = `0.375241004400694`, **abs_diff = 0.0**
  (agree to all 16 printed digits; ≥10 required). ✅
- **Genuineness control:** forcing the canonical path to `mside_sin_coeff = 8.0` yields
  `0.375241003792822`, which differs from the independent value by `6.08e-10` — proving the
  independent constant is truly **4**, not 8, and (again) that coeff 4 gives the **higher**
  objective (tighter-but-valid).

A final repo-wide grep across `path_b_*.py`, `white_full_convex*.py`, `symmetric_push.py`,
`_run_lasserre3_test.py` for any `8 / (m*np.pi)` sine site returns **zero** matches — no stale
uncorrected copy survives. (Caveat: `path_b_independent` shares `white_full_convex`'s cell-bound
*helpers* by import, so it is independent at the constraint-assembly/constant level, not a
from-scratch re-derivation of the cell bounds — but for the specific 5.6/5.7 coefficient under
test it is a genuinely separate, separately-edited copy, and it matches.)

---

## 7. Verdict

- **correct_coeff = 4** (confirmed by independent derivation, cross-check, and White's email).
- **verdict = NEUTRAL** at the binding center (`cde_n30_iter3`, near row4): the imaginary/sine
  cell-consistency constraint is effectively **slack** there, so the 8 → 4 correction does not
  materially change the bound. Where the constraint *is* binding (row7, off the binding row) the
  old 8 was **conservative** (coeff 4 improves the local bound by +1.6e-3) — never an overclaim.
- **max_abs_delta** at the headline-governing binding center ≈ **2.74e-8** (full-precision
  `prob.value`). The largest |delta| *anywhere measured* is row7's **+1.60e-3**, but that is in
  the **conservative (improving)** direction and at a **non-binding** center, so it does not
  reduce the published bound. **No center exhibits a material negative delta** ⇒ **no overclaim**.
- **material = false.** `|max delta|` at the binding center (~2.7e-8) is ~four orders of magnitude
  below the 6th decimal place of the published headline, so it cannot change the published 6-dp
  core-region bound.

### Implication for the core-region headline

The current core-region (region (5.16)) headline **µ ≥ 0.380284 (conservative) / 0.3802973
(corrected-tail)** is **unaffected** by the 8 → 4 correction. The change is a
**correctness/provenance fix**, not a numeric lever. At the binding center the effect is
immaterial (~2.7e-8); the only sizeable effect is *positive* (improving) and at a non-binding
center, so it cannot lower the bound and is too small/peripheral to raise the 6-dp headline.

**Caveat on exactness:** these light-N solves rigorously establish the **direction** (corrected
coeff 4 ≥ coeff 8 everywhere; immaterial at the binding center) but only **estimate the
magnitude**. The exact corrected headline requires a **production-config full-cover recompute**
with the corrected coefficient (`N=20000, T=4000, bochner_n=40, pm_k_max=20`) across all 7 rows
plus the `path_b_*` ellipse-extension cover. This memo is also a **core-region** statement; it is
NOT yet a certified full-space µ bound (full-space promotion over White's far regions
R6/R7/R9/R16/R17 is a separate follow-up).

---

## 8. Recommended next steps

1. **Production-config full-cover recompute with corrected coeff** (`mside_sin_coeff = 4.0`,
   the new default): `N=20000, T=4000, bochner_n=40, pm_k_max=20`, all 7 rows, then run the
   `path_b_*` ellipse-extension cover to confirm the residual region (5.16) is still covered.
   This converts the *direction* established here into the *exact* corrected 6-dp headline.
   Expectation (from §5): unchanged at 6 dp; row4/cde binding region moves by ~1e-8, off-binding
   row7 gains ~1e-3 but does not govern the MIN.
2. **Re-run the `path_b` ellipse-extension at production scale** with the corrected coeff and
   re-verify all three implementations (`path_b_rigorous`, `path_b_independent`,
   `path_b_analytical`) agree to ≥10 digits — the constant changed in two of them and the third
   inherits transitively, so the cross-check should be re-confirmed at production config (not just
   the light Config A/B above).
3. **(Separate follow-up, out of scope here)** Full-space promotion over White's far regions
   **R6/R7/R9/R16/R17** to turn the core-region (5.16) bound into a certified unconditional µ
   lower bound. This is independent of the 8 → 4 correction.
4. **Reply to White** confirming: Bochner-PSD validated; 5.6/5.7 corrected 8 → 4 (now the
   default, parametrized for reproducibility); 5.8/5.9 already used `2m−1`; and that we
   empirically reproduce his "not material" finding (binding-center Δ ≈ 2.7e-8), with the
   byproduct that at the off-binding row7 center the correction is a clear +1.6e-3 improvement.

---

*Scratch scripts: `code/_white_corr_row4.py`, `code/_white_corr_row7.py` (and a
cde_n30_iter3 solve). Do not edit `findings.md` from here — the main loop prepends the
one-liner.*
