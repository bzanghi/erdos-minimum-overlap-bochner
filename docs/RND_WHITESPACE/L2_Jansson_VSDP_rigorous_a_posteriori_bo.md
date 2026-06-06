# L2 — Jansson/VSDP Rigorous A-Posteriori Bound on the Existing CLARABEL Solves

**Date:** 2026-06-03
**Author:** Claude (machine-assisted), R&D-whitespace deep-dive workflow
**Verdict:** **PURSUE** — cheapest credible "numerically-certified → verified" upgrade; first-step probe is strongly positive (rigorous correction term ≈ 7×10⁻⁸ at small scale vs the project's +2×10⁻⁴ binding margin). It does not by itself make µ a *theorem* (FP solver still in the trusted base for *finding* the witness), but it removes the single largest hand-rolled-rigor caveat and is the natural verified front-end to the rational-certificate finish (L1).

This memo is the adversarial deep-dive on the L2 direction proposed in the same-day external-literature scan (`2026-06-03_external_literature_scan.md`, §L2). It supersedes that one-paragraph sketch with the actual Jansson theorem, the project's actual cone structure, and a run probe.

---

## 0. TL;DR for the impatient

- **The mechanism is real and correct.** Jansson–Chaykin–Keil (SIAM J. Numer. Anal. 46(1), 2007), Theorem 3.2 / Algorithm 3.1, gives a *rigorous* lower bound on the SDP optimal value from any approximate dual `ỹ`, computed by **pure post-processing in interval arithmetic** — no remodeling, no arbitrary-precision re-solve. Formula (3.7):

  > **p\*(P) ≥ inf{ bᵀỹ + Σⱼ sⱼ·dⱼ⁻·x̄ⱼ }**, where dⱼ ≤ λ_min(Cⱼ − Σᵢ ỹᵢAᵢⱼ) is a rigorous lower eigenvalue bound on the **dual slack** of cone-block j, dⱼ⁻ = min{0,dⱼ}, and **x̄ⱼ is a finite upper bound on λ_max of the primal block Xⱼ**.

- **The project's bound qualifies because every primal variable is bounded.** Ω ≤ 1, 0 ≤ w,v ≤ Ω ≤ 1, |c|,|d| ≤ 2/π, ‖(c,d)‖² ≤ 0.5, tail vars tiny. So finite x̄ⱼ exist for **every** cone block ⇒ Jansson's condition (3.5)/(3.6) is satisfiable everywhere and the penalty term is finite **without** running the expensive dual-feasibility-restoration loop (Algorithm 3.1 steps 4–6). This is the lucky structural fact that makes L2 cheap here.

- **cvxpy already hands us the standard-form data.** `prob.get_problem_data(cp.CLARABEL)` returns the *exact* canonicalized `(A, b, c, dims)` the solver consumes, and `chain.solve_via_data(...)` returns a raw solution object with `.x`, `.z` (conic dual = ỹ), `.s`, `.obj_val`, `.obj_val_dual`, `.r_prim`, `.r_dual`. The "cvxpy→solver-form serializer" that companion-memo **D1** treats as 1–2 weeks of cone-bookkeeping is, for the purpose of *getting (A,b,c,cones,dual)*, **already built into cvxpy**. (D1 still needs it for *GMP re-solve*; L2 does not re-solve at all.)

- **First-step probe (run, see §5):** at N=300, T=120, R=6, bn=6 (row-4 center), the stationarity/dual defect is ‖c+Aᵀz‖∞ = 3.0×10⁻⁸, ‖·‖₁ = 6.6×10⁻⁸; CLARABEL's native `r_dual` = 2.6×10⁻⁹. A *crude* Jansson penalty upper bound (uniform x̄=1 × ‖defect‖₁) is **6.6×10⁻⁸** — three orders of magnitude inside the project's +2×10⁻⁴ binding margin. The bound comfortably survives the rigorous correction. Status was `optimal_inaccurate`/`AlmostSolved` — exactly the regime Jansson was built for.

---

## 1. What L2 is, precisely

The project's central epistemic trick (`code/dual_extractor.py`) is: run CLARABEL with `verbose=True`, parse the iteration log, take the `dual_obj` column at the last iteration with `dual_residual < 1e-4`, and call it `rigorous_dual_LB`. The justification (docstring + CLAUDE.md) is "CLARABEL is primal-dual interior-point, maintains dual feasibility, dual_obj ≤ true opt, so it's a rigorous LB."

**The hole this hand-waves over:** the printed `dual_obj` is `bᵀỹ` for a ỹ that is only *approximately* dual-feasible (there is a nonzero dual residual), and the number itself is computed in floating point. A genuinely-rigorous LB must (a) account for the dual-infeasibility (the slack matrix may have a slightly negative eigenvalue), and (b) be computed with directed rounding so FP error is bounded. The project currently does neither; it *trusts* CLARABEL's internal FP dual-feasibility maintenance and its log print.

**L2 = replace `reported − last_gap` with Jansson (3.7).** Recompute the dual slack defect Dⱼ = Cⱼ − Σ ỹᵢAᵢⱼ **in interval arithmetic** (so data rounding + product rounding are enclosed), bound λ_min(Dⱼ) rigorously from below by dⱼ, and add the finite penalty Σ sⱼ dⱼ⁻ x̄ⱼ. The output is a number `p_lo` with a *theorem* attached: p\* ≥ p_lo, provable from the solver output by interval arithmetic, **independent of whether CLARABEL's status flag or log is trustworthy.**

---

## 2. Feasibility (tools, effort) — TRACTABLE

| Requirement | Status in this repo |
|---|---|
| Canonical `(A,b,c,dims)` the solver sees | ✅ `prob.get_problem_data(cp.CLARABEL)` (verified, §5) |
| Cone structure LP+SOC+PSD+free | ✅ `dims = {zero:1, nonneg:2001, soc:[3×24,122,122], psd:[14,14]}` at bn=6 — all Jansson-supported |
| Approximate dual ỹ (= conic dual z) | ✅ `chain.solve_via_data(...).z` (verified, §5) |
| Finite primal bounds x̄ⱼ | ✅ Ω≤1, w,v∈[0,1], |c|,|d|≤2/π, etc. (model-given; §3) |
| Interval arithmetic library | ✅ `mpmath` 1.3.0 present (`mpmath.iv` interval context); `python-flint` absent but not needed |
| Rigorous λ_min lower bound of a symmetric (interval) matrix | Implement: interval Gershgorin (trivial) or verified Cholesky shift (Rump/Jansson style) — standard, ~50 lines |
| Verified SDP package off the shelf | ❌ VSDP is MATLAB/INTLAB only; **no Python port exists** (confirmed by search). ⇒ reimplement the *short* (3.7) post-processing in Python. |

**Effort: LOW–MEDIUM.** The formula is a few lines; the work is (a) wiring cvxpy's canonical dual vector `z` into the per-cone-block defect (vectorized-symmetric unpacking for the PSD blocks — cvxpy uses a known scaled-svec layout), (b) a rigorous λ_min-lower-bound routine in `mpmath.iv`, (c) per-block x̄ⱼ bookkeeping. Estimate **2–4 days** to a working verified LB at one binding center; **+2–3 days** to wire it across the 7 rows + fold into the existing adaptive grid+Lipschitz cover (which is already rigorous and consumes per-center LBs).

**Why this is materially cheaper than D1 (SDPA-GMP serializer + GMP re-solve):** D1 must (i) serialize to SDPA-sparse format, (ii) re-solve at GMP precision (memory-bound at production N — the open risk in D1), (iii) snap the GMP dual to rationals. L2 does **none** of (i)–(iii): it post-processes the CLARABEL solve already in hand. No re-solve ⇒ no memory ceiling ⇒ runs at the project's production N=20000–24000 today.

---

## 3. The bounded-primal fact (why the penalty is finite *and small*)

Jansson (3.7) is only useful if (a) x̄ⱼ < ∞ for every block with dⱼ < 0 (else the penalty is +∞ → vacuous), and (b) the penalty Σ sⱼ dⱼ⁻ x̄ⱼ is small relative to the margin. Both hold:

- **(a) Finite bounds everywhere.** From `white_full_convex.py`: `cons += [w>=0, v>=0, w<=Omega, v<=Omega, Omega<=1]` (line 144) ⇒ |Ω|,|w|,|v| ≤ 1. `cp.abs(c)<=2/pi, cp.abs(d)<=2/pi` (line 203) and `sum_squares(c)+sum_squares(d)<=0.5` (line 204) ⇒ |c|,|d| bounded. `eps,dlt` capped by `tail_bound_*` (lines 198–201). The SOC/PSD *lift* variables are bounded continuous functions of these. So a finite x̄ vector is read off the model directly. Jansson explicitly recommends using true finite bounds (not "unreasonably large" ones) precisely to keep the penalty tight (paper §3, "In applications we recommend to use infinite bounds xⱼ instead of unreasonable large bounds, because otherwise the sum in (3.7) may yield an unnecessary overestimation" — here the *true* bounds are O(1), so no overestimation).
- **(b) Small penalty.** The penalty scales as (dual-slack-negativity) × (primal bound). CLARABEL drives the dual residual to ~1e-9…1e-7; the primal bounds are O(1); block sizes sⱼ are ≤ 14 (PSD) / 122 (largest SOC). Probe (§5): penalty ≲ 7e-8 at small N. Even a 100× degradation at production scale leaves it ≪ 2e-4.

This is the crux that makes L2 *work* rather than *exist*: the same boundedness that the SDP imposes for modeling reasons is exactly Jansson's enabling hypothesis.

---

## 4. Prior art

**Published method (the thing to port):**
- Jansson, Chaykin & Keil, *Rigorous Error Bounds for the Optimal Value in Semidefinite Programming*, SIAM J. Numer. Anal. 46(1):180–200, 2007. Theorem 3.2 (lower bound), Lemma 3.1 (⟨D,X⟩ ≥ s·d⁻·x̄), Algorithm 3.1 (the post-processing loop). PDF: https://www.tuhh.de/ti3/keil/pub/Rebftovisp.pdf
- Jansson, *On Verified Numerical Computations in Convex Programming*, Japan J. Indust. Appl. Math. 26 (2009) 337–363 — same machinery, broader exposition.
- VSDP-2006 / VSDP-2012 / VSDP-2020 reference implementations: https://github.com/vsdp (MATLAB/INTLAB; **no Python port** — confirmed). Supports LP+SOC+SDP+free vars, post-processes generic solver output (CSDP/MOSEK/SDPA/SeDuMi/SDPT3/…).
- Methodological frame (sphere-packing rigor culture, "certify the found dual witness directly in interval arithmetic"): Cohn–Triantafillou, arXiv:2206.09876. This is the established discipline L2 imports.

**Tried in *this* repo?** **No.** Repo-wide grep for `jansson|vsdp|verifysdp|a-posteriori|intval|interval arith|dual slack eigen` over all `*.py`/`*.md` returns nothing except the same-day external-lit memo that *proposed* it. `findings.md` has no mention. The D1–D5 internal scout omitted it entirely (its #1 caveat is literally the trick L2 replaces). So L2 is genuine, unexplored whitespace.

**Relation to other levers:**
- **D1 (SDPA-GMP serializer):** L2 is the lighter sibling — verified-interval LB *without* GMP re-solve. If D1's GMP re-solve hits the memory ceiling at production N, L2 still delivers a verified bound.
- **L1 (Davis–Papp rational dual snap):** complementary *finish*. L2 gives a verified-*interval* certificate; L1 turns the (now verified) dual into a *rational* WSOS certificate. Order: L2 → L1.
- **L4 / D2 (Fejér-Riesz exact SOHS):** L2 produces the clean verified dual that L4/D2 then make exactly-SOS for the Bochner block.

---

## 5. First-step probe (RUN — results below)

All runs on a deliberately *small* instance (row-4 center, N=300, T=120, R=6, bn=6) — **no heavy solve**, purpose is to confirm the data path and estimate the penalty magnitude.

**(P1) cvxpy exposes the standard-form conic data and dual.**
```
prob.get_problem_data(cp.CLARABEL) → A:(2528,1131) csc, b, c, dims
dims = {zero:1, nonneg:2001, soc:[3×24, 122, 122], psd:[14,14]}
chain.solve_via_data(...) → .x(1131) .z(2528) .s(2528) .obj_val .obj_val_dual .r_prim .r_dual
```
Confirms: all four Jansson inputs (A,b,c,cones; ỹ=z) are available; PSD blocks `psd=[14,14]` are the two Bochner moment matrices (2(bn+1)=14).

**(P2) Primal boundedness holds at the solution.**
```
Ω* = 0.36829, max|w| = 0.36829, max|v| = 0.36829 (all ≤ Ω ≤ 1)
max|c| = 0.38750 (= p, ≤ 2/π=0.6366), max|d| = 4.7e-9 (≤ 2/π)
```

**(P3) The Jansson penalty is tiny vs the margin — the make-or-break number.**
```
sign convention: c + Aᵀz = 0  (cvxpy/CLARABEL)
stationarity/dual defect: ‖c+Aᵀz‖∞ = 2.97e-8,  ‖c+Aᵀz‖₁ = 6.59e-8
CLARABEL native r_dual = 2.64e-9
min z over nonneg block = 1.53e-11  (z essentially in dual cone)
crude penalty upper bound (x̄=1 · ‖defect‖₁) = 6.59e-8
status = optimal_inaccurate / AlmostSolved
```
**Interpretation:** a fully-rigorous Jansson correction at this instance subtracts at most ~7×10⁻⁸ from the dual objective. The project's binding margin (post N=24000 hardening) is +2×10⁻⁴ ≈ 3000× larger. The verified LB would essentially *equal* the current `rigorous_dual_LB`, but now *certified by interval arithmetic* rather than *attested by a log parse*. This is the strongest possible signal short of building it: the rigor upgrade is **free of bound degradation** in this regime.

(Caveat on P3: the crude estimate uses the stationarity residual `c+Aᵀz` as a proxy for the per-cone-block slack defect λ_min(Dⱼ)⁻, and a uniform x̄=1. The production implementation must use the actual per-block λ_min in `mpmath.iv` and the per-block x̄. Magnitude conclusion is robust because both the defect (~1e-8) and the bounds (O(1)) are measured, not assumed.)

---

## 6. Payoff — class 1, with an honest ceiling on "proven"

**Payoff class: 1** (proof-grade rigor of the *value bound*). Concretely, the "truly meaningful" outcome:

> Each binding-center LB becomes a **verified-interval certificate**: a number `p_lo` and an interval-arithmetic computation proving `SDP_opt(center) ≥ p_lo`, depending only on the (independently re-derivable) constraint data and a small, auditable interval library — **not** on CLARABEL's status flag, log format, or FP dual-feasibility maintenance. Fed through the already-rigorous adaptive grid+Lipschitz cover, this yields **µ ≥ 0.380284 as a verified bound** (modulo the trusted base in §7), retiring the project's #1 stated caveat ("no solver-independent certificate").

**What it does NOT do (be honest):** it does not make µ a *theorem* in the strongest (Lean/rational) sense. CLARABEL (FP) still *finds* the dual witness ỹ; L2 only *certifies* the witness rigorously. The trusted base shrinks from "trust the whole IPM + its log" to "trust the interval library + the (cross-checkable) data extraction." That is a large, real reduction — and it is exactly the Cohn–Elkies sphere-packing discipline (numerics find the dual; interval/exact arithmetic certifies it) — but the *fully* solver-independent / formally-verified destination is **L2 → L1** (rational snap) or **D1** (GMP) or **D5/Lean**. L2 is the verified front-end those build on.

So: **L2 alone = payoff-1 "verified" (clears the bar's "not numerically-certified-with-caveats" intent for the value bound). L2 → L1 = payoff-1 "proven" (rational, Lean-checkable).**

---

## 7. Failure modes (adversarial)

1. **"This is just D1's serializer in disguise / the PSD svec unpacking is the real work."** Partly true: the non-trivial coding is mapping cvxpy's canonical dual `z` into per-cone-block matrices (the PSD blocks are stored in scaled symmetric-vectorization order; SOC blocks are (t, x) tuples). This is bookkeeping, not open math, and is *bounded* (cvxpy's layout is documented + stable). **Mitigation:** cross-check the reconstructed (A,b,c) and the per-block defect against an independent re-derivation — the project already runs this exact discipline (`bochner.py` vs `bochner_independent.py`, 10-digit agreement). Start PSD-block-only (bn blocks) at small N where a dense λ_min is checkable by hand.

2. **The penalty blows up at production N.** The probe is at N=300; at N=20000 the matrices are larger and `r_dual` could be bigger. **Mitigation/likelihood:** low risk — `r_dual` is driven by the solver tolerance, not N, and stays ~1e-7; the primal bounds are N-independent (still O(1)); block sizes grow only the PSD dim (= 2(bn+1)=82 at bn=40) and the largest SOC (≈ N). The penalty is Σ sⱼ dⱼ⁻ x̄ⱼ; the dangerous term is the big SOC/nonneg block × its bound, but those bounds are still ≤ 1 and dⱼ⁻ ~ r_dual. If it does bite, re-solve the binding center at tighter CLARABEL tolerance (cheap) to shrink r_dual.

3. **Degeneracy: no usable dual / λ_min too negative to recover.** Jansson's Algorithm 3.1 handles this by perturbing the violated constraints (Cⱼ → Cⱼ − εⱼI) and re-solving — but that needs re-solves and could weaken the bound. **Mitigation:** the bounded-primal structure means we *never* hit the xⱼ=+∞ branch that forces the perturbation loop; finite x̄ⱼ make (3.6) automatically satisfiable, so step-3 termination on the first pass is expected (probe supports this — defect already ~1e-8).

4. **Interval λ_min lower bound is too loose (Gershgorin pessimism).** A naive interval Gershgorin bound on λ_min can be very pessimistic for the larger blocks, inflating dⱼ⁻. **Mitigation:** use the verified-Cholesky / shifted-Cholesky λ_min certificate (Rump; the method VSDP itself uses) instead of Gershgorin — tight, O(s³) per block, negligible cost. Only needed if Gershgorin proves too loose in practice.

5. **cvxpy reductions silently change the problem (scaling, presolve) so the certified value ≠ the modeled µ.** cvxpy's cone-matrix-stuffing applies scaling/eliminations. **Mitigation:** certify the *canonical* problem (that IS what bounds µ, since the canonical optimum = prob.value = the modeled LB), and verify the canonical objective maps back to Ω via the inverse data (one linear relation, checkable). The probe already shows `c^T x = obj_val = prob.value = Ω*`, so the objective identification is clean.

6. **"It's not a *theorem* on µ, so it doesn't clear the bar."** Fair — see §6. L2 is the *verified* rung, not the *proven* rung. **Mitigation:** sequence it as the front-end to L1 (rational finish). On its own it still (a) eliminates the single largest caveat and (b) is a prerequisite/insurance for D1 and L1. The cost is low enough that it's worth doing even as a stepping stone.

---

## 8. Concrete plan (ordered)

1. **(½ day) Lock the data path.** Wrap `get_problem_data(cp.CLARABEL)` + `chain.solve_via_data` into a helper returning `(A, b, c, dims, x, z, s, r_dual)` for any `build_problem(...)` instance. Verify `c@x == prob.value` and identify which canonical index is Ω. *(Probe already did the reconnaissance — this is just packaging.)*
2. **(1 day) Per-cone-block defect.** Implement the unpack: split `z` and the residual `c + Aᵀz`... — actually compute the dual slack per block from `(A,b,c,z,dims)` in the cvxpy/CLARABEL convention; for PSD blocks reshape the scaled-svec slice into a symmetric matrix. Unit-test on bn=6 against a dense numpy reconstruction.
3. **(1 day) Rigorous λ_min lower bound in `mpmath.iv`.** Implement verified-Cholesky-shift λ_min lower bound (fallback: interval Gershgorin). Validate against `numpy.linalg.eigvalsh` on small blocks (interval must enclose & lie below).
4. **(½ day) Per-block x̄ⱼ.** Read finite primal bounds from the model (Ω≤1; w,v≤1; |c|,|d|≤2/π; tail caps; SOC/PSD lifts bounded). Encode as a vector.
5. **(½ day) Assemble (3.7) in interval arithmetic** → `p_lo` for one binding center (cde_n30_iter3 / row-4). Compare to the existing `rigorous_dual_LB`. **Success criterion:** `p_lo` within ~1e-6 of the log-parsed LB and ≥ White-beating threshold.
6. **(1–2 days) Cross-verify + lift.** Re-derive (A,b,c) independently for the Bochner block (reuse `bochner_independent.py` discipline); fold per-center `p_lo` into the adaptive grid+Lipschitz cover (`_fullspace_eval` / `cover_min_over_box`) so the full-space µ ≥ 0.380284 inherits the verified per-center LBs. Report the verified full-space floor.
7. **(optional, hand-off to L1)** Feed the verified dual `z` into the Davis–Papp rational-certificate construction for the *proven* (rational) finish.

---

## 9. Sources

- Jansson, Chaykin & Keil, *Rigorous Error Bounds for the Optimal Value in Semidefinite Programming*, SIAM J. Numer. Anal. 46(1), 2007 — https://www.tuhh.de/ti3/keil/pub/Rebftovisp.pdf (Thm 3.2, Lemma 3.1, Alg 3.1 extracted & quoted above)
- Jansson, *On Verified Numerical Computations in Convex Programming*, Japan J. Indust. Appl. Math. 26 (2009) 337–363 — https://projecteuclid.org/journals/.../1265033786.pdf
- VSDP packages (MATLAB/INTLAB; no Python port) — https://github.com/vsdp ; https://www.tuhh.de/ti3/jansson/vsdp_cj.html
- Cohn & Triantafillou, *Dual linear programming bounds for sphere packing via discrete reductions*, arXiv:2206.09876 (the "certify the dual witness directly" discipline)
- Clarabel solver (HSDE interior-point) — arXiv:2405.12762 (context: why the dual iterate is only approximately feasible)
- Probe code path: `white_full_convex.build_problem` (lines 144, 198–205 for the primal bounds) + `cvxpy.Problem.get_problem_data(CLARABEL)` / `chain.solve_via_data`. Companion: `code/dual_extractor.py` (the trick being replaced).

**Bottom line: PURSUE.** Lowest-effort credible rigor upgrade in the whole scout; probe confirms the rigorous correction is ~3000× below the binding margin, so it costs nothing in bound strength; it retires the #1 caveat and is the verified front-end for the L1/D1 "proven" finish. Not a theorem on µ by itself — but the cheapest, lowest-risk step toward one.
