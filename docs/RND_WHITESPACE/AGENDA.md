# R&D Whitespace Agenda — toward a PROVEN µ bound

**Date:** 2026-06-03
**Status:** SCOUT/PLAN of record. Nothing here is claimed proven. This file ranks the *validated* directions (Phase-1 scout → Phase-2 deep-dive → Phase-3 skeptic) and fixes the recommended thrust.
**Bar (Ben's standard):** a *truly meaningful* result — (1) a proof-grade / solver-independent / formally-verified lower-bound certificate making µ ≥ ~0.3803 a **theorem**; (2) a clean analytic structural theorem; or (3) a substantially stronger bound (past the ~0.380558 framework ceiling, or closing much of the gap from the UB side). Incremental +1e-5 numeric gains do **not** clear the bar.

**Current state being upgraded:** µ ≥ 0.380284 over White's full parameter space, NUMERICALLY certified with load-bearing caveats. Open gap to the Together UB 0.380871 is 5.87e-4.

---

## TL;DR

- **#1 thrust: L2 (Jansson/VSDP verified a-posteriori bound) + interval recomputation of the problem data (`tb`, cell integrals, Fourier coeffs).** This is the cheapest credible route from *numerically-certified* to *verified*, and — only when paired with the interval data recomputation — it is the **minimal FULL verified chain** for the existing µ ≥ 0.380284. Independently reproduced this session: the dual defect `‖c + Aᵀz‖∞ ≈ 2.97e-8` at N=300 with the poly-moment cuts **in the cone**, ~3000× below the +2e-4 binding margin → the rigorous correction costs ~zero bound strength.
- **#2 thrust (parallel hedge): D5-Phase-B** — a small-N, Bochner-only (NO poly-moment) exact-arithmetic rational certificate, the only route yielding a clean **PROVEN constant** (µ ≥ ~0.3793, still White-beating) **decoupled from the poly-moment data caveat**.
- **#3: D1 narrowed** — SDPA-GMP used only on the small binding instance to crush the dual residual where CLARABEL's ~1e-7 is too coarse for the thin Bochner-block margin. Production GMP re-solve is memory-blocked (52–73 GB dense Schur); the production path *is* L2.
- **DROP:** D3 (infeasibility-corner Farkas — seam certified non-load-bearing), D4 (saturation ceiling — target likely ill-posed; Σ m·λ grows with N). **COLLAPSE:** D2/L1/L4 are one thing — the Fejér–Riesz/SOHS finish on the Bochner block, a *sub-step* of L2/D1, not standalone certificates.

---

## The one finding that reframes the whole scout

Every rigor lever scouted (Jansson, GMP re-solve, rational dual snap) certifies **dual feasibility** for the SDP *with the given float64 matrix entries*. None of them certifies that those entries (`tb` poly-moment tail bounds, cell-min integrals, Fourier coefficients) are **valid discretization data**. The chain to µ needs both. The 2026-05-22 trap (tail bound `tb` ~20% too small → invalid bound, corrected the headline) proves this surface is real and biting, and `poly_moment.py` confirms `tb` is a precomputed **float64** fed as the RHS of LP rows `m_k ≥ −tb`.

**Consequence:** since the headline 0.380284 is load-bearing on poly-moment cuts, *no single proposed lever* makes the full 0.380284 a theorem. The two honest destinations are:

- **(a) FULL verified 0.380284:** L2 (verified dual) **+** an `mpmath.iv` interval recomputation of `tb` / cell-integrals / Fourier-coeffs (verified data). Both halves are required; both are tractable in days.
- **(b) Clean PROVEN weaker constant:** D5-Phase-B Bochner-only, which omits poly-moment entirely and so sidesteps the data surface, at the cost of a weaker (still White-beating) constant.

This reframing is why L2-alone is ranked #1 *with the explicit data-recomputation rider*, and why D5-Phase-B is a co-equal hedge rather than a fallback.

---

## Ranked agenda

| Rank | Direction | Payoff class | Effort | Verdict |
|------|-----------|:---:|:---:|---------|
| 1 | **L2 — Jansson/VSDP verified a-posteriori bound + interval data recomputation** | 1 | Low–Med (~1–2 wk) | PURSUE (keystone) |
| 2 | **D5-Phase-B — small-N Bochner-only exact rational certificate** | 1 | Medium | PURSUE (parallel hedge) |
| 3 | **D1-narrowed — SDPA-GMP residual-crusher on the small binding instance** | 1 | Medium | PURSUE only re-scoped |
| 4 | **D2/L1/L4 — Fejér–Riesz/SOHS Bochner-block finish (Magron alg.)** | 1 | Med (subsumed) | MAYBE — only as a sub-step of #1/#2 |
| — | D3 — Farkas/analytic certs for infeasibility corners | 1 (claimed) / ~0 (real) | Low | **DROP** (seam non-load-bearing) |
| — | D4 — saturation/complementarity analytic ceiling C∞ | 2 | High | **DROP** (target likely ill-posed) |

### Rank 1 — L2 + interval data recomputation
**Why:** Cheapest credible "certified → verified" upgrade and the only lever whose make-or-break number was *independently reproduced* (this session: `‖c+Aᵀz‖∞ ≈ 2.97e-8`, `‖·‖₁ ≈ 6.59e-8` at N=300 with poly-moment cuts in-cone). cvxpy already exposes everything Jansson needs — `get_problem_data(cp.CLARABEL)` → `(A,b,c,dims)`; `chain.solve_via_data(...)` → `.x, .z, .r_dual` — so the "serializer" D1 treats as 1–2 weeks of bookkeeping is, for *extracting the dual and cone data*, already built into cvxpy. No re-solve ⇒ no memory ceiling ⇒ runs at production N=20000–24000 today. The bounded-primal hypothesis is real (Ω≤1; 0≤w,v≤1; |c|,|d|≤2/π; tail caps) so Jansson Alg 3.1 terminates in one pass (no perturbation loop). **The rider that makes it a FULL certificate:** add an `mpmath.iv` recomputation of `tb` (the IBP tail bound), the cell integrals, and the Fourier coefficients, so the certified statement is "µ ≥ p_lo for the SDP whose data is *verified-valid*", not merely "for the SDP with these float entries".
**Honest ceiling:** L2 alone retires the dual-feasibility caveat but **not** the poly-moment data caveat. CLARABEL (FP) still *finds* the witness; L2 *certifies* it (the Cohn–Elkies discipline). Verified ≠ rational/Lean-proven, but the trusted base shrinks from "the whole IPM + its log parse" to "a small interval library + cross-checkable data extraction".

### Rank 2 — D5-Phase-B (Bochner-only exact rational certificate)
**Why:** The only pick producing a clean PROVEN constant *decoupled from the load-bearing poly-moment caveat*. Bochner-only (exact-(W.1) cell integral, no poly-moment cuts) first beats White 0.379005 at N≈900–1000 (rigorous dual_LB ≈ 0.37925, +2.45e-4). Because it omits poly-moment, it sidesteps the entire tail-bound data surface that makes 0.380284 fragile. The exact-arithmetic rational route (de Laat–Dostert SDP-rounding precedent, arXiv:2001.00256) is medium effort.
**Honest ceiling:** the certificate is **large** (~9164 conic rows at N=1000, 7089 LP) — the "verify one PSD matrix" framing understates the target by ~3 orders of magnitude, so **DROP the standalone-Lean pitch** (Mathlib has `Matrix.PosSemidef` but ZERO conic/SDP weak-duality and no Bochner–Herglotz; that is a Mathlib-contribution-sized, months-long build). Stop at the exact-arithmetic rational certificate — still proof-grade, still payoff-class-1.

### Rank 3 — D1-narrowed (GMP residual-crusher)
**Why:** Validated *only* in re-scoped form. The literal production GMP re-solve is memory-blocked (N=10000 → 36,083 vars / 74,060 eq rows → dense GMP Schur 52–73 GB; only N≤2000 fits 4 GB), so the scalable production tool is L2. GMP retains a *narrow, real* role: the binding Bochner-block Jansson correction (≈ 42·dⱼ⁻·21) needs λ_min(D_Bochner) ≳ −1e-8, but CLARABEL's float residual is ~1e-7 — so for the *binding* row, GMP-on-the-small-instance (or N≥24000 margin-widening) is **required, not optional**, to keep the correction inside the margin. Also requires fixing a stale SDPA-GMP dylib rpath (use `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib`, GMP 6.3.0 — 10-min fix; spin off separately).

### Rank 4 — D2/L1/L4 (Fejér–Riesz / SOHS Bochner-block finish)
**Why MAYBE-only:** These are **one** direction, not three — the exact trig-SOHS finish on the Bochner block (use **Magron arXiv:2202.06544**, which has a concrete algorithm + bit-complexity bound; prefer it over Davis–Papp). It is the *finishing move* of a rational-certificate pipeline, **not a standalone certificate**: the binding-center dual is a mixed LP+SOC+PSD point (~1890 scalar components at N=300; the Bochner block is ~390 of them and the `1−f≥0` block is identically inactive), and PSD-ness of the dual Bochner block is *automatic* from any interior-point solver. So "rational Fejér–Riesz of the Bochner block ⇒ certificate" is a non sequitur — the certificate is the **whole** mixed-cone dual. The favorable facts (Bochner dual is complex-rank ≤5; one block) are real but the work is dominated by the surrounding cone + the data caveat. **L1 (Davis–Papp) as stated is dropped:** wrong cone class (it certifies one univariate-interval WSOS cone; White's dual is LP-dominated mixed-cone), boundary obstruction (binding bound is a saturated optimum; Scheiderer — rational polys on the boundary may have no rational SOS), and it needs a re-solve anyway (dual_extractor returns only the scalar). Its only unique block (Bochner) is done better by Magron.

---

## DROPPED — with reasons (do not re-propose)

- **D3 (Farkas/analytic certs for infeasibility-excluded corners) — DROP.** Self-defeating: `lp_research_state/FULLSPACE_VERIFICATION.md` (the "Infeasibility-exclusion rigor — NOT load-bearing anywhere" section) states every region's certified floor is set by the **feasible** part where the cover clears on grid+Lipschitz geometry; excluded corners have cover values 0.47–1.6 ≫ target. A Farkas cert changes neither the value (0.3802838) nor its validity. The clean analytic handle (|c[0]|, |d[0]| ≤ 2/π) covers only the *far* slabs, not the empirically-infeasible annulus (q≈0.35–0.40 at p=0.37, well below 2/π=0.637), which is nonlinear-constraint-driven with no closed form. The solver exposes 0/82 nonzero duals on infeasible corners, so the proposed "ask the solver for the dual ray" first step fails. **Keep only** two one-line slab lemmas (|c[0]|, |d[0]| ≤ 2/π) as preprint hygiene, plus the monotone-cover lemma (cover_min ≤ V because V=+∞ on infeasible points) which *eliminates* the need for any infeasibility cert.

- **D4 (saturation/complementarity analytic ceiling C∞≈0.380558) — DROP.** Target plausibly ill-posed: `LEVER_I_PRIME_FINAL.md` records Σ m·λ **growing +73%** from N=30K→40K (5.97→10.32), with the (π/2N) decay barely outpacing it; the author flags the N-uniform Σ m·λ bound — the real open lemma — as possibly **FALSE**, so there may be no finite C∞ to prove. The proposed lever (a-priori |ξ| ≤ C·Ω) bounds the **wrong** quantity (Theorem 2 bounds Σ λ, unweighted; the 2R bridge to Σ m·λ loses 20× → vacuous ceiling 0.381209 > UB even with the exact empirical |ξ|=1.46Ω). The scale-invariance first step is blocked (PRO-14: the c[0] anchor, sum_squares(c,d)≤0.5, and the 2/π box all fail to scale). Even full success is class-2 and bounds only the **framework**, not µ.

---

## Recommended thrust plan (ordered, multi-step path to the meaningful result)

**Track L2 (keystone, production-scalable):**
1. **Data path:** wrap `get_problem_data(cp.CLARABEL)` + `chain.solve_via_data` into a helper returning `(A,b,c,dims,x,z,r_dual)`; verify `c@x == prob.value` and locate Ω's canonical index. *(Reconnaissance done this session — works.)*
2. **Per-cone-block dual-slack defect** in CLARABEL's convention (`c + Aᵀz`): for PSD blocks unpack the scaled symmetric-vectorization slice into a symmetric matrix; unit-test vs a dense numpy reconstruction at bn=6.
3. **Rigorous λ_min lower bound in `mpmath.iv`** (verified-Cholesky-shift; fallback interval Gershgorin); validate the interval encloses and lies below `numpy.linalg.eigvalsh` on small blocks.
4. **Per-block x̄ⱼ** finite primal bounds (Ω≤1; w,v≤1; |c|,|d|≤2/π; tail caps; SOC/PSD lifts bounded) → assemble Jansson (3.7) `p_lo = inf{bᵀỹ + Σ sⱼ·dⱼ⁻·x̄ⱼ}` in interval arithmetic for one binding center; success = within ~1e-6 of the current rigorous_dual_LB and ≥ White-beating threshold.
5. **DATA RIDER (the half that makes it FULL):** recompute `tb` (IBP tail bound), cell-min integrals, and Fourier coefficients in `mpmath.iv`; confirm the float64 values used by the solve lie inside the verified intervals, with margin. *(This closes the long-pole the 2026-05-22 trap exposed.)*
6. **Lift to µ:** fold per-center verified `p_lo` into the existing adaptive grid+Lipschitz cover (`path_b_*` / `_fullspace_eval`); cross-check via `path_b_independent.py` to 10+ digits; report the verified full-space floor.

**Track D5-Phase-B (parallel hedge, clean PROVEN constant):**
7. Pin the smallest White-beating Bochner-only exact-(W.1) instance (≈N=1000, T=500, R=10, bochner_n=14; push to N=1500–2000 to buy headroom if ellipse extension consumes the +2.45e-4).
8. Re-solve that center at GMP precision (SDPA-GMP, after the dylib fix) → high-precision dual; back off to strictly-interior c < Ω; snap the dual to rationals (de Laat–Dostert); verify dual feasibility + objective ≥ target in exact arithmetic (sympy/python-flint). Use **Magron (D2/L4)** for the Bochner-block SOHS piece if a clean rational PSD factorization is needed.
9. Lift to µ via `path_b_analytical` re-pointed at the exact program; confirm the 7 ellipses still cover White's region (5.16) with positive margin after extension.

**Convergence:** Track L2(+rider) → FULL verified µ ≥ 0.380284. Track D5-B → standalone PROVEN µ ≥ ~0.3793 with no poly-moment caveat. Run in parallel; they share the data-extraction and GMP infra and hedge each other.

---

## First action (cheapest, highest-information)

**Implement and run the L2 Jansson per-block defect + rigorous `mpmath.iv` λ_min on ONE binding center, and assemble `p_lo` (Track L2, steps 1–4).** The dual defect and cone/data extraction were already reproduced this session (`‖c+Aᵀz‖∞ ≈ 2.97e-8` at N=300 with poly-moment cuts in-cone); the only unproven-in-practice piece is whether the *rigorous interval* λ_min lower bound on the two Bochner PSD blocks (plus the LP/SOC correction) keeps `p_lo` above the White-beating threshold with the margin intact. That single number is make-or-break for the entire #1 thrust and is a ~1-day mpmath exercise on the already-extracted data — no new solves, no GMP, no serializer.

---

## Expected payoff

- **Realistic (Track L2 + data rider):** µ ≥ 0.380284 promoted from *numerically-certified* to a **verified** bound — a number `p_lo` plus an interval-arithmetic computation proving SDP_opt(center) ≥ p_lo from independently re-derivable, interval-validated data, with no dependence on CLARABEL's status flag or log format. Clears bar item (1) in its *verified* sense; the strongest (rational/Lean) sense is the L2→Magron rational finish (Track D5/D2).
- **Realistic (Track D5-Phase-B):** a clean, defensible, publishable **first proof-grade improvement over White** — µ ≥ ~0.3793 as a theorem with a checkable rational dual certificate, solver-independent, decoupled from the poly-moment caveat. Modest in constant, unimpeachable in rigor; clears bar item (1).
- **Does NOT deliver:** a stronger *value* (neither track pushes past the conjectured C∞≈0.380558 ceiling — this is a RIGOR upgrade, not a strength upgrade), nor any UB-side gap closing (entirely unscouted and in-scope for the bar, but not on this agenda).

---

## Caveats (honest)

1. **Data-certification is the long pole.** L2 (and GMP, and rational snap) certify the program *as written*; the chain to µ needs `tb`/cell-integrals/Fourier-coeffs to be valid data. The full verified 0.380284 requires the `mpmath.iv` data-recomputation rider (step 5), not L2 alone. The 2026-05-22 tail-bound trap proves this surface is live.
2. **Verified ≠ proven-in-the-strongest-sense.** L2 yields a *verified* (interval) certificate; CLARABEL (FP) still *finds* the witness. The fully solver-independent / rational / Lean destination is the Magron rational finish or D5-Phase-B.
3. **Bochner-block residual vs margin (Track L2/D1).** The binding-row correction needs λ_min(D_Bochner) ≳ −1e-8; CLARABEL's ~1e-7 float residual may be insufficient for the binding row → GMP-on-small-instance or N≥24000 margin-widening may be required. If neither closes it, only a *weaker* White-beating constant (≈0.3795/0.380) is re-certifiable at full rigor, not the full 0.380284 (still payoff-class-1).
4. **D5 ellipse-extension headroom (risk F1).** The +2.45e-4 single-center margin can be consumed by the off-center Ω drop in the cover; may force N=1500–2000, inflating the certificate. The fully-exact µ-bound also needs the cover/Lipschitz step in interval/rational arithmetic (Cohn–Triantafillou discipline) — real extra work.
5. **Certificate size (Track D5).** ~9164 conic rows at N=1000 — large but checkable in exact arithmetic; **not** Lean-tractable as scoped (Mathlib lacks conic duality). Stop at the rational certificate.
6. **The serializer is NOT the keystone.** The internal D1-D5 memo framed the cvxpy→SDPA-S serializer as the keystone; the validated finding is that L2 reaches the verified goal *without* it (cvxpy already exposes the dual + cone data), and the production GMP re-solve the serializer enables is memory-blocked anyway. Build the serializer only for the small D5/D1 GMP instance.

---

## Provenance

- Phase-1 internal scout: `docs/RND_WHITESPACE/2026-06-03_whitespace_scout.md` (D1–D5).
- Phase-1 external scout: `docs/RND_WHITESPACE/2026-06-03_external_literature_scan.md` (L1–L4).
- Phase-2 deep-dives: `docs/RND_WHITESPACE/{D1..D5,L1,L2}_*.md`.
- Phase-3 skeptic: megamemory node `whitespace-stress-test-2026-06-03-data-certification-is-the-missed-long-pole`.
- This-session independent reproductions: Jansson dual defect `‖c+Aᵀz‖∞ ≈ 2.97e-8` / `‖·‖₁ ≈ 6.59e-8` at N=300 with poly-moment cuts compiled into the cone (`build_problem(300,120,6,0.004,0.004,0.3875,0.3875,-0.02,0.02,bochner_n=6)`; cone dims = 1 eq / 2001 ineq / 26 SOC / two PSD-14 Bochner blocks); `dual_extractor.py` confirmed log-parse-scalar at 1e-4 residual threshold (no dual vector); `poly_moment.py` confirmed float64 `tb` as LP-row RHS; `FULLSPACE_VERIFICATION.md` confirmed infeasibility exclusions non-load-bearing.
