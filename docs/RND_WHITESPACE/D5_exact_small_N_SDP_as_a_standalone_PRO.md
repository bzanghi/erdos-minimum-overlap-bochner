# D5 — Exact small-N SDP as a standalone PROVEN theorem at a weaker White-beating constant

**Author:** Claude (machine-assisted), whitespace deep-dive workflow
**Date:** 2026-06-03
**Status:** Assessment only. SCOUTS and PLANS; certifies nothing.
**Companion memos:** `2026-06-03_whitespace_scout.md` (D1–D5 internal track; D5 first stated there as a hedge), `2026-06-03_external_literature_scan.md` (L1–L4 citation track).

---

## 0. Restatement of the direction (as proposed in Phase 1)

Instead of chasing the full thin-margin µ ≥ 0.3803 certificate over White's entire parameter
space, prove a **fully-exact, solver-independent** theorem at whatever constant the *smallest
tractable* exact SDP yields — still beating White (2023)'s µ ≥ 0.379005. Concretely:

1. Use `white_full_convex_exact.py`'s exact-(W.1) program (cell **integral**, not cell-min relaxation)
   at small N with **only the Bochner block** — no poly-moment cuts (sidesteps the load-bearing
   poly-moment caveat #2), single binding center + ellipse extension.
2. Serialize to SDPA-GMP, snap the dual to rationals, verify dual feasibility in exact arithmetic.
3. If scoped tightly, formalize the **finite** PSD/feasibility check in Lean
   (`Matrix.PosSemidef` + rational linear algebra).

**Payoff class claimed:** 1 (proof-grade / plausibly machine-checked).

---

## 1. FEASIBILITY

### 1.1 The exact program exists and already beats White at small N (NEW first-step data)

`white_full_convex_exact.py` (ERD-9, written 2026-05-31) already builds exactly the program D5
needs: the EXACT cell-envelope integral `I_m(j) = ∫ cos(πmx/2) dx` (so (W.1) is an *identity* in the
Fourier variables, not a relaxation), with an optional Bochner block and **no poly-moment**. So the
primal-side construction is done; D5 inherits it.

**First-step probe (run this session, binding row-4 center h=0.004, p=0.3875, q=±0.02, Bochner-only):**

| N | T | R | bochner_n | Ω (CLARABEL primal) | rigorous dual_LB | beats White? | wall |
|---|---|---|-----------|---------------------|------------------|--------------|------|
| 600 | 300 | 8 | 10 | 0.378481 | — | no (−5.2e-4) | 4.4 s |
| 800 | 400 | 9 | 12 | 0.378942 | — | no (−6.3e-5) | 2.9 s |
| **1000** | **500** | **10** | **14** | **0.379251** | **0.37925** | **YES (+2.45e-4)** | 4 s |

(For comparison, the archived ERD-9 `exact_integral_verification.json` shows the exact program at the
same center: N=500→0.378342, N=1000→0.379350 — consistent; the exact program converges *faster* than
the relaxed one, e.g. Ω_relaxed at N=1000 was only 0.376249.)

**Finding:** the smallest single-binding-center exact+Bochner instance that clears White's 0.379005
is **N≈900–1000, T≈500, bochner_n≈14**, with a rigorous-dual margin of **only +2.45e-4** *before
ellipse extension*. This pins the certificate size to formalize — and that is exactly where the
direction runs into trouble (§1.3).

### 1.2 The arbitrary-precision / rational route is genuinely feasible (medium effort)

- **SDPA-GMP is built + smoke-tested to 1e-75** (`lp_research_state/bin/sdpa_gmp`). The **only** missing
  piece is the cvxpy→SDPA-sparse serializer (`sdpa_gmp_wrapper.py` docstring flags it as deferred).
  This is shared with D1/L1 — *not* D5-specific.
- The cone is mixed and fully serializable: at N=1000 the compiled CLARABEL cone is
  **zero=21, LP-nonneg=7089, 40×SOC(size 3) + 2×SOC(size 502), 2×PSD(side 30)**; A is **9164×4083 with
  121,122 nonzeros**. SDPA-S handles LP+SOC+PSD blocks natively, so serialization is mechanical
  bookkeeping, not new math.
- **Exact rounding of an SDP dual is a published, named technique with working precedent:**
  de Laat–Dostert–Oliveira, *Exact semidefinite programming bounds for packing problems*
  (arXiv:2001.00256) — "round the floating-point SDP output to rationals (or a quadratic extension),
  without requiring strict feasibility, working for large problems." So step 2's *rounding* primitive
  is real and de-risked by literature.
- **TWO inherited obstructions from the L1 deep-dive (megamemory node
  `l1-davis-papp-rational-dual-snap-adversarial-assessment`), both of which hit D5 directly:**
  - **(O1) Boundary obstruction.** White's binding bound is an *optimum* with active constraints /
    complementary slackness (PRO-6 saturation) — i.e. the certified point sits on the **boundary** of
    the SOS/PSD cone, where a rational SOS/PSD decomposition **may not exist** (Scheiderer; Davis–Papp
    verbatim: "polynomials on the boundary may not have a rational SOS decomposition... interior ones
    do"). To get an *exact* certificate you must **back off to a strictly-interior, slightly-weaker
    constant** (certify µ ≥ c with c < c*). This is *fine for D5's purpose* (D5 already accepts a weaker
    constant) and the bit-size cost is only logarithmic in the back-off ε (Davis–Papp), so the back-off
    is affordable — **but it confirms D5's certified number is strictly below even the small-N Ω**, and
    must stay above White's 0.379005 after *both* the back-off and the ellipse extension. At N=1000 the
    pre-extension headroom is only +2.45e-4; this is the binding budget.
  - **(O2) No dual vector is extracted today.** `dual_extractor.py` parses only the dual *objective
    scalar* (5 sig figs) from CLARABEL's log — it does **not** emit a dual *vector*. So "snap the dual
    to rationals" has nothing to snap from the current pipeline; D5 needs the full dual vector from
    SDPA-GMP (⇒ the serializer must be built) or a re-solve. This re-confirms D5 is **not** independent
    of the D1 serializer / GMP-resolve infrastructure.

**Verified-bound shortcut (L2):** Jansson/VSDP a-posteriori interval bounds would upgrade the *current*
CLARABEL solve to a *verified* LB with no serializer and no GMP re-solve. That is strictly cheaper than
D5's GMP path for the "solver-independent" half of the payoff (gives interval-rational, not symbolic-
rational, but is proof-grade in the sense that the bound is *proven from the output*).

### 1.3 The Lean step — D5's distinguishing claim — is the weak link

D5's novelty over D1/D2/L1/L4 is "formalize the *finite* PSD/feasibility check in Lean." Adversarial reading:

- **Mathlib HAS `Matrix.PosSemidef`** (`Mathlib.LinearAlgebra.Matrix.PosDef`) and rational linear
  algebra. But `Matrix.PosSemidef` is defined as `IsHermitian ∧ ∀ x, 0 ≤ xᴴ M x` — it is **not
  decidable by `rfl`/`decide`**. To *prove* a specific 30×30 rational matrix is PSD you must exhibit a
  witness: a rational `LDLᵀ`/Cholesky with nonneg diagonal, *and* prove it reconstructs M (rational
  matrix-multiply equality — feasible but verbose) *and* invoke a PSD-from-LDLᵀ lemma. The two 30×30
  Bochner blocks are the *easy* part.
- **The hard part is everything around the PSD blocks.** A lower-bound certificate is **weak SDP
  duality**: a dual-feasible point whose objective is the bound. PRO-27 (`LEAN_LEMMA_INVENTORY.md`)
  found Mathlib has **zero** conic/SDP duality and **no** Bochner–Herglotz / Toeplitz-PSD theorem. So a
  Lean proof must either (a) re-derive weak duality for this specific mixed LP+SOC+PSD program from
  scratch, or (b) hand-verify the scalar inequality `⟨b, y⟩ ≥ target` together with dual feasibility
  `Aᵀy ⪯_K* c` across **all 9164 cone rows** — including the **7089 LP rows** that scale with N. That is
  a 9164-row rational feasibility certification, *not* "verify one explicit PSD matrix."
- **Then the ellipse extension** (the step that converts a single-center value into a bound on **µ**,
  not just on Ω at one (h,p,q)) must *also* be formalized: 7 concave-quadratic ellipses covering
  White's residual region (5.16), each with its own per-center dual. Mathlib can do quadratic
  inequalities, but this is a second, independent formalization of comparable size, and it is
  load-bearing for the word "µ".

**Net:** the finite certificate is real and exact-arithmetic-checkable (de Laat–Dostert precedent), but
it is **large (9164 rows), not "medium-size."** Magron's Coq nonlinear-optimization work
(`jfr.unibo.it/article/.../4319`) explicitly scopes formal SOS verification to *medium-size*
semialgebraic certificates; Flyspeck (HOL Light + Isabelle, person-years) is the precedent for pushing
large SDP/nonlinear certificates through a proof assistant. **The Lean half of D5 is months of
foundational work (build conic weak-duality in Mathlib) + a large bespoke certificate — not the small,
self-contained target the thesis implies.** The exact-arithmetic-rational half (skip Lean) is the
realistic deliverable.

---

## 2. PRIOR ART

### 2.1 Inside this project

- **`white_full_convex_exact.py` (ERD-9, 2026-05-31)** — the exact-(W.1) primal program D5 builds on.
  Already exists; was written to falsification-test the Step-E saturation theorem, *not* for rigor.
- **`exact_integral_verification.json`** — exact program already solved at the binding center for
  N∈{200,500,1000}; the N=1000 value 0.379350 already (numerically) beats White. D5's "find the smallest
  White-beating N" first step is **partially already done** (this memo completes it with dual-LB +
  cone-size data).
- **`path_b_*.py`** — three independent implementations of White's §5.1 ellipse extension, already
  validated to 10+ digits on the *relaxed* program. D5 would re-point them at the exact program (small
  change; same `find_ellipse_h_p(..., target=0.379005)` machinery).
- **Companion D1/D2 (`2026-06-03_whitespace_scout.md`)** — D5 is explicitly stated there as the
  "fast hedge if production-scale exact certification proves intractable," consuming D1's serializer
  and D2's Fejér-Riesz dual. **D5 is not standalone whitespace** — it is a *scoping choice* (small N,
  Bochner-only, weaker constant) layered on the D1/D2/L1/L4 rigor machinery.
- **PRO-27 `LEAN_LEMMA_INVENTORY.md`** — already assessed the Lean landscape: Mathlib has the Fourier
  primitives and `Matrix.PosSemidef` but **no SDP duality, no Bochner–Herglotz**; recommended *against*
  formalizing the general saturation theorem (12+ mo). D5's claim that a *specific finite* certificate
  dodges this is **partially true** (the PSD blocks are finite) **but partially false** (weak duality
  for the conic program and the ellipse extension are still missing/large).

### 2.2 External (citation-grounded)

- **de Laat, Dostert, Oliveira, arXiv:2001.00256** — exact SDP rounding for packing bounds; the
  canonical "snap the dual to rationals, verify exactly, no strict feasibility needed, scales" precedent.
  This is the realistic backbone of D5 step 2.
- **Davis & Papp, arXiv:2105.11369 / 2305.19039** — rational dual certificate from the *dual* with a
  *provable bit-size bound*, univariate-on-interval case (= White's dual after x=cos θ). Lighter than
  SDPA-GMP; the bit-size bound is exactly the "is the certificate checkable" guarantee D5 needs.
- **Magron et al., arXiv:2202.06544** — exact SOHS for trigonometric univariate polynomials (the
  Fejér-Riesz finish for the Bochner block, companion-L4/D2).
- **Magron, *Formal Proofs for Nonlinear Optimization* (JFR)** — Coq verification of SOS certificates,
  scoped to *medium-size*; the honest ceiling on "formally checked" for this kind of certificate.
- **Flyspeck (Hales et al., arXiv:1501.02155)** — formal proof of Kepler via HOL Light + Isabelle;
  precedent + warning for the person-year cost of large nonlinear/SDP certificates in a proof assistant.
- **No competing *proven* LB on µ exists.** Confirmed again (external scan #4): post-White movement is
  UB-side only (AlphaEvolve 0.380924, Together 0.380871). The LB-as-theorem whitespace is genuinely open.

---

## 3. CONCRETE PLAN (ordered)

**Phase A — fix the certificate (cheap, mostly done).**
1. ✅ (this memo) Pin the smallest White-beating exact+Bochner instance: **N≈1000, T=500, R=10,
   bochner_n=14**, single binding center, rigorous dual_LB 0.37925 (+2.45e-4). Cone: 9164 rows.
2. Re-point `path_b_analytical.find_ellipse_h_p` at the exact program; confirm the 7 ellipses at
   the exact-program duals still cover region (5.16) at target 0.379005 with **positive** margin
   (the +2.45e-4 single-center headroom must survive extension — see Failure Mode F1). Decide the
   *honest* theorem constant (likely µ ≥ 0.379x, x < 5, after extension shrinkage), or push N to
   1500–2000 to buy headroom.

**Phase B — solver-independent rational certificate (medium; the realistic payoff).**
3. Build the cvxpy→SDPA-S serializer (shared with D1; LP+SOC+PSD bookkeeping) **or** take the L1
   route (snap CLARABEL's dual directly via Davis–Papp, no serializer).
4. Re-solve the binding center at GMP precision (SDPA-GMP) → high-precision dual.
5. Snap the dual to rationals (de Laat–Dostert) and verify dual feasibility + objective ≥ target in
   **exact arithmetic** (sympy/`python-flint`). **Deliverable: a checkable rational certificate of
   µ ≥ (constant) > 0.379005, independent of any FP solver.** This alone is payoff-class-1.
6. (Optional, cheapest validity upgrade) Run L2 Jansson/VSDP on the existing CLARABEL solve as an
   independent verified cross-check of step 5's constant.

**Phase C — Lean formalization (HIGH effort; optional, the only part that yields "machine-checked").**
7. Formalize the two 30×30 rational Bochner PSD blocks via rational `LDLᵀ` (tractable, ~days).
8. Formalize weak duality for the specific mixed LP+SOC+PSD program **OR** the explicit 9164-row
   dual-feasibility scalar inequality chain. **This is the wall**: build conic weak-duality in Mathlib
   (absent) or a very large bespoke `decide`-style certificate. Realistic estimate: **months**.
9. Formalize the 7-ellipse covering of (5.16). Second large piece.
10. Only if 7–9 land: a *machine-checked* µ ≥ (constant > 0.379005) theorem.

**Recommended stopping point: end of Phase B.** That delivers a genuine payoff-class-1 result
(solver-independent rational certificate beating White) at medium effort. Phase C is where D5's
distinguishing claim lives, and it is not bounded/small.

---

## 4. PAYOFF

- **Bar class: 1** (proof-grade / solver-independent), with an *aspirational* "machine-checked" upgrade
  if Phase C completes.
- **"Truly meaningful" outcome (Phase B):** µ ≥ ~0.3793 (or higher if N pushed) as a **theorem with a
  checkable rational dual certificate**, no SDP solver in the trusted base, **and decoupled from the
  poly-moment caveat** (Bochner-only) and the thin-margin production machinery. This is a clean,
  defensible, *publishable* "first proof-grade improvement over White" — modest in constant but
  unimpeachable in rigor. That genuinely clears Ben's bar item (1).
- **Phase C upside:** the first *machine-checked* lower bound on the Erdős minimum-overlap constant.
  High prestige, but gated on building conic SDP duality in Mathlib (a Mathlib-contribution-sized
  project in its own right).

**Caveat on the constant:** D5 deliberately yields a *weaker* number than the 0.380284 production
bound. Its value is rigor, not size. If the production-scale rational certification (D1/L1 at N=24000
with poly-moment) turns out tractable, D5's weaker-constant theorem is partly redundant. D5's real role
is the **hedge**: guaranteed proof-grade output even if full-scale exact certification fails.

---

## 5. FIRST-STEP PROBE — RESULT (run this session)

Already reported in §1.1. Summary of the finding:

> The exact-(W.1)+Bochner program (no poly-moment), single binding center, **first beats White's
> 0.379005 at N≈900–1000 / T=500 / bochner_n=14**, with rigorous dual_LB **0.37925 (+2.45e-4)** *before*
> ellipse extension. The certificate at that size is **9164 cone rows / 4083 vars / 121k nnz**
> (7089 LP + 42 SOC + 2 PSD(30)). So the exact-arithmetic certificate is **feasible but large**; the
> "verify one PSD matrix" framing understates it by ~3 orders of magnitude in row count.

This is a real, decision-relevant finding: it *confirms feasibility of the exact-arithmetic rational
route* (the constant is reachable at tractable N) while *refuting the implied small-Lean-target* (the
certificate is large and the surrounding weak-duality + ellipse-extension formalization is missing).

---

## 6. FAILURE MODES (adversarial)

- **F1 — Ellipse-extension eats the margin (most likely).** The +2.45e-4 is a *single-center* value;
  White's ellipse extension propagates it outward and the dual objective *decreases* off-center. The
  binding-row margin after extension could fall below 0 at N=1000, forcing larger N (1500–2000) just to
  state a White-beating *full-µ* theorem — which inflates the certificate further. `findings.md`
  repeatedly warns the "MIN-over-rows-AND-ranges" rigorous Δ is structurally near 0 without headroom.
  *Mitigation:* push N to buy headroom; accept a thinner constant; or certify only at the binding center
  and lean on White's *published* ellipse geometry (but that re-imports a White number into the bound,
  contradicting the "standalone" goal).
- **F2 — Lean weak-duality is a Mathlib project, not a lemma.** §1.3: conic SDP duality is absent from
  Mathlib (PRO-27). Phase C step 8 is the true cost and is unbounded at the "specific finite certificate"
  framing. *Mitigation:* stop at Phase B (rational certificate, exact-arithmetic-verified) — still
  proof-grade, skip "machine-checked."
- **F3 — Large rational certificate has high bit-height.** 9164 rows snapped to rationals near a thin
  margin can yield large numerators/denominators (bit size ∝ 1/distance-to-boundary). *Mitigation:*
  Davis–Papp's bit-size bound + re-solving off the margin (larger N) keep height controlled; de
  Laat–Dostert explicitly handle large, non-strictly-feasible instances.
- **F4 — `optimal_inaccurate` at the binding center.** N=1000 already returns `optimal_inaccurate` from
  CLARABEL; the rigorous dual_LB (0.37925) is the project's log-parse trick, *not yet* a certificate.
  The whole point is to replace it — but if SDPA-GMP also struggles at this conditioning, the GMP dual
  may need a slightly perturbed (interior) center, marginally weakening the constant. *Mitigation:* L2
  Jansson/VSDP is purpose-built for exactly the `optimal_inaccurate`/degenerate regime.
- **F5 — Redundancy.** If D1/L1 certify the full 0.380284 at production scale, D5's weaker constant is
  subsumed. D5 only "wins" as the hedge. *Mitigation:* run D5's Phase B in parallel as cheap insurance;
  it shares all infrastructure (serializer/dual-snap) with D1/L1, so marginal cost is low.
- **F6 — Boundary obstruction stacks with thin margin (inherited from L1, see §1.2 O1).** The certified
  point is a cone-boundary optimum, so an exact rational certificate requires backing off to a strictly
  interior c < Ω. The certified µ must then clear White's 0.379005 after BOTH the interior back-off AND
  the ellipse-extension shrinkage — two debits against a pre-extension headroom of only +2.45e-4 at
  N=1000. *Mitigation:* push N to 1500–2000 (cheap; the exact program converges fast — N=1000 already
  0.379251) to widen the headroom before spending it on back-off + extension. This is the single
  most important quantitative check before committing to a specific theorem constant.

---

## 7. VERDICT

**MAYBE (lean pursue, but DEMOTE the Lean claim and FOLD into D1/L1).**

The exact-arithmetic rational-certificate half (Phase B) is a genuine, feasible, payoff-class-1 result
and a sound *hedge* — and the first-step probe confirms a White-beating exact+Bochner instance is
reachable at tractable N≈1000. But D5 is **not independent whitespace**: it shares its entire rigor
backbone (SDPA-GMP serializer / Davis–Papp dual-snap / exact rounding) with D1, L1, L4, and its only
distinguishing feature — the "small finite Lean check" — is **overstated**: the certificate is 9164
rows (not "one PSD matrix"), and the surrounding conic weak-duality + 7-ellipse extension are missing
from Mathlib and are months of work (F2). **Recommendation:** pursue **Phase B as a deliberately-small,
Bochner-only, single-binding-center instance of D1/L1** (the cleanest, poly-moment-free, exact rational
certificate beating White), treat it as the safe payoff-class-1 floor of the rigor program, and **defer
Lean (Phase C) indefinitely** unless someone independently lands conic SDP duality in Mathlib. As a
standalone "machine-checked theorem" pitch, **drop**; as "the smallest clean exact rational certificate
that beats White," **pursue (inside D1/L1).**
