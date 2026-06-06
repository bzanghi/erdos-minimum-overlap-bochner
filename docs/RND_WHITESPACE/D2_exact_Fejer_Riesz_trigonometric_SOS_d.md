# D2 — Exact Fejér–Riesz / Trigonometric-SOS Dual Certificate for the Binding Center

**Author:** Claude (machine-assisted), whitespace deep-dive workflow — 2026-06-03
**Direction under assessment (Phase-1 proposal D2):** At the binding center, extract the SDP dual on the Bochner/Toeplitz PSD block (a nonnegative trigonometric polynomial = sum of Hermitian squares by Fejér–Riesz), then prove its nonnegativity EXACTLY via a rational LDLᵀ / rational Fejér–Riesz factorization — yielding a closed, human/Lean-checkable nonnegativity certificate with no SDP solver in the trusted base.
**Payoff class claimed:** 1 (proof-grade / solver-independent / formally-checkable LB), partial 2.
**This memo SCOUTS and PLANS. It certifies nothing. No new bound is claimed proven.**

---

## TL;DR verdict: **MAYBE → pursue, but ONLY as a sub-step of D1/L-track, and re-scoped.**

D2's *literal* thesis — "the certificate is a rational Fejér–Riesz factorization of the Bochner dual block" — is **necessary but far from sufficient**, and the way it is phrased is misleading about where the work and the risk live. Three findings from this deep-dive:

1. **The Bochner block is only ONE cone among many.** The binding-center dual lives over a mixed LP + SOC + (two small) PSD cones. A genuine solver-independent certificate requires the **entire dual vector** to satisfy dual-feasibility/stationarity **exactly** — not just PSD-ness of the Bochner block (which the solver already guarantees). Fejér–Riesz handles ~390 of the dual's components; it does nothing for the other ~1500 (at a tiny N=300 instance) / ~100k+ (at production N). **D2 in isolation is not a certificate.**
2. **The favorable structural facts are real and were verified here.** The Bochner dual is genuinely low-rank: at production `bochner_n=20` the f≥0 block has complex-Hermitian rank **5** (real-form rank 10 of dim 42), and the 1−f≥0 block is **identically inactive** (rank 0). Low rank ⇒ the SOHS has few terms ⇒ the Fejér–Riesz piece itself is cheap. **And the published machinery exists** (Magron exact-SOHS; Davis–Papp rational dual certificates with bit-size bounds) — D2 is not inventing an algorithm.
3. **The thin-margin failure mode is real and biting NOW.** A naive rational snap of the bn=20 Bochner dual gives `min_eig = −2.1×10⁻⁴` at denominator 10⁴ — **larger than the ~5×10⁻⁵ binding margin.** Davis–Papp's theory says exactly why: rational dual certificates exist for cone-*interior* points with bit size ∝ 1/distance-to-boundary; the binding center is, by definition, on the active boundary. So D2-as-stated is at material risk unless the center is first hardened (N≥24000) AND a *principled* certificate construction (not naive rounding) is used.

**Net:** D2 is best understood as the *Fejér–Riesz finishing move* inside the D1 / external-L-track rigor program, not a standalone direction. The cleaner, lower-risk path to the SAME payoff-class-1 destination is **L1 (Davis–Papp rational dual certificate from the dual, no Gram round, no serializer)** + **L4 (Magron exact-SOHS for the Bochner block)** + **L2 (Jansson/VSDP verified a-posteriori bound) as the fast first rung.** D2's distinctive value is the explicit, checkable Fejér–Riesz/SOHS form for the (small, low-rank) PSD block — worth producing, but only after the surrounding cone is handled.

---

## 1. FEASIBILITY (tractable with available tools? at what effort?)

### 1.1 The cone structure — verified by direct probe

Built the binding-center program (`white_full_convex.build_problem`, row4 region center `h=0.004, p=0.3875, q∈[−0.02,0.02]`) and inspected the constraint cones:

| Instance | PSD blocks | Inequality (LP box + cell-env SOC) | Equality | Total scalar dual components |
|---|---|---|---|---|
| N=300, T=120, R=6, bn=6 | 2 (dim 14 each) | 63 | 1 | **1890** |

At production scale (N≈10000–24000, T=4000, R=10, bn=20) the LP part alone is ≈ 2N box + 2N nonneg ≈ 40k–96k scalar duals, plus 2R cosine-SOC + 4R sine two-sided rows, plus the two 42×42 Bochner blocks. **The Bochner blocks are a vanishingly small fraction of the dual.**

**Consequence (the central scoping correction to D2):** "prove the Bochner block's nonnegativity by rational Fejér–Riesz ⇒ a closed certificate" is a *non sequitur* as stated. PSD-ness of the dual Bochner block is **automatic** from any interior-point solver (the dual iterate is in the cone). The thing that actually certifies the bound is the **dual-feasibility system**: a nonnegative combination of ALL constraints whose aggregate is ≤ the objective. Fejér–Riesz/SOHS is the tool for re-expressing the PSD-block contribution in a checkable form, but the certificate is the *whole assembled dual*, verified exactly. This is exactly the gap that D1 (serialize + GMP re-solve + rational snap of the full dual) or L1 (Davis–Papp full rational dual certificate) close — and D2 cannot close it alone.

### 1.2 The favorable facts — verified by direct probe

**Bochner dual rank (binding center):**

| bochner_n | f≥0 block: dim / real-form rank / complex-Hermitian rank | 1−f≥0 block |
|---|---|---|
| 6  | 14 / 2 / **1** | rank 0 (inactive) |
| 10 | 22 / 2 / **1** | rank 0 (inactive) |
| 20 | 42 / 10 / **5** | rank 0 (inactive) |

The dual eigenvalues come in equal pairs (`4.769e-3, 4.769e-3, 1.855e-3, 1.855e-3, …`) — the signature of the real embedding `[[Re,−Im],[Im,Re]]` of a complex Hermitian matrix, so the *complex* rank is half the real-form rank. **Complex-Hermitian rank 5 at bn=20 is genuinely low.** A rank-5 Toeplitz-dual Gram ⇒ a sum of ≤ 5 Hermitian squares ⇒ a small Fejér–Riesz object. The 1−f block being identically inactive means only ONE trigonometric SOHS is needed, not two. **This part of D2 is easy.**

### 1.3 The thin-margin obstruction — verified by direct probe

Snapping the bn=20 f≥0 dual matrix `Z` to rationals (naive `round(Z·den)/den`) and re-checking PSD-ness:

| denominator | min eigenvalue after rounding | ‖Z − Z_rounded‖_F |
|---|---|---|
| 10    | 0.0       | 7.4e-3 |
| 100   | 0.0       | 7.4e-3 |
| 1000  | **−1.0e-3** | 7.0e-3 |
| 10000 | **−2.1e-4** | 9.6e-3 |

The non-monotonicity is a rounding artifact, but the message is robust: **naive rational rounding pushes the dual matrix out of the PSD cone by O(1e-4..1e-3)**, which is *larger* than the binding margin (R6 +2.5e-5, R17 +5.1e-5, core exactly 0.3802838). A negative-eigenvalue correction of 2e-4 added back as slack would **eat the entire LB improvement and then some.** This is precisely D2's own stated failure mode, now quantified and confirmed to bite at the current margin.

The fix is NOT to round harder; it is to use a **principled construction** (Davis–Papp produces a rational certificate that is provably IN the cone, not a rounded-then-corrected one; Magron uses perturbation+exact-compensation so the output is exactly PSD by construction) AND to **harden the center to N≥24000 first** to move off the active boundary (Davis–Papp bit size ∝ 1/dist-to-boundary; findings.md already recommends this hardening).

### 1.4 Tooling adequacy

| Tool | Adequate for D2? | Note |
|---|---|---|
| **mpmath / sympy** | YES | exact LDLᵀ / Fejér–Riesz factorization of a rank-5 Toeplitz dual; rational arithmetic for dual-feasibility verification — both well within mpmath/sympy. |
| **`white_full_convex_exact.py`** | YES (important) | already builds the exact-(W.1)-integral primal, removing the only non-exact primal ingredient (cell-min approx), so an exact dual is *meaningful* — otherwise the primal slack contaminates the certificate. |
| **SDPA-GMP** | YES, but only via the missing serializer (D1) | needed to get a *high-precision dual* to snap; CLARABEL's FP dual is too coarse (see 1.3). The cvxpy→SDPA-S serializer is the one engineering gap (D1). |
| **CLARABEL dual** | Marginal | usable to *prototype* (as done here) but its FP precision (`optimal_inaccurate`, residual ~1e-7) is the reason naive snapping fails; not adequate for the final certificate. |
| **Lean/Mathlib** | Partial (hedge only) | Mathlib HAS `Matrix.PosSemidef` and rational linear algebra, so verifying ONE explicit rational PSD Gram + ONE finite dual-feasibility inequality system is in-scope (this is the D5 dovetail). It LACKS the general Bochner-Herglotz / SDP-duality theory (PRO-27 / LEAN_LEMMA_INVENTORY), so only a *specific finite certificate check* is formalizable, not the framework. |

**Effort estimate:** **Medium-high, and dominated by D1 (the serializer + full-dual rational snap), not by the Fejér–Riesz step itself.** The Fejér–Riesz/SOHS step on a rank-5 block is ~1–3 days with mpmath or a from-scratch port of Magron's root-isolation variant. The surrounding work — serialize the mixed cone to SDPA-S, GMP-re-solve, rational-snap and exactly verify the **entire** dual (LP + SOC + PSD), handle the thin margin — is the real cost (≈ 2–4 weeks), and it is **shared with D1/L1**, not unique to D2.

---

## 2. PRIOR ART

### 2.1 Inside this repo — has any piece been tried?

- **NO direct attempt.** `grep -iE "fejer|riesz|sum.of.hermitian|rational.*gram|exact.*dual|SOS"` over `findings.md` returns **zero** Fejér-Riesz / SOS / rational-certificate hits. The exact-SOS-dual idea is genuinely unexplored here. ✅ whitespace.
- **The companion scout memos (same day) already name this direction.** `docs/RND_WHITESPACE/2026-06-03_whitespace_scout.md` lists it as **D2** verbatim and flags that it **consumes D1's GMP dual** ("D2/D3/D5 all consume its serializer/GMP dual; build D1 first"). The companion **external** memo (`2026-06-03_external_literature_scan.md`) supplies the missing algorithms (L4 = Magron, L1 = Davis–Papp) and the lighter alternatives (L2 = Jansson/VSDP). **This deep-dive's net-new contribution over those memos is the empirical cone/rank/snap probe in §1, which converts D2's hand-waved "Fejér–Riesz applies cleanly" into a quantified, partly-cautionary picture.**
- **Lasserre (RETRACTED)** is the cautionary sibling: it truncated a moment expansion with no tail bound. D2 is the *opposite* case (the Bochner block is finite-dimensional and exactly PSD-representable, so no tail issue) — this asymmetry is correctly stated in the proposal and confirmed here.
- **`white_full_convex_exact.py` exists** and is the right primal substrate (exact (W.1) integral). ✅ ready.
- **PRO-6 / LEVER_I_PRIME_THEOREM**: the dual/KKT structure has been studied (3 scalar cell-envelope shadow prices ξ, τ, ν₃; the strict-form complementarity with the Bochner dual matrix Z is *open*). Relevant because the full dual-feasibility system D2 must verify is parameterized partly by these.

### 2.2 External literature — the load-bearing papers (verified to exist & be on-point)

- **Magron, Safey El Din, Schweighofer, Vu — "Exact SOHS decompositions of trigonometric univariate polynomials with Gaussian coefficients," ISSAC 2022, [arXiv:2202.06544](https://arxiv.org/abs/2202.06544).** Confirmed: three hybrid numeric→symbolic algorithms computing weighted-sum-of-Hermitian-squares decompositions for trigonometric univariate polynomials positive on the unit circle, **built on Riesz–Fejér**, with Gaussian-integer coefficients and **polynomial bit-complexity**. This is *literally* the exact Fejér–Riesz engine D2 needs for the Bochner block. **De-risks D2's "is the rational Gram low-height enough" worry with a proven output-size bound.**
- **Davis & Papp — "Rational Dual Certificates for WSOS Polynomials with Boundable Bit Size," J. Symbolic Comput. 2023, [arXiv:2305.19039](https://arxiv.org/abs/2305.19039)** (foundations: [arXiv:2105.11369](https://arxiv.org/abs/2105.11369), SIAM J. Opt. 2022). Confirmed: rational nonnegativity certificates **from the dual cone**, with bit size bounded by degree, #vars, and **distance from the cone boundary**; featured special case is **univariate polynomials nonnegative over a bounded interval in several bases** — exactly White's cosine-substituted test family. **This is a strictly better route than D2's naive "snap the Gram": it certifies from the dual the project already extracts, with no rounding/projection, and the bit-size-∝-1/dist-to-boundary result is precisely the §1.3 obstruction made rigorous.**
- **Jansson, Chaykin & Keil — verified SDP bounds (SIAM J. Numer. Anal. 2007) + VSDP-2012.** The *fast* a-posteriori verified-LB route (handles LP+SOC+SDP cones natively, built for the `optimal_inaccurate` regime). Not Fejér–Riesz, but the natural **first rung** under D2 (verify before you symbolically certify).
- **Cohn–Elkies / Cohn–Triantafillou (sphere packing, [arXiv:2206.09876](https://arxiv.org/abs/2206.09876))** — the established "certify the numerically-found dual witness directly in exact/interval arithmetic" culture; the methodological template for the whole D2/L-track.

**Verdict on prior art:** D2 is genuine whitespace *for this project* (never tried), the enabling algorithms are **published, peer-reviewed, and confirmed on-point**, and the companion memos already correctly position it as a sub-step of D1. The deep-dive adds the quantitative cone/rank/snap evidence.

---

## 3. CONCRETE PLAN (ordered)

**Pre-req (this is really D1 / L-track; D2 cannot start without a precise dual):**
0. Build the cvxpy→SDPA-S serializer for the mixed LP+SOC+PSD cone (D1 keystone), OR implement the Jansson/VSDP a-posteriori verified-LB post-processor (L2) as the lighter first rung. Either gives a *trustworthy* dual / verified LB to work from. (CLARABEL's FP dual is demonstrably too coarse — §1.3.)

**D2 proper (the Fejér–Riesz finishing move):**
1. **Harden the binding center to N≥24000** (and use `white_full_convex_exact.py` for exact (W.1)) so the active trigonometric polynomial sits in the cone *interior* with a margin > the rational-certificate slack. Re-extract the dual at GMP precision (via step 0). *Rationale: Davis–Papp bit size ∝ 1/dist-to-boundary; §1.3 shows the FP-margin snap fails by ~2e-4.*
2. **Isolate the active face.** Confirm (as §1.2 found at small N) that the 1−f≥0 Bochner block is inactive and the f≥0 block is low-rank (complex rank ~5 at bn=20). Reduce the certificate to that active face — drop inactive cones from the exact-verification system to keep bit size small. *Mitigation for D2's "too high-dimensional for low-height rational Gram" failure mode.*
3. **Exact SOHS of the Bochner block.** Run Magron's algorithm (arXiv:2202.06544, root-isolation variant — ~1 file, mpmath/sympy) on the GMP dual block to get an **exact** weighted-sum-of-Hermitian-squares with Gaussian-integer coefficients; check the bit size against the paper's predicted bound. *This is the literal Fejér–Riesz certificate for the PSD piece.*
4. **Rational dual certificate for the FULL dual.** Apply Davis–Papp (arXiv:2305.19039) to assemble a rational dual-feasible point over the **entire** cone (LP multipliers + SOC multipliers + the SOHS Bochner contribution). Verify dual feasibility / weak-duality inequality in **exact rational arithmetic** (sympy/`python-flint`). The objective of this exact dual point is the certified LB. *This — not step 3 alone — is the certificate.*
5. **Lift to a bound on µ (do NOT skip).** A single-center exact dual certifies the SDP optimum *at that center*, which is NOT yet a bound on µ. Run the `path_b_*` ellipse-extension (or the adaptive grid+Lipschitz cover) with the exact per-center duals to certify coverage of the residual region — the CLAUDE.md "critical caveat." For a *fully* exact µ-bound the cover/Lipschitz step must also be done in interval/rational arithmetic (this is the Cohn–Triantafillou "discrete reduction" discipline; real additional work).
6. **(Hedge / D5 dovetail) Bochner-only, small-N, Lean.** If the full poly-moment-augmented certificate is too large, certify the **Bochner-only** program (no poly-moment cuts — sidesteps load-bearing caveat #2) at the smallest N that still beats White's 0.379005 (augmented binding-center solves already sit ~0.3796–0.3803, well above White — `findings.md` "0.379653" recipe). Formalize THAT finite rational PSD + dual-feasibility check in Lean (`Matrix.PosSemidef` + rational linear algebra are in Mathlib). Output: a *machine-checked* theorem µ ≥ (something > White), decoupled from the thin-margin production machinery.

---

## 4. PAYOFF

**Bar class: 1 (proof-grade / solver-independent / formally-checkable certificate of a lower bound), with a partial-2 upside.**

- **If steps 1–5 succeed** at the production binding center: µ ≥ ~0.3803 becomes a **theorem** whose trusted base is exact rational arithmetic (and Fejér–Riesz / SOHS), with **no SDP solver and no floating point in the trusted base** — directly the stated #1 payoff ("making µ ≥ ~0.3803 a theorem"). This is the single most valuable outcome the workflow is chasing.
- **If only the hedge (step 6) succeeds:** a clean, possibly **Lean-machine-checked** theorem µ ≥ (a constant > White's 0.379005), Bochner-only, no poly-moment caveat. Strictly weaker constant, but unambiguously "truly meaningful" by the bar (the word *theorem*, even *machine-checked theorem*, on a White-beating bound).
- **Partial class 2:** the exact SOHS decomposition of the binding-center Bochner dual could expose a clean analytic structure (e.g. a small fixed set of Hermitian squares with recognizable coefficients), feeding the saturation/KKT-with-Z program (D4 / PRO-6). The complex rank being only ~5 at bn=20 makes this plausible — a rank-5 Toeplitz dual is a candidate for a closed-form SOHS.

**What it does NOT do:** it does not push past the conjectured framework ceiling C_∞ ≈ 0.380558, and does not touch the open gap from the UB side. D2 is a *rigor-upgrade* of the existing bound, not a *stronger* bound (class 3). That is appropriate — the bar explicitly counts proof-grade rigor as a high payoff.

---

## 5. CHEAP FIRST-STEP PROBE — RUN, with findings

Three probes were run (all small, no heavy solves):

**(a) Cone structure** (N=300,T=120,R=6,bn=6): 2 PSD blocks (dim 14), 63 inequality cones, 1 equality, **1890 scalar dual components total**. ⇒ *the Bochner block is a tiny fraction of the dual; Fejér–Riesz alone cannot be the certificate.*

**(b) Bochner dual rank at the binding center:** f≥0 block complex-Hermitian rank **1** (bn≤10) → **5** (bn=20); 1−f≥0 block **identically inactive** (rank 0) at all levels. ⇒ *the Fejér–Riesz/SOHS piece is genuinely small and low-rank — the easy part of D2 is confirmed easy, and only one trigonometric SOHS is needed.*

**(c) Rational-snap stress test (the crux):** naive rounding of the bn=20 f≥0 dual to denominator 10⁴ gives `min_eig = −2.1×10⁻⁴`, **exceeding the ~5×10⁻⁵ binding margin.** ⇒ *D2's own stated failure mode is real and biting now; naive snapping is insufficient, confirming the need for (i) N≥24000 hardening to widen the margin and (ii) a principled rational-certificate construction (Davis–Papp / Magron) rather than rounding.*

**First-step finding (one line):** the Fejér–Riesz piece is small and low-rank (✅), but the binding-center dual is a 1890+-component mixed-cone object whose naive rational snap violates PSD by ~2e-4 ≫ the 5e-5 margin (⚠) — so D2 is viable only as the *finishing move* of a full-dual rational-certificate pipeline (D1/L1) on a *hardened* (N≥24000) center, not as a standalone snap of the Bochner block.

---

## 6. FAILURE MODES (honest)

1. **Scope error (most likely "soft fail").** Treated literally — "rational Fejér–Riesz of the Bochner block ⇒ certificate" — D2 **does not produce a certificate**, because the bound is certified by the *full* mixed-cone dual, not the PSD block. *Mitigation:* re-scope as the SOHS finishing step inside D1/L1; that is what this memo recommends.
2. **Thin-margin / boundary blow-up (quantified, real).** The binding center is on the active boundary; Davis–Papp bit size ∝ 1/dist-to-boundary, and §1.3 shows naive snapping fails by ~2e-4 ≫ 5e-5. *Mitigation:* harden to N≥24000 first (widens margin; already recommended in findings.md), and use principled construction not rounding. *Residual risk:* even at N=24000 the margin may stay thin enough that the rational certificate has large bit size — checkable but ugly, and possibly not Lean-tractable.
3. **The serializer (D1) is the real gate.** D2 cannot start from CLARABEL's FP dual (§1.3). If D1's cvxpy→SDPA-S serializer for the mixed cone proves harder than estimated (SOC + two PSD blocks + 100k LP rows bookkeeping), D2 slips with it. *Mitigation:* L2 (Jansson/VSDP) gives a verified LB from the existing solve as a cheaper first rung, and the Bochner-only hedge (step 6) needs a much smaller serialization.
4. **µ-lift still required and still partly non-exact.** Even a perfect exact single-center certificate is not a bound on µ without the ellipse/cover step, which is itself currently Lipschitz-numeric. A *fully* exact µ-bound needs that step in interval/rational arithmetic too (real extra work; the Cohn–Triantafillou discipline). *Mitigation:* the hedge (Bochner-only small-N + ellipse extension done in interval arithmetic at small N) is the contained version.
5. **Poly-moment entanglement.** The production 0.380284 bound is load-bearing on poly-moment cuts. Certifying the *full* augmented dual means certifying those too (a Hausdorff/IBP-tail argument), which is more than Fejér–Riesz. *Mitigation:* the Bochner-only hedge sidesteps this entirely (at a weaker, still-White-beating constant).
6. **SOHS analytic structure may be uninformative (kills the partial-class-2 upside).** The rank-5 dual might have no recognizable closed coefficients (consistent with the NEGATIVE wide-basis PSLQ result that µ has no clean closed form). *Mitigation:* none needed — the class-1 rigor payoff stands regardless; class-2 was only a bonus.

---

## 7. Verdict

**MAYBE (lean pursue) — but re-scoped: pursue D2 as the Fejér–Riesz/SOHS finishing move of the D1 / external-L-track rational-dual-certificate pipeline, on an N≥24000-hardened center, NOT as a standalone snap of the Bochner block.** The enabling algorithms are published and on-point (Magron exact-SOHS; Davis–Papp rational dual certificates), the favorable structure is verified (Bochner dual complex-rank ~5, 1−f block inactive), but the standalone framing is a non sequitur and the thin-margin snap fails by ~2e-4 ≫ 5e-5 today. Sequence: **L2 (verified a-posteriori) → D1/L1 (full rational dual certificate) → D2/L4 (exact SOHS of the low-rank Bochner block) → ellipse-lift in interval arithmetic**, with the **Bochner-only small-N Lean hedge** as the contained fallback that still beats White.

---

### Provenance
- Probes: `white_full_convex.build_problem` at row4-region center, CLARABEL, bn∈{6,10,20} (cone tally, dual rank, rational-snap stress) — run 2026-06-03, light configs, no heavy solves.
- Repo inputs: `lp_research_state/code/white_full_convex.py`, `white_full_convex_exact.py`, `dual_extractor.py`, `sdpa_gmp_wrapper.py`; `docs/archive/PRO6_COMPLEMENTARITY_PROOF.md`, `LEAN_LEMMA_INVENTORY.md`, `FULLSPACE_RIGOR_MEMO.md`; `docs/RND_WHITESPACE/2026-06-03_whitespace_scout.md` (D1–D5) and `2026-06-03_external_literature_scan.md` (L1–L4); `findings.md`; megamemory graph.
- External (verified to exist & on-point via web search): arXiv:2202.06544 (Magron, ISSAC 2022); arXiv:2305.19039 + 2105.11369 (Davis–Papp, J. Symb. Comput. 2023 / SIAM J. Opt. 2022); Jansson/VSDP (SIAM J. Numer. Anal. 2007); arXiv:2206.09876 (Cohn–Triantafillou).
