# R&D Whitespace Scout — 2026-06-03

**Author:** Claude (machine-assisted), whitespace-scout workflow
**Scope:** Scout genuinely unexplored directions toward THE BAR (proof-grade rigor / clean analytic theorem / substantially stronger bound). This file SCOUTS and PLANS; it does NOT certify any new bound.

**Current state (verbatim from session brief + findings.md + megamemory):**
- Full-space LB µ ≥ 0.380284 (0.3802838) over White's entire (E(M), c1, d1) space, binding = core, NO White number in the bound. Numerically certified (CLARABEL + dual extraction + adaptive grid+Lipschitz cover + Bochner-PSD + poly-moment Hausdorff cuts).
- UB µ ≤ 0.380871 (Together, March 2026). Open gap 5.87e-4.
- Conjectured SDP framework asymptotic ceiling C_∞ ≈ 0.380558.

**The defining problem for THE BAR:** the existing bound is *numerically certified with load-bearing caveats*, not *proven*. The caveats (from megamemory + FULLSPACE_VERIFICATION.md + FULLSPACE_RIGOR_MEMO.md) are, verbatim in spirit:
1. CLARABEL is a floating-point interior-point solver; the LB is `reported − last_gap` parsed from its log. No solver-independent certificate.
2. Load-bearing on poly-moment cuts (rigorous only since the 2026-05-22 tail-bound fix) and on a handful of fresh dual centers (pm_k_max=20, N=20000).
3. Deep-q/high-p infeasibility exclusions are **solver-attested only** — "CLARABEL 'infeasible' at multiple interior points with NO Farkas/dual-ray certificate (rigor = empirical-robust-multipoint)." They are argued not-load-bearing, but they are the weakest rigor seam.
4. Thin margins: core exactly 0.3802838; thinnest regions R6 +2.5e-5, R17 +5.1e-5, R16 +7.1e-5.

The whitespace that clears THE BAR is therefore overwhelmingly on the **rigor-upgrade** axis (convert the existing bound into a theorem) and the **analytic-structure** axis (saturation/KKT theorems), not on chasing +1e-5 of numeric.

---

## Ruled-out / tried ledger (compact — do NOT re-propose as new)

| Lever | Outcome | Why (1-line) |
|---|---|---|
| Bochner-PSD moment matrix M_n(f), M_n(1−f) ⪰ 0 | **WORKED** (core lever, +~2e-3 at scale; White-validated by email) | necessary cond. for f∈[0,1]; valid SDP tightening |
| Poly-moment Hausdorff cuts m_2k ≥ −tail_bound_k | **WORKED** (after 2026-05-22 tail fix; load-bearing for 0.380284) | exact IBP tail bound `|α_j^(k)|≤4k/(π²j²)`; pre-fix was an overclaim |
| CDE cover-refinement / adaptive grid+Lipschitz dual cover | **WORKED** (+7.6e-5 from placement alone; the wide-region certifier) | max-of-more-ellipses tightens; adaptive subdivision recovers grid_min |
| Lasserre level-2 (truncated moment lift) | **RETRACTED** | truncated `(f²)̂(m)` with no Fejér-Riesz tail bound → overclaim |
| Wide-basis PSLQ closed-form hunt (50-digit, 33 consts) | **NEGATIVE (near-definitive)** | µ has no clean closed form; FP rate ~1e-28..1e-41 |
| Rechnitzer high-precision ansatz transfer | **PARTIAL — UB-side only** | Hölder-Plancherel is L²-specific; `a/√k` ansatz mismatches h's 1/k decay |
| Spectral / translation-operator reformulation | **NEGATIVE (4× loose)** | clean duality µ=1−sup_h inf_t⟨h,T_t h⟩ but naive Rayleigh bound loose |
| Beurling-Selberg extremal majorants | **CANCELLED** | no fit to the function class; absent from Mathlib too |
| Barnard-Steinerberger / Madrid-Ramos autocorrelation transfer | **NEGATIVE** | different function class (their bracket 0.37–0.411) |
| Direct sup_t SDP (drop cell-envelope) | **INVALID bound** | cell-envelope is necessary for VALIDITY (decouples (a,b)↔(c,d)); reconstructed f∈[−1.5,1.7] |
| M-side Bochner (SOC-relaxed AND exact Schur) | **EMPIRICALLY DEAD (≤1e-7)** | one-sided slack U_m≥|f̂(m)|² absorbs all M-Toeplitz PSD content |
| ILP/SAT exact M(n) to n≤20 | **UNINFORMATIVE** | all known M(n)/n for n≤18 are >0.40; Haugland reached n≤43 |
| Lean/Mathlib mining for the saturation theorem | **BLOCKED (12+ mo)** | Mathlib lacks SDP duality, Bochner-Herglotz, Beurling-Selberg, Lasserre |
| Alternative basis (wavelet/Chebyshev to kill Gibbs) | **REFUTED at dual level** | 78% of binding dual mass at Fourier lags 0–5, not in the high-freq tail |
| bochner_n 30→40 joint-scaling | **REFUTED** | only +3.3e-6; complementarity benefit plateaus at bn=30 |
| Together UB local refinement (SLP/basin-hop, n≤1200) | **STALLED** | h* equioscillates over ~230 active shifts = true minimax local opt |
| KKT functional equation at Together's h* | **PARTIAL (byproduct: µ<0.380871 strict)** | h* fails the KKT eqn to 1e-2 ≫ 1e-9 → not a tight optimum, can't seed analytic solve |

**Key structural facts to BUILD ON (not ruled out — these are leads):**
- **(W.1) is an EXACT identity, not a relaxation.** White's cell-envelope encodes `∫₀² M(x)cos(πmx/2)dx = (4sin(πm/2)/πm)·a_m − 2(a_m²+b_m²)`. The ONLY relaxation is approximating the cell integral `I_m(j)` by its cell-min `α_m^-(j)`. `white_full_convex_exact.py` already builds the exact-integral version. This means a *finite-dimensional, exact-in-the-Fourier-variables* SDP exists.
- **Per-m relaxation residual is provably tiny and N-shrinking.** `δ_m(j) ≤ (πm)²L³/24` (Case B, numerically tight) and `πmL²/4` (Case A). Cumulative `ResidualGain ≤ (π/2N)·Σ m·λ_m`, with `Σ m·λ_m ≤ 6.1` rigorously stable across rows (LEVER_I_PRIME_THEOREM §1). Break-even N for the saturation theorem to go non-vacuous is ~12,750 (cosine) / ~16,378 (cosine+sine).
- **The cone structure is mixed LP + SOC + small PSD blocks** (`white_full_convex.py`: `sum_squares`, `quad_form(psd_wrap)`, two `bmat(...) >> 0` Bochner blocks). Fully serializable to SDPA-S; cleanly amenable to rational/exact dual certification.
- **SDPA-GMP is built + smoke-tested to 1e-75.** The ONLY missing piece is a cvxpy→SDPA-S serializer (`sdpa_gmp_wrapper.py` docstring flags it explicitly as deferred).

---

## Whitespace directions (aimed at THE BAR)

### D1 — cvxpy→SDPA-S serializer + arbitrary-precision re-certification (payoff class 1)
**Thesis.** Build the one missing engineering piece (cvxpy→SDPA-sparse serializer), feed the binding centers (core row4-region + the thin gate centers R16/R17) into SDPA-GMP at GMP precision, and extract a dual point whose objective certifies µ ≥ ~0.3803 *independently of CLARABEL's floating-point log-parse trick*. Then snap the GMP-precision dual to exact rationals and verify dual feasibility in exact arithmetic (sympy/mpmath) — yielding a checkable rational certificate.
**Why whitespace.** The serializer is explicitly deferred and never built; SDPA-GMP is otherwise idle. Every current bound rests on `reported − last_gap` from CLARABEL — the project's "central epistemic trick," not a certificate. This is the single highest-leverage rigor upgrade and is *engineering, not open math*.
**First step.** Serialize ONE small instance (e.g. the exact-(W.1) row4 program at modest N, T, bochner_n) to `.dat-s`; confirm SDPA-GMP's primal/dual match CLARABEL to the digits CLARABEL prints; diff the dual objective.
**Payoff class:** 1 (proof-grade / solver-independent certificate). **Effort:** medium (serializer ~1–2 wk incl. SOC+PSD cone bookkeeping; rational snap-and-verify on top). **Fails if:** SDPA-GMP can't fit production N in memory (mitigate: certify a small-N *valid* LB — still a theorem, just weaker; or certify only the binding sub-instance); or the GMP dual doesn't snap cleanly to low-height rationals (mitigate: keep the certificate at high-precision-interval rather than exact-rational — still solver-independent).

### D2 — exact Fejér-Riesz / trigonometric-SOS dual certificate for the binding center (payoff class 1, partial 2)
**Thesis.** At the binding center the SDP dual is a nonnegative trigonometric polynomial (the Bochner/Toeplitz PSD block is exactly a sum-of-Hermitian-squares by Fejér-Riesz). Extract the dual SOS/Gram decomposition at GMP precision, then prove its nonnegativity *exactly* via a rational LDLᵀ (or rational Fejér-Riesz factorization). This converts the binding-row LB into a closed, human/Lean-checkable nonnegativity certificate — no SDP solver in the trusted base.
**Why whitespace.** Never attempted. The retracted Lasserre work failed for lack of a tail bound; this is the *opposite* — the Bochner block is already finite-dimensional and exactly PSD-representable, so Fejér-Riesz applies cleanly. The exact-(W.1) variant (`white_full_convex_exact.py`) removes the only non-exact ingredient on the primal side, making an exact dual meaningful.
**First step.** On the smallest binding instance, take CLARABEL's (or D1's GMP) dual matrix on the Bochner block, compute its eigendecomposition, and test whether rounding the Gram factor to rationals + a small rational slack keeps it PSD and dual-feasible (certifies a slightly-weaker-but-exact LB).
**Payoff class:** 1 (formally checkable certificate); 2 if the SOS structure reveals a clean analytic form. **Effort:** medium-high (couples to D1 for the GMP dual). **Fails if:** the dual is too high-dimensional for a low-height rational Gram (mitigate: certify a reduced face / only the active constraints); or the exact slack needed to absorb rounding eats the thin 5e-5 margin (mitigate: combine with N≥24000 hardening recommended in findings.md to widen the margin first).

### D3 — Farkas/dual-ray certificates for the infeasibility-excluded corners (payoff class 1)
**Thesis.** The deep-q / high-p corners of White's wide outside regions (R6/R7/R8/R9) are currently excluded as "CLARABEL says infeasible at several interior points" — explicitly flagged as having NO Farkas certificate. For each such corner, produce an explicit dual-ray / Farkas certificate (a nonnegative combination of constraints that is infeasible), ideally in exact arithmetic. There is a clean analytic handle already found: the program enforces `|d|≤2/π` and `sum_squares(c,d)≤0.5`, so `d1=d[0]` is bounded and feasibility provably collapses for `|d1|` beyond ~0.5 (megamemory R8 node). That bound is a candidate *analytic* infeasibility proof for the bulk of the excluded mass.
**Why whitespace.** These exclusions are the single named "not certificate-grade" seam in the entire full-space argument. Upgrading them from "solver-attested" to "Farkas-certified" (or analytic) removes the last empirical link.
**First step.** Take one cleanly-infeasible corner (e.g. R8 h=0.016, p=0.10, q=−0.60); ask the solver for the infeasibility certificate / dual ray; cross-check against the analytic `sum_squares(c,d)≤0.5` ∧ `d1≥q_lo` contradiction; write it as an exact inequality chain.
**Payoff class:** 1 (closes a rigor seam → contributes to making the full-space bound a theorem). **Effort:** low-medium (per-corner; the analytic ball-bound may cover most corners in one shot). **Fails if:** some infeasible corners are infeasible only at production N (discretization), not analytically (mitigate: those are exactly where the cover already clears target on the feasible part, so they can be dropped without loss — verify each).

### D4 — close the saturation/complementarity theorem to a clean analytic statement (payoff class 2)
**Thesis.** Finish PRO-6 + LEVER-I-prime into an unconditional, publishable analytic theorem of the form: *"White's SDP program augmented with the full cell-envelope + Bochner-PSD + poly-moment stack cannot exceed C_∞ ≈ 0.380558 in the N→∞ limit"* — OR refute the ceiling. Two concrete sub-targets: (a) prove the conjectural a-priori shadow-price bounds `|ξ|≤Ξ, τ≤T, ν_3≤V` that currently make Theorem 2 only *verifiable-per-solve* rather than *unconditional* (LEVER_I_PRIME_THEOREM §2.2 — note PRO-14 already found `|ξ|≤Ω` is FALSE, the replacement is `|ξ|≤~1.5Ω`, an open constant); (b) prove the strict-form complementarity `r_CB ≤ max(r_C, r_B)` by extending the KKT identity to include the Bochner dual matrix Z.
**Why whitespace.** The tautological identity `(‡)` and the corrected residual are done; the *unconditional* version is blocked only on these named lemmas, none of which has been seriously attacked since the corrections. A clean ceiling theorem is itself a meaningful result (it tells you exactly how much of the gap is beyond this framework) and directly informs whether D1/D2 can ever reach the UB.
**First step.** Attack the shadow-price bound `|ξ|≤C·Ω`: ξ is the normalization multiplier; PRO-14's empirical `|ξ|/Ω = 1.46 ± 0.02` is strikingly stable across disparate rows — chase the scale-invariance argument (under `∫M=c` rescaling) that would pin C.
**Payoff class:** 2 (clean analytic structural theorem). **Effort:** high (real math; weeks). **Fails if:** the shadow prices have no uniform a-priori bound (then the theorem stays conditional/verifiable-only — still publishable as "verifiable per-instance," which is the current honest state).

### D5 — exact small-N SDP as a *standalone proven* theorem at a weaker constant (payoff class 1, hedge)
**Thesis.** Rather than chase the full 0.3803 certificate, prove a clean, fully-exact, solver-independent theorem at whatever constant the *smallest tractable* exact SDP yields (e.g. µ ≥ 0.3795 or 0.380, beating White 0.379005). Use the exact-(W.1) program at small N with only the Bochner block (no poly-moment, sidestepping caveat #2), serialize to SDPA-GMP, snap the dual to rationals, verify in exact arithmetic, and — if scoped tightly enough — formalize the *finite* PSD/feasibility check in Lean (Mathlib HAS `Matrix.PosSemidef` and rational linear algebra; the gap is only the general theory, not a specific finite certificate check).
**Why whitespace.** The Lean inventory ruled out formalizing the *general* saturation theorem, but a *specific finite rational certificate* is a different, much smaller target — verifying one explicit PSD matrix and one explicit feasibility inequality system is finite linear algebra Mathlib can already do. This is the most direct path to the word "theorem" (even "machine-checked theorem") on a bound that still beats White, decoupled from the thin-margin production machinery.
**First step.** Determine the smallest (N, T, bochner_n) for which the exact-(W.1) + Bochner SDP (no poly-moment, single binding center, with ellipse-extension) gives a rigorous LB strictly above White's 0.379005; that fixes the certificate size to formalize.
**Payoff class:** 1 (proof-grade, plausibly Lean-formalizable). **Effort:** medium (builds on D1/D2; Lean step is additional but bounded). **Fails if:** even the smallest White-beating instance is too large to snap to checkable rationals or to formalize (mitigate: stop at the exact-arithmetic-verified rational certificate, skip Lean — still solver-independent and proof-grade).

---

## Notes on prioritization
- **D1 is the keystone** — D2, D3, D5 all consume its GMP dual / serializer. Build it first.
- D1+D2+D3 together would convert the *existing* 0.380284 into a genuine theorem (the stated #1 payoff). D5 is the fast hedge if production-scale exact certification proves intractable. D4 is the independent analytic-theorem track (stands alone, and tells you the ceiling).
- All five avoid the ruled-out levers and target the rigor/analytic seams the verification memos themselves flag as the weak points.
- Honest framing throughout: none of these is claimed proven here. They are scoped routes to *make* it proven.
