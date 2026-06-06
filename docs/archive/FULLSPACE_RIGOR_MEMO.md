# Full-Space Promotion — STAGE 2 RIGOR VERIFIER (adversarial)

**Date:** 2026-05-25
**Role:** Independent, adversarial certification of the two load-bearing Stage-1 claims, with **no prior session memory** and **no assumption that Stage 1 is correct**.
**Code (additive, throwaway):** `lp_research_state/code/_fullspace_rigor.py`
**Numbers:** `lp_research_state/parallel_results/fullspace_rigor.json`
**Inputs read:** `FULLSPACE_PROMOTION_STAGE1.md`, `…/fullspace_stage1.json`, `…/phase5_N20K_bn40_dualext.json`, `CLAUDE.md`, `findings.md`, White (2023) arXiv:2201.05704 (full text fetched from ar5iv).

---

## Bottom line

| Claim under test | Verdict |
|---|---|
| **(2) core anchor ≈ 0.380284** | **CONFIRMED** — independent from-scratch recompute = **0.3802837847**, agrees with the recorded value to **1.3 × 10⁻¹¹**. |
| **(1) µ ≥ 0.380000 FULL-SPACE** | **NOT INDEPENDENTLY CERTIFIED** — rests on reading White's rounded "0.38" as ≥ 0.380000 for 5 regions; my independent unaugmented recompute **cannot reach 0.380000 there** at tractable N (it plateaus ~0.375–0.379, ~0.4 % below White at White's *own* N), and the augmented Φ cover is far below 0.380 there too. This is **not a refutation of White**, but the literal-floor assumption is **unverified and mildly contraindicated**. |

The defensible, fully-independent headline from this audit is the **core**: **µ ≥ 0.380284** over region (5.16). The jump to **full-space µ ≥ 0.380000** depends on a White rounding convention I could not corroborate.

---

## TASK 1a — the rounding-convention question (read White's text)

White's Table 2 column header is literally **"Optimum lower bound"**, and the surrounding text says:

> *"In Table 2 we display **the value of the objective function of the dual program for a verified feasible point** in the dual space for several choices of input."*

A dual-feasible point's objective is, by weak duality, a **valid lower bound** on µ, and White's whole verification apparatus (Appendix II §8.2: constraints "strictly satisfied by a margin that exceeds the worst-case-scenario for floating-point rounding errors") is designed to make each entry a *certified* LB. So the entries are genuine lower bounds — **in the right direction**.

**But the adversarial question is whether the displayed "0.38" guarantees ≥ 0.380000.** On this, White's text is silent:

- The "0.38" entries are shown to **2 decimal places** with **no statement that they are rounded *down*** (truncated). Nothing in the paper rules out "0.38" being a round-to-nearest of e.g. 0.3798.
- The one place White needs more precision — the binding strip (region 18, his true global min) — he prints **0.37925 (5 sig figs)**. This is *weak evidence* that he shows extra digits exactly when a value is near a decision boundary, hence that the "0.38" entries may be comfortably above 0.380. It is **not** a proof.
- The decisive combination sentence uses the literal values: *"Combining all data from Table 2 shows that **either µ ≥ 0.37925 or** [the parameters lie in (5.16)]."* — note he draws **0.37925**, not 0.380, as the global figure.

**Finding 1a:** White's Table-2 numbers are valid lower bounds, but **the rounding *direction* is not stated**. Reading "0.38" as "≥ 0.380000" is an **assumption**, not something the paper establishes. This is exactly the kind of convention a headline must not silently lean on — hence Task 1b.

---

## TASK 1b — independent recompute of White's *own* unaugmented §5 bound

I re-solved White's **unaugmented** §5 program (`build_problem_with_dual_handles(..., bochner_n=0)` — **no Bochner, no poly-moment, no T3/T5/T5'**), dual-extracted to a conservative LB (`rigorous_dual_LB − 1e-5`), the project's rigor convention. **The "value" is never used.**

### Box-validity and the mechanism correction (the central subtlety)

White obtains **one** number per Table-2 region. The mechanism matters:

- White states the full-range program optimum *"ends up close to 0.25"* and that he uses a **"divide and conquer" strategy of breaking up the ranges into small chunks; the minimum of all the optimums over these smaller intervals will be our lower bound."* For the **core** he uses a **7-ellipse feasible-point cover** (Table 3).
- My first pass solved each region as a **single full-range box**. That is **box-valid** (parameters live only in RHS of (5.3),(5.4),(5.12),(5.13), so one full-range solve lower-bounds µ over the whole box — White Lemma 9/10), **but it is a far weaker mechanism** than White's divide-and-conquer.

Concrete proof that single-box ≠ White's mechanism:

| quantity (unaugmented, N=5000) | value | White |
|---|---|---|
| single-box solve of the **core** (5.16) | 0.37510 | — |
| **7-ellipse cover** of the core (my recompute) | 0.378335 | **0.379005** (his cover, N=25000) |

So a single wide-box solve **understates** the true regional bound by ~4 × 10⁻³; my numbers below must be read with that in mind.

**Validity of my program** (so the numbers are trustworthy LBs, not artifacts):
- At White's **Table-3 core centers** my unaugmented program **reproduces or EXCEEDS** White's published "initial objective" values: row4 **0.37925** vs White 0.37905; row7 **0.38085** vs 0.3794; strip **0.37943** vs 0.37925. My program is therefore **not weaker** than White's (if anything slightly stronger, likely exact vs Lipschitz cell bounds).
- **Cross-solver:** R16 via **SCS = 0.37460** vs CLARABEL 0.37453 (agree to 7 × 10⁻⁵) — not a CLARABEL artifact.
- **Sanity:** augmented row4 (bn=40, pmk=20, N=10000) = **0.38016**, matching the known headline ~0.3803.
- **N-monotonicity confirmed** (White: "the optimum increases with N"); every value is a valid LB at every N.

### Per-region independent unaugmented box LB (single-box, N=5000, conservative)

| # | E(M) | c1 | d1 | White disp. | our box LB | ≥ 0.380000? |
|---:|---|---|---|---|---|:--:|
| 1 | (0.75,2) | (0,1) | (−1,1) | 0.38 | 0.43942 | YES |
| 2 | (0.4,0.75) | (0,1) | (−1,1) | 0.38 | 0.40854 | YES |
| 3 | (0.2,0.4) | (0,1) | (−1,1) | 0.38 | 0.39689 | YES |
| 4 | (0.1,0.2) | (0,1) | (−1,1) | 0.38 | 0.38325 | YES |
| 5 | (0.08,0.1) | (0,1) | (−1,1) | 0.38 | 0.38279 | YES |
| 6 | (0,0.08) | (0,1) | (−1,−0.05) | 0.38 | **0.37650** | **NO** |
| 7 | (0,0.08) | (0,1) | (−0.05,−0.025) | 0.38 | **0.37349** | **NO** |
| 8 | (0,0.08) | (0,1) | (0.05,1) | 0.38 | 0.38266 | YES |
| 9 | (0,0.08) | (0,1) | (0.025,0.05) | 0.38 | **0.37361** | **NO** |
| 10 | (0,0.08) | (0,0.25) | (−0.025,0.025) | 0.38 | 0.39894 | YES |
| 11 | (0,0.08) | (0.25,0.3) | (−0.025,0.025) | 0.38 | 0.38955 | YES |
| 12 | (0,0.08) | (0.3,0.33) | (−0.025,0.025) | 0.38 | 0.38332 | YES |
| 13 | (0,0.08) | (0.5,1) | (−0.025,0.025) | 0.38 | 0.39148 | YES |
| 14 | (0,0.08) | (0.45,0.5) | (−0.025,0.025) | 0.38 | 0.38116 | YES |
| 15 | (0.06,0.08) | (0.33,0.45) | (−0.025,0.025) | 0.38 | 0.38297 | YES |
| 16 | (0,0.06) | (0.33,0.45) | (−0.025,−0.02) | 0.38 | **0.37452** | **NO** |
| 17 | (0,0.06) | (0.33,0.45) | (0.02,0.025) | 0.38 | **0.37453** | **NO** |
| 18 | (0,0.06) | (0.33,0.35) | (−0.02,0.02) | **0.37925** | 0.37942 | NO (expected) |

**12 of 17 "0.38" regions clear 0.380000** from the (weak) single-box recompute. **Five do not**: R6, R7, R9, R16, R17. Region 18 (strip) gives 0.37942 — *above* White's 0.37925 (validating White's strip number and the variable map) but below 0.380000, as expected; the strip is meant to be lifted only by the augmented cover.

### Are the 5 "below" regions an N-undershoot, or real?

I pushed them up in N (each unaugmented solve is light: N=10000 ≈ 6 s, ~1.1 GB; N=20000 ≈ 10–15 s, ~2.0 GB — run sequentially):

| region | N=5000 | N=10000 | N=20000 | White (his N) |
|---|---|---|---|---|
| R6 | 0.37650 | 0.37811 | — | 0.38 (N=10000) |
| R7 | 0.37349 | 0.37428 | **0.37469** | 0.38 (N=10000) |
| R9 | 0.37361 | 0.37550 | — | 0.38 (N=10000) |
| R16 | 0.37452 | 0.37507 | **0.37539** | 0.38 (N=20000) |
| R17 | 0.37453 | 0.37532 | **0.37612** | 0.38 (N=20000) |

They rise monotonically but **plateau ~0.375–0.376, NOT 0.380** — and crucially, **at White's own N (R16/R17: N=20000) they sit ~0.4 % below his reported "0.38".** Fine 3-D subdivision and even the **augmented** full-box program do not lift them (augmented R7 = 0.3732, R16 = 0.3747, R6 = 0.3759).

This is **not** explained by my program being weaker than White's (it reproduces/exceeds him at the core centers). The remaining explanations are: White obtained "0.38" for these regions via a finer divide-and-conquer than I replicated (the right mechanism — but then his min over sub-boxes would also have to clear 0.380, which my sub-box probes near the binding c1≈0.39 manifold do **not**), **or** "0.38" is a round-to-nearest of a value in [0.375, 0.380) for these specific regions. **I cannot independently confirm ≥ 0.380000 for R6, R7, R9, R16, R17.**

### What the µ ≥ 0.380000 claim actually rests on

For R6, R7, R9, R16, R17 the augmented **Stage-1 12-anchor Φ cover** is also below 0.380000:

| region | Stage-1 augmented Φ-cover min | source of any ≥0.380000 floor |
|---|---|---|
| R6 | 0.2775 | **White's 0.38 only** |
| R7 | 0.3659 | **White's 0.38 only** |
| R9 | 0.3659 | **White's 0.38 only** |
| R16 | 0.38026 | White's 0.38 (Φ is 1.8e-4 short) |
| R17 | 0.38026 | White's 0.38 (Φ is 1.8e-4 short) |

So the full-space **µ ≥ 0.380000** claim leans **entirely on White's literal "0.38" floor** for these five regions — the exact assumption Task 1a found to be unstated and Task 1b found to be unreachable by independent computation at tractable N. **Verdict 1: not independently certified.**

---

## TASK 2 — independent re-derivation of the core anchor (≈ 0.380284)

Re-implemented from scratch in `_fullspace_rigor.py:task2` (does **not** call `_verify_cover_dualext.py` or reuse `_fullspace_eval.py`): load the 12 conservative dual anchors (`primal − 1e-5`) from `phase5_N20K_bn40_dualext.json`; rebuild each Φ_c(h,p,q) directly from `dual_objective_shift`; cover = max over centers; min over the **core** (5.16) on a **4001 × 4001** (h,p) grid at q ∈ {−0.02, 0, +0.02}; subtract the Lipschitz cell-error `eps_grid = L_max · ½ · diag(cell)`.

```
grid_min     = 0.3802859582   at (h = 0.00399, p = 0.39227)
eps_grid     = 2.174e-06       (L_max = 0.1491)
rigorous_LB  = 0.3802837847
```

Agreement:
- vs recorded Stage-1 `core_eval` **0.3802837846529683**: **+1.27 × 10⁻¹¹**
- vs canonical `find_ellipse_h_p` cross-check path (a second, independent code path): **+1.27 × 10⁻¹¹**
- vs the rounded conservative headline 0.380284: −2.15 × 10⁻⁷ (just the grid/Lipschitz margin)

**Verdict 2: CONFIRMED.** Two independent code paths agree to 10⁻¹¹, far inside the 10⁻⁶ bar — no bug. The conservative core anchor **0.3802838** (headline 0.380284) is solid.

---

## Discrepancies / bugs found

1. **No bug in the core-anchor evaluation** (Task 2 reproduces to 1e-11).
2. **Mechanism mismatch in the full-space argument:** certifying a Table-2 region requires White's **divide-and-conquer / feasible-point cover**, not a single full-range box solve (a single box understates by ~4e-3; demonstrated on the core: 0.37510 single-box vs 0.379005 cover). Any future per-region certification must use the cover/subdivision, at production N.
3. **Unverified rounding assumption with mild contra-evidence:** the µ ≥ 0.380000 full-space claim depends on White's "0.38" meaning ≥ 0.380000 for R6, R7, R9, R16, R17. My independent unaugmented recompute plateaus ~0.4 % below 0.38 for these even at White's own N, and my program is *not* weaker than White's (it matches/exceeds him at the core centers). This does **not** prove White wrong, but it means **0.380000 is not independently established** by this audit.

## Recommendation

- **Quote the core result (µ ≥ 0.380284 over (5.16)) as the rigorously-independent headline.** It is confirmed to 1e-11.
- **Do not headline full-space µ ≥ 0.380000** until one of:
  (a) White's **exact** per-region dual objectives are obtained (his paper: "available upon request") and shown ≥ 0.380000 for R6, R7, R9, R16, R17; or
  (b) those five regions are re-certified with the **correct mechanism** — divide-and-conquer min over fine sub-boxes **plus** the augmented Φ cover, at production N (N≈20000–25000) — and shown ≥ 0.380000. My probes suggest this is **not guaranteed** to succeed near the binding c1≈0.39 manifold and must be actually computed, not assumed.

## Rigor guardrails honored

Conservative dual-extracted anchors only (never raw "value"); no Lasserre values; "center solve" vs "single-box bound" vs "cover/divide-and-conquer bound" vs "µ bound" kept crisp throughout; failures reported honestly (no forced agreement); **no git commands**; **additive-only** files (`_fullspace_rigor.py`, this memo, `fullspace_rigor.json`); SDP solves run **sequentially** at light configs on the shared machine.
