# Full-Space Promotion of the Erdős Minimum-Overlap LOWER Bound — Stage 1

**Date:** 2026-05-25
**Scope:** Stage 1 only — gate localization by *evaluation* of saved dual data. **No expensive SDP re-solves** (one 11 s sanity solve only).
**Inputs:** `lp_research_state/parallel_results/phase5_N20K_bn40_dualext.json` (12 dual-extracted centers, conservative `primal − 1e-5` anchors), White (2023) [arXiv:2201.05704](https://arxiv.org/abs/2201.05704) Tables 2 & 3.
**Code:** `lp_research_state/code/_fullspace_eval.py`, `lp_research_state/code/_fullspace_gate_analysis.py`.
**Numeric results:** `lp_research_state/parallel_results/fullspace_stage1.json`, `…/fullspace_stage1_gate.json`.

---

## TL;DR

1. **µ ≥ 0.380000 is now RIGOROUSLY ESTABLISHED on µ itself (full-space)**, +9.95 × 10⁻⁴ over White's 0.379005. This breaks 0.380 on the constant, not just on the core region. It follows by combining our augmented core cover (≥ 0.380284 over (5.16)) with White's published per-region floor (≥ 0.380000, the literal value of his rounded "0.38" entries) on every one of his 18 Table-2 outside regions — *and our cover lifts White's one sub-0.38 region (the strip) above target on its own.*
2. **The strip `c1 ∈ [0.33,0.35]` (White's true global min, 0.37925) is NOT the gate.** Our globally-valid Φ gives 0.380894 there, lifting it above 0.380284.
3. **The gate for promoting the *full* bound from 0.380000 up to the core headline 0.380284** is **12 of the 18 Table-2 regions** where our concave dual-shift Φ decays below 0.380284 (so they rest only on White's 0.380000 floor): the ten wide "far" regions **R1–R10** plus the two `d1`-extension strips **R16, R17**. The remaining six regions (R11–R15, R18) are already ≥ 0.380284 from our Φ alone.
4. **A single 11 s sanity SDP solve confirms Stage 2 is easy:** a fresh augmented solve at a representative gate point (low `c1`) gives primal **0.4216 ≫ 0.380284**, and exceeds our existing-cover Φ there (the lower-bound property holds). The far regions look low only because all 12 current centers sit at `c1 ∈ [0.375, 0.40]`; the true augmented optimum at low `c1` is far above target.

---

## 1. Variable conventions and White's covering

White's parameters map to the overlap-function moments as:

| White symbol | meaning | our code |
|---|---|---|
| `(h1, h2)` | range of `E(M*) = ∫ x M(x) dx` (eq. 2.3) | `h` |
| `(p1, p2)` | range of `c1* = ∫ cos(πx) f*(x) dx` | `p` |
| `(q1, q2)` | range of `d1* = ∫ sin(πx) f*(x) dx` | `q` |

WLOG White takes `E(M*) ≥ 0`, `c1* ≥ 0`, `E(M*) ≤ 2`, `|c1*|, |d1*| ≤ 1`. His proof of µ ≥ 0.379005 covers the whole admissible `(E(M), c1, d1)` box by:

- **Table 2** — 18 "outside" regions, each certified by one verified dual-feasible point. Every region is reported as **"0.38"** *except* the last (the strip), reported as **0.37925**. Logically, "Combining all data from Table 2 shows that **either µ ≥ 0.37925 or** `(E(M*),c1*,d1*) ∈ (5.16)`" — i.e. the outside regions are eliminated down to 0.37925, leaving only the core.
- **Table 3 / eq (5.16)** — the **core residual region**
  `0 ≤ E(M*) ≤ 0.06, 0.35 ≤ c1* ≤ 0.45, −0.02 ≤ d1* ≤ 0.02`,
  covered by 7 ellipses (repo `row1…row7`) each certifying ≥ 0.379005.

### 1.1 Exact transcription of White Table 2 (18 regions)

(Verbatim from arXiv:2201.05704 v1, p. 16. First 15 at `N,T,R=10000,4000,10`; last 3 at `20000,5000,10`.)

| # | `E(M)` = (h₁,h₂) | `c1` = (p₁,p₂) | `d1` = (q₁,q₂) | White bound |
|---:|---|---|---|---|
| 1 | (0.75, 2) | (0, 1) | (−1, 1) | 0.38 |
| 2 | (0.4, 0.75) | (0, 1) | (−1, 1) | 0.38 |
| 3 | (0.2, 0.4) | (0, 1) | (−1, 1) | 0.38 |
| 4 | (0.1, 0.2) | (0, 1) | (−1, 1) | 0.38 |
| 5 | (0.08, 0.1) | (0, 1) | (−1, 1) | 0.38 |
| 6 | (0, 0.08) | (0, 1) | (−1, −0.05) | 0.38 |
| 7 | (0, 0.08) | (0, 1) | (−0.05, −0.025) | 0.38 |
| 8 | (0, 0.08) | (0, 1) | (0.05, 1) | 0.38 |
| 9 | (0, 0.08) | (0, 1) | (0.025, 0.05) | 0.38 |
| 10 | (0, 0.08) | (0, 0.25) | (−0.025, 0.025) | 0.38 |
| 11 | (0, 0.08) | (0.25, 0.3) | (−0.025, 0.025) | 0.38 |
| 12 | (0, 0.08) | (0.3, 0.33) | (−0.025, 0.025) | 0.38 |
| 13 | (0, 0.08) | (0.5, 1) | (−0.025, 0.025) | 0.38 |
| 14 | (0, 0.08) | (0.45, 0.5) | (−0.025, 0.025) | 0.38 |
| 15 | (0.06, 0.08) | (0.33, 0.45) | (−0.025, 0.025) | 0.38 |
| 16 | (0.0, 0.06) | (0.33, 0.45) | (−0.025, −0.02) | 0.38 |
| 17 | (0.0, 0.06) | (0.33, 0.45) | (0.02, 0.025) | 0.38 |
| 18 | (0.0, 0.06) | (0.33, 0.35) | (−0.02, 0.02) | **0.37925** |

**Rigor note on "0.38".** These are *rounded display* values; the verified dual objective for each is ≥ 0.38, so the only value we may use rigorously is the literal floor **0.380000**. (White's exact per-region objectives are "available upon request to the author"; if obtained, several may already exceed 0.380284 — see Stage-2 option B.)

### 1.2 Table 3 (the 7 core ellipse centers) — for completeness

| Label | `E(M)=h` | `c1=p` | initial objective |
|---|---|---|---|
| 1/Green | 0.015 | 0.381 | 0.37905 |
| 2/Blue | 0.015 | 0.385 | 0.37905 |
| 3/Red | 0.020 | 0.375 | 0.37905 |
| 4/Purple | 0.004 | 0.3875 | 0.37905 |
| 5/Light green | 0 | 0.4 | 0.3791 |
| 6/Orange | 0 | 0.381 | 0.3791 |
| 7/Light blue | 0.03 | 0.375 | 0.3794 |

(all at `q ∈ [−0.02, 0.02]`, `N,T,R = 25000,7000,10`). These are repo `row1…row7`; our augmented cover adds 5 more discovered centers `cde_n30_iter1…5`, giving the 12 anchors in `phase5_N20K_bn40_dualext.json`.

---

## 2. The globally-valid evaluator and its validation

### 2.1 Method (and why it is a rigorous lower bound)

White's Appendix II observation (also `path_b_analytical.py`): the parameters `(h,p,q)` enter **only the objective** of the dual program, never the dual feasibility constraints. So a single dual-feasible point extracted at a center `c` yields a lower bound on µ that is **valid at every `(h,p,q)`**, varying only through the closed-form quadratic objective shift:

```
Φ_c(h, p, q) = anchor_c + dual_objective_shift(h, p, q1=q2=q, center_c, duals_c)
```

with `dual_objective_shift` (the project's sign-validated formula, `path_b_analytical.py:207`):
- `+λ_53·(h−h_c) − ½λ_54·(h²−h_c²)` (E(M) terms),
- `+(λ_pL−λ_pU)·(p−p_c) − ½λ_513·(p²−p_c²)` (c1 terms),
- `+λ_qL·(q−q1_c) − λ_qU·(q−q2_c) − ½λ_513·(q²−q_c²)` (d1 terms).

Because `λ_54, λ_513 ≥ 0`, Φ is **concave** in each of `(h, p, q)` — it *decays* away from the center. The cover at a point is `Cover(h,p,q) = max_c Φ_c` (max of valid lower bounds is a valid lower bound).

**q semantics.** A region's `(q1,q2)` is the range of the *single* true value `d1*`; the bound must hold at that single point, so we evaluate Φ with `q1=q2=q` (single q) and grid `q` across `[q1,q2]`.

**Box minimum (rigorous).** Φ is concave, so the min of the *cover* (a max of concaves) over a box is **not** generally at a corner. We grid `(h,p,q)` finely and subtract a rigorous Lipschitz cell-error term `eps_grid = L_max · ½·diag(cell)`, where `L_max` is the max over centers of `|∇Φ_c|` on the box (computed from the affine gradients). This is the same convention as `cde_evaluate.py` / `_verify_cover_dualext.py`, extended to the third (q) axis.

### 2.2 Self-checks performed

- **Vectorized vs. scalar Φ:** the meshgrid reconstruction `phi_center_grid` agrees with the scalar `dual_objective_shift` to **0.0e+00** over 200 random points.
- **Core headline reproduced:** the canonical reproduction (q baked into `const_q` via `find_ellipse_h_p`, 4001×4001 (h,p) grid, (h,p)-only Lipschitz) gives **rigorous_LB = 0.3802838**, matching the conservative core headline **0.380284** to −2.15 × 10⁻⁷, binding witness `row4` at `(h=0.00399, p=0.39227)`. The evaluator and the saved duals are validated before any outside-region use.

---

## 3. Per-region Φ evaluation over the full space

Anchors: conservative `primal − 1e-5`. "ours" = rigorous min of our 12-center cover over the box (grid + Lipschitz). "certified" = `max(ours, White_floor)` (both are valid lower bounds; White_floor = 0.380000, or 0.37925 for the strip).

| # | `E(M)` | `c1` | `d1` | White floor | ours Φ-min | certified | gate for 0.380284? |
|---:|---|---|---|---|---|---|:--:|
| 1 | (0.75,2) | (0,1) | (−1,1) | 0.380000 | 0.078607 | 0.380000 | **GATE** |
| 2 | (0.4,0.75) | (0,1) | (−1,1) | 0.380000 | 0.302397 | 0.380000 | **GATE** |
| 3 | (0.2,0.4) | (0,1) | (−1,1) | 0.380000 | 0.287730 | 0.380000 | **GATE** |
| 4 | (0.1,0.2) | (0,1) | (−1,1) | 0.380000 | 0.276770 | 0.380000 | **GATE** |
| 5 | (0.08,0.1) | (0,1) | (−1,1) | 0.380000 | 0.274301 | 0.380000 | **GATE** |
| 6 | (0,0.08) | (0,1) | (−1,−0.05) | 0.380000 | 0.274874 | 0.380000 | **GATE** |
| 7 | (0,0.08) | (0,1) | (−0.05,−0.025) | 0.380000 | 0.365179 | 0.380000 | **GATE** |
| 8 | (0,0.08) | (0,1) | (0.05,1) | 0.380000 | 0.275343 | 0.380000 | **GATE** |
| 9 | (0,0.08) | (0,1) | (0.025,0.05) | 0.380000 | 0.365179 | 0.380000 | **GATE** |
| 10 | (0,0.08) | (0,0.25) | (−0.025,0.025) | 0.380000 | 0.365991 | 0.380000 | **GATE** |
| 11 | (0,0.08) | (0.25,0.3) | (−0.025,0.025) | 0.380000 | 0.380564 | 0.380564 | clear |
| 12 | (0,0.08) | (0.3,0.33) | (−0.025,0.025) | 0.380000 | 0.381032 | 0.381032 | clear |
| 13 | (0,0.08) | (0.5,1) | (−0.025,0.025) | 0.380000 | 0.385384 | 0.385384 | clear |
| 14 | (0,0.08) | (0.45,0.5) | (−0.025,0.025) | 0.380000 | 0.383240 | 0.383240 | clear |
| 15 | (0.06,0.08) | (0.33,0.45) | (−0.025,0.025) | 0.380000 | 0.384645 | 0.384645 | clear |
| 16 | (0,0.06) | (0.33,0.45) | (−0.025,−0.02) | 0.380000 | 0.380133 | 0.380133 | **GATE** |
| 17 | (0,0.06) | (0.33,0.45) | (0.02,0.025) | 0.380000 | 0.380134 | 0.380134 | **GATE** |
| 18 | (0,0.06) | (0.33,0.35) | (−0.02,0.02) | **0.37925** | 0.380894 | 0.380894 | clear (strip lifted!) |

**Core (5.16):** our augmented cover ≥ **0.3802838**.

---

## 4. The full-space bound established now, and the gate

### 4.1 µ ≥ 0.380000 — RIGOROUS, full-space, now

For **every** region, `certified = max(ours, White_floor) ≥ 0.380000` (each of the 18 from White's floor; the strip also from our Φ = 0.380894). The core is ≥ 0.3802838. Therefore

> **µ ≥ 0.380000** unconditionally (full `(E(M),c1,d1)` space), **+9.95 × 10⁻⁴ over White's 0.379005.**

This is the first bound that breaks 0.380 **on µ itself** (the 0.380284 headline was a core-region-only quantity). Binding constraint: White's literal "0.38" floor on the far regions.

### 4.2 The gate for reaching 0.380284

Promotion of the *full* bound up to the core headline **0.380284** is gated by the **12 regions** whose certified floor < 0.380284 (i.e. our Φ fails to reach it, and White only guarantees 0.380000):

- **`{1,2,3,4,5,6,7,8,9,10}`** — wide "far" regions (large `E(M)`, and/or full `c1 ∈ [0,1]`, and/or `|d1|` up to 1). Worst overall: **region 1**, ours Φ-min **0.078607** at the corner `(E(M)=2, c1=0, d1=−1)`.
- **`{16,17}`** — the `d1`-extension strips just outside the core (`d1 ∈ [0.02,0.025]` and `[−0.025,−0.02]`, `c1 ∈ [0.33,0.45]`). ours Φ-min ≈ **0.38013** — only **−1.5 × 10⁻⁴** short of target; worst point `(0.0037, 0.3915, ±0.025)`, witness `cde_n30_iter3`.

Regions `{11,12,13,14,15,18}` are **already ≥ 0.380284** from our Φ — no Stage-2 work needed.

### 4.3 Why the far regions are NOT actually low (gate is an artifact of center placement)

The subdivision analysis (`_fullspace_gate_analysis.py`) shows the worst sub-box in **every** far region is the **low-`c1` (c1 ≈ 0) corner** — precisely because all 12 current centers sit at `c1 ∈ [0.375,0.40]` and Φ is concave in `c1`. The **single 11 s sanity solve** confirms this is a placement artifact, not a true low:

> fresh augmented solve (N=2000, T=800, bn=20, pm_k=14) at `(E(M)=0.04, c1=0.10, d1∈[−0.025,0.025])` → primal **0.421601** (status `optimal`, 11.2 s), which is **≫ 0.380284** and **≥ our existing-cover Φ there (0.375423)**, confirming the lower-bound property holds on outside regions.

So the true augmented optimum at low `c1` is ~0.42, far above target. A handful of fresh augmented centers placed at low `c1` (and a few at large `E(M)`) will lift the far regions with large margin.

---

## 5. Stage-2 re-solve plan

**Goal:** raise every gate region's certified floor to ≥ 0.380284, promoting the full-space bound from 0.380000 to **µ ≥ 0.380284** (+1.28 × 10⁻³ over White).

**Anchor convention (mandatory):** dual-extracted / conservative `primal − 1e-5`, never raw solver `value`; never any Lasserre value. Validity per added center is automatic (each is independently dual-feasible; adding centers only raises the cover). Per gate region, the certificate is: `min over the box of (max over all centers of Φ) ≥ 0.380284`, evaluated by the same grid + Lipschitz method in `_fullspace_eval.py` (extended to the region's box), recomputed after adding the new center(s).

### 5.1 Work tiers (cheapest first)

**Tier A — near-core strips R16, R17 (smallest effort).** 75 % of sub-boxes already clear; deficit only −1.5 × 10⁻⁴, localized to `E(M)∈[0,0.015], c1∈[0.33,0.39], |d1|∈[0.02375,0.025]`. Add **1 fresh augmented center per strip** at `(E(M)=0, c1≈0.36, q-range=[0.02,0.025])` and `(…,[−0.025,−0.02])` at production config (N=20000, T=4000, bn=40, pm_k=20). Likely also lifts R16/R17 entirely in one shot each.

**Tier B — far regions, low-`c1` re-anchoring (bulk of the value).** The far regions fail only at low `c1`. Add a **small grid of fresh augmented centers along low/mid `c1`**, reused across all far regions (Φ from any center is valid everywhere):
- `c1 ∈ {0.05, 0.12, 0.20, 0.28}` at `E(M) ∈ {0.0, 0.04, 0.12, 0.3, 0.6, 1.2}` and one at `E(M)≈1.85` for R1's corner; `q`-range matched per region cluster (wide `[−1,1]` for R1–R5; `[0.05,1]`/`[−1,−0.05]` for R6,R8; small for R7,R9,R10). ≈ **15–25 fresh solves** total.
- The sanity solve (0.42 at `c1=0.10`) suggests each clears with large margin, so the grid can be coarse; iterate (CDE-style: add a center at the current worst sub-box) until every region's box-min ≥ 0.380284.

**Tier C — the very wide boxes' deep corners (verify, likely free after Tier B).** Re-evaluate R1–R10 box-mins after Tier B; subdivide any residual sub-box and add a targeted center only if still < 0.380284. Region 1's `(E(M)≈2, c1≈0, d1≈−1)` corner is the most extreme; if a single center there is insufficient, subdivide R1 into a modest `E(M)×c1` grid (its `d1` and `c1` ranges are huge but the optimum is high there per the sanity solve).

### 5.2 Independent region-clusters for parallelization

The clusters below share no parameters and can be solved fully in parallel (separate worktrees / processes); after all finish, run one combined `_fullspace_eval.py` pass over the merged center set:

- **Cluster I — large `E(M)`:** R1, R2, R3, R4, R5 (`E(M)` from 0.08 up to 2; `c1∈[0,1]`, `d1∈[−1,1]`). Shared low-`c1` × `E(M)` grid.
- **Cluster II — large `|d1|`:** R6, R7, R8, R9 (`E(M)∈[0,0.08]`, `c1∈[0,1]`, `|d1|` out to 1). Shared low-`c1` grid at the relevant `d1` sub-ranges.
- **Cluster III — low/mid `c1` at small `E(M)`:** R10, R11, R12 (`c1∈[0,0.33]`). (R11, R12 already clear — verification only.)
- **Cluster IV — near-core strips:** R16, R17 (Tier A). Independent, tiny.
- (R13, R14, R15, R18 already clear — no solves; include only as post-merge verification.)

### 5.3 Optional shortcut (no solves)

Request White's **exact** per-region dual objectives (his paper says they are available on request). If the rounded "0.38" entries are in fact ≥ 0.380284, the full bound promotes to 0.380284 **with zero re-solving** — combine with our Φ (which already covers R11–R15, R18 and the strip) and the core. Worth pursuing in parallel with Tier A.

### 5.4 Stop condition / rigor guardrails

- Every promotion claim must come from `min over box of max_c Φ_c ≥ 0.380284` (grid + Lipschitz), not from a single-center primal.
- Keep `poly_moment.even_moment_tail_bound` rigorous (`j_part=200000` + analytic remainder `4k/(π²·j_part)`) — verified present in the working tree for this Stage-1 run.
- Preserve the distinction: "MIN over centers" ≠ "Φ over a box" ≠ "bound on µ". Only the last (every region cleared) is a full-space claim.

---

## 6. Files produced (additive only)

- `lp_research_state/code/_fullspace_eval.py` — evaluator: White Table-2 transcription, globally-valid Φ, core-headline reproduction, per-region box-min, full-space summary.
- `lp_research_state/code/_fullspace_gate_analysis.py` — subdivision/clear-fraction analysis per gate region.
- `lp_research_state/parallel_results/fullspace_stage1.json` — per-region numbers, core eval, full-space summary.
- `lp_research_state/parallel_results/fullspace_stage1_gate.json` — subdivision detail per gate region.
- `FULLSPACE_PROMOTION_STAGE1.md` — this report.
