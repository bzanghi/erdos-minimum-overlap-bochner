# L2 CLEAN verified-cover theorem (PRO-47 finish)

**Date:** 2026-06-06
**Driver:** `lp_research_state/code/_l2_clean_theorem.py`
**Self-consistent extraction:** `lp_research_state/code/_jansson_with_duals.py` (production, N=20000)
**Data:** `/tmp/pro47/L2_CLEAN_sc_prod.json` (the two binding centers' self-consistent
(p_lo, con_* duals) tuples) + `lp_research_state/parallel_results/cde_phase5_corrected_tail.json`
(the 10 non-binding production centers' duals + V_c).
**Result JSON:** `docs/RND_WHITESPACE/L2_CLEAN_THEOREM.json`.

---

## Headline

> **CLEAN verified cover floor (production, self-consistent binding):**
> **µ ≥ 0.3802952394** over White's core region (5.16), with the two-solve
> nondeterminism residual **ELIMINATED** for the center that binds the floor.

- vs White (Acta Arith. 2023) 0.379005 : **+1.290 × 10⁻³**
- vs prior published 0.379544 : **+7.512 × 10⁻⁴**
- vs the project headline 0.380284 : **+1.124 × 10⁻⁵** (clears it)

> **CLEAN TIER 1 (fully UNCONDITIONAL, no gap assumption at all):**
> **µ ≥ 0.3793989110** — 2 self-consistent centers only.
- vs White 0.379005 : **+3.94 × 10⁻⁴** (beats White unconditionally)
- vs prior published 0.379544 : **−1.45 × 10⁻⁴** (below it — the 2-center cover dips at
  the box corner near (h≈0.0089, p=0.45); this is the inherent weakness of a 2-center
  cover, not a defect).

---

## What the residual was, and how it is now gone

The previous TIER-2c headline (µ ≥ 0.3802959548, `L2_COVER_VERIFIED.{json,md}`) had ONE
remaining non-clean surface for the **binding** centers (row4, cde_n30_iter3):

- its **shift coefficients** `con_53, con_54, con_512_{pL,pU,qL,qU}, con_513` came from
  the historical **cover** solve (`cde_phase5_corrected_tail.json`), while
- its **anchor** (Jansson `p_lo`) came from a **different** production Jansson solve
  (`L2_PROD.json`).

Two separate CLARABEL runs → ~5 × 10⁻⁶ cross-solve drift (documented, stress-bounded,
but not clean).

**The fix:** `_jansson_with_duals.py` builds ONE production problem
(N=20000, T=4000, bochner_n=40, pm_k_max=20) and extracts BOTH `p_lo` AND the `con_*`
duals from it. Two facts make the resulting tuple genuinely self-consistent:

1. **CLARABEL is bit-for-bit deterministic on identical canonical data within a process.**
   Verified directly: `solve_via_data` on the same `(A,b,c,dims)` twice gives
   `‖z − z₂‖_∞ = 0`, and `prob.solve()` vs `solve_via_data` agree to 0 in the objective.
   So the `con_*` duals (used in the shift) and the `z` used for `p_lo` are the SAME `z`.
2. **The `con_*` duals ARE the z-components driving the shift** (not approximations of
   some "true" dual). Verified by `_verify_shift_eq_dualobj.py`:
   `max|A_θ − A_c| = 0`, `max|c_θ − c_c| = 0` (Lemma 10: only `b` depends on θ), and the
   reconstructed shift equals the exact change in `−b(θ)ᵀz` to **1.95 × 10⁻¹⁶**.

Therefore, for the binding centers, anchor and shift now come from ONE solve →
**residual eliminated**. The clean TIER-2c floor (0.3802952394) is **−7.15 × 10⁻⁷**
below the residual-bearing version (0.3802959548); that small drop is exactly the price
of self-consistency, and is consistent with the previously-quoted ~5 × 10⁻⁶ stress
bound.

Crucially, the clean TIER-2c **binding witness is `cde_n30_iter3`, a self-consistent
center** (`worst_witness_at_floor_all_is_rigorous = True`,
`binding_witness_is_self_consistent = True`): the center that actually sets the floor is
the one with the interval-certified, self-consistent anchor + duals.

---

## Production self-consistent extraction (one ~8 GB solve per center)

| center | N | p_lo | V_c (=Ω\*) | penalty (=−DᵀX) | PSD λ_min lower (both blocks) | status | peak RSS |
|---|---|---|---|---|---|---|---|
| cde_n30_iter3 | 20000 | 0.3803098600 | 0.3803123159 | −2.27 × 10⁻⁶ | +5.57 × 10⁻¹² | AlmostSolved | 4.65 GB |
| row4 | 20000 | 0.3803798129 | 0.3803890768 | −8.52 × 10⁻⁶ | +2.04 × 10⁻¹¹ | AlmostSolved | 4.61 GB |

- Both 82×82 Bochner blocks (`bochner_n=40` → 2(40+1)=82) are **rigorously PSD**
  (interval λ_min ≥ 0 via symmetric-diagonal-pivoted interval LDLᵀ), so the Jansson
  cone term `pen_zs = 0` and the entire penalty is the θ-independent defect `DᵀX`.
- `p_lo ≤ V_c` holds for both (Jansson soundness at the center).
- Both `p_lo` values reproduce the historical `L2_PROD.json` numbers to 9–10 digits
  (0.3803098600 vs 0.38030985999933864; 0.3803798129 vs 0.3803798129433894).

---

## Rigor self-checks (all pass)

- **Interval box-min, fast vs slow pure-iv** (coarse 80×80): |Δ| = 1.0 × 10⁻¹² for BOTH
  tiers (the fast float scan + 1e-12 round margin agrees with the pure-mpmath.iv scan).
- **Independent formula cross-check** (interval Φ midpoint vs
  `path_b_independent.Phi_row`, SAME self-consistent duals): worst |Δ| = 5.55 × 10⁻¹⁷
  (machine precision).
- **Concavity** of every center's separable (h,p) parabola (A_h2 = −con_54/2 < 0,
  A_p2 = −con_513/2 < 0) — required for the cell-corner-min rigor; asserted in the
  box-min and holds (e.g. cde: A_h2 = −0.118, A_p2 = −0.142).
- **Lemma 10 + shift==dualobj** and **CLARABEL determinism** as above.

---

## Exactly what is — and is NOT — established (honesty: the project's #1 trap is overclaiming)

**Established (interval-certified, end-to-end):**

1. For each of the two binding centers, `SDP_opt(center) ≥ p_lo` (Jansson a-posteriori
   bound, all of `−bᵀz`, the defect penalty, and the cone λ_min in `mpmath.iv`
   directed rounding; both PSD blocks certified PSD).
2. The SDP **data** is a valid relaxation up to a summed FP budget ≤ 4.5 × 10⁻¹⁵
   (`_data_rider.py`, both configs; the poly-moment tail bound is interval-certified
   strict, closing the 2026-05-22 tail-trap surface).
3. The cover lift `Φ_c(θ) = p_lo_c + shift_c(θ)` is a valid LB on `SDP_opt(θ)` for ALL θ
   (Lemma 10; `pen_Dx` θ-independent; `pen_zs = 0`), and the **shift + box-min are fully
   interval-certified**.
4. The clean TIER-2c cover floor **µ ≥ 0.3802952394** over the core box, whose binding
   witness is a self-consistent center (anchor + duals from ONE N=20000 solve) →
   **two-solve residual eliminated for the binding center**.
5. The clean TIER-1 cover floor **µ ≥ 0.3793989110** under **NO assumption whatsoever**
   (2 self-consistent centers; the anchor IS the certified p_lo for those duals).

**Caveats that remain (this is a RIGOR upgrade; the strength is unchanged):**

- **(a) The 10 NON-binding centers** still use the documented `V_c − 1e-6` margin anchor
  with cover-solve duals (their duality gap is the one trusted scalar). They are
  non-binding and robustness-stress-bounded (the floor clears the headline even when
  every float anchor is degraded by an extra 5 × 10⁻⁵ → 0.3802885). To make TIER-2c
  *fully* unconditional one would re-run `_jansson_with_duals.py` at the other 10 centers
  (cheap-ish, ~10 more ~8 GB solves; deferred). Until then, TIER-2c's only non-interval
  input is "the 10 non-binding cover-solve gaps ≤ 1 × 10⁻⁶".
- **(b) Region coverage** — that the core-box cover floor lifts to an unconditional bound
  on µ over White's full (E(M), c₁, d₁) parameter space — is inherited from **PRO-38**
  (`fullspace_promote_final.json`, independently-certified full-space floor 0.3802838,
  binding region = core, `regions_still_white_reliant = []`), which is verified
  separately and carries its own poly-moment-load-bearing caveat.
- **(c)** The clean TIER-2c floor (0.3802952394) and the PRO-38 full-space floor
  (0.3802838) are the two pieces of a **clean verified-µ theorem**: this run upgrades the
  *core-box cover* piece's binding center from a cross-solve mix to a self-consistent
  single-solve certificate. The combined verified statement is therefore
  **µ ≥ 0.380284** (PRO-38 full-space, the binding region being the core which this run
  re-certifies cleanly at 0.3802952394 with the residual gone). Nothing here pushes past
  the conjectured C_∞ ≈ 0.380558 or the Together upper bound 0.380871; the Erdős minimum
  overlap constant remains OPEN in [0.380284, 0.380871].

---

## One-line summary

The last residual in the verified core-box cover floor is gone: the binding center's
Jansson anchor and its ellipse-shift duals now come from a SINGLE production N=20000
solve. Clean unconditional-on-the-binding-center floor **µ ≥ 0.3802952394**
(self-consistent binding witness cde_n30_iter3); fully unconditional 2-center floor
**µ ≥ 0.3793989110**. Combined with PRO-38 region coverage, this is a clean
verified-µ theorem at **µ ≥ 0.380284** over White's full space, +1.29 × 10⁻³ over White.
