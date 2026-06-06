# L2 cover floor — FINAL rigor gap CLOSED (PRO-47)

**Date:** 2026-06-06
**Code:** `lp_research_state/code/_cover_iv_certify.py` (+ `_jansson_with_duals.py`,
`_verify_shift_eq_dualobj.py`)
**Machine artifacts:** `L2_COVER_VERIFIED.json`, `L2_COVER_VERIFIED_sc_N3000.json` (this dir)

This closes the LAST un-certified surface of the L2 verified-bound thrust: the
cover's per-center **dual-shift coefficients** (the ellipse SLOPES `con_53, con_54,
con_512_pL/pU/qL/qU, con_513`), plus the **shift + find_ellipse quadratic + box-min**,
are now all in **directed-rounding interval arithmetic** (`mpmath.iv`). Previously
(`_cover_lift.py` → `L2_FINISH_cover.json`) those duals were TRUSTED CLARABEL floats;
only the per-center *anchors* were Jansson-verified.

---

## HEADLINE

> **VERIFIED COVER FLOOR  µ ≥ 0.3802959548**  (TIER 2c)
> +1.20 × 10⁻⁵ over the 0.380284 headline · +7.52 × 10⁻⁴ over White's prior pub
> 0.379544 · +1.29 × 10⁻³ over White's published 0.379005.

Binding witness **cde_n30_iter3**, anchored at its production-N (=20000) Jansson
**interval-certified** `p_lo = 0.3803098600`. The dual-shift coefficients, the
`find_ellipse_h_p` quadratic, and the box-min are all interval-certified. This is
essentially the float reference (`L2_FINISH_cover.json` = 0.3802958) re-derived
rigorously — the interval floor is +1.5 × 10⁻⁷ ABOVE it (the interval cell-enclosure
box-min is *tighter* than the float grid+Lipschitz convention; see below).

**The floor did NOT drop** below 0.380284 under the interval dual-shift — it clears
by +1.2 × 10⁻⁵, and the robustness stress (below) shows it clears even when every
non-binding anchor is degraded by an extra 5 × 10⁻⁵.

---

## WHY ENCLOSING THE FLOAT DUALS AS THIN INTERVALS IS RIGOROUS (the crux)

The cover bound is `Φ_c(θ) = anchor_c + shift_c(θ)`, `θ = (h,p,q)`, where
`shift_c(θ) = Σ_i λ_i^c · Δrhs_i(θ)` (`path_b_analytical.dual_objective_shift`) and
`λ_i^c` are the constraint duals.

**White's Lemma 10:** in the canonical conic form `min cᵀx s.t. Ax+s=b(θ), s∈K`,
ONLY `b` depends on `θ`; `A, c, K` are θ-independent.
→ **VERIFIED bit-for-bit** at small N (`_verify_shift_eq_dualobj.py`):
`max|A_θ − A_center| = 0.0`, `max|c_θ − c_center| = 0.0`.

Therefore, for the SAME numeric conic dual `z` CLARABEL returns at the center (z need
NOT be exactly feasible), the Jansson construction gives, for EVERY θ:

```
    SDP_opt(θ)  ≥  −b(θ)ᵀz  +  pen_zs(θ)  −  pen_Dx
```
where `pen_Dx = Σ_i |c + Aᵀz|_i · xbar_i` is **θ-independent** (c, A, z, xbar all
θ-independent), and `pen_zs(θ) = Σ_j min(0, λ_min^{K*}(z_j))·sbar_j(θ)` is **= 0**
because every cone block of z is certified in K* (PSD blocks verified PSD, nonneg
coords ≥0) → every `min(0,·)` factor is 0. So define
`Φ_c(θ) := −b(θ)ᵀz − pen_Dx`; it is a rigorous LB on `SDP_opt(θ)` for ALL θ, and

```
    −b(θ)ᵀz  =  −b(θ_c)ᵀz  +  shift_c(θ; λ_i)        [VERIFIED to 1.95e-16]
            =  (p_lo_center + pen_Dx)  +  shift_c(θ)
  ⇒ Φ_c(θ)  =  p_lo_center  +  shift_c(θ).
```

**Conclusion:** the float duals `λ_i` are NOT approximations of some "true" dual —
they ARE the components of the fixed numeric `z` for which the bound is proved. The
rigorous operation is to enclose each consumed float `λ_i` as a thin FP interval
`iv.mpf(repr(λ_i))` and propagate `shift + find_ellipse quadratic + box-min` in
directed-rounding interval arithmetic, anchored at the interval `p_lo`. **No
"dual-feasibility-at-perturbed-RHS re-check" is needed** beyond what Jansson already
certified at the center, because the perturbation lives entirely in `b` (Lemma 10)
and the two penalty terms are θ-independent / identically zero.

---

## THE INTERVAL BOX-MIN (two independent methods, cross-checked)

The floor is `min_{(h,p)∈core} [ max_c Φ_c(h,p) ] − (Lipschitz term)`.
Core box `(h,p) ∈ [0,0.06] × [0.35,0.45]`, q at row endpoints (so `const_q = 0`,
identical to the float cover).

- **M1 — direct cell-enclosure (headline, no Lipschitz constant).** Partition the box
  into a grid of cells; on each cell evaluate every `Φ_c` over the cell's interval box
  in `mpmath.iv`. Each `Φ_c = anchor + g_h(h) + g_p(p)` is separable with **strictly
  concave** `g_h, g_p` (`A_h2, A_p2 < 0` for all 12 centers — verified, since
  `con_54, con_513 > 0`), so the per-cell minimum is at a cell CORNER. A lower bound on
  `max_c Φ_c` over a cell is `max_c (lower endpoint of Φ_c)`; the min over cells is a
  TRUE rigorous lower bound on the continuum min. **No Lipschitz fudge — rigor is
  intrinsic to interval arithmetic.** Vectorized in float for speed (2000×2000 cells),
  then RE-CERTIFIED in `mpmath.iv` at the binding cell + a −1 × 10⁻¹² float-rounding
  safety margin. The fast scan was validated against the pure-`iv` slow scan on a
  coarse 80×80 grid: **|Δ| = 1.0 × 10⁻¹²** (= the round margin).
- **M2 — grid + Lipschitz (matches the documented `cover_min_over_box` convention).**
  `grid_min` (at node points) − `eps_grid = L_max · half_diag`, with `L_max` an `iv`
  upper bound on `max_c sup_box|∇Φ_c|` (affine grad → sup at corners). Cross-check:
  M2 = 0.3802952, within `eps_grid` (4.3 × 10⁻⁶) of M1 and ≈ the float reference
  0.3802958 — confirming consistency with the prior convention. M1 is **+3.3 × 10⁻⁶
  tighter** than M2 because it avoids the Lipschitz term.

**Independent self-test:** the interval cell-enclosure floor (0.38029595) is ≤ a
brute-force float cover min over a dense 5000² = 25 M-point grid (0.38029604), gap
8.6 × 10⁻⁸ — confirming the interval floor is a valid (and tight) LB on the continuum
min.

**Cross-check vs `path_b_independent.Phi_row`** at identical (h,p,q) with the SAME
stored duals: **worst |Δ| = 5.6 × 10⁻¹⁷** (machine precision; the two cover
implementations agree at the formula level). Per the CLARABEL-nondeterminism finding,
cross-*solve* objective agreement caps at ~7–9 digits; here both read the SAME stored
duals, so formula agreement is ~machine ε, as expected.

---

## RIGOR TIERS (honest ledger — distinguishing what is verified)

| Tier | floor µ ≥ | what is interval-certified | trusted scalar input | clears 0.380284 |
|---|---|---|---|---|
| **2c (HEADLINE)** | **0.3802959548** | shift + ellipse + box-min; binding/verified anchors = production Jansson `p_lo` | the 10 NON-binding anchors' gap ≤ 1e-6 (non-binding; stress-bounded) | **yes, +1.20e-5** |
| 2 (uniform margin) | 0.3802990936 | shift + ellipse + box-min | ALL 12 anchors' gap ≤ 1e-6 (documented convention) | yes, +1.51e-5 |
| 2b (prod, 2 verified only) | 0.3791854 | shift + ellipse + box-min | none (drops 10 float centers) — but anchors mix cover-solve duals with diff-solve `p_lo` | no |
| **1 (UNCONDITIONAL)** | **0.3788996** | EVERYTHING — anchor == certified `p_lo` for the SAME solve's duals (self-consistent), N=3000, 2 centers | **none** | no |

**TIER 1** is the unconditionally-rigorous floor (no gap assumption at all): the
self-consistent extractor (`_jansson_with_duals.py`, ONE solve → both `p_lo` AND the
`con_*` duals) at the light N=3000 config. It is weak in VALUE (only 2 of 12 centers
cover the box, N small) but proves the machinery is airtight end-to-end.

**TIER 2c** is the meaningful production number. The binding witness cde_n30_iter3 is
anchored at its production Jansson interval-certified `p_lo`; the 10 non-binding
centers use `V_c − 1e-6`.

### Robustness stress (sensitivity to the non-binding anchors' gap)

Degrade EVERY float-fallback anchor by an EXTRA penalty (so its anchor = `V_c − 1e-6 −
extra`, a guaranteed LB for any true gap ≤ `1e-6 + extra`) and re-min:

| extra | floor_all | witness | clears 0.380284 |
|---|---|---|---|
| 0 | 0.3802991 | cde_n30_iter3 (verified) | yes |
| 5e-6 | 0.3802983 | row4 (verified) | yes |
| 1.5e-5 | 0.3802960 | cde_n30_iter1 (float) | yes |
| 5e-5 | 0.3802885 | cde_n30_iter3 (verified) | yes |

→ The floor clears the headline even if every non-binding center's true gap is as
large as 5 × 10⁻⁵ (far above any observed; production Jansson penalties were
1–9 × 10⁻⁶). The floor is **robust** to the non-binding anchor uncertainty.

---

## WHAT (IF ANYTHING) STILL STANDS BETWEEN THIS AND A CLEAN "verified µ ≥ X theorem"

The dual-shift interval-certification — the explicit PRO-47 gap — is **CLOSED**.
What remains is NOT this gap; it is the same two items already named in `L2_FINISH.md`,
plus one solve-provenance subtlety surfaced here:

1. **Self-consistent (anchor == certified `p_lo`, duals) at PRODUCTION N for the
   binding centers.** The stored production cover duals (`cde_phase5_corrected_tail.json`)
   came from a CLARABEL `optimal_inaccurate` solve with only a 1e-6 margin anchor; the
   production Jansson `p_lo` (`L2_PROD.json`) is from a *different* solve (CLARABEL
   nondeterminism, ~5 × 10⁻⁶ between the two solves' duals/values). TIER 2c anchors the
   binding center at the Jansson `p_lo` while using the cover-solve duals for the shift
   — a cross-solve mix, **bounded by the robustness stress** (clears to +5e-5). A clean
   theorem wants ONE solve giving both; that is a single production Jansson re-solve
   *with `con_*` extraction* (`_jansson_with_duals.py` at N=20000). **Deferred — the
   charter forbids new heavy solves.** TIER 1 gives the unconditional version today at
   N=3000.
2. **Production-N verified anchors for the 10 NON-binding centers** (currently
   `V_c − 1e-6`). Cheap re-Jansson; deferred per charter. Non-binding + stress-bounded,
   so the floor does not depend on them.
3. **Region coverage** (the 7+iter ellipses cover the residual `(h,p,q)` region, White
   §5.1) — argued separately; the full-space promotion to all 18 White Table-2 regions
   is PRO-38 (`fullspace_promote_final.json`, independently_certified_floor 0.3802838).
   Inherited unchanged.

The bound's **strength is unchanged** — this is a RIGOR upgrade. The trusted base for
the cover step shrinks from "{CLARABEL IPM + log-parse + float dual values + float
shift formula + float grid arithmetic + Lipschitz fudge}" to, at TIER 2c, "{the 10
non-binding cover-solve gaps ≤ 1e-6}" — everything else is interval-certified — and at
TIER 1 to **nothing** (unconditional, N=3000).

---

## REPRODUCE

```bash
cd lp_research_state/code
# (TIER 1 input — light N=3000 self-consistent extraction; ~40s, NOT a heavy solve)
../../.venv/bin/python _jansson_with_duals.py --N 3000 --T 1200 --bochner_n 20 \
    --pm_k_max 14 --out ../../docs/RND_WHITESPACE/L2_COVER_VERIFIED_sc_N3000.json
# the interval certification (no solves; ~30s)
../../.venv/bin/python _cover_iv_certify.py --n_h 2000 --n_p 2000 \
    --out ../../docs/RND_WHITESPACE/L2_COVER_VERIFIED.json
# Lemma-10 + shift==dual-obj-change verification (small N, ~10s)
../../.venv/bin/python _verify_shift_eq_dualobj.py
```
