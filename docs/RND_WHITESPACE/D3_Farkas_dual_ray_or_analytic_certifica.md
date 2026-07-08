# D3 — Farkas / dual-ray (or analytic) certificates for the infeasibility-excluded corners

**Scout date:** 2026-06-03. **Status:** SCOUT ONLY — nothing claimed proven.
**Verdict (one line):** DROP as a path to "make 0.380284 a theorem" — the seam it closes is already certified non-load-bearing; KEEP a 1-page analytic lemma as cheap preprint hygiene.

---

## 0. The direction restated

Replace the *solver-attested* infeasibility of the deep-q / high-p corners of White's wide
outside regions (R6/R7/R8/R9, and the high-E(M) regions R1–R5) with explicit
Farkas / dual-ray certificates, ideally exact. Proposed analytic handle: the program
enforces `|d| ≤ 2/π` and `sum_squares(c) + sum_squares(d) ≤ 1/2`, so `d1 = d[0]` is bounded
and feasibility "provably collapses for `|d1|` beyond ~0.5", giving a single analytic
infeasibility proof covering the bulk of the excluded mass.

Payoff class claimed: **1** (close the last empirical seam → contribute to promoting the
numerically-certified `µ ≥ 0.380284` to a theorem).

---

## 1. FEASIBILITY — tractable, but the easy part is trivial and the hard part has no exposed ray

### 1a. The analytic handle is REAL and even cleaner than stated (verified)

The relevant constraints in `lp_research_state/code/white_full_convex.py`
(and the dual-handle build `path_b_analytical.build_problem_with_dual_handles`):

```
203:  cons += [cp.abs(c) <= 2.0/np.pi, cp.abs(d) <= 2.0/np.pi]    #  |c_k|,|d_k| ≤ 2/π
204:  cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)   #  Σ c² + Σ d² ≤ 1/2
205:  cons += [c[0] >= p1, c[0] <= p2, d[0] >= q1, d[0] <= q2]    #  c[0]=p∈[p1,p2], d[0]=q∈[q1,q2]
```

The *binding* analytic bound on `d1 = d[0]` is **`|d[0]| ≤ 2/π ≈ 0.63662`** (line 203), which
is *tighter* than the `sum_squares ≤ 1/2 ⇒ |d[0]| ≤ √0.5 ≈ 0.7071` the direction names. So:

- **Any sub-box with `q1 > 2/π` (≈ 0.637) is infeasible by a one-line box-vs-box contradiction:**
  the parameter box demands `d[0] ≥ q1 > 2/π`, but line 203 demands `d[0] ≤ 2/π`. Disjoint.
  Exact, rational, no solver. This is a *legitimate analytic infeasibility certificate*.
- **Symmetric free bonus the direction missed:** `c[0] = p` and `|c[0]| ≤ 2/π` ⇒ **any sub-box with
  `p1 > 2/π` is infeasible by the identical contradiction.** This covers the high-p tails of
  R1–R5, R8, R13 (p∈[0.5,1.0]), etc. (Verified: the feasmap's `p=0.70,q=0.05` and `p=1.0,q=0.05`
  infeasibilities are exactly this slab — see §First-step probe.)

These two slabs (`p > 2/π` ∪ `|q| > 2/π`) are trivially, exactly certifiable. No tooling needed
beyond writing the inequality chain by hand.

### 1b. …but the trivial slabs cover only the FAR corners; the real infeasibility starts much earlier

Empirically (existing artifact `code/_verify_R8_feasmap.py`, reproduced this session at small N):

| point (h,p,q)        | status     | note |
|----------------------|------------|------|
| (0, 0.37, 0.30)      | optimal    | feasible, val 0.515 |
| (0, 0.37, 0.35)      | solver-fail (near boundary) | |
| (0, 0.37, **0.40**)  | **infeasible** | **q = 0.40 ≪ 2/π = 0.637** |
| (0, 0.37, 0.50)      | infeasible | |
| (0, 0.0, 0.40)       | optimal    | feasible at p=0 |
| (0, 0.0, 0.60)       | infeasible | q=0.60 < 2/π |

So the infeasibility boundary in the `(p,q)` plane is a **curve around `q ≈ 0.35–0.60` that
depends on p**, well *inside* the trivial slab `q > 0.637`. The region

> `{ |q| ≤ 2/π } ∩ { p ≤ 2/π } ∩ { infeasible }`  (e.g. p≈0.37, 0.40 ≤ q ≤ 0.63)

is genuinely infeasible but **NOT** reachable by either box contradiction. It is driven by the
*interaction* of `con_513` (line 209, the cosine-moment lower bound `(L/2)(a₊₂·(w+v)) ≥ rhs_513`),
the quadratic cell-consistency constraints (lines 185, 193), and the `sum_squares` ball — a
nonlinear, multi-constraint infeasibility. There is **no one-line analytic handle** for this slice.
A Farkas certificate here is a nonnegative combination of ~80 constraints (including PSD/SOC
Bochner blocks) — derivable in principle, but not "an exact inequality chain" you write by hand.

**Quantitative split (verified, §First-step probe):** of White's 18 Table-2 regions, **zero are
fully covered by the trivial slabs.** Every wide region (R1–R9) has a feasible low-p / low-|q|
part that needs the cover regardless, plus a "hard" infeasible annulus between the empirical
boundary (~0.4) and the trivial slab (0.637) that the analytic handle does not reach.

### 1c. The solver does NOT expose a Farkas ray (verified — kills the proposed first step)

The direction's first step is "ask the solver for the infeasibility certificate / dual ray."
**This does not work with the current stack.** On the infeasible corner (0, 0.37, 0.50) at N=1200:

- **CLARABEL:** `status=infeasible`, and **0 of 82 constraints** carry a nonzero `dual_value`.
- **SCS:** `status=infeasible`, likewise **0 of 82** nonzero duals.

cvxpy does not populate dual variables with the infeasibility certificate for these conic
solvers here; there is no ready-made dual ray to read off. Extracting a usable Farkas ray would
require either (i) going under cvxpy to the raw solver certificate (CLARABEL exposes a primal/dual
infeasibility certificate in its native Rust/Python API, not surfaced through `Constraint.dual_value`),
or (ii) re-solving a dedicated Farkas LP/SDP (`find y ≥ 0 : Aᵀy = 0, bᵀy < 0` over the cone),
then rounding `y` to rationals and verifying `Aᵀy = 0 ∧ bᵀy < 0` in exact arithmetic. Both are
real engineering, not "read the dual_value." This is the same missing-infrastructure problem as
the cvxpy→SDPA-S serializer (D1): the certificate is *in there*, but not extractable as-is.

### 1d. Tooling verdict

- Trivial slabs (`p,|q| > 2/π`): trivially tractable, exact, by hand. ✅
- Hard annulus Farkas certificate: tractable only via a dedicated Farkas-LP re-solve + exact
  rounding (mpmath/python-flint) + verification; ~1–2 files; SDPA-GMP could verify the rounded
  ray at high precision but the *extraction* path must be built first. Medium effort.
- Solver-emitted ray via cvxpy: ❌ not available.

---

## 2. PRIOR ART — already analyzed and explicitly retired as non-load-bearing

This is the decisive finding for payoff. The infeasibility exclusions have **already been
investigated** in this repo, and the canonical verification memo
[`lp_research_state/FULLSPACE_VERIFICATION.md`](../../lp_research_state/FULLSPACE_VERIFICATION.md)
(§"Infeasibility-exclusion rigor — NOT load-bearing anywhere", lines 123–133) states verbatim:

> "Several region reports characterize deep-q / high-p corners as SDP-infeasible (solver-attested
> CLARABEL 'infeasible' at multiple interior points; **NO Farkas/dual-ray certificate extracted —
> rigor level = empirical-robust-multipoint, NOT certificate-grade**). **This rigor gap does not
> affect the verdict.** In every region the certified floor is set by the FEASIBLE part of the box
> where the cover already clears the target on pure grid+Lipschitz geometry; the excluded regions
> have cover values far above target. **No region's floor depends on excluding an infeasible
> sub-box.**"

Concretely (memo per-region table, lines 142–152):

| region | repro grid_min (true cover inf) | infeas rigor |
|--------|---------------------------------|--------------|
| R5 | 0.385981 | solver-attested, **not load-bearing** |
| R6 | 0.381652 (strip 0.3816) | solver-attested, **not load-bearing** |
| R7 | 0.380594 | solver-attested, **not load-bearing** |
| R8 | 0.382416 | solver-attested, **not load-bearing** |
| R9 | 0.380714 | none (whole box feasible) |
| R16/R17 | 0.3804–0.3805 | none (whole box feasible) |

And §"R6 / R8" (lines 114–121): the deep-q infeasibility-exclusion is "a **RED HERRING** for the
floor in both regions (the excluded deep-q points have cover/SDP values 0.47–1.6 ≫ target …)."

The megamemory verification node (`pro-38-full-space-candidate-…-r16-refute-overruled`) says the
same: *"INFEASIBILITY EXCLUSIONS ARE NOT LOAD-BEARING anywhere (so verdict is CONFIRMED, not
PARTIAL)."* And `_promote_R8_centers.py` already frames the R8 strategy as "cover the feasible
sub-box; the complementary part is cleanly SDP-infeasible (vacuous)."

**Why "not load-bearing" is logically airtight (the cover is an upper-bounding envelope):**
The per-region floor is `min over the box of  Φ_cover(h,p,q)`, where `Φ_cover = max over dual
centers` is a **valid lower bound on the true SDP value `V(h,p,q)` at every feasible point**, and
at infeasible points `V = +∞ ≥ Φ_cover`. Taking the min of `Φ_cover` over the *whole* box
(including would-be-infeasible corners) is therefore a *conservative* lower bound on
`min over feasible points of V` — the quantity that actually bounds µ. Including infeasible
points can only *lower* the reported floor (never inflate it), so if the whole-box `grid_min`
already clears target (it does: 0.380594–0.385981 ≫ 0.3802838), the bound is valid **without any
infeasibility argument at all.** The exclusions were a presentational convenience, not a logical
load.

**Net prior-art status:** D3's target seam was identified, characterized as "not certificate-grade,"
and *deliberately retired as immaterial* by the verification that established µ ≥ 0.380284. There is
no Farkas-extraction code in the repo (`grep Farkas/dual_ray` → only prose mentions), so the *work*
is genuinely undone — but the *reason it was left undone* is that it doesn't matter.

**External lit (light):** Farkas/infeasibility certificates for conic programs are standard
(Permenter–Parrilo facial-reduction; Klep–Schweighofer Positivstellensatz infeasibility; VSDP /
Jansson rigorous a-posteriori bounds, already on the L2 lever list in the external scan). None of
this is novel, and none changes the load-bearing analysis. arXiv:2405.13625 (Kolmogorov–Naldi–Zapata,
degenerate-SDP feasibility certs) is the closest modern reference, already filed as "context only"
in the same-day external scan.

---

## 3. CONCRETE PLAN (if pursued anyway — scoped as preprint hygiene, not as a theorem-maker)

Ordered, smallest-first:

1. **(½ day, exact, by hand) Write the two trivial-slab lemmas.**
   - Lemma A (high-|q|): On any parameter box with `q1 > 2/π`, the program is infeasible:
     `d[0] ≥ q1 > 2/π` contradicts `|d[0]| ≤ 2/π`. ∎
   - Lemma B (high-p): symmetric with `c[0] = p`, `p1 > 2/π`. ∎
   - These are publishable one-liners; they exactly certify the *far* corners of R1–R8 / R13.
2. **(1 day) Enumerate which Table-2 sub-boxes are fully covered by Lemma A/B** (done in §1; none
   are *fully* covered, but the far-corner sub-boxes are). State, per region, the residual
   "hard annulus" `{ p ≤ 2/π, |q| ≤ 2/π, infeasible }`.
3. **(1 day, the only real value-add) Re-frame the hard annulus via the monotone-cover argument**
   instead of a Farkas certificate: prove the clean lemma "the whole-box `grid_min` of the
   conservative cover `Φ_cover` lower-bounds `min over feasible points of V`, because `V = +∞` on
   infeasible points and `Φ_cover ≤ V` on feasible ones." This *eliminates the need for any
   infeasibility certificate at all* and is the honest, minimal closure of the seam. (It is exactly
   the argument FULLSPACE_VERIFICATION.md already uses — formalize it as a one-paragraph lemma.)
4. **(only if a reviewer insists on a literal Farkas ray) Build the Farkas-LP extractor:** for a
   chosen hard corner, solve `min 0 s.t. Aᵀy = 0, bᵀy = -1, y ∈ K*` (the alternative system), round
   `y` to rationals, verify `Aᵀy = 0 ∧ bᵀy < 0` in exact/interval arithmetic (mpmath or
   python-flint), optionally re-verify the PSD/SOC dual blocks at GMP precision via SDPA-GMP. This
   is the same extraction infrastructure D1 needs; reuse it. ~1–2 files, medium effort, *and it
   buys nothing the §3.3 lemma doesn't already give.*

---

## 4. PAYOFF — class 1 in name only; near-zero marginal value toward "proven"

- **Stated payoff:** class 1, "removes the last empirical link in promoting 0.380284 to a theorem."
- **Actual payoff:** the link it removes was *already established as not part of the chain*. The
  full-space floor is set by the FEASIBLE-part cover; the infeasible corners contribute `+∞` and
  are conservatively included by the monotone-cover argument. A Farkas certificate (or the trivial
  slab lemmas) changes neither the value (0.3802838) nor its validity.
- **What "truly meaningful" outcome would it produce?** Only the §3.3 one-paragraph monotone-cover
  lemma is meaningful, and even that is *replacing* an informal sentence the verification memo
  already states. It tidies the preprint's rigor narrative; it does not make a
  numerically-certified bound into a theorem. The things that *would* make it a theorem are the
  cover's own machinery (proof-grade SDP duals via D1/D2, exact arithmetic on the *feasible*-part
  dual cover) — D3 is orthogonal to all of them.
- **Honest reclassification:** payoff class **1-cosmetic** (or "0" against THE BAR). It is preprint
  hygiene, not a step toward proof-grade rigor of the bound.

---

## 5. FIRST-STEP PROBE — RUN (this session)

Two cheap probes were executed at small N (no heavy solves):

**Probe 1 — trivial box contradiction & coverage map (`white_full_convex` constraints):**
- Confirmed `2/π = 0.63662`, `√0.5 = 0.70711`; binding bound on `|d[0]|` and `|c[0]|` is `2/π`.
- Mapped all 18 Table-2 regions: **none is fully covered** by the slabs `p>2/π` ∪ `|q|>2/π`; the
  slabs carve off only the far corners. Each wide region retains a feasible low-p/low-|q| part.

**Probe 2 — solver Farkas-ray availability (N=1200, bn=12, corner (0, 0.37, 0.50)):**
- CLARABEL: `infeasible`, **0/82** nonzero `dual_value`.
- SCS: `infeasible`, **0/82** nonzero `dual_value`.
- High-p check: (0,0.70,0.05) and (0,1.0,0.05) both `infeasible` — these *are* exactly the
  trivial `p > 2/π` slab (Lemma B), confirming the free high-p handle.
- Empirical infeasibility boundary at p=0.37 sits at `q ≈ 0.35–0.40` (val 0.515 at q=0.30, then
  infeasible by q=0.40) — i.e. **0.24 below the trivial slab `q > 0.637`**, confirming the
  "hard annulus" with no analytic handle.

**Finding:** the analytic handle exists but is partial (far corners only); the solver gives no
free ray; and — decisively — the seam is already certified non-load-bearing, so closing it does
not move the bound toward "proven."

---

## 6. FAILURE MODES (why it does not pan out)

1. **(realized) The seam is not load-bearing.** FULLSPACE_VERIFICATION.md + megamemory both
   establish that no region's floor depends on the infeasibility exclusion. A certificate here is
   immaterial to the bound. This is the dominant failure mode and it is *already true*, not a risk.
2. **(realized) The clean analytic handle is partial.** `|d|,|c| ≤ 2/π` certifies only the
   `>2/π` slabs; the empirically-infeasible annulus `0.4 < |q| < 0.637` (moderate p) is driven by
   nonlinear constraint interaction and has no closed-form handle.
3. **(realized) The solver exposes no dual ray** through cvxpy for CLARABEL/SCS — the proposed
   first step ("ask the solver for the dual ray") fails; a dedicated Farkas-LP extractor must be
   built, which is D1-grade infrastructure work for zero marginal payoff.
4. **Exact-rounding fragility** (if §3.4 pursued): rounding an interior-point infeasibility
   certificate to rationals and re-verifying `Aᵀy = 0` exactly can fail on near-degenerate corners
   (the same conditioning that makes them "optimal_inaccurate" elsewhere); mitigated by SDPA-GMP
   but, again, for no gain.
5. **Opportunity cost.** Effort spent here is effort not spent on D1 (cvxpy→SDPA-S serializer) or
   D2/L1/L4 (exact SOS/Fejér-Riesz dual on the *feasible* binding center), which are the actual
   theorem-makers.

---

## 7. VERDICT

**DROP** as a route to promoting `µ ≥ 0.380284` to a theorem: the infeasibility-exclusion seam is
already certified non-load-bearing by FULLSPACE_VERIFICATION.md (the cover's whole-box `grid_min`
clears target with `V=+∞` at infeasible points), so a Farkas/analytic certificate there changes
neither the value nor its validity. **KEEP** only the cheap byproducts for preprint hygiene: the
two one-line trivial-slab lemmas (`|d[0]| ≤ 2/π`, `|c[0]| = p ≤ 2/π`) and the one-paragraph
monotone-cover lemma that *replaces* the need for any infeasibility certificate. Route real rigor
effort to D1/D2/L1/L4 (proof-grade duals on the feasible binding center) instead.
