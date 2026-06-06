# Adversarial Assessment — Three-point / k-point (triple-correlation / bispectrum) SDP lift

**Date:** 2026-06-06
**Workflow:** goal-trial-1, deep-dive vetting of ranked approach #2 (MERGED across 3 lenses).
**Status:** PROPOSAL ASSESSMENT — not a result. No bound is claimed here.
**Verdict:** **promising** (genuine, principled, non-redundant lever — but with one specific failure
mode that is *more* dangerous here than the proposal's own risk section admits, and a tail-bound
rigor trap that is structurally identical to the retracted Lasserre-2).

---

## 0. What is actually being proposed (precise object)

Lift White's purely **2-point** Fourier program to a **3-point** relaxation. The genuine
new decision object is the **bispectrum / triple correlation** of the single density `f`:

```
R3_f(s,t) = ∫ f(x) f(x+s) f(x+t) dx      (real-space triple correlation)
B(m1,m2)  = f̂(m1) f̂(m2) conj(f̂(m1+m2)) (Fourier-side bispectrum; the 3rd-order moment)
```

Build a **symmetrized 3-point Gram / moment block** over a finite inner-shift grid, impose
PSD **alongside** (never replacing) the existing 2-point Bochner block + the load-bearing
cell-envelope. Run the **finite small-N exact** version first (no tail bound), detect gain,
only then add the Henrion–Rudi periodic-Hausdorff truncation (arXiv:2401.07734) for rigor and
scale via symmetry reduction.

**Refs the proposal cites (all real and correctly attributed):** Bachoc–Vallentin *JAMS* (2008)
/ arXiv:math/0608426; de Laat–Vallentin *Math. Program.* (2015) arXiv:1311.3789;
Cohn–de Laat–Salmon arXiv:2206.09876; three-point bounds arXiv:2206.15373;
ClusteredLowRankSolver arXiv:2202.12077 / SDPB; Henrion–Rudi arXiv:2401.07734.

---

## 1. FEASIBILITY with available tools

**Verdict: feasible to *prototype* (torus, small T₀) in cvxpy/CLARABEL in <1 day; the *rigorous*
version is hard and is the real cost.**

- **The lifted SDP fits cvxpy.** Introduce a complex array `b[m1,m2]` (the bispectrum) as new
  variables; impose (a) the **linking/localizing PSD tie** that connects `b` to the program's
  existing `(a_m,b_m) = (2 Re f̂(m), −2 Im f̂(m))`, and (b) PSD of a Gram block built from `b`.
  This is the standard moment-lift pattern; the project already wired an analogous lift in
  [`lp_research_state/code/lasserre3_toeplitz.py`](../../lp_research_state/code/lasserre3_toeplitz.py)
  (a trilinear/quadrilinear `W_m, Q_m` lift with bordered LMIs), so the cvxpy plumbing is
  known-good and re-usable.
- **Size.** A bispectrum block indexed by `|mᵢ| ≤ T₀` is `O(T₀²) × O(T₀²)`. At `T₀ = 6–8`
  (~36–64 index pairs) this is a ~64-dim PSD block — trivial for CLARABEL, even joint with the
  existing N=300–3000 program. Tractable as a *probe*.
- **GMP / SDPA-GMP available** for the precision cross-check the project standardizes on
  ([`lp_research_state/bin/sdpa_gmp`](../../lp_research_state/bin/sdpa_gmp)). The kissing-number
  3-point literature *required* SDPA-GMP/QD (Mittelman–Vallentin; Machado–Oliveira) — so GMP is
  not optional at scale, and the repo already has it built.
- **Symmetry reduction** (the merged Approach-2 enabler) is pure linear algebra
  (de Klerk–Pasechnik–Schrijver regular-\*-representation); feasible but unbuilt — see
  [`NA_optimization_hierarchies.md`](NA_optimization_hierarchies.md) Approach 2.
- **The expensive, research-grade piece** is the **rigorous tail bound on the truncated cubic
  moment sum** (Henrion–Rudi periodic Hausdorff). The interval/periodic-boundary work (µ lives on
  `[0,2]`, Henrion–Rudi is on the torus) is real reformulation effort the paper itself flags.

**Effort:** **medium for the tail-free torus probe (~1 day to a go/no-go signal); high–to–research
for a *rigorous* `µ`-bound** (the periodic-boundary reformulation + a provable cubic tail bound +
GMP-scale solves + ellipse re-cover are each non-trivial, and the last two are the same class of
work as the rigor track).

---

## 2. PRIOR ART

**The technique is real and powerful — but its power is documented in a setting whose source of
strength does NOT obviously transfer to this abelian 1-D problem.**

- **Sphere packing / kissing numbers (the cited precedent).** Bachoc–Vallentin (*JAMS* 2008,
  arXiv:math/0608426) introduced the 3-point SDP bound; it gives many best-known kissing-number
  upper bounds. Improved/scaled by Mittelman–Vallentin (SDPA-QD/GMP, d≤14) and Machado–Oliveira
  (symmetry reduction + SDPA-GMP, d≤16; "Improving the SDP Bound for the Kissing Number by
  Exploiting Polynomial Symmetry"). de Laat–Vallentin (*Math. Program.* 2015, arXiv:1311.3789)
  generalized to a t-point hierarchy for packing graphs. **Crucial structural caveat (see §5):**
  in these problems the 3-point gain is driven by the **non-abelian** symmetry group `SO(n)` and
  the geometry of point configurations on a sphere — the 3-point function lives on a *stabilizer
  subgroup* and is analyzed via `SO(n−1)` representation theory (Gegenbauer → matrix-valued
  zonal kernels). That machinery is what makes the 3-point relaxation *strictly* and *usefully*
  larger than the 2-point Delsarte LP.
- **Cohn–Elkies → Bachoc–Vallentin "broke the plateau" / non-sharpness (Cohn–de Laat–Salmon
  arXiv:2206.09876, *Adv. Math.* 2024).** Correctly cited; this is the strongest piece of
  evidence that 2-point ceilings *can* be pierced by higher-arity relaxations. Already flagged in
  [`docs/archive/OUT_OF_BOX_CROSS_DOMAIN.md`](../archive/OUT_OF_BOX_CROSS_DOMAIN.md) §2a as
  "never-built" here, though that memo mis-routed it into a non-sharpness *diagnostic* rather than
  a value attack.
- **1-D autocorrelation extremal constants (the genuinely adjacent class — and the telling gap).**
  A targeted literature sweep (merit factor / LABS / flat Littlewood polynomials: Borwein–Choi–
  Jedwab arXiv:1205.0626; Borwein–Ferguson–Knauer; Günther–Schmidt; BBMST *Ann. Math.* 2020
  arXiv:1907.09464) surfaced **no instance where a bispectrum / triple-correlation SDP delivered
  an improved LOWER bound for a 1-D autocorrelation extremal constant.** The triple
  correlation / bispectrum is a standard *signal-processing* object (Bartelt–Lohmann phase
  retrieval; translation-invariant), but it has **not** been turned into a working
  *lower-bound SDP* for this problem family. That absence is meaningful negative evidence: if
  the 3-point lift were an easy win for flat-autocorrelation problems, the merit-factor community
  (which has chased these constants for 40 years) would plausibly have used it.

**Has any piece been tried *here*?** Partially, and the partial attempt is a warning:
- [`lasserre3_toeplitz.py`](../../lp_research_state/code/lasserre3_toeplitz.py) built a
  **trilinear/quadrilinear moment lift** — but on the program's **auxiliary inner variables**
  `(c_k, d_k)` (the Fourier coeffs), *not* on the genuine bispectrum `B(m1,m2)` of `f`. Its own
  docstring records the key fact: by **Fejér–Riesz**, climbing Lasserre levels "on the circle"
  gives **no** new constraints — the localizing matrix of order `T_loc` already captures all
  nonneg degree-`2T_loc` trig test functions. This is exactly why the naive 3-point version fails
  (§3). The genuine-bispectrum lift was **never built** (no `bispectr` / `R3` / `B(m1,m2)` traces
  in `lp_research_state/code/`; only the inner-variable `lasserre3` lift exists).
- The **Lasserre-2** density-side attempt is **RETRACTED** for a bad tail bound
  (`communications/lasserre_tail_bound.md`; megamemory "Tail-bound rigor trap"). The proposed
  3-point lift is **distinct in object** (a genuine cubic/3-point moment, not a degree-lift of the
  pairwise object) — but it inherits the **same rigor trap** (a truncated cubic sum needs an
  analytic remainder; see §5.3).

**Near-duplicate within this very workflow:** this approach is the *same object* as Approach 1 of
[`NA_optimization_hierarchies.md`](NA_optimization_hierarchies.md) (de Laat–Vallentin 3-point
measure hierarchy + Henrion–Rudi tail, with the identical `B(m1,m2)` first experiment) and HA-2's
companion certificate in [`harmonic_analysis_lens.md`](harmonic_analysis_lens.md). The
[`RANKED_SHORTLIST_2026-06-06.md`](RANKED_SHORTLIST_2026-06-06.md) already merged all three into its
#2. This memo is the adversarial deep-dive on that merged item; it does not introduce a fourth copy.

---

## 3. FIRST EXPERIMENT (cheapest signal) — and a result from running the cheap analytic check

The proposal's nominated first experiment (a small-T₀ torus SDP solve at row 4) is the right
*eventual* probe, but it is a **solve**, not a literature/analytic check, so per task discipline I
did **not** run it. Instead I ran the **tail-free analytic precursor** that must pass *before* the
solve is even worth wiring — and it returns a sharp, decision-relevant signal.

**Analytic probe (run, tail-free, no SDP):** *Is a "3-point moment matrix" an independent
constraint, or does it collapse onto the 2-point Bochner/Toeplitz block already in the program?*

**Finding (load-bearing — this is the crux of the whole approach):**
- A moment matrix indexed by **pairs** `(i,j)` with entries `E[χ_{i+j} \overline{χ_{k+l}}] =
  f̂((i+j)−(k+l))` has entries that are `f̂` at index **differences** ⇒ it is **Hermitian–Toeplitz
  ⇒ it IS exactly the Bochner block** already in
  [`white_full_convex.py`](../../lp_research_state/code/white_full_convex.py). Numerically its
  min-eigenvalue tracks the existing 2-point Bochner matrix (probe: min eig `≈ −1e-16`, i.e. the
  pair-indexed Gram is the same PSD object, adding nothing). **The "easy" 3-point matrix gives
  ZERO gain** — and this is precisely the Fejér–Riesz collapse the `lasserre3` docstring warns of.
- The **genuine** 3-point content is the **cubic** moment `B(m1,m2) =
  f̂(m1) f̂(m2) \overline{f̂(m1+m2)}`. This is **NOT a linear functional of the program's decision
  variables** `(a_m,b_m)`. To impose its PSD-ness you must **introduce `B` as NEW lifted
  variables**, add a **linking PSD tie** to `(a_m,b_m)`, and — for any *rigorous* bound — supply a
  **tail bound on the truncated cubic sum** `Σ_{|m|>T₀} (cubic)`.

**Consequence for sequencing.** The decisive *free* questions are therefore answerable *before*
any heavy solve, and they sharply gate the whole approach:
1. **(Done, above)** Does the naive pair-indexed block add anything? **No** (it is the existing
   Bochner block). ⇒ Any real test must lift the genuine cubic `B`.
2. **(Cheap, ~½ day, recommended next, still tail-free)** Wire the **genuine** lifted-`B` block on
   the **torus** at `T₀ = 6` jointly with the existing Bochner level-1, solve at row 4 small-N,
   and check the two go/no-go conditions the proposal itself names: **(i)** does the dual objective
   **exceed** the row-4 2-point value? **(ii)** is the `B`-block **active** (nonzero dual)? Plus the
   **mandatory PRO-22 validity gate**: reconstruct `f`, confirm `f ∈ [0,1]` and `Ω ≤ sup_t (f*f)`.
   This is the cheapest computation that could give a *positive* signal; it is a solve, so it is
   out of scope for *this* memo but is the correct immediate next action.

**Reported finding from the cheap check:** the approach's entire viability rests on the *genuine
cubic* lift, because the obvious "3-point Gram" silently reduces to the 2-point block. This both
(a) confirms the approach is **not** a bare re-tread of the existing Bochner block, and (b) routes
it straight into the cubic-tail-bound rigor trap that retired Lasserre-2.

---

## 4. WHY IT COULD BEAT THE ~0.380558 SATURATION (the genuine upside)

This is the **single most principled candidate to pass the ceiling**, and the reasoning is sound:

- **The ceiling is provably a 2-point phenomenon.** PRO-6
  ([`docs/archive/PRO6_COMPLEMENTARITY_PROOF.md`](../archive/PRO6_COMPLEMENTARITY_PROOF.md))
  proves the joint Bochner+cell-envelope stack saturates at `C_∞ ≈ 0.380558` because **both**
  tightenings are pairwise: the full augmentation only ever constrains the **degree-2 moments
  `|f̂(m)|²`** of a single density. PRO-22
  ([`docs/archive/LEVER_SUPT_DIRECT.md`]) shows the cell-envelope is the load-bearing **2-point
  Parseval link** `M̂(m) = a_m f̂(m) − 4|f̂(m)|²` and that dropping it is invalid (off by ~5× the
  gap). So the ceiling is a property of the **2-point moment cone**, not a tractability artifact.
- **A genuine 3-point relaxation lives in a strictly larger cone.** The bispectrum `B(m1,m2)`
  constrains the **phases** `arg f̂(m)` *jointly* (`arg(m1)+arg(m2)−arg(m1+m2)`), which the power
  spectrum `|f̂(m)|²` — and hence the entire current program — leaves **free**. The
  Bartelt–Lohmann theorem (bispectrum phase + magnitudes determine the signal up to translation)
  confirms `B` carries strictly more information than `|f̂|²`. This is the *necessary* condition
  for a gain: there genuinely exist `(a_m,b_m)` the 2-point program admits whose phase
  configuration a bispectrum-PSD constraint could exclude.
- **The sphere-packing precedent is exactly this move.** Cohn–Elkies (2-point LP) → Bachoc–
  Vallentin (3-point SDP) is the canonical case of higher arity breaking a 2-point plateau and
  *proving non-sharpness* (d=12,16 via Cohn–de Laat–Salmon). The structural analogy —
  `R_h(k)=⟨h,T_k h⟩` is intrinsically 2-point; triples encode *which inner shifts can
  simultaneously be low* — is the right intuition, and it is the **axis the retracted Lasserre-2
  did NOT move** (Lasserre-2 lifted the *degree* of the *same pairwise* object; the Fejér–Riesz
  collapse means that gives nothing new on the circle, as `lasserre3`'s own docstring states).
- **Convergence is a theorem.** Unlike the retracted Lasserre attempt, the de Laat–Vallentin
  t-point hierarchy *provably converges* to the true extremal value (the Euclidean analogue of
  flag-algebra / theon convergence). So the question is *at which level* the gain appears, not
  *whether* the hierarchy is sound.

**In one line:** the upside is real and well-founded — `µ`'s `min-over-translations` / `L∞`
structure is exactly where higher-arity moments are *expected* to bite, precisely because the
`L²` sibling (White's autoconvolution, arXiv:2210.16437) is "almost tight" at degree 2 while the
`L∞` object is not.

---

## 5. RISKS / why it might fail (adversarial — including risks the proposal under-weights)

**5.1 The transfer risk is bigger than "the cut may be implied" (the proposal's stated risk).**
The sphere-packing 3-point gain is powered by the **non-abelian `SO(n)`** symmetry: the 3-point
function decomposes via `SO(n−1)` zonal matrix kernels, and *that* representation theory is what
makes the 3-point cone usefully larger than the 2-point LP. The Erdős problem's symmetry is the
**abelian** translation group on an interval (plus a `Z/2` reflection). Over an abelian group the
"3-point" structure is much flatter — there is no rich stabilizer-subgroup representation theory to
exploit, and the bispectrum's extra content is "only" the joint phase relations. **Whether joint
phase constraints actually move a *lower bound* for a translation-invariant problem is exactly the
question with zero positive precedent in the 1-D autocorrelation literature (§2).** This is the
deepest risk: the cited precedent's *mechanism* may not be present here.

**5.2 The naive version is a guaranteed null (now established analytically, §3).** The obvious
"3-point Gram" is the existing Bochner block (Fejér–Riesz collapse). So the approach *cannot* be
tested cheaply by the naive matrix — it *must* lift the genuine cubic `B`, which immediately raises
the cost and the rigor bar. Any positive signal requires the harder object from the start.

**5.3 The cubic-tail rigor trap — structurally identical to the retracted Lasserre-2.** A *finite*
torus probe is tail-free and can detect a gain, but a **rigorous bound on `µ`** needs the truncated
cubic moment sum `Σ_{|m|>T₀}` to carry an **analytic remainder**. This is the *same* failure mode
that retired Lasserre-2 (megamemory "Tail-bound rigor trap"; `communications/lasserre_tail_bound.md`,
which derived the Fejér–Riesz tail bound and showed it *kills the gain* at tractable `T_max`). The
Henrion–Rudi periodic-Hausdorff bound (arXiv:2401.07734 Prop. 7) is the proposed principled
replacement — but (i) it is *periodic* and µ lives on `[0,2]` with a boundary (real reformulation
work), and (ii) **there is a live danger the rigorous cubic tail eats the gain at the same `T₀`
where the gain first appears**, exactly as the Fejér–Riesz tail did for Lasserre. The probe must
therefore report not just "gain at finite T₀" but "gain *net of the Henrion–Rudi tail* at a
tractable T₀" — and the Lasserre precedent says that net is the thing most likely to vanish.

**5.4 Validity (PRO-22 gate) must be re-derived for the lifted encoding.** PRO-22 proved a naive
higher-moment SDP that drops the cell-envelope is **invalid** (it under-counts `sup_t (f*f)` by ~5×
the gap). The lifted-`B` block bolts *onto* the cell-envelope (correct by construction in the
proposal), but the **linking tie** between `B` and `(a_m,b_m)` is a new place validity can leak;
the reconstruct-`f`-and-check-`Ω ≤ sup_t(f*f)` gate is mandatory and non-negotiable.

**5.5 Scale.** Even if a gain survives at small `T₀`, the bispectrum block is `O(T₀²)×O(T₀²)`; the
honest rigorous version may not show net gain until `T₀` is past tractable size — "same wall, one
level up" (the proposal's own risk #1, correctly stated). Symmetry reduction (merged Approach 2)
is the mitigation but is itself unbuilt, and over an abelian group buys less than it does in the
`SO(n)` setting where it was developed.

**5.6 Bound-neutral negative is a likely outcome.** The most probable result is one of: (a) the
finite probe shows the `B`-block is **inactive** at row 4 (no gain — stop, ~½ day spent); or (b)
it shows a finite-T₀ gain that the rigorous tail erases (informative but bound-neutral, and a
*second* instance of the Lasserre pattern). A genuine net LB advance is the high-variance tail of
the distribution, not the mode.

---

## 6. Bottom line

| Dimension | Assessment |
|---|---|
| **Verdict** | **promising** — most principled ceiling-piercing candidate; real precedent; non-redundant with the ledger; but transfer mechanism is unproven for abelian 1-D and the rigorous version inherits the Lasserre tail trap. |
| **Feasibility** | Tail-free torus probe: medium, <1 day. Rigorous µ-bound: high → research (periodic-boundary reformulation + provable cubic tail + GMP scale + ellipse re-cover). |
| **Prior art** | Genuine and strong in *sphere packing/kissing* (Bachoc–Vallentin JAMS 2008; Machado–Oliveira; de Laat–Vallentin 2015; Cohn–de Laat–Salmon 2024). **Zero positive precedent** for a bispectrum SDP improving a *lower* bound on a *1-D autocorrelation* extremal constant. Here: only the inner-variable `lasserre3_toeplitz.py` lift exists; the genuine bispectrum lift is unbuilt; Lasserre-2 (adjacent, retracted) flags the tail trap. |
| **Why it could beat saturation** | The ceiling is a *proven* 2-point phenomenon (PRO-6/PRO-22). The bispectrum constrains joint phases the power spectrum leaves free, lives in a strictly larger cone, has a *convergence theorem* (de Laat–Vallentin), and follows the exact LP→SDP arity jump that broke the sphere-packing 2-point plateau. |
| **Effort** | Medium for go/no-go; high–to–research for a rigorous bound. |
| **Decisive next step** | The cheap analytic check is **done** (naive 3-point Gram = existing Bochner block; must lift genuine cubic `B`). Next: the **genuine-`B` torus probe** at `T₀=6`, row 4, small-N — check (i) dual > 2-point value, (ii) `B`-block active, (iii) PRO-22 validity gate — *before* any tail-bound or scaling work. |

**Honest framing for the orchestrator:** this is the right #2. It is the only proposal in the
slate that can *transcend* (not merely approach) `C_∞ ≈ 0.380558` with theory behind it. But it is
**not** a quick engineering win: the naive version is a proven null, the genuine version is a cubic
lift that walks straight into the same tail-bound rigor trap that retracted Lasserre-2, and its
cited precedent draws its power from non-abelian symmetry that this abelian 1-D problem lacks. Fund
the ~½-day genuine-`B` torus probe as a cheap, decisive go/no-go; treat a positive net-of-tail
signal as the high-variance prize, and a bound-neutral "B-block inactive / tail eats the gain"
outcome as the modal (and still useful) result.
