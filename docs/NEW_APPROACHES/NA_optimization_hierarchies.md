# New Approaches — Lens: Optimization Hierarchies (beyond Lasserre/Bochner)

**Date:** 2026-06-06
**Author:** workflow subagent (lens = optimization-hierarchies)
**Scope:** Relaxation/proof hierarchies that are NOT in the ruled-out ledger: flag
algebras / continuous-combinatorics limit objects, the de Laat–Vallentin
moment(-SOS) hierarchy on the *measure* formulation with rigorous truncation,
and representation-theoretic *symmetry reduction* of the existing Fourier SDP to
break the stated large-`T_max` tractability barrier. Copositive/CP and NPA are
discussed and (honestly) demoted.

These are PROPOSALS, vetted but not claimed. Read alongside `findings.md`,
`docs/archive/LEVER_F3_FULL_SATURATION.md` (the saturation diagnostic),
`docs/archive/LEVER_SUPT_DIRECT.md` (= PRO-22, the load-bearing cell-envelope),
and `lp_research_state/code/white_full_convex.py`.

---

## Why this lens, sharply

Every lever in the ledger is a **two-point** relaxation. The objective is
`R_h(k) = ⟨h, T_k h⟩` — intrinsically a pairwise/2-point functional — and White's
program plus all augmentations (Bochner-PSD on `f, 1−f`; poly-moment Hausdorff;
Lasserre-2 truncated; ellipse cover) only ever constrain the **degree-2 moments**
`|f̂(m)|²` of a single density. PRO-22 (`LEVER_SUPT_DIRECT.md`) proved the
cell-envelope is the load-bearing object: it is the Parseval LINK
`∫₀² M(x)cos(πmx/2)dx = (4 sin(πm/2)/πm) a_m − 2(a_m²+b_m²)` tying `M`'s Fourier
coefficients to `|f̂(m)|²`. Drop it → the bound goes invalid (off by 5× the gap).
The framework ceiling `C_∞ ≈ 0.380558` is therefore a property of the **degree-2,
two-point** moment cone, not a tractability artifact.

The two structurally distinct ways to beat a degree-2/2-point ceiling are:
(A) **go to higher-degree / higher-arity moments** of the SAME measure — a genuine
Lasserre/flag-algebra step the project never took (Lasserre-2 was on the
*wrong* object and retracted for a bad tail bound); or
(B) **keep the degree but make the SDP vastly bigger** by exploiting the symmetry
group so the effective bandwidth `T_max` jumps an order of magnitude — directly
hitting the documented "large-T intractability" barrier. Both are below.

---

## Approach 1 — de Laat–Vallentin moment(-SOS) hierarchy on the MEASURE formulation, with rigorous Hausdorff truncation (the Euclidean analogue of flag algebras)

**Technique.** Recast `µ = inf_h sup_k M(h,k)` as an infinite-dimensional
*generalized moment problem* (GMP) over a measure on the *graph of the density*,
and apply the **de Laat–Vallentin SDP hierarchy for packing/energy problems**
(an infinite-dimensional generalization of Lasserre's moment method to continuous
domains; "topological packing graphs", moment matrices `M_t(y) ⪰ 0`). Level `t`
sees `t`-point interactions of the configuration, so level `t ≥ 3` constrains
joint statistics of triples `(x, x+k₁, x+k₂)` — moments the two-point cell-envelope
**cannot represent**. Use the **Henrion–Rudi periodic moment-SOS construction
(arXiv:2401.07734)** to get the *rigorous* inner/outer truncation: it gives moment-
SOS for nonnegative measures on the unit ball of a Sobolev space of periodic
trigonometric functions with **explicit finite-truncation Hausdorff error bounds
(their Prop. 7)** — i.e. a principled, rate-equipped replacement for the
hand-rolled poly-moment tail that the retracted Lasserre attempt botched.

**Citations.**
- de Laat, Vallentin, *A semidefinite programming hierarchy for packing problems
  in discrete geometry*, Math. Program. (2015); Vallentin, *Conic optimization for
  extremal geometry* (survey), arXiv:2510.06960 (2025) — the t-point moment
  hierarchy + symmetry machinery in one place.
- Cohn, Salmon, on the Euclidean limit of the Lasserre bound (level-1 = Cohn–Elkies
  LP); Cohn–de Laat–Salmon, arXiv:2206.09876 — discrete-reduction duals.
- Henrion, Rudi (or Korda et al.), arXiv:2401.07734 (2024) — moment-SOS on Sobolev
  of periodic trig functions WITH explicit Hausdorff truncation bounds.
- Lasserre, *A semidefinite programming approach to the generalized problem of
  moments*; Tacchi et al., arXiv:2011.08139 (general-case convergence).
- (Framing) Coregliano, Razborov, *Semantic limits of dense combinatorial objects*
  (theons), arXiv:1910.08797 — the analogue statement that a hierarchy of SDPs
  converges to the true extremal limit of a continuous limit object.

**Why untested here (vs ledger).** The ledger's "Lasserre level-2" was (i) applied
to `(f²)̂(m)` of the **density** with a truncation that **had no analytic tail
bound** (RETRACTED, see `communications/lasserre_tail_bound.md`), and (ii) still a
degree-lifting of the *two-point* object. This proposal is different on both axes:
it is the **measure/GMP (de Laat–Vallentin) hierarchy**, whose level `t` encodes
genuine **t-point** configuration moments (not just higher powers of the same
2-point quantity), and it imports a **published rigorous Hausdorff truncation**
(Henrion–Rudi Prop. 7) so the tail is bounded by construction. The megamemory
rigor scout filed Henrion–Rudi only as a *certification footnote (L3)* for the
existing poly-moment cuts — never as a new-value attack on a measure hierarchy.
No `de_laat/vallentin/moment.hierarchy/3.point/theon` traces exist in
`lp_research_state/code/`. Honest cousin: `symmetric_push.py` adds a `lasserre_T_max`
knob, but that drives the same retracted truncated-`(f²)̂` object, not a t-point
measure moment matrix.

**First experiment (cheap).** Build the **3-point level** on the *torus* version
first (avoids the interval-boundary work): variables = the bispectrum-like moments
`B(m₁,m₂) = f̂(m₁) f̂(m₂) \overline{f̂(m₁+m₂)}` for `|mᵢ| ≤ T₀` small (`T₀ = 6–8`),
the constraint = the 3-point moment matrix indexed by `{1, e^{iπm x}}` is PSD, plus
the existing Bochner level-1 blocks. Solve at one row center (row 4) at small `N`.
Three decision outcomes: (a) the level-3 dual objective at row 4 **exceeds**
`0.380558` → the 2-point ceiling is genuinely broken, scale up + add Henrion–Rudi
tail; (b) it equals the 2-point value → 3-point adds nothing here, stop; (c)
infeasible/invalid → the torus↔interval boundary is essential (informative; pivot
to the interval encoding). Critically validate by reconstructing `f` and checking
`f ∈ [0,1]` and `(f*f)` sup — the PRO-22 validity gate is mandatory.

**Why it could beat saturation.** The ceiling is a 2-point/degree-2 phenomenon
(PRO-22). A 3-point moment matrix constrains the SDP's chosen `(a_m,b_m)` against
*triple* correlations of the SAME `f`, which the cell-envelope's pairwise Parseval
link leaves free. The autoconvolution sibling (White 2210.16437) is "almost tight"
in `L²` precisely because the `L²` object is fully captured at degree 2; `µ`'s
`min-over-translations`/`L∞` structure is exactly where higher-arity moments are
expected to bite (the external scout's own "sister-problem signal"). And unlike the
retracted Lasserre, convergence-to-`µ` of the de Laat–Vallentin hierarchy is a
*theorem* (the analogue of flag-algebra/theon convergence), so the only question is
the level at which the gain appears, not whether the hierarchy is sound.

**Risk.** (1) The 3-point moment matrix is `O(T²) × O(T²)`; even `T₀ ≈ 8` is a
~64-dim block and the honest interval (non-torus) version with rigorous truncation
may not show a gain until `T₀` is past tractable size — same wall, one level up.
(2) The interval-boundary (Henrion–Rudi is periodic; `µ` lives on `[0,2]` with
boundary) is real work the paper flags as a reformulation. (3) It might simply
*confirm* the 2-point ceiling is sharp in the limit (a meaningful but
bound-neutral negative). (4) Validity (`f ∈ [0,1]`) must be re-derived for the
3-point encoding — PRO-22 shows naive higher-moment SDPs can be invalid.

---

## Approach 2 — Representation-theoretic SYMMETRY REDUCTION of White's SDP (regular *-representation block-diagonalization) to push `T_max` an order of magnitude

**Technique.** Block-diagonalize White's Fourier SDP using the **regular
*-representation / matrix-*-algebra reduction (de Klerk–Pasechnik–Schrijver,
Math. Program. 2007)** and the **Bachoc–Vallentin / Gatermann–Parrilo** symmetry-
in-SOS machinery. The problem carries an explicit symmetry: White proves *all even
cosine coefficients of `M` are nonpositive* and the extremal `f` satisfies the
reflection `f(x) ↔ f(2−x)`, so the relevant group is (at least) the `Z/2` reflection
times the parity/sign structure of the cosine/sine basis. Under this group the
Bochner moment matrices and the cosine/sine cell-envelope constraints are
*invariant kernels*; isotypic decomposition splits each PSD block into a direct sum
of much smaller blocks (the imaginary/odd part decouples from the real/even part),
so an SDP that today costs `bochner_n = 20–30` at `N = 20k–40k` can be solved at
**far larger effective `T_max`/Bochner level for the same flops** — directly
attacking the documented "large-T intractability" barrier that gates every
higher-level scan. This is exactly how Cohn–Elkies-type Fourier SDPs are pushed to
huge bandwidth in extremal geometry (Vallentin, *Symmetry in SDPs*, arXiv:0706.4233;
survey arXiv:2510.06960).

**Citations.**
- de Klerk, Pasechnik, Schrijver, *Reduction of symmetric SDPs using the regular
  *-representation*, Math. Program. 109 (2007) 613–624.
- Bachoc, Vallentin, *New upper bounds for kissing numbers from SDP*, JAMS (2008);
  Vallentin, *Symmetry in semidefinite programs*, arXiv:0706.4233.
- Gatermann, Parrilo, *Symmetry groups, semidefinite programs, and sums of squares*,
  J. Pure Appl. Algebra (2004).
- de Klerk, Sotirov, group-symmetry SDP for QAP, Math. Program. (2010);
  de Klerk, *Numerical block diagonalization of matrix *-algebras*, Math. Program.
  (2011) — the algorithmic block-diagonalizer.

**Why untested here (vs ledger).** The barrier is named in CLAUDE.md and the
saturation memos as *tractability at large `T_max`* — but the project's only use of
symmetry is `lp_research_state/code/symmetric_push.py`, which makes the *ad-hoc,
CONDITIONAL* assumption `f` is **even** (`d_k = 0`), trivially zeroing the imaginary
Bochner block. That (i) only covers rows 5 and 6 (rows with `h>0` contradict
even-`f`, so it is NOT an unconditional bound on `µ`), and (ii) is a manual variable
drop, not a representation-theoretic isotypic decomposition of the full invariant
SDP. The systematic regular-*-representation reduction — which block-diagonalizes
*unconditionally* (it exploits the symmetry of the *constraints*, not an assumption
on the optimizer) and applies to ALL 7 rows — is genuinely new. No
`representation/regular.representation/block.diag/isotypic/irrep` machinery exists in
the code. This is an *enabler* (lets existing rigorous levers run far larger), so it
composes with — does not re-tread — the Bochner-PSD + ellipse base.

**First experiment (cheap).** Pure linear algebra, no new theory. Take the existing
`build_problem(...)` PSD/SOC blocks at a small config (`N=300, bochner_n=12`). Form
the symmetry group `G = ⟨reflection x↦2−x⟩` action on the cosine/sine basis and on
the Bochner index set; compute the **commutant algebra** (matrices commuting with the
`G`-action) numerically via `de Klerk`'s block-diagonalization routine (or a Schur-
decomposition of the averaging projector `P_G = (1/|G|)Σ_g ρ(g)`). Verify (a) each
PSD block splits into the predicted even/odd sub-blocks, and (b) re-solving the
*reduced* SDP reproduces the unreduced dual objective to 10+ digits (the project's
standard cross-check). If it matches, measure the flop/memory ratio; if the ratio is
≥ 3–4×, immediately re-run the binding row 4 at `bochner_n = 60–80` (today infeasible
in 4 GB) and read the new `rigorous_dual_LB`.

**Why it could beat saturation.** F3 (`LEVER_F3_FULL_SATURATION.md`) and the
F3 §5.3 "practical takeaway" explicitly project that pushing `bochner_n ≥ 40–50`
plus tighter cell-envelope yields `+4–6×10⁻⁴`, plausibly `µ ≥ 0.38058` — but flagged
it as *infeasible at tractable scale*. Symmetry reduction is the standard tool that
**makes that scale tractable**: it does not change the relaxation, it changes the
cost of solving it, so the projected gain that was "research, not engineering" becomes
an engineering computation. It is the single highest-leverage move against the stated
barrier because it multiplies the reach of *every* rigorous lever already validated.

**Risk.** (1) The honest worst case: if the only symmetry is `Z/2`, the block split is
~2× — useful but not the order-of-magnitude the barrier needs; the payoff hinges on a
*larger* effective group (e.g. a dihedral/translation structure in the cell index set
or in `Z_{2R}`), which must be checked, not assumed. (2) It is an enabler, so even a
perfect 10× still only buys `µ ≈ 0.38058` (the F3 projection) — it cannot pass the
`C_∞ ≈ 0.380558` 2-point ceiling, only approach it; combining with Approach 1 is what
would transcend the ceiling. (3) Block-diagonalization bookkeeping (svec layout of the
scaled PSD blocks, matching cvxpy's canonicalization) is finicky — same class of cost
as the rigor track's PSD-unpack work; mitigate with an independent re-derivation
cross-check per the repo's `_independent` discipline.

---

## Approach 3 (lower-promise, logged for completeness) — Completely-positive / set-copositive reformulation of the discrete `M(n)` envelope

**Technique.** The *finite-`n`* min-overlap problem (partition `{1..2n}`, minimize
max-shift overlap) is a `{0,1}` quadratic min-max. Standard-quadratic / mixed-binary
quadratic programs admit **exact** completely-positive reformulations
(Burer 2009; Bomze–de Klerk), and the CP cone is approximated by Parrilo/de
Klerk–Pasechnik SOS-based outer cones giving *lower* bounds on the optimum. Apply the
set-copositive outer hierarchy to the discrete `M(n)` to get certified lower bounds
that, unlike the SAT/ILP point values (ruled out as "uninformative"), come with a
*tunable hierarchy* and *dual certificates*, then study the `n→∞` limit.

**Citations.** Burer, *On the copositive representation of binary and continuous
nonconvex quadratic programs*, Math. Program. 120 (2009); Bomze, de Klerk, *Solving
standard quadratic optimization problems via LP and copositive programming*, J. Glob.
Opt. (2002); Dür, *Copositive programming — a survey* (2010); Yıldırım et al.,
set-copositive outer approximations, Math. Program. (2025).

**Why untested here (vs ledger).** The ledger rules out *ILP/SAT for exact `M(n)`*
(point values, no hierarchy) and *Lasserre-2 on the density* (retracted). Copositive
programming on the **discrete** `M(n)` is neither: it is a different cone (CP, not
moment-SOS-on-density) attacking the **combinatorial** object directly with a
convergent outer hierarchy. No copositive/CP code exists beyond the false-positive
grep hits (`symmetric_push.py`, `sweep_T5p.py` matched the substring "positive").

**First experiment.** For small `n` (say `n ≤ 40`), write `M(n)` as a standard
quadratic min-max over the simplex of indicator weights, form Burer's CP
reformulation, solve the **level-0/1 Parrilo outer relaxation** (doubly-nonnegative
+ first SOS layer) for a certified LB, and compare to the known exact `M(n)/n` and to
White's `0.379005`. Decision: if the level-1 LB tracks `µ` and *improves with the
level* at fixed `n`, study the `n→∞` scaling; if it is loose at small `n` (likely),
stop.

**Why it could (weakly) work.** It bounds the *actual discrete* quantity with
duals, sidestepping the continuous relaxation's 2-point ceiling entirely. **Honest
assessment: lowest promise of the three.** CP/copositive hierarchies are notoriously
slow to converge and the `n→∞` limit of a per-`n` CP bound has no obvious uniform
control — this is most likely an *uninformative-at-tractable-`n`* outcome, logged so
the lens is complete and the idea is not silently re-proposed later.

---

## One-line ranking for the orchestrator

1. **Approach 2 (symmetry reduction)** — highest leverage, lowest theory risk; a pure
   linear-algebra enabler that directly attacks the *named* barrier and multiplies
   every existing rigorous lever. Cheapest decisive first experiment.
2. **Approach 1 (de Laat–Vallentin 3-point measure hierarchy + Henrion–Rudi tail)** —
   the only proposal that can *transcend* (not just approach) the 2-point
   `C_∞ ≈ 0.380558` ceiling, with a convergence theorem behind it; higher risk/cost.
3. **Approach 3 (copositive on discrete `M(n)`)** — completeness hedge; honestly
   low promise.

**Combination worth flagging:** 2 ⊕ 1 — use symmetry reduction to make the 3-point
measure moment matrix tractable at a `T₀` large enough to show a gain. That pairing is
the realistic route to a substantive LB advance.
