# Adversarial assessment — NPA / Helton–McCullough noncommutative-moment certificate for sup_t ⟨h, T_t h⟩

**Date:** 2026-06-06
**Status:** Adversarial deep-dive of shortlist candidate **#4** (also Approach B of
`CROSS_DOMAIN_WILDCARD_2026-06-06.md`). **This is a VET, not a claim.**
**Verdict: WEAK** (demoted from its #4 shortlist slot). The proposal's two
load-bearing novelty claims — (i) it certifies a *new* object whose relaxation gap is
"unrelated to White's C_∞", and (ii) word-length truncation is "tail-free, dodging the
rigor trap" — **both fail under scrutiny**, verified computationally below. The
genuinely-true positives (NPA has never been pointed at this problem class; the
operator-localizer framing of the L∞ bound is elegant) do not rescue it: the relaxation
**collapses to the already-proposed multi-point-moment lane (shortlist #2)** via a more
exotic and heavier encoding, and the continuum→Z_m step is in the **wrong direction** for
a valid lower bound on µ.

---

## The proposal (recap)

Certify `µ = 1 − sup_h inf_t ⟨h, T_t h⟩` (the PRO-29 duality, `docs/archive/PRO29_SPECTRAL.md`)
as a **noncommutative polynomial eigenvalue-optimization** problem. Variables: the
one-parameter unitary shift group `{T_t}` (`T_s T_t = T_{s+t}`, `T_t* = T_{−t}`) and
multiplication-by-`h`. Build the NPA moment matrix `Γ` over words in `{T_{t_i}, h}` up to
length 2(–4); impose `Γ ⪰ 0`, unitarity localizers, and `0 ⪯ h ⪯ 1`, `⟨h⟩=1` as **operator
localizing inequalities** `[h]Γ ⪰ 0`, `[1−h]Γ ⪰ 0`. Truncate by **word length** (claimed
finite-dim/tail-free). First experiment: discretize shifts to cyclic `Z_m` (m=8–16),
solve level-2 NPA (~80×80) in cvxpy.
Refs: Navascués–Pironio–Acín NJP (2008); Pironio–Navascués–Acín SIAM J. Optim. (2010);
Helton–McCullough Positivstellensatz; arXiv:2402.02126 (2024); arXiv:2510.08427 (2025).

---

## (1) FEASIBILITY with available tools — MODERATE, no blocker

- `cvxpy 1.8.2` + **MOSEK** (licensed, confirmed `installed_solvers()`) + CLARABEL/SCS.
  An ~80×80 PSD with a handful of localizing blocks is trivial for any of these.
- `ncpol2sdpa` is **NOT installed** (`importlib.util.find_spec` → None). It is pure-Python
  on top of cvxpy/MOSEK/SDPA, so `pip install ncpol2sdpa` would work; or one can hand-roll
  the moment matrix (≈ the toy below). `TSSOS` is Julia-only (not in this stack).
- `sympy`, `mpmath`, `numpy`, `scipy` all present → symbolic word algebra + GMP-precision
  re-solve of any final small SDP are available. SDPA-GMP at `lp_research_state/bin/` could
  certify a final certificate at GMP precision (the proposal's "rational-roundable" upside is
  real *if* a useful number ever appears).

**Effort to a first honest number: LOW–MEDIUM (~3–5 days),** dominated by getting the
crossed-product moment bookkeeping and the continuum→Z_m argument right, not by solver time.
But see (3)/(5): the cheap first number is very likely *uninformative or invalid-direction*.

---

## (2) PRIOR ART

**The technique itself (NPA / Helton–McCullough nc-moment / SOHS hierarchy):** mature and
widely used, but **exclusively in quantum information / operator-algebra eigenvalue problems**:
- Navascués–Pironio–Acín, *A convergent hierarchy of SDPs characterizing the set of quantum
  correlations*, NJP 10 (2008) 073013 — the origin; Bell-inequality / quantum-correlation sets.
- Pironio–Navascués–Acín, *Convergent relaxations of polynomial optimization with noncommuting
  variables*, SIAM J. Optim. 20 (2010) 2157 — eigenvalue/trace optimization, the general engine.
- arXiv:2510.08427 (2025) "spectral-gap certificates for **qubit Hamiltonians**";
  arXiv:2402.02126 (2024) "upper-bound hierarchies for nc-poly optimization" — both *eigenvalue*
  bounds, but for spin/qubit Hamiltonians, not group-action extremal constants.
- NC-symmetry reduction (arXiv:2112.10803), term-sparsity (arXiv:2010.06956),
  state/trace/moment polynomials (arXiv:2412.12342) — all quantum-info-flavoured.

**Web/literature search (June 2026) found NO application of NPA to an additive-combinatorics
extremal constant of the overlap/autocorrelation type.** So pointing NPA at `µ` is *genuinely
novel territory for the technique* — this is the proposal's one robust merit. **But novelty of
venue ≠ a new cone:** the literature also offers no precedent that NPA *helps* on a
group-action correlation extremal problem, and the closest cited successes are qubit
Hamiltonians where the algebra is intrinsically noncommutative in a way ours is not (below).

**Tried here?** No NPA/nc code exists in-repo (grep over `lp_research_state/code/` for
`npa|noncommut|helton|mccullough` → only name-collisions in `build_then_solve.py`,
`_jansson_verify.py`; confirmed in `RANKED_SHORTLIST_2026-06-06.md` "all 19 techniques are
genuinely unbuilt"). The *adjacent* lever **is** tried-and-documented: `PRO29_SPECTRAL.md`
(the naive Rayleigh quotient that was 4× loose) and PRO-22 `LEVER_SUPT_DIRECT.md` (the direct
sup_t SDP that gave an INVALID bound because it dropped the cell-envelope/validity link).
This proposal is positioned as the fix to both — see (5) for why the fix does not land.

---

## (3) FIRST EXPERIMENT — run as cheap structural/analytic checks (DONE, no heavy solves)

Rather than fire a 5-day Z_m SDP, I ran the *decisive* cheap structural tests, because the
proposal's value rests entirely on two algebraic claims that can be settled analytically.

### 3a. What object does NPA-on-{T,h} actually certify? — VERIFIED, and it COLLAPSES

The relevant von Neumann algebra is generated by (a) the abelian multiplication algebra
`L^∞(R/2Z)` (mult-by-`h` operators **commute among themselves** — `h` is a *function*, not a
generic operator), and (b) the shift group acting on it by automorphisms. This is the
**crossed product** `L^∞(X) ⋊ Z` (resp. `⋊ R`). A moment of a word against the cyclic
vector `ψ ≡ 1` is

```
⟨ψ| h T_{t1} h T_{t2} h ⋯ |ψ⟩  =  ∫_X h(x) h(x+t1) h(x+t1+t2) ⋯ dx
```

i.e. a **multi-point correlation of h**: a length-2k alternating word `h T h T … h` gives a
**k-point correlation**. I verified this on 2000 random words in `{H, T_j}` on `Z_10`: every
word-moment `⟨ψ|w|ψ⟩` equals the corresponding product-of-translates integral to `<1e-9`.

**Consequence (the kill):** NPA level-L on `{T,h}` constrains *exactly the ≤L-point
correlations of h*. That is the **same object class** as shortlist **#2** (the
de Laat–Vallentin / triple-correlation `R₃(s,t)` / bispectrum `B(m₁,m₂)` lift). The
"noncommutativity" of `T` and mult-`h` is **real but is precisely the mechanism that
generates multi-point correlations** — it is not an independent new certified object. So the
claim "the relaxation gap is unrelated to White's C_∞" is only half-true: it is a
higher-arity-moment relaxation (which *can* exceed the 2-point ceiling, like #2), but it is
**not a different lever from #2** — it is the same lever in heavier clothing. The genuine
3-point cone of #2 is the *cleaner, lighter* way to reach the identical constraints.

Also confirmed: `T_k ψ = ψ` (the constant function is shift-invariant), so the **shift
generators alone do nothing at level 1**, and the objective `⟨h, T_k h⟩` first appears at a
**level-2 word** `h T_k h`. The minimal informative relaxation is level 2 — matching the
proposal, but underscoring that there is no "free" structure from the group beyond what the
correlations already carry.

### 3b. Is the "tail-free" claim valid? — NO, the direction is BACKWARDS

Word-length truncation is finite-dimensional/tail-free **for a FIXED finite `Z_m`** — true.
But the bound is then a statement about `Z_m`, while `µ` lives on the continuum `R/2Z`. The
continuum→`Z_m` reduction is the load-bearing step, and it is in the **invalid direction**:

```
inf over the COARSER Z_m shift set  ≥  inf over the continuum shift set
  ⇒ Z_m  OVER-estimates  µ_dual = sup_h inf_t ⟨h,T_t h⟩
  ⇒ Z_m  UNDER-estimates µ = 1 − µ_dual.
```

A coarse `Z_m` has only `m−1` shifts to minimize over and generically **misses the binding
continuum shift** `t*` (Together's optimizer binds at shift ≈ 33/n; a small-`m` grid skips it),
so `inf_t` over the grid is *too large* → `µ_dual` too large → `µ` too **small**. A naive `Z_m`
NPA number is therefore **not a valid lower bound on µ** — it is the *same failure shape as
PRO-22*, where dropping the validity link produced an invalid-direction bound. Recovering rigor
requires exactly the **Naimark/quadrature/character-tail argument** the proposal lists as
"risk (1)". So the headline "tail-free, structurally dodges the rigor trap" is **false**: the
truncation simply **relocates** the rigor trap from a Fourier tail to a shift-discretization
tail. (The proposal's hedge — "controls a modulus-1 character, cleaner than an f² tail" — may
make the argument *doable*, but it is still a tail bound that must be supplied; nothing is
dodged.) This is the project's single most-repeated overclaim pattern
(`memory/project_tail_bound_rigor_trap.md`), and the proposal walks back into it.

### 3c. Quantitative toy (illustrative, NOT the real geometry)

A faithful tiny relaxation on `Z_6` was built (full crossed-product moment mechanics:
`T_k ψ = ψ`, objective at level 2, 2-pt-correlation Bochner block). Random search gives
`sup_h min_{k≥1} R₂[k] ≈ 0.90` on `Z_6`; the 2-point Bochner SDP relaxes to `1.0`. Both are
**far from the continuum** `µ_dual ≈ 0.619`, because at `m=6` near-constant `h` keeps every
overlap large and the continuum interval geometry is absent. This toy is *not* evidence either
way on tightness — it only confirms (3a)/(3b) mechanics and demonstrates that **small `m` does
not see the continuum constant at all**, so the cheap first number is uninformative until `m`
is large enough to resolve `t*` — at which point the encoding is no longer "trivial 80×80".

---

## (4) WHY IT COULD (in principle) BEAT THE ~0.380558 SATURATION

To be fair, the honest upside, stated precisely:
- It is a **higher-arity-moment relaxation** (3a), and PRO-6/PRO-22 establish that the
  `C_∞ ≈ 0.380558` ceiling is a **2-point / degree-2 phenomenon** of White's cell-envelope +
  Bochner cone. A ≥3-point relaxation lives in a strictly larger cone with a *different*
  ceiling — so it *can*, in principle, exceed `C_∞`. (This is exactly the shortlist-#2
  argument; the NPA route inherits it.)
- The **L∞ bound enters as an operator localizer** `[h]Γ, [1−h]Γ ⪰ 0` — a clean way to carry
  `0 ≤ h ≤ 1` *inside* the moment matrix, which is the constraint whose omission made
  PRO-29 4× loose and killed the B-S/M-R transfers. This framing is genuinely nice and is the
  proposal's best idea.
- NPA is exact at finite level for some highly-symmetric eigenvalue problems, and the shift
  group is a large symmetry — *if* (big if) finite-level sharpness held, a clean
  GMP-rational certificate would be a fundamentally different, rigorous proof route.

**But:** every one of these upsides is *also* available — more cheaply and without the
crossed-product overhead — to shortlist **#2** (direct 3-point/bispectrum block, with the L∞
box carried by the same Bochner moment-matrix `M ⪰ 0` machinery already in the repo). The NPA
route adds no reach over #2; it adds encoding complexity and a discretization tail.

---

## (5) RISKS / why it will most likely fail

1. **It is not a new lever — it reduces to shortlist #2 (3a).** The whole novelty pitch
   ("noncommutative ⇒ new certified object, gap unrelated to C_∞") collapses: the certified
   object is the ≤L-point correlations of `h`, identical to #2's cone class. So the *correct*
   experiment is just to run #2 (lighter, repo-native validity gate). If #2 shows no gain
   (its own #1 risk: the triple cut may be implied by 2-point + cell-envelope for this pairwise
   objective), NPA inherits that negative; if #2 shows a gain, NPA is a strictly worse way to
   capture it. **NPA is dominated by #2 either way.**
2. **The "tail-free" rigor claim is backwards (3b).** A naive `Z_m` NPA number
   **under-estimates µ** (invalid direction), reproducing PRO-22's failure shape. A valid LB
   requires a Naimark/quadrature/character-tail bound — the rigor trap is relocated, not dodged.
3. **Validity link (the cell-envelope lesson) is unaddressed.** PRO-22 proved that any new
   attack must *either* keep White's `M↔f` Parseval link *or* supply an equivalent validity
   certificate. The operator localizers `0⪯h⪯1` constrain `h`'s correlations but do **not** by
   themselves reconstruct the `M(t) = (f⋆f)(t)` ↔ `|f̂|²` link at the SDP optimum; the same
   `f ∉ [0,1]` reconstruction failure that voided PRO-22 must be re-checked for the NPA optimum
   (the localizers help but are not obviously sufficient at level 2 — "unknown a priori",
   proposal's own risk (2)).
4. **Continuum shift group `t ∈ R/2Z` vs `Z_m`** (3b/3c): small `m` literally does not resolve
   the binding shift; large enough `m` to resolve `t*` blows up the moment matrix and the
   "trivial 80×80" tractability claim. The genuinely hard, useful regime is not cheap.
5. **Tooling:** `ncpol2sdpa` absent (installable) / `TSSOS` Julia-only; the encoding is
   more exotic than the repo's cvxpy-native SDPs, so even the *negative* result costs more
   engineering than re-confirming #2's negative.

---

## Recommendation

**Do NOT pursue NPA as a standalone attack.** It is dominated by shortlist **#2** (3-point /
bispectrum lift): same certified cone (proven in 3a), lighter encoding, repo-native validity
gate, and no extra shift-discretization tail. The one transferable idea — **carrying `0≤h≤1`
as an operator/moment-matrix localizer** — should simply be **folded into #2** (it already is,
via the Bochner `M ⪰ 0` block). Run **#2 first**; only if #2 reveals a real, active 3-point
gain *and* one wants a specifically noncommutative re-derivation for a different proof aesthetic
would the NPA framing be worth revisiting — and even then the continuum→Z_m Naimark argument
(risk 2) is mandatory non-optional work, not a free lunch.

**Net:** genuinely-novel *venue* for NPA, elegant L∞-localizer framing, but the central
novelty/rigor claims do not survive (verified 3a/3b), and it is dominated by an already-ranked,
cheaper proposal. **Verdict: weak.**

---

## Artifacts / verification

- Crossed-product = multi-point-correlation identity: verified on 2000 random words in
  `{H, T_j}` on `Z_10` (all word-moments matched product-of-translates to `<1e-9`).
- `T_k ψ = ψ`, objective first appears at level-2 word `h T_k h`: confirmed numerically.
- Continuum→`Z_m` direction analysis: `inf over coarser grid ≥ inf over continuum` ⇒ `Z_m`
  over-estimates `µ_dual` ⇒ under-estimates `µ` (invalid LB direction).
- Tooling: `cvxpy 1.8.2` + MOSEK/CLARABEL/SCS present; `ncpol2sdpa` absent (pip-installable);
  `sympy`/`mpmath`/SDPA-GMP available for a hypothetical final certificate.
- (Toy SDPs were ad-hoc REPL scripts, not saved as repo modules per the `_`-throwaway convention.)
