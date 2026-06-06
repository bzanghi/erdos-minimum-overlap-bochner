# NEW APPROACHES — Probabilistic / Information-Theoretic Lens

**Date:** 2026-06-06
**Lens:** probabilistic & information-theoretic (entropy / second-moment / large-deviation / convolution-information inequalities)
**Status:** PROPOSALS ONLY — not results. Each is vetted against the ruled-out ledger.

## Setup recap (the object every approach below acts on)

White's analytic form: `µ = inf_f sup_t A(t)`, where `f : [0,2] → [0,1]`, `∫f = 1`,
and `A = f ⋆ f̃` is the autocorrelation, `Â(ξ) = |f̂(ξ)|² ≥ 0`.

**The structural fact the convex-optimization lane never exploited as an *analytic* tool:**
`A` is nonnegative, even, integrates to `(∫f)² = 1`, and has a *nonnegative Fourier transform*.
Hence `A` is (Wiener–Khintchine) the **autocovariance of a stationary process** and `A/A(0)`
is a **positive-definite probability density** supported in `[-2,2]`. The whole min-overlap
problem is: *how flat can a positive-definite probability density on a bounded interval be in
sup-norm, given its value structure?* This is a question that entropy, Fisher information, and
the Turán extremal-PD theory speak to directly — and none of those has been used here as a
*lower-bound generator* (LEVER_B tried `cp.entr` only as a **DCP constraint inside the same SDP**,
which is a different and strictly weaker move — see "why_untested" in Approach 1).

---

## Approach 1 — Entropy/uncertainty floor on the autocorrelation, with the Turán extremal-PD program as its analytic dual

**Technique / object to compute.**
Treat `ρ := A/A(0)` as a positive-definite probability density on `[-2,2]`. For a density whose
support has diameter `D` and whose sup is `S = ‖ρ‖_∞`, differential entropy obeys
`h(ρ) ≤ log(1/S)` (uniform-on-the-superlevel bound) **and** `h(ρ) ≥` an explicit floor coming from
*positive-definiteness + bounded support*. The key never-used inequality: for a **positive-definite**
density, the Fourier transform `ρ̂ ≥ 0` is itself (up to scale) a density, so the entropic
uncertainty principle (Beckner–Hirschman, *Ann. Math.* 1975; Bialynicki-Birula–Mycielski 1975)
gives `h(ρ) + h(ρ̂) ≥ log(eπ)`. Because `ρ` is **supported in `[-2,2]`**, `h(ρ) ≤ log 4`, which
*forces* `h(ρ̂) ≥ log(eπ/4)` — a nontrivial lower floor on the *spread of the spectrum* `ρ̂ = |f̂|²/A(0)`.
Feed that spectral-spread floor back through the Turán/Boas–Kac extremal characterization of
positive-definite functions of bounded support (Boas–Kac 1945; Gorbachev arXiv:math/0211291;
Kolountzakis–Révész arXiv:math/0204086; recent dual: arXiv:2510.10172) to lower-bound
`A(0) = ‖f‖_2²` and the achievable `sup_t A(t)`.

Concretely the chain is:
`sup_t A(t) ≥ A(0)·(achievable flatness ratio)`, where `A(0)` is bounded **below** by an
entropy-constrained Turán dual, and the flatness ratio is bounded below by the
spectral-entropy floor (you cannot make `A` flat in real space without spreading `|f̂|²`,
which entropy forbids past a computable point).

**Citation.**
- Beckner, *Inequalities in Fourier analysis*, Ann. Math. 102 (1975) — sharp Hausdorff–Young + entropic UP.
- Hirschman 1957 / Bialynicki-Birula–Mycielski 1975 — entropic uncertainty `h(ρ)+h(ρ̂) ≥ log(eπ)`.
- Turán extremal-PD problem: Boas–Kac 1945; Gorbachev–Tikhonov; Kolountzakis–Révész
  (arXiv:math/0204086); Efimov–Gorbachev–Martyanov dual (arXiv:2510.10172, 2025);
  Springer JFAA 2024 (10.1007/s00041-024-10068-0).

**Why untested here (vs ledger).**
The ledger rules out "Bochner moment-matrix PSD" and a *cancelled* "Beurling–Selberg majorant"
attempt, and LEVER_B records that an entropy term was tried — but **only as `cp.entr(f) ≥ τ`
inside CVXPY**, which DCP rejected (concave ≤ const). That is fundamentally different from what
is proposed here in three ways: (i) the entropy is applied to the **autocorrelation `A`, not to
`f`** — `A` is the positive-definite object, `f` is not; (ii) it is used as an **analytic
one-line inequality producing a numeric floor**, not as a solver constraint, so DCP-encodability
is irrelevant; (iii) it is **coupled to the Turán extremal-PD program**, which is a *different
extremal SDP than White's* (Turán maximizes `∫(PD function)` for fixed support and `f(0)=1`;
White minimizes `sup A`). The Turán/Boas–Kac dual has **never appeared in this repo** — the
inequality scavenger hunt (OUT_OF_BOX_INEQUALITIES.md) scanned Plancherel, Hausdorff–Young,
Bombieri, Selberg, Beckner-as-HY, Wiener, Ingham, Heisenberg, Logan, etc., but **not the
entropic uncertainty principle and not the Turán extremal-PD characterization** — those two are
the missing pieces.

**First experiment (cheap).**
Numerically, take the Together UB minimizer `h*` (`lp_research_state/data/together_f_star.json`),
form `A = h* ⋆ h̃*` by FFT, normalize `ρ = A/A(0)`, and **measure all four quantities**:
`h(ρ)`, `h(ρ̂)`, `‖ρ‖_∞`, `A(0)`. Check empirically how slack the chain `‖ρ‖_∞ ≥ exp(-h(ρ))`
and `h(ρ)+h(ρ̂) ≥ log(eπ)` are at the conjectured optimum. If the entropic UP is within ~10⁻³ of
tight at the optimum, the analytic floor it produces will be near `µ` and the approach is live;
if it is slack by 30%+, the floor saturates well below 0.3803 and the approach is dead. ~20 lines,
no solver. **This single measurement is decisive and costs minutes.**

**Why it could beat SDP saturation.**
White's program is the terminal point of *linear* Fourier-duality reasoning (the scavenger hunt
proved this). Entropy is a **strictly concave / nonlinear** functional — it sees information that
no finite collection of linear moment constraints can, which is exactly why the Bochner-PSD
hierarchy plateaued. The entropic UP is *dimension-free and tight for Gaussians*; min-overlap's
near-optimizer is famously **not** Gaussian-shaped (it's a flat-topped near-box), so the
**entropy deficit relative to the Gaussian is large and computable**, and that deficit is precisely
a quantitative obstruction to flatness — i.e., a lower bound on `sup A`. This is a different
mathematical resource than PSD-ness of moment matrices.

**Risk.**
The entropic UP is tight only for Gaussians, so for the box-like optimizer the *uncertainty*
inequality itself may be loose, even though the *entropy deficit* argument runs the other way.
The real risk is that the two halves (entropy→spectral-spread, spectral-spread→Turán→`sup A`)
each lose a constant, and the composed constant lands below 0.293 (the existing partition-of-unity
floor). The first experiment is designed to detect exactly this slack before any theorem is attempted.

---

## Approach 2 — Fisher-information / log-Sobolev variational surrogate for the L∞ functional (a *convex, convolution-superadditive* relaxation the entropy attempt couldn't reach)

**Technique / object to compute.**
The reason `sup_t A(t)` resists a clean convex-analytic lower bound is that `sup` is not smooth
and not polynomial in `f̂` (this is the documented structural blocker that killed the Lever-H
L²-autoconvolution transfer). **Replace `sup_t A` by a Fisher-information control on the
autocorrelation.** Define `J(A) = ∫ (A')²/A` (Fisher information of the density `A`). Two facts make
this a *new* and *rigorous* lever:
1. **Roughness ⇒ Fisher cost.** A function that is small everywhere but integrates to 1 on a
   bounded interval must be *wiggly*, and wiggliness is exactly large `∫(A')²`. There is already a
   crude version of this in OUT_OF_BOX_INEQUALITIES.md item (the `∫(A−1/8)² ≤ (8/π)²∫(A')²`
   Poincaré line) — but that used **Poincaré (L² roughness)**, which is the weak cousin. Fisher
   information `∫(A')²/A` is **scale-correct** and ties directly to entropy via de Bruijn's identity
   and the **log-Sobolev inequality** (Gross 1975; Carlen, *J. Funct. Anal.* 1991 showing LSI ⇔ Stam).
2. **Convolution superadditivity.** `A = f ⋆ f̃` is itself a self-convolution, and Fisher
   information is *superadditive under convolution* (Stam/Blachman: `1/J(X+Y) ≥ 1/J(X)+1/J(Y)`).
   So `J(A)` is controlled by `J(f)`, giving a **closed inequality chain entirely in terms of `f`**
   that lower-bounds the achievable `sup A` through `‖A‖_∞ ≥ (something)/√J(A)` (a 1-D
   Gagliardo–Nirenberg / Faber–Krahn-type sup bound: `‖g‖_∞ ≥ c·‖g‖_2²/‖g'‖_? ` adapted to densities).

The object to compute: the variational problem `min_f [ Gagliardo–Nirenberg sup-floor of A ]`
subject to `∫f=1, 0≤f≤1, supp f ⊆ [0,2]` — which, unlike the raw `sup`, **is a smooth convex
program** in the appropriate variables (Fisher information `∫(A')²/A` is *jointly convex* in `(A, A')`).

**Citation.**
- Stam, *Information and Control* 2 (1959); Blachman, IEEE IT 1965 — Fisher-info convolution ineq.
- Carlen, *J. Funct. Anal.* 101 (1991) — superadditivity of Fisher info; LSI ⇔ Stam.
- Gross, *Amer. J. Math.* 97 (1975) — logarithmic Sobolev inequality.
- de Bruijn identity / Costa EPI (Costa, IEEE IT 1985) — entropy–Fisher bridge.
- arXiv:1608.05431 (links between LSI and convolution inequalities for entropy & Fisher info).

**Why untested here (vs ledger).**
LEVER_B explicitly flags **"log-Sobolev surrogates"** as a *genuinely-new direction that was
NOT done* because it "cannot be encoded in CVXPY's DCP framework without lifting." That is a
statement about the **old (entropy-as-constraint) framing failing in CVXPY** — it is not a
statement that the **Fisher-information variational floor fails mathematically.** The crucial
distinction: `∫(A')²/A` is the *perspective function* of `(u,v) ↦ u²/v`, which **IS jointly
convex and IS DCP-representable** (`cp.quad_over_lin`), unlike `cp.entr` of `f`. So this approach
is *both* (a) outside the Fourier-PSD class the ledger exhausted, *and* (b) actually encodable —
resolving the exact blocker LEVER_B hit. The Poincaré-roughness line in
OUT_OF_BOX_INEQUALITIES.md is the only neighbor; it is strictly weaker (L² not information-scale,
no convolution superadditivity, gave nothing past 1/8).

**First experiment (cheap).**
Two parts, both minutes:
(i) On the Together `h*`, compute `J(A) = ∫(A')²/A` and the Gagliardo–Nirenberg constant `c` such
that `‖A‖_∞ ≥ c·A(0)²·J(A)^{-1/2}` (calibrate `c` empirically, then prove the sharp `c` later).
See how close the inequality is to tight at the optimizer.
(ii) Build the *small* convex program `min_f c·A(0)²·J(A)^{-1/2}` in CVXPY using `quad_over_lin`
for `∫(A')²/A`, at modest `N` (e.g. N=500), and read the optimal value. If it exceeds 0.293
(partition-of-unity) it is already informative; if it approaches 0.380 it is a contender.

**Why it could beat SDP saturation.**
The Bochner hierarchy controls `f` through *linear* moment data and got stuck because `sup A` has
no polynomial expansion in `f̂`. Fisher information sidesteps `sup` entirely: it bounds flatness
through a **differential (gradient) functional** that the moment hierarchy is blind to, and its
**convolution superadditivity is a genuinely nonlinear structural law** (Stam) with no
finite-moment proxy. This is the most likely of the three to give a *clean closed-form* improvement
because every link (de Bruijn, Stam, LSI, Gagliardo–Nirenberg) has a known sharp constant.

**Risk.**
The sup-norm lower bound `‖A‖_∞ ≥ c·‖A‖_2²/√J(A)` is the load-bearing and least-standard link;
the sharp 1-D constant for *positive-definite* densities may be weaker than the generic
Gagliardo–Nirenberg constant, and if the box-optimizer has *low* Fisher information (flat top →
small `A'` on the plateau, large only at the shoulders) the floor could be soft. Mitigant: the
shoulders of the near-box optimizer are steep, so `J(A)` is plausibly large there — exactly the
regime where the bound bites. Experiment (i) settles it.

---

## Approach 3 — Second-moment / planted-vs-extremal interpolation on the *discrete* `M(n)`, transferring the Overlap-Gap-Property moment machinery (a fundamentally different proof, bypassing the continuous SDP)

**Technique / object to compute.**
Abandon the continuous `f`-functional entirely and attack the discrete `M(n)` by the
**first/second-moment method**, the way the random-number-partitioning literature does. For a
balanced partition `[2n] = A ⊔ B` and shift `k`, let `X_k(A) = |A ∩ (B+k)|` be the overlap.
The min-overlap constant is `µ = lim (1/n) min_A max_k X_k`. The **interpolation move**: introduce
a temperature-`β` Gibbs measure over shifts, `Z_β(A) = Σ_k e^{β X_k(A)}`, so that
`max_k X_k = lim_{β→∞} (1/β) log Z_β`. Then **lower-bound the worst-case `min_A max_k` by the
first-moment / annealed free energy** `(1/β) log E_A[Z_β(A)]` *combined with a concentration
(second-moment) certificate* that `log Z_β` does not fall far below its mean for **any** `A`.
The object to compute is the **annealed free energy curve** `F(β) = (1/n) log E_A[Σ_k e^{β X_k}]`
and its Legendre transform — the **large-deviation rate function `I(x)`** for a single overlap
`X_k/n → x`, which is explicit (it's a `2×2` relative-entropy / hypergeometric rate, because
`X_k` is a sum of indicator overlaps with a hypergeometric/binomial law).

The lower-bound logic: if for **every** partition there are `~n` essentially-independent shifts
`k`, and each `X_k/n` has large-deviation rate `I(x)` of *being below* a target `x*`, then
`P(all shifts below x*) ≤ exp(-n·(effective count)·I(x*))`; a union/second-moment bound over the
`exp(n·H)` partitions shows that for `x* < µ_candidate` **no** partition can hold all shifts below
`x*`, giving `min_A max_k ≥ x*·n`. This is **exactly the Gamarnik–Kızıldağ OGP/first-moment
template** repurposed from "typical" to "extremal."

**Citation.**
- Gamarnik–Kızıldağ, *Algorithmic obstructions in the random number partitioning problem*,
  Ann. Appl. Probab. 33 (2023) (also ISIT 2022, DOI 10.1109/ISIT50566.2022.9834647) — first/second
  moment + OGP for number partitioning.
- Gamarnik, *The overlap gap property: a topological barrier…*, PNAS 118 (2021) — the moment-method
  framework.
- Achlioptas–Naor–Peres style **second-moment method for sharp thresholds** (Nature 2005;
  *Ann. Math.* 2005) — the canonical "second moment closes the first-moment gap" technology.
- Mézard–Parisi–Zecchina interpolation / Guerra–Toninelli interpolation (CMP 2002) — the
  free-energy interpolation that makes annealed↔quenched rigorous.

**Why untested here (vs ledger).**
The ledger's only discrete attack is "ILP/SAT for exact `M(n)` at small `n` (uninformative)" —
that is **exhaustive small-`n` enumeration**, the polar opposite of an **asymptotic
moment-method bound** that produces a *constant* directly. No second-moment, large-deviation, or
interpolation argument appears anywhere in findings.md or the archive (the one second-moment use
in `erdos_lower_density_research.md` is for a **different** problem — Sidon/`B_h` additive-basis
density, not the overlap functional). The OGP/number-partitioning literature has **never been
cited in this repo**. This is a genuinely orthogonal proof route: it never touches White's `f`,
never solves an SDP, and produces `µ ≥ x*` from a *counting* inequality.

**First experiment (cheap).**
Pure numerics, no solver:
(i) Write the single-shift overlap law: for a uniformly random balanced `A` and fixed `k`, `X_k`
is (asymptotically) `Binomial(n, 1/2)`-like with mean `n/2`... **but min-overlap is about the
*minimax*, so compute the joint law of `(X_k)_k` correlations** — first measure
`Cov(X_k, X_{k'})` empirically by sampling random balanced partitions of `[2n]` for `n = 50..500`
and tabulating `max_k X_k/n`. (ii) Compute the explicit rate function `I(x)` for one overlap and
the annealed `F(β)`, then evaluate the first-moment threshold `x*` where
`H(partitions) = effective_count · I(x*)`. **If `x*` from the crude annealed bound already exceeds
~0.30, the second-moment refinement is worth the heavy lifting; if `x* ≈ 1/4` (Erdős' original
trivial bound) the method only re-derives 1955-era results and is dead.** Sampling + rate-function
eval is ~50 lines.

**Why it could beat SDP saturation.**
The SDP framework is *analytic and continuous*; it cannot see the **integrality and counting
entropy** of the discrete partition lattice — there are only `C(2n,n)` partitions, and the
second-moment method exploits that finite count directly. The OGP literature repeatedly shows that
moment methods give **sharp constants** for partition-type problems where continuous relaxations
are loose (the random-NPP statistical-to-algorithmic gap is the canonical example). If min-overlap's
extremal value is governed by a *clustering/overlap-gap* geometry (plausible: near-optimal
partitions are rigid), the moment method could pin `µ` from below at a value the continuous SDP
*provably cannot reach*, since it is a bound on a different (combinatorial) object.

**Risk — stated honestly.**
This is the highest-variance proposal. Two real dangers: (1) **Worst-case vs typical.** OGP/moment
methods natively bound the *typical* (random-instance) value; Erdős min-overlap is a *worst-case
minimum over all partitions*. The first-moment bound on `min_A max_k` is valid (union bound over
partitions is exactly a worst-case statement), but the **second-moment step that usually tightens
it is calibrated for typical instances** and may not transfer, leaving only a weak annealed bound.
(2) **The overlaps `X_k` across shifts `k` are strongly dependent** (they share the same `A`), so
the "effective count" of independent shifts could be `O(1)` rather than `O(n)`, collapsing the
large-deviation gain. The first experiment measures `Cov(X_k,X_{k'})` precisely to quantify this
before committing. If the effective count is small, the method reproduces only the elementary
`µ ≥ 1/4` and should be abandoned.

---

## Cross-cutting note

All three share one philosophical bet: **the min-overlap optimum is information-theoretically
special** (a maximally-flat positive-definite density / a rigid extremal partition), and the
SDP plateaued precisely because *linear moment duality is blind to the nonlinear information
content* (entropy, Fisher information, counting entropy). Approaches 1–2 stay on White's
continuous side but swap the *resource* from PSD-linear to information-nonlinear; Approach 3
leaves the continuous side entirely for the discrete moment method. The three first-experiments
are mutually independent and each costs minutes, so they can be run in parallel as a triage before
any theorem is invested in.
