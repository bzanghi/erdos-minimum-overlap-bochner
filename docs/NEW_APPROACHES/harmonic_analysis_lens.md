# NEW APPROACHES — Harmonic-Analysis Lens

**Workflow:** goal-trial-1 novel-approach generation. **Lens:** functional/harmonic analysis.
**Date:** 2026-06-06. **Status:** PROPOSALS ONLY — not results. Vet before claiming.

Setup recap (Together/§5.1 analytic form):
`µ = inf_h sup_t A_h(t)`, `A_h(t) = ∫ h(x) h(x+t) dx = (h * h̃)(t)`, `h:[0,2]→[0,1]`, `∫h=1`.
So `Â_h(ξ) = |ĥ(ξ)|² ≥ 0` (A is a **positive-definite** function), `A_h ≥ 0`, even, `∫A_h = 1`,
supported in `[-2,2]`. Current: µ ≥ 0.380284 (SDP, rigorous), µ ≤ 0.380871 (Together).
SDP-framework ceiling conjectured ≈ 0.380558 (PRO-6).

These three proposals are deliberately OUTSIDE the convex-optimization lane that saturated.
Each names the technique, a grounded citation, why it is untested here (vs the ruled-out ledger),
the concrete first computation, why it could beat the SDP saturation, and the honest risk.

---

## HA-1. Direct L⁴ / autoconvolution-sup inequality WITH the L∞ constraint built in (Martin–O'Bryant + Yu, *strengthened* by `h ≤ 1`)

**Technique.** Martin–O'Bryant (arXiv:0807.5121, *Illinois J. Math.* 2009) prove, for nonneg
`f` supported on an interval of length `I`: `sup_x (f*f)(x) ≥ 0.631 ‖f‖₁²/I`, via an adaptation
of **Yu's number-theoretic averaging technique** (NOT an LP). Crucially their proof uses ONLY
`f ≥ 0`, `f ∈ L¹∩L²` — it does **not** use any pointwise upper bound on `f`. The Erdős problem
is the autocorrelation (reflection) analogue: `A_h = h * h̃`. For our `h` (`I=2`, `‖h‖₁=1`) the
raw Martin–O'Bryant constant gives only `sup A_h ≥ 0.631/2 ≈ 0.3155` — below the current bound —
**precisely because it throws away `h ≤ 1`.** The proposal: re-run their averaging argument with
the extra hypothesis `0 ≤ h ≤ 1` (equivalently `h² ≤ h`, i.e. `‖h‖₂² ≤ ‖h‖₁ = 1`) inserted into
the energy step. Concretely, Yu/Martin–O'Bryant bound `sup A_h` from below by relating it to
`∫A_h² = ∫|ĥ|⁴ = ‖ĥ‖₄⁴` (the **additive energy / L⁴ norm** of ĥ) and to the support length; the
`h ≤ 1` constraint is exactly a sharp *upper* bound on `‖ĥ‖₂² = A_h(0) = ‖h‖₂² ≤ 1` and lets you
trade a flat autocorrelation against concentration. The object to derive is a closed-form
`sup A_h ≥ Φ(‖h‖₂²)` that is monotone increasing as `‖h‖₂²↓` and is then minimized over the
admissible range `‖h‖₂² ∈ [1/2, 1]` (1/2 = the indicator `1_{[0,1]}`).

**Citation.** Martin–O'Bryant, arXiv:0807.5121; Yu, "An analogue of the large sieve" (the
averaging trick they adapt); Cloninger–Steinerberger (sup-autoconv ≥ 1.28 on `(-1/4,1/4)`,
shows the constant is support-and-class dependent — so re-deriving for our class is the right move).
Contrast with White's own *L²* autoconvolution work (arXiv:2210.16437) which solves the
`‖f*f‖₂` problem to 0.0014% — the `sup` (L∞) problem is the open sibling and the one we need.

**Why untested here (vs ledger).** The ruled-out "Barnard–Steinerberger / Madrid–Ramos transfer"
(PRO-32) failed for THREE specific reasons: (i) missing L∞ bound, (ii) t-range direction wrong,
(iii) support mismatch. This proposal is **the opposite move**: it does not transfer a finished
constant; it RE-OPENS the Martin–O'Bryant proof and inserts the L∞ bound as a NEW hypothesis — i.e.
it attacks exactly the gap (missing L∞) that killed the transfer. The OUT_OF_BOX_INEQUALITIES memo
considered Plancherel/Cauchy–Schwarz as *one-liners* and got 1/8; it never re-derived the
Martin–O'Bryant/Yu averaging argument with the box constraint. The "additive-energy/L⁴" lane was
named in scoping but never computed. No SDP, no Bochner, no ellipse — a genuinely different proof
engine.

**First experiment (cheap, <1 day).** (a) Numerically verify the target: on Together's near-optimal
`h*` and on `1_{[0,1]}`, compute `‖h‖₂²`, `‖ĥ‖₄⁴ = ∫A_h²`, and `sup A_h`; tabulate the ratio
`sup A_h · I / (‖h‖₁²)` to see how far above 0.631 the *constrained* optimum sits (expect ≈ 0.76,
since 0.380871·2 ≈ 0.7617) — this is the empirical headroom the L∞ bound buys. (b) Symbolically
(sympy/mpmath) reproduce Martin–O'Bryant's key averaging inequality for a smoothing kernel, then
add the `‖h‖₂² ≤ 1` constraint via Lagrange and see whether the floor rises from 0.3155 toward 0.38.
(c) PSLQ-free sanity: check whether the *known* answer 0.380871 is consistent with a clean
`Φ(‖h‖₂²)` form at the conjectured optimal `‖h‖₂²` (read off Together's h*).

**Why it could beat the SDP saturation.** The SDP ceiling (≈0.380558) is a property of the
*finite Fourier truncation + Bochner-PSD relaxation* — a specific dual cone. An averaging/L⁴
inequality is a **primal analytic** bound living in a different cone (real-space second moments of
A), so it is not subject to the same truncation saturation; if the L∞-strengthened constant lands
above 0.380558 it pierces the conjectured framework ceiling with a one-page proof rather than a
larger SDP. It also could yield a *closed-form* lower bound (the holy grail the PSLQ hunt sought
from the wrong end).

**Risk.** The honest danger: with `I=2` fixed and `‖h‖₁=1`, the L∞ bound may only lift the constant
from 0.3155 to somewhere still below 0.38 (the autoconvolution-sup floor for the *full* class is
genuinely lower than µ because µ also exploits the inf-over-t and the partition structure). If the
strengthened constant tops out < 0.380284 it is a weaker LB and only of expository value. Medium-high
risk on the headline; low risk on producing a clean, citable, framework-independent bound that at
least *reproves* a substantial fraction of µ by new means.

---

## HA-2. Sign / Fourier-concentration uncertainty for the positive-definite autocorrelation (Bourgain–Clozel–Kahane → Gonçalves–Oliveira e Silva–Steinerberger), via a Poisson-summation / modular-form LP on `A_h`

**Technique.** `A_h` is positive-definite (`Â_h = |ĥ|² ≥ 0`), nonneg, even, `∫A_h = 1`, band-in-space
to `[-2,2]`. The **sign-uncertainty** program (Bourgain–Clozel–Kahane 2010; Cohn–Gonçalves 2019;
Gonçalves–Oliveira e Silva–Steinerberger, "New Sign Uncertainty Principles", arXiv:2003.10771,
*Discrete Anal.*) studies exactly functions whose value AND whose transform are sign-constrained,
and bounds the geometry of the sign-change / mass-concentration trade-off. The minimum-overlap
functional is a **two-sided sign-uncertainty instance with the positive-definite (+1) constraint
globally on the transform**: we want `A_h ≥ 0` with `Â_h ≥ 0`, `A_h(0)=‖h‖₂²≤1`, `A_h` supported in
`[-2,2]`, and we minimize `sup_{|t|≥t₀} A_h(t)` — i.e. force the autocorrelation to *spread* its
unit mass flatly out to the support edge rather than spike near 0. The concrete object: build the
**Poisson-summation / Cohn–Elkies-style certificate** specialized to this constraint set — choose an
auxiliary even function `g` with `g ≤ A_h` pointwise outside `[-t₀,t₀]` and `ĝ` controllable, and
optimize `g` to maximize the forced floor. This is the sphere-packing-LP *dual on A directly*, not
on f — and the sign-uncertainty refinements give sharper certificates than the naive Logan first-zero
bound (which the OUT_OF_BOX memo correctly noted collapses to 1/8).

**Citation.** Bourgain–Clozel–Kahane, *Ann. Inst. Fourier* 2010 ("Principe d'Heisenberg et fonctions
positives"); Cohn–Gonçalves, *Invent. Math.* 2019 (the −1 / +12-dim eigenfunction); Gonçalves–
Oliveira e Silva–Steinerberger arXiv:2003.10771 (best 1-D numerics: `A₊(1) < 0.555` conjecturally,
lower bound `1/√(2πe)`); Carneiro–Quesada-Herrera "Generalized sign Fourier uncertainty"
(arXiv:2006.00959) for the weighted/derivative variants that match the autocorrelation's even structure.

**Why untested here (vs ledger).** The ledger rules out **Beurling–Selberg extremal majorants
(cancelled)** — but Beurling–Selberg is the *majorant of an indicator* tool; the SIGN-uncertainty
program is a different and newer object (it bounds the last sign change / concentration of a
sign-constrained pair, not a one-sided majorant). The ledger's "direct sup_t SDP (INVALID — cell
envelope load-bearing)" is the *primal* SDP on the discretized A; this proposal is an *analytic dual
certificate* on the *continuous* A using a globally-valid Poisson-summation identity (no cell
envelope, hence no validity hole). The spectral/translation-operator reformulation (PRO-29) is
naive-Rayleigh and 4× loose; this is not a Rayleigh quotient — it is an LP-certificate on A with the
sign-uncertainty geometry that the Rayleigh bound ignores. Cohn–de Laat–Salmon non-sharpness
machinery (flagged in OUT_OF_BOX_CROSS_DOMAIN) is the *companion* tool for the SAME certificate, never
built here.

**First experiment (cheap, ~1–2 days).** (a) Pose the finite-dimensional certificate LP: discretize
`g` on `[-2,2]` in a cosine basis (degree ~40), enforce `ĝ ≤ 0` on the relevant band and
`g(t) ≤ 0` for `|t| ≤ t₀`, maximize `ĝ(0)/g(0)`-type ratio that lower-bounds `sup A_h`. This is a
SMALL LP (different variables than White's — it is the LP dual *on A*, not on f) — run it for a sweep
of `t₀` (the effective min-shift) and read the certificate floor. (b) Cross-check the floor against
the sign-uncertainty constant `A₊(1)` after the affine change of variables that maps "last sign change
of `A_h − c`" to the GOSS normalization; if the specialization is tight it should reproduce a number in
the 0.30–0.38 band. (c) Use Poisson summation over the lattice `2ℤ` (matching the support) to turn the
certificate into an *exact* inequality, à la Cohn–Elkies, then evaluate at Together's h* to confirm no
contradiction with 0.380871.

**Why it could beat the SDP saturation.** The certificate is constructed on the *autocorrelation A
directly* with its true global positive-definiteness, bypassing the cell-envelope relaxation whose
load-bearing looseness is exactly what caps the current program (the "direct sup_t SDP is invalid"
result shows the cell envelope is *necessary for validity but lossy*). A Poisson-summation certificate
is valid WITHOUT the envelope, so its ceiling is set by certificate quality, not by the
Bochner-truncation saturation. The sign-uncertainty refinements (modular-form-built `g`, as in the
12-D sharp case) are dramatically stronger than the generic LP and have repeatedly hit *sharp*
constants in sibling problems.

**Risk.** Two real risks: (i) the autocorrelation problem may map to a sign-uncertainty instance whose
optimal certificate is itself only ≈ the SDP value (the two LPs could be dual descriptions of the same
bound, giving no gain — this is the Cohn–de Laat–Salmon "is it the same cone?" question, which is
*itself* worth answering rigorously as a no-go theorem). (ii) Constructing a near-sharp `g` analytically
(the modular-form step) is hard and may need a Viazovska-style ansatz that does not exist off-the-shelf
in dim 1 for this constraint. Medium risk; high informational value either way (a sharp certificate
*or* a proved equivalence/no-go).

---

## HA-3. Validated-numerics (interval Newton / Krawczyk) solution of the PRO-23 KKT functional equation — a computer-assisted PROOF of µ, not an SDP relaxation

**Technique.** PRO-23 already derived the Euler–Lagrange / KKT functional equation at the optimum:
on the active shift set `S = argmax_t A_h(t)` there exist weights `γ_t ≥ 0` (`Σγ=1`) and a scalar `κ`
with `Σ_{t∈S} γ_t [h*(x+t)+h*(x−t)] = κ` on the interior `{0<h*<1}`, `≥κ` on `{h*=0}`, `≤κ` on
`{h*=1}`, and `µ = sup_t A_{h*}(t)`. PRO-23 *failed* only because it fed in Together's h* (not a tight
optimum, residual ~10⁻²) and tried to read off γ — it never attempted to **solve the system from
scratch with a validated solver**. The proposal: posit a finite active set `|S|=m` (bang-bang
structure ⇒ h* is piecewise determined by the equation on a finite partition of `[0,2]`), reduce the
functional equation + the m tie-conditions `A_{h*}(t_i)=µ` + the constraints `∫h=1`, `t_i∈S` to a
finite nonlinear system `F(z)=0` in `z=(breakpoints, γ, {t_i}, κ, µ)`, and apply the **interval
Krawczyk operator** `K(F, ž, Z)` (Rump/INTLAB; Moore; Plum's computer-assisted PDE program): if
`K(...) ⊂ int(Z)` then F has a UNIQUE zero in `Z` and the printed interval for the `µ`-component is a
**rigorous two-sided enclosure of µ** — simultaneously a lower AND upper bound, i.e. it could close
the gap.

**Citation.** Krawczyk operator + interval Newton: Moore, *Interval Analysis* (1966); Rump, INTLAB /
"Verification methods" (*Acta Numerica* 2010); Plum, "Computer-assisted proofs for semilinear elliptic
BVPs" — the canonical template for proving existence+enclosure of solutions of nonlinear functional
equations. Directly analogous *and freshly precedented in this exact problem family*: **Rechnitzer,
"The first 128 digits of an autoconvolution inequality" (arXiv:2602.07292, 2026)** — rigorous
ball-arithmetic enclosure of the *L²* autoconvolution constant ν₂² to 128 digits. That paper PROVES a
sibling extremal constant by validated numerics; HA-3 is the autocorrelation-sup analogue applied to the
KKT system. Also: van den Berg–Lessard rigorous-numerics-in-dynamics program (interval Newton in
Fourier/Chebyshev coefficient space) — the standard machinery for validating coefficient-space fixed
points, exactly our setting since h* is naturally a finite step/Fourier object.

**Why untested here (vs ledger).** The ledger's "rigor track (Jansson/VSDP, SDPA-GMP, rational SOS)" is
explicitly about **re-certifying the EXISTING SDP bound 0.380284** — a different goal. HA-3 does NOT
certify the SDP; it solves the *exact optimality system* and produces µ itself with a guaranteed
enclosure, independent of any relaxation. PRO-23 stopped at "Together's h* isn't tight ⇒ blocked"; it
never set up `F(z)=0` and never ran a Krawczyk test — the validated-numerics solve is genuinely
unexecuted. No interval/Krawczyk/INTLAB tooling exists in the repo (`grep` for krawczyk/interval/newton
= empty; python-flint/arb not installed — only mpmath/sympy). This is a new tool and a new object.

**First experiment (cheap, ~2–3 days).** (a) Pin the combinatorial structure: from Together's h* read
the *qualitative* active-set size and the bang-bang breakpoint pattern (how many interior arcs, how many
0/1 plateaus). PRO-23 already has this data (`together_active_set.json`). (b) Take the SMALLEST credible
ansatz — e.g. the symmetric `m`-active-shift reduction at low `m` (start `m=2..4`) — and write `F(z)=0`
in mpmath. (c) Implement a bare interval Krawczyk step by hand in mpmath (interval = midpoint±radius;
`K(ž,Z)=ž−Y·F(ž)+(I−Y·F'(Z))(Z−ž)`, `Y≈F'(ž)⁻¹`); seed `ž` from a float Newton solve started at the SDP
binding point `(h≈0.004, p≈0.3875)`. (d) Check the inclusion `K⊂int(Z)`; even a *non-rigorous* float
Newton that converges to a self-consistent `(h*,γ,κ,µ)` with KKT residual ≪10⁻⁹ (vs PRO-23's 10⁻²) is
an immediate win — it would be the first *tight* KKT point and would pin µ far better than the open gap.

**Why it could beat the SDP saturation.** It sidesteps the framework entirely. The SDP ceiling is a
relaxation artifact; the KKT system is the *exact* condition for the true optimum, so a validated
enclosure of its unique root is a rigorous enclosure of the *true* µ — it can in principle deliver
`µ ∈ [a,b]` with `b−a` at machine-or-better precision, closing the 5.87×10⁻⁴ gap outright rather than
nibbling the LB. Rechnitzer's 128-digit precedent on the sibling L² problem is direct evidence the
approach is feasible at extreme precision for this family.

**Risk.** The biggest risk is the **combinatorial structure of the optimizer**: if the true active set
is large/degenerate (PRO-23 saw 437 near-ties at 10⁻⁹ on the *discrete* surrogate) or genuinely
infinite/Cantor-like, the finite `F(z)=0` reduction is wrong and Krawczyk has no isolated root to
enclose — the inclusion test simply fails to certify (honest, non-misleading failure). Also the
Jacobian `F'` can be ill-conditioned near the degenerate optimum, inflating interval radii past the
inclusion threshold. Medium-high risk on a *full* gap-closing enclosure; LOWER risk on the intermediate
deliverable (a tight float-Newton KKT point that sharply localizes µ and unblocks PRO-23 Step 4). The
right-sized first bet is the small-`m` symmetric reduction, which is cheap to falsify.

---

## Cross-cutting notes

- HA-1 and HA-2 are LOWER-BOUND / different-proof attacks that live OUTSIDE the Bochner-truncation cone
  (real-space L⁴ energy; Poisson-summation certificate on A) — either could in principle exceed the
  conjectured 0.380558 framework ceiling because neither is subject to its saturation mechanism.
- HA-3 is the only one of the three that could close the gap from BOTH sides at once (an enclosure),
  and it is the highest-variance / highest-payoff bet; its small-`m` first experiment is cheap.
- All three are non-redundant with the ledger and with each other: distinct objects (constrained
  averaging constant; sign-uncertainty LP certificate; validated KKT root) and distinct tools
  (Yu/Martin–O'Bryant; Cohn–Elkies/sign-uncertainty; Krawczyk/INTLAB).
- Shared cheap precursor for HA-1 and HA-3: a clean high-precision computation of `‖h*‖₂² = A_{h*}(0)`
  and the active-set geometry of Together's h* (mpmath) — feeds the L∞-strengthened constant (HA-1)
  and the KKT ansatz (HA-3).
