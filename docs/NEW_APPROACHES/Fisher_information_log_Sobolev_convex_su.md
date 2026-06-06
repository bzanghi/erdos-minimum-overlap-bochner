# Approach assessment — Fisher-information / log-Sobolev convex surrogate for the L∞ functional

**Date:** 2026-06-06
**Assessor:** adversarial-vetting subagent (workflow: "new attacks on µ")
**Source proposal:** RANKED_SHORTLIST_2026-06-06.md #5; `probabilistic_info_theory.md` Approach 2.
**Verdict: WEAK.** The approach rests on a *misidentification of the min-overlap functional*
that, once corrected, voids the load-bearing inequality; and the central feasibility claim
("`∫(A')²/A` is DCP-encodable via `cp.quad_over_lin`, unlike `cp.entr(f)`") is **false as applied
to this problem** — empirically confirmed: CVXPY returns curvature `UNKNOWN`. Two independent,
each-fatal defects. Documented below with the cheap experiments that settle them.

---

## TL;DR (the two kill-shots)

1. **Wrong target functional.** The min-overlap value `µ ≈ 0.3809` is **NOT** the global sup
   `‖A‖_∞` of the autocorrelation. For *every* admissible `h`, Cauchy–Schwarz forces
   `A_h(t) = ⟨h, h(·+t)⟩ ≤ ‖h‖₂² = A_h(0)`, so `‖A_h‖_∞ = A_h(0)` is **always attained at the
   origin**. Numerically on Together's optimizer: `A_h(0) = 0.7746749 = ‖h‖₂²`, while `µ = 0.380871`.
   The actual functional is the **constrained** sup `sup_{|t| ≥ t₀} A_h(t)` (shifts bounded *away*
   from the origin; effectively `t₀ ≈ 2/3` — see numbers below). The proposed sup-floor
   `‖A‖_∞ ≥ c·A(0)²·J(A)^{-1/2}` lower-bounds the **global** sup, which equals `A(0)`. Substituting
   `‖A‖_∞ = A(0)` collapses the inequality to `J(A) ≥ c²·A(0)²` — a Stam/Cramér–Rao bound on the
   autocorrelation's Fisher information, **decoupled from µ entirely**. It says nothing about the
   *restricted* sup that defines the problem.

2. **The DCP claim is false here.** `cp.quad_over_lin(u,v) = u²/v` is jointly convex in `(u,v)` —
   true. But it is a convex function *of the optimization variable* only when `u` and `v` are
   **affine** in that variable (convex-composition theorem). Here `u = A'(t)` and `v = A(t)` are
   both **quadratic** (bilinear) in `f` (`A(t) = ∫f(x)f(x+t)dx`). `quad_over_lin(quadratic,
   quadratic)` is **not** convex in `f`. CVXPY agrees: curvature `UNKNOWN`, `is_dcp() == False`.
   The proposal's framing — "unlike `cp.entr(f)` which DCP correctly rejected, this IS encodable" —
   is exactly the trap it claims to escape, one composition-rule deeper. Moving to `A` as the
   variable does not save it: the objective `A(0)²·J(A)^{-1/2}` is non-convex in `A` too
   (`J^{-1/2}` is convex-*decreasing*, composed with convex `J`, then multiplied by convex `A(0)²` —
   product of convex is not convex; CVXPY: `UNKNOWN`), and the constraint "`A = f⋆f̃` for some
   admissible `f`" is precisely the nonconvex membership that the SDP relaxes via Bochner. Drop it →
   the already-ruled-out invalid direct-sup SDP. Keep it (Bochner) → back inside the saturated cone.

Either defect alone sinks the approach. They are independent.

---

## (1) Feasibility with available tools

| Step claimed | Tool | Reality |
|---|---|---|
| `J(A)=∫(A')²/A` as `cp.quad_over_lin` over `f` | cvxpy | **Not DCP** in `f` (args quadratic). Verified `UNKNOWN`. |
| `min_f c·A(0)²·J(A)^{-1/2}` "smooth convex" | cvxpy | **Not DCP** even with `A` as variable. Verified `UNKNOWN`. |
| Calibrate `c` on Together `h*` | numpy/FFT | Trivial, done below (minutes). |
| Sharp constants (de Bruijn, Stam, LSI, GN) | mpmath/sympy | Available, but the *sup-lower-bound* link they'd need does not exist as a standard inequality (GN/Stam are upper-bound / additivity tools). |

The "small CVXPY solve at N=500" the proposal promises **cannot be written** in DCP form. One could
solve the non-convex surrogate with a generic NLP (scipy `minimize`, IPOPT) and *hope* for a global
optimum, but (a) that forfeits the rigor that is the entire point of a lower bound — a non-convex
local min is not a certified bound — and (b) it would be lower-bounding the wrong quantity (defect 1).
**Feasibility: the literal proposal is not implementable as stated; a non-convex relaxation is
implementable but non-rigorous and aimed at the wrong functional.**

## (2) Prior art

- **Stam (Inf. & Control 1959), Blachman (IEEE-IT 1965):** Fisher-information convolution
  inequality `1/J(X+Y) ≥ 1/J(X)+1/J(Y)`. **Superadditivity / upper-bounds-on-roughness tool.**
- **Carlen (JFA 101, 1991):** superadditivity sharpening; LSI ⇔ Stam.
- **Gross (Amer. J. Math. 97, 1975):** logarithmic Sobolev inequality. **Smoothness/entropy tool.**
- **Costa (IEEE-IT 1985), de Bruijn identity:** entropy–Fisher bridge.
- **Gagliardo–Nirenberg:** in 1-D, `‖g‖_∞ ≤ C‖g‖₂^{1/2}‖g'‖₂^{1/2}` — an **UPPER** bound on the sup.

**Crucial prior-art gap (honest):** every cited tool produces *upper* bounds on a sup or
*super/sub-additivity* of an information functional. **None produces a lower bound on a sup-norm.**
The proposal needs `‖A‖_∞ ≥ (…)/√J` — a *lower* bound on a sup — which is not a standard inequality.
The only honest lower bound on a global sup is mass/width: `‖A‖_∞ ≥ (∫A)/|supp A| = 1/4` (the
Erdős 1955 trivial bound), which has *no* `J`-dependence to raise it. arXiv search for
"Gagliardo-Nirenberg lower bound sup norm", "Stam inequality extremal Fourier sup-norm",
"Fisher information autoconvolution extremal" returned **no relevant hits**. In-repo grep
(per RANKED_SHORTLIST) confirms `fisher/stam/log-sobolev` appear only as name-collisions — never
built. **No piece of this has been tried here, and the literature does not supply the missing link.**

*Adjacent in-repo:* `OUT_OF_BOX_INEQUALITIES.md` has the Poincaré line `∫(A−⅛)² ≤ (8/π)²∫(A')²`
(weaker L²-roughness cousin); it gave nothing past `1/8`. That is the closest neighbor and it died.

## (3) First experiment — RUN (cheap, decisive)

Reconstructed Together's `h*` on `[0,2]` (600-cell step fn, `∫h=1`), refined ×16, FFT autocorrelation.

```
A_h(0) = ‖h‖₂²                                  = 0.7746749   ← the GLOBAL sup ‖A_h‖_∞
sup_t A_h(t)                                    = 0.7746749   at t = 0   (== A(0), as proven)
µ (Together)                                    = 0.380871
sup_{|t|≥0.667} A_h(t)                          = 0.390093    at t ≈ 0.667   ← THIS tracks µ
A_h(t) first crosses 0.3809 near                t ≈ 0.68
J(A) = ∫(A')²/A                                 = 11.07
c needed for ‖A‖_∞ = c·A(0)²·J^{-1/2}            = 4.296   (this just re-expresses A(0)=‖h‖₂²)
```

**Reading:** the global sup the inequality controls (0.7747) is the wrong number by a factor of ~2.
The quantity that *is* µ is a sup over `|t| ≥ t₀`, which a global Fisher-information floor does not
see. The "calibrated `c`" is not a slack measurement of a bound on µ — it is just the identity
`‖A‖_∞ = A(0)` rearranged. **The decisive measurement the proposal asked for, performed, shows the
inequality targets the wrong object.**

*Scaling sub-check (in the proposal's favor, for fairness):* under mass-preserving stretch
`A_λ(t)=(1/λ)A(t/λ)`, `‖A‖_∞~λ⁻¹`, `A(0)²~λ⁻²`, `J~λ⁻²` ⇒ RHS `~λ⁻¹`. The inequality **is**
scale-consistent, and on Gaussians the `A(0)²` form is *exact* with `c=√(2π)=2.5066` (constant
across all scales). So the proposal's instinct about a clean constant is not wrong — it is just a
clean constant for the **wrong functional** (`‖A‖_∞=A(0)`), giving a Stam-type bound `J(A)≥2π‖h‖₂⁴`
that is orthogonal to the min-overlap question.

*DCP sub-check (CVXPY, venv):*
```
A_k = Σ f_i f_{i+1}            (one autocorr lag)   curvature UNKNOWN, is_dcp False
quad_over_lin(affine, affine)                       CONVEX        (the only valid use)
quad_over_lin(quad_in_f, quad_in_f)                 UNKNOWN, is_dcp False   ← the claimed J(A) encoding
A(0)²·J^{-1/2}                                       UNKNOWN, is_dcp False
minimize A(0)²·J^{-1/2}                              is_dcp False
```

## (4) Why it could (in principle) beat the 0.380558 saturation — and why that hope doesn't cash

The *strategic* premise is sound and worth preserving for other approaches:
- **The ceiling `C_∞ ≈ 0.380558` is provably a property of White's cell-envelope + Bochner cone
  specifically** (PRO-6 equality `f_CB=f_C`; PRO-22 shows dropping the cell-envelope makes the bound
  invalid). A bound built from a *different* analytic functional has a *different* ceiling.
- **The Lever-H blocker** (`sup A` has no polynomial expansion in `f̂`) is real, and a *differential*
  functional like `∫(A')²/A` genuinely sees information the moment hierarchy is blind to.
- **Convolution superadditivity** (Stam) is a genuinely nonlinear structural law with no
  finite-moment proxy.

But none of this cashes, because:
- The functional that needs bounding is a **constrained sup** `sup_{|t|≥t₀}A`, not the global sup.
  Fisher information of `A` controls global flatness; it has no native grip on "how large must `A`
  be *away from the origin*." The Stam/GN machinery would have to be re-derived as a *local* sup
  lower bound — a tool that does not exist and is not obviously true (a function can be globally
  smooth yet small on `|t|≥t₀`, e.g. a bump concentrated near 0).
- Even granting such a tool, the convex program to optimize it over admissible `f` is non-convex
  (defect 2); any rigorous convexification re-imposes Bochner ⇒ the same `C_∞`.

So the approach does not transcend the ceiling; it either (a) bounds the wrong quantity, or (b) to
bound the right quantity rigorously, falls back into the saturated cone.

## (5) Risks / why it fails (summary, ranked)

1. **[FATAL, certain] Wrong functional.** `‖A‖_∞ = A(0) = ‖h‖₂² ≠ µ`. The min-overlap sup is
   restricted to `|t|≥t₀`. A global-sup floor is the wrong instrument. *Verified numerically.*
2. **[FATAL, certain] Non-convex program.** `quad_over_lin(quad,quad)` and `A(0)²·J^{-1/2}` are not
   DCP in `f` or in `A`. The "small CVXPY solve" cannot be written. *Verified in CVXPY.*
3. **[FATAL, structural] No sup-lower-bound inequality exists.** GN/Stam/LSI are upper-bound/
   additivity tools. The needed `‖A‖_∞ ≥ c·‖A‖₂²/√J` for a *local* sup is not in the literature and
   is not generically true.
4. **[if forced through anyway] Rigor loss.** A non-convex NLP local minimum is not a certified
   lower bound — defeats the purpose.
5. **[soft-floor, as the proposal itself flagged]** The near-box optimizer has a flat top (small
   `A'`, small `J` contribution there) so even a valid floor would be soft where it matters.

---

## Recommendation

**Do not pursue as a lower-bound generator on µ.** Re-rank below the SDP-adjacent new-cone attacks
(#1 symmetry reduction, #2 three-point/bispectrum lift, #4 NPA) and the analytic-certificate
attacks (#6 sign-uncertainty), all of which act on the *correct* restricted functional and have a
coherent rigor path.

**Salvageable fragment (low priority):** the *entropy-deficit* idea in the sibling Approach 1
(`probabilistic_info_theory.md`) is a different mechanism — it works on the spectral side
`ρ̂=|f̂|²/A(0)` via the Beckner–Hirschman entropic UP, not on a global real-space sup, so it is not
killed by defect 1. It still must clear its own composed-constant slack (its first experiment is the
20-line entropic-UP-slack measurement). If any thread of the information-theoretic lens survives, it
is that one, not the Fisher-information sup-floor. Worth a separate cheap triage, not part of this
approach.

**One genuinely-new sub-idea this assessment surfaced** (offered, not claimed): the correct object
is `sup_{|t|≥t₀} A_h(t)` with `Â_h=|f̂|²≥0`. That is *literally a one-sided / band-restricted
sign-uncertainty problem on a positive-definite function* — which is exactly what shortlist **#6**
(sign-uncertainty / Poisson-summation certificate on `A`) targets. The Fisher detour is subsumed by,
and strictly weaker than, #6. Fold any residual interest there.

---

### Reproduction

```bash
source .venv/bin/activate
cd lp_research_state/code
python3 - <<'PY'
import json, numpy as np
d=json.load(open('../data/together_f_star.json'))['together']
bp=np.array(d['breakpoints']); vals=np.array(d['values']); dx_cell=bp[1]-bp[0]
N=len(vals)*16; x=(np.arange(N)+0.5)*(2/N)
h=vals[np.minimum((x/dx_cell).astype(int),len(vals)-1)]; dx=2/N
A=lambda t:(lambda sh: np.sum(h*np.r_[h[sh:],np.zeros(sh)])*dx if sh>=0 else np.sum(h*np.r_[np.zeros(-sh),h[:N+sh]])*dx)(int(round(t/dx)))
print("A(0)=",A(0)); ts=np.linspace(0,2,801); Av=np.array([A(t) for t in ts])
for t0 in (0,0.667): m=ts>=t0; print(f"sup_t>={t0}",Av[m].max(),"@",ts[m][Av[m].argmax()])
PY
```

Files referenced: `lp_research_state/data/together_f_star.json`,
`lp_research_state/data/together_diagnostic_results.json` (`exact_parseval.f_norm_sq=0.7747`),
`docs/NEW_APPROACHES/harmonic_analysis_lens.md` line 86 (`sup_{|t|≥t₀}` definition),
`docs/archive/PRO6_COMPLEMENTARITY_PROOF.md` (C_∞ ceiling),
`docs/archive/LEVER_B_DISCOVERY.md` (the `cp.entr` DCP rejection this proposal claimed to escape).
