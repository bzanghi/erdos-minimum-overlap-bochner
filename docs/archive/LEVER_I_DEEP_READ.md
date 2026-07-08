# Lever I — Deep Read: Discrete-Reduction Dual Bounds (arXiv:2206.09876) vs White's Bochner-PSD Hierarchy

**Date:** 2026-05-10
**Paper:** Rupert Li, *Dual Linear Programming Bounds for Sphere Packing via Discrete Reductions*, **Adv. Math.** 460 (2024). [arXiv:2206.09876v3, 6 Nov 2022].
**Note on attribution:** The task brief (and `OUT_OF_BOX_CROSS_DOMAIN.md` §2a) labels this paper "Cohn-de Laat-Salmon (CdLS)". This is incorrect. arXiv:2206.09876 is **Rupert Li, sole author**, with Henry Cohn as mentor (UMN Duluth REU 2022). The actual Cohn-de Laat-Salmon paper is arXiv:2206.15373 (*Three-point bounds for sphere packing*), which produces better upper bounds via SDP, not a non-sharpness theorem. The deep read below is of Li, since that is the paper whose technique was requested for evaluation.

---

## 1. TL;DR

**Verdict: DOES NOT APPLY** (to the goal of a finite ceiling on White's Bochner-PSD hierarchy).

Li's discrete-reduction technique is **a method for *constructing dual feasible points*** in the Cohn-Elkies LP — i.e., for producing **numerical upper bounds on the LP's optimal value**. It is *not* a non-sharpness theorem in the sense the task is hoping for. Li uses these dual points to show that for sphere packing in d ∈ {3,4,5}, the Cohn-Elkies LP's optimal value strictly exceeds the conjectured packing density — but only because **a tighter, independent upper bound on the packing density is already known** (e.g. Hales for d=3, or Cohn-de Laat-Salmon's three-point SDP for d=4,5). The discrete reduction does **not** itself bound the LP from above without that external input. Worse, Li's Theorem 3.2 only **lower-bounds** the (continuous) Cohn-Elkies bound by a (discrete) LP bound — so the discrete LP is only useful in the direction: discrete-dual feasible point ≤ discrete LP ≤ continuous Cohn-Elkies LP, and we conclude continuous ≥ discrete-dual-value. To prove "continuous LP ≤ C*" one would need the *opposite* inequality, which Li does not have. **There is no machinery in Li that converts to "Bochner-PSD hierarchy ≤ C*" for the Erdős program.**

A genuine non-sharpness theorem for sphere packing in d=12,16 exists (Cohn-Triantafillou, ref [10] in Li, using modular forms to construct *primal* feasible dual auxiliary functions over R^d directly), but that machinery is dimension-specific (requires modular forms of weight d/2 for even d, half-integral weight for odd d) and has no analog in White's program.

---

## 2. Li's Technique Summary (≈300 words)

The Cohn-Elkies LP (Problem 2.2) is:
> minimize (r/2)^d f(0) over radial smooth f: R^d → R with f̂(0) ≥ 1, f(x) ≤ 0 for |x| ≥ r, f̂(y) ≥ 0 for all y.

This is **infinite-dimensional**; its dual is over measures on radii, also infinite. Li introduces the **discrete Cohn-Elkies LP** (Problem 3.1) over functions g_m: Z_m^d → R with analogous constraints:
> minimize (r/2)^d g_m(0) s.t. ĝ_m(0) ≥ m^{-d/2}, g_m(x) ≤ 0 for |x| ≥ r in Z_m^d, ĝ_m(y) ≥ 0.

**Theorem 3.2** is the load-bearing reduction: given any feasible f for the continuous LP, *restrict* f to Z^d (call it g), then *periodize* g to Z_m^d (call it g_m). Provided m ≥ 2r, g_m is feasible for Problem 3.1 with the same objective bound (in fact ≤). This implies **Corollary 3.3**: continuous-LP-optimum ≥ discrete-LP-optimum on Z_m^d. The discrete LP is finite-dimensional, so its standard LP dual (Problem 4.1) has no gap, and any feasible point of the discrete dual yields a *numerical lower bound* on the continuous LP. Li then constructs dual feasible points by computer search over G_d-orbits (G_d = Z_2^d ⋊ S_d acts on Z_m^d), rounding to rationals, and certifying with interval arithmetic; results in §§6-9.

**Crucial logical structure:** Li's chain is

    discrete dual feasible λ ≤ discrete LP value ≤ continuous Cohn-Elkies LP value (= continuous primal).

He never bounds the continuous LP value from *above* by this method (he can't — the discrete LP is a *lower bound* on the continuous one). The "Cohn-Elkies bound is not sharp in d=3,4,5" conclusion uses an external upper bound on δ_d: Hales' Kepler proof for d=3, and the three-point SDP of Cohn-de Laat-Salmon for d=4,5. Li shows discrete-LP value > δ_d, which forces continuous-LP value > δ_d, which means continuous-LP ≠ δ_d.

---

## 3. Structural Matching to White's Program

| Required structure                  | Cohn-Elkies LP             | White's SDP                                       | Match? |
|-------------------------------------|----------------------------|---------------------------------------------------|--------|
| Translation-invariant on R^d        | Yes, R^d → Z^d → Z_m^d     | No — problem is on [0,2], not a group             | **NO** |
| LP (no PSD, only sign/positivity)   | Yes — pure LP              | **SDP with Bochner-PSD matrices**                 | **NO** |
| Self-dual Fourier domain            | R^d (Pontryagin self-dual) | Trigonometric Fourier on [0,2]; OK                | OK     |
| Poisson summation gives reduction   | Yes — core of Thm 3.2      | Has Poisson-summation analogs, but…              | partial|
| Continuous → discrete inclusion     | f|_{Z^d} → g_m              | What does "restrict to Z_m" even mean for f∈[0,1]?| ?     |
| Direction of inequality (LB on LP)  | Yes — gives LB on LP value | We need **upper bound on LP value**               | **NO** |

The two showstoppers (in priority order):

**(a) Wrong direction of the inequality.** Li's Corollary 3.3 says continuous-LP ≥ discrete-LP. For the Erdős program we want continuous-LP ≤ C*. Li's technique, even if it could be ported, would give us a *lower bound on the LP's optimal value*, which is **what we already have via dual extraction in `dual_extractor.py`**. We are not after a better lower bound on the LP value; we are after a proof that the LP value is bounded *above* by some C* < µ. Li doesn't address that direction at all. (Cohn-Triantafillou via modular forms *does* — they construct a primal-feasible f over R^d directly with high f(0), which is an *upper bound* on the LP value. That technique is fundamentally different.)

**(b) Bochner-PSD constraints are not preserved under restriction-and-periodize.** Li's reduction (Thm 3.2) crucially uses that `g_m(x) ≤ 0 for |x| ≥ r` is preserved by periodization (sum over m·Z^d only adds non-positive terms when m ≥ 2r). The Bochner-PSD constraints in White's program are *Hermitian Toeplitz matrices in the Fourier coefficients c_k, d_k*. Periodization on the physical-side variable corresponds to *sampling* on the Fourier side, which is a linear operation — fine for the LP constraints (`f̂(y) ≥ 0`), but **for the Bochner moment matrix `M_n(f) = [f̂(j-k)]_{j,k=0..n} ⪰ 0`, the sampled version is a different PSD constraint, generally neither implied by nor implying the original**. Li's restriction map is constructed precisely for the LP setting where the only "PSD" content is sign of Fourier coefficients (a 1×1 Toeplitz, trivially compatible with sampling). Once you have a level-n Toeplitz PSD constraint, the discrete reduction breaks.

(There is an extension: SDP discrete reductions in coding theory exist via Bachoc-Vallentin / Schrijver-style equivariant decompositions, but those *enable computation*, not non-sharpness proofs.)

**(c) Problem-structural mismatch.** Cohn-Elkies's domain (R^d) is a group, with a clean Fourier inversion, Poisson summation, and a known near-optimal *configuration* (root lattice). The Erdős min-overlap problem is on [0,2] — not a group — with an objective `sup_k ∫ h(x) h(x+k) dx` that is itself a `sup` over an auxiliary parameter; White's program is a *Lagrangian relaxation* of this sup, not a direct LP. The "ellipse-extension argument" (§5.1 of White) is needed to convert any per-(h,p,q)-point LP value into a global bound. None of Li's machinery survives that conversion.

What White's program **does** share with Cohn-Elkies:
- Fourier-analytic LP/SDP structure ✓
- A near-optimal primal construction (Together's h*) ✓
- Convergence as truncation parameters T, R → ∞ ✓

But the analog of "f(x) ≤ 0 for |x| ≥ r" (the pointwise negativity constraint that drives Li's reduction) is **not present** in White's program — White's cell constraints are bidirectional bounds on `w_j, v_j ∈ [0, Ω]` and the global moment constraints, not sign constraints on a free function.

---

## 4. Verdict — DOES NOT APPLY

The task asks specifically: "Could CdLS prove a finite ceiling on White's Bochner-PSD hierarchy?" The answer is **no**, for three independent reasons:

1. **Wrong direction of the inequality.** Li's technique produces lower bounds on the LP optimal value, not upper bounds. We need the upper bound to prove a "saturation ceiling".
2. **PSD constraints don't survive the discrete reduction.** Li's restriction-and-periodize map only preserves LP-type sign constraints, not Toeplitz PSD constraints. Bochner-PSD on level n=20 (our setting) is not preserved.
3. **Cohn-Elkies's "f ≤ 0 for |x| ≥ r" has no analog in White's program.** The reduction critically uses this. White uses cell-bound LP constraints + Fourier moment constraints, not a free function with pointwise sign constraints.

Even if a "discrete-reduction-style" analysis were attempted on a per-constraint basis, the **fundamental epistemic blocker** is point (1): Li's method, by construction, gives the same type of bound as our `dual_extractor.py` already gives — a rigorous LB on the SDP value. We do not lack rigorous LB on the SDP value; we have it. The diagnostic finding is that the SDP value itself has plateaued empirically. To prove a *finite ceiling* on what the hierarchy can prove, we'd need to construct a **primal-feasible f** for the limiting (continuous, N,T,R → ∞) program with f(0)/f̂(0) **close to** the saturation value — that is the Cohn-Triantafillou / modular-form direction, not the discrete-reduction direction.

For completeness: the closest analog of "non-sharpness theorem" in our context would be:

> **Open problem.** Construct, by some analytic technique, a feasible point of the *limit* of White's SDP (Bochner-PSD level → ∞, cell refinement → ∞) with objective value ≥ µ_best_upper. This would prove the hierarchy cannot close the gap.

This is a hard direct construction problem, structurally close to what Cohn-Triantafillou does for sphere packing, and **completely different from Li's discrete reduction**.

---

## 5. Recommended Action

**Do not pursue Li-style discrete reduction.** It cannot, even with substantial adaptation, produce the saturation-ceiling theorem the project is hoping for.

The two productive directions, in priority order:

**5a. Cohn-Triantafillou-style primal construction (high value, high effort).** Cohn-Triantafillou (ref [10] in Li, *Math. Comp.* 91 (2021)) constructs *primal-feasible auxiliary functions* over R^d using **modular forms** to give a strict upper bound on the Cohn-Elkies LP optimum (proving non-sharpness in d=12,16). The analog for White's program would be: construct an explicit f satisfying all of White's SDP constraints (cell bounds, Fourier moments, **Bochner-PSD at every level n**) with objective Ω ≥ 0.380X. The challenge is the Bochner-PSD constraints at *every* level n simultaneously, which in the limit is equivalent to f being a probability density on [0,2] (Bochner's theorem). So the construction reduces to: **find an explicit density h* on [0,2] with mass 1, h* ∈ [0,1] pointwise, and `inf_k R_{h*}(k) ≥ 1 − 0.380X` exactly.** This is essentially the same problem as constructing Together's upper bound, but with **exact arithmetic / certified bounds** instead of step-function search. Probably **months of work**; the right experts to involve are White and the Together team.

**5b. Empirical-to-rigorous saturation via Lasserre-tail-bound style.** The Lasserre tail-bound analysis (already in `communications/lasserre_tail_bound.md`) shows that one specific augmentation cannot push past the current bound at tractable T_max. A more systematic version: for each augmentation family (Bochner level n, Lasserre level k, M-side Bochner, Hankel, iterated covers), derive an *unconditional* tail-bound that quantifies the gain at level → ∞ in terms of an explicit residual. If every residual is small, this gives a saturation-style theorem (modulo finite levels). This is **weeks of work**, with a real chance of yielding a paper-worthy theorem of the form "no PSD augmentation in the family F can yield µ ≥ 0.3801279 + δ for δ > δ*". Less clean than 5a but tractable in-session.

**Concrete next action:** Pivot from 5a/5b debate to action 5b, since it is in the session's reach. Specifically: enumerate the augmentation family currently in `white_full_convex.py` (`bochner_n`, `mside_bochner_n`, `lasserre_T_max`, `use_T5p`, etc.) and for each, write down the analog of the Lasserre tail-bound argument. If all of them admit a finite, computable "ceiling gain" as level → ∞, the cumulative bound is the saturation theorem. **This does not require Li's paper at all** — it is a self-contained tail-bound analysis specific to White's program.

---

## Appendix — Specific Citations Used

- **Problem 2.2** (p. 5): Cohn-Elkies LP statement.
- **Problem 3.1** (p. 5): Discrete Cohn-Elkies LP on Z_m^d.
- **Theorem 3.2** (p. 6): The discrete reduction. Hypothesis `m ≥ 2r`. Proof restricts to Z^d and periodizes to Z_m^d; the periodization step preserves `g_m(x) ≤ 0 for |x| ≥ r` only because m ≥ 2r ensures all shifts `x + n` for nonzero `n ∈ mZ^d` have `|x+n| ≥ r` in Z^d.
- **Corollary 3.3** (p. 8): Cohn-Elkies LP ≥ discrete Cohn-Elkies LP.
- **Problem 4.1** (p. 9): Discrete dual LP. Five constraint families: λ(y) ≥ 0, λ = µ̂, µ(x) ≥ 0 for |x| ≥ r, µ(x) = 0 for 0 < |x| < r, µ(0) = (r/2)^d. Objective maximize m^{-d/2} λ(0).
- **§§6-8** (pp. 11-15): Numerical dual feasible points in d = 3, 4, 5, 6-12. Each is a rational vector certified by interval arithmetic to be feasible for Problem 4.1, with objective exceeding Hales' / Cohn-de Laat-Salmon's known upper bound on the packing density.
- **§9** (pp. 15-18): General odd-d construction, but the bound 1/(2(d+1)) decays faster than known packing densities for d > 13, so the construction doesn't extend the non-sharpness conclusion beyond what is achievable by the §§6-8 ad hoc method.
- **§10** (p. 19): Open problem — rigorous convergence of discrete LP to continuous LP as m, r → ∞ is left unresolved. (This is also a problem for any Erdős-side adaptation.)

---

## Commit instructions

```
cd /Users/benzanghi/Documents/Claude/Projects/Erdos
git add LEVER_I_DEEP_READ.md
git commit -m "Lever I deep read: arXiv:2206.09876 applicability to White's Bochner-PSD hierarchy (verdict: DOES_NOT_APPLY)"
```
