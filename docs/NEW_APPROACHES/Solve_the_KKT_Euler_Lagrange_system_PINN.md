# Approach assessment — Solve the KKT / Euler–Lagrange system: PINN proposer → interval-Krawczyk certifier (MERGED ml + harmonic)

**Date:** 2026-06-06
**Assessor:** Claude (goal-workflow subagent, adversarial-vetting pass)
**Verdict:** **PROMISING (conditional)** — strong *intermediate* win available (a KKT-tight float-Newton `h*` that tightens the UB and unblocks PRO-23 Step 4 and feeds the dual-γ lane), but the headline claim (a *rigorous two-sided enclosure of the true μ* via one interval-Krawczyk step) is **weak/at-risk**, blocked by a structural obstruction this repo has already measured: the optimizer's active set is a near-degenerate ~222-atom near-continuum, not the small finite set the certifier needs.

This is **PROPOSAL + VETTING ONLY**. No bound is claimed. The cheap checks below are analytic/data inspections (no heavy solves).

---

## 0. What the approach is, precisely

Two stages bolted to PRO-23's exact stationarity system (`LEVER_FUNCTIONAL_EQUATION.md`):

- **EL / KKT free-boundary system.** For the true optimizer `h*:[0,2]→[0,1]`, `∫h*=1`, with active shift set `S = argmax_t M(t)`, dual weights `γ_t ≥ 0` (`Σγ=1`) and scalar `κ`:
  ```
  Σ_{t∈S} γ_t [ h*(x+t) + h*(x−t) ]  =  κ     on the interior {0<h*<1}   (Euler–Lagrange)
                                      ≥  κ     on the lower set {h*=0}
                                      ≤  κ     on the upper set {h*=1}
  ```
- **Stage (a) — PINN proposer.** Two networks: `h_θ(x)` with a sigmoid `[0,1]` head, plus learnable free-boundaries and a dual-atom parametrization `γ_φ`. Train on the coupled residual (interior EL + complementarity + the max-condition `M(t)=κ_M` on `S`, `≤` off it). Output: a structurally-resolved candidate `h*` and its active set `S`. `max_t M(h_θ)` evaluated on a verified grid is a true UB (the net is only a *proposer*).
- **Stage (b) — interval-Krawczyk certifier.** Reduce to a finite square system `F(z)=0`, `z = (breakpoints x_i, weights γ_t, shifts t_i, κ, μ)`. Float-Newton to a candidate `ẑ`; then one **interval-Krawczyk** step `K(Z) = ẑ − C·F(ẑ) + (I − C·DF(Z))(Z − ẑ)` on a box `Z ∋ ẑ`. If `K(Z) ⊂ int(Z)` ⇒ unique root in `Z` ⇒ a **rigorous two-sided enclosure of μ**.

The PINN-proposer half is essentially **Approach 2 of `ml_experimental_lens.md`** (already on file). The genuinely new contribution of *this* merged proposal is the **interval-Krawczyk certifier on the finite KKT reduction** — the part that promises to *close the gap two-sidedly* rather than nibble the LB.

---

## 1. FEASIBILITY with available tools

**Tooling present (verified this session):**
- `mpmath 1.3.0` with `mpmath.iv` interval arithmetic (`from mpmath import iv`) — works to arbitrary `dps`. This is a usable, if slow, backend for an interval-Krawczyk step: build `DF` as an interval matrix, do interval linear algebra by hand. There is **no `python-flint`/`arb`** (so no fast rigorous Bessel/ball-arithmetic à la Rechnitzer) and **no INTLAB** (MATLAB `verifynlss` not available). Krawczyk must be hand-rolled on `mpmath.iv`.
- **No `torch` / `jax`** in `.venv` (both import-fail). A PINN therefore needs either an install or a pure-numpy/`mpmath` autodiff. Given the EL residual is an explicit, differentiable functional of a *piecewise* `h` with `O(K)` breakpoints, a "PINN" is overkill: a **plain parametric least-squares / float-Newton on the finite reduction** (numpy + `scipy.optimize`, both available) is the right tool and far cheaper. The "neural" framing buys little here; the cheap intermediate win does not need a net.
- `cvxpy`/CLARABEL, SDPA-GMP present (for cross-checks / the LP feasibility re-test of any candidate `S,γ`).

**Feasibility verdict by stage:**
- **Stage (a) intermediate win — float-Newton/LS KKT-tight point:** **HIGH feasibility, low effort.** Solving `min_{h,γ,κ} ‖interior-EL residual‖² + complementarity penalties` with the active set *fixed* is exactly PRO-23's Step 3 LP **plus** making `h` a free variable (PRO-23 held `h` = Together's and only solved for `γ,κ`). Letting `h` move is a modest extension of existing code. Reaching residual `≪ 1e-9` for a *self-consistent* (h, γ, S) is plausible — but see §5: it may converge to a *different, smaller-support* optimum than Together's, which is the interesting outcome.
- **Stage (b) Krawczyk certificate — LOW feasibility at the true active set.** The finite reduction is square only if `|S|` (number of γ unknowns) and the number of active/complementarity constraints balance the breakpoints. With the measured `|S| ≈ 222` distinct shift magnitudes (444 with signs), `z` has several hundred components and `DF(Z)` is a several-hundred-square **near-singular** interval matrix (degeneracy quantified in §3). `K(Z) ⊂ int(Z)` requires `‖I − C·DF‖ < 1` with `C ≈ DF(ẑ)⁻¹`; near-singular `DF` makes `C` huge and the interval inflation `(I−C·DF(Z))(Z−ẑ)` explode — the test fails. The proposal's own mitigation (m = 2–4 atom reduction) is **falsified by the data** (§3): a 2–4 atom `S` cannot satisfy KKT (PRO-23 got residual 7.6e-3 with a *richer* set), so a Krawczyk certificate on a tiny reduction would certify a root of the *wrong* system.

**Bottom line on feasibility:** the *proposer*/intermediate-win half is feasible and cheap with on-hand tools (no net needed). The *certifier* half — the only part that delivers the advertised two-sided μ enclosure — is **not feasible at reasonable effort** with the actual active-set structure, for structural (degeneracy) reasons, not merely tooling reasons.

---

## 2. PRIOR ART

**The technique is real and mature** (citations verified this session via web search):

- **Interval Newton / Krawczyk for verified roots & first-order-optimality systems:** Krawczyk (1969); Moore (*Interval Analysis*, 1966); **Rump** added the uniqueness test and `verifynlss` in INTLAB (Rump, *Acta Numerica* 2010 survey; INTLAB toolbox). Modern usage for *systems of first-order optimality equations*: Krawczyk verification of KKT zeros is explicitly an application area. Recent certified-root pipelines using the interval Krawczyk test: "An improved verification algorithm for nonlinear systems based on Krawczyk operator" (*J. Comput. Appl. Math.* 2019); "Certifying zeros of polynomial systems using interval arithmetic" (arXiv:2011.05000); "Certified surface approximations using the interval Krawczyk test" (arXiv:2602.07718, 2026); `NumericalCertification` in Macaulay2 (arXiv:2208.01784). These precedents are **low-dimensional polynomial systems with isolated, non-degenerate roots** — the regime where Krawczyk shines.

- **Verified solutions of extremal-constant / variational problems:** the closest *successful* sibling is **Rechnitzer, "The first 128 digits of an autoconvolution inequality," arXiv:2602.07292 (2026)** — confirmed real. But note carefully (this is decisive for prior-art transfer): Rechnitzer does **NOT** certify a *root of a KKT system* by Krawczyk. He (i) posits a **smooth** inverse-square-root real-space ansatz `Σ a_j (1−4x²)^{j−1/2}`, (ii) evaluates the L² objective rigorously via **Arb ball arithmetic + Bessel asymptotics** with a Hurwitz-zeta tail bound, and (iii) gets the LB from a **Hölder–Plancherel** dual, not from a Newton/Krawczyk root. PRO-26 already analyzed this and concluded the smooth ansatz and the Hölder LB **do not transfer** to our bang-bang min-max problem (`PRO26_RECHNITZER_ANALYSIS.md` §3). So Rechnitzer is *feasibility evidence that high-precision rigorous certification of a sibling constant is possible*, but it is **not** evidence that the *Krawczyk-on-KKT-root* mechanism works for *this* problem class — different optimizer geometry (smooth singularity vs bang-bang free boundary), different LB mechanism.

- **PINN free-boundary / obstacle / bang-bang control:** obstacle-PINN (arXiv:2304.03552), single-loop bilevel control of obstacle problems (arXiv:2601.04120), Pontryagin-informed nets for bang-bang. Real and standard, but these are *proposers* (approximate solvers), never the rigorous part.

**What has been tried *here*:**
- **PRO-23 / `LEVER_FUNCTIONAL_EQUATION.md`:** derived the exact EL/KKT system (Step 2) and ran the LP feasibility test (Step 3) — but **only evaluated Together's fixed `h*`** against it (best residual **7.6e-3**, vs `1e-9` needed). It never forward-solved for the unknown `(h*, S, γ)`, and never attempted any interval certificate. So the *forward solve* is genuinely untested.
- **`ml_experimental_lens.md` Approach 2** already proposes the **PINN proposer** half (same EL system, same two-network idea, same "PINN is only a proposer; verify `max_t M` on a grid" discipline). **This merged proposal's PINN stage is therefore a near-duplicate of an existing on-file proposal.** Its *new* element over Approach 2 is the **interval-Krawczyk certifier**.
- **`computer_assisted_formal_global_optimization.md` Approach C** already proposes **rigorous interval branch-and-bound over a bang-bang finite reduction** with the SDP as the box oracle — i.e. the *rigorous-interval-over-finite-reduction* idea, but as exhaustive B&B (limit = μ by exhaustion) rather than a single Krawczyk uniqueness step. The present proposal is a **lighter-weight cousin of Approach C** (Krawczyk on one Newton root, instead of B&B over the whole box).

**Net prior-art finding:** Neither half is fully novel to the repo — the proposer ≈ existing Approach 2, the rigorous-interval idea ≈ existing Approach C. The *specific composition* (PINN-seed → float-Newton → single Krawczyk uniqueness step on the KKT reduction, for a two-sided μ enclosure) is **a new combination not on file**, and the interval-Krawczyk-on-KKT-root mechanism specifically has **no precedent for this problem's degenerate-free-boundary structure**.

---

## 3. FIRST EXPERIMENT — done (data inspection, no heavy solve)

The cheapest decisive signal is *structural*: does the finite reduction `F(z)=0` even exist at small, well-conditioned dimension? I inspected `lp_research_state/data/together_active_set.json` (the saved PRO-23 active set). **Findings (run this session):**

| Active tolerance (gap ≤) | # shifts |
|---|---|
| 1e-12 | **2**  (just literal max ±33) |
| 1e-10 | **12** |
| 1e-9  | **437** |
| 1e-6  | 437 |
| 1e-5  | 444 |

- **Distinct active shift magnitudes |j|: 222**, of which **217 are ±-symmetric** — so the active set is essentially `{±j}` over 222 magnitudes.
- **Active |j| spans 22 … 256** — a *wide band*, NOT a tight cluster near the max ±33.
- The gap profile is a **near-continuum of near-ties**: jumping the tolerance from 1e-10 to 1e-9 takes `|S|` from 12 → 437. There is no clean small plateau.

**Interpretation (the make-or-break read):**
1. The proposal's stated cheap-falsification mitigation — *"smallest symmetric m = 2–4 reduction (cheap to falsify)"* — is **already falsified by this data**. A 2–4-atom `S` (e.g. `{±33}`, or `{±33,±30}`) is not the true active set; PRO-23 showed even the full 437-set only reaches KKT residual 7.6e-3, and the 12-set residual is *worse* (0.3). So a finite reduction built on a tiny `S` would certify a root of a *misspecified* system — an invalid μ enclosure.
2. The genuine reduction needs `|S| ≈ 222` magnitudes ⇒ several-hundred-dimensional `z`, with a `DF` that is **structurally near-singular** (437 shifts agreeing to ~1e-9 ⇒ near-linearly-dependent EL rows ⇒ near-zero singular values). This is exactly the regime where the Krawczyk inclusion `K(Z)⊂int(Z)` **cannot hold** (it requires a well-conditioned `DF` so that `‖I−C·DF(Z)‖<1`).
3. Therefore the *cheapest* analytic check already discriminates the verdict: **the certifier half is structurally obstructed**, while the **proposer/intermediate-win half** (forward-solve for a self-consistent small-support `(h,γ,S)` and re-test KKT residual) is the part that could yield signal — and it might find that the *true* optimal support is genuinely smaller than Together's 437 (if Together is at a degenerate non-global stall, as PRO-23 hypothesized), in which case the picture improves. That is the next experiment, but it is a *solve*, not a free check, so it is deferred per task scope.

**Signal delivered:** the proposal's own kill-switch (tiny-m reduction) fires negative on inspection; the certifier's success hinges on a non-degenerate small `S` that the data says does not exist at Together's optimum. Whether a *re-optimized* `h*` has a small clean support is the open question that gates everything.

---

## 4. WHY IT COULD BEAT THE ~0.380558 SATURATION (or close the gap / give a different proof)

The argument is genuinely different in *kind* from the SDP lane, and that is its appeal:

- **It targets the EXACT optimality condition, not a relaxation.** PRO-6's `C_∞ ≈ 0.380558` is a property of the **2-point cell-envelope + Bochner** relaxation (`computer_assisted_formal_global_optimization.md` framing; PRO-22 shows the envelope is load-bearing *for that encoding*). The KKT/EL system is the *exact* stationarity condition of the original `inf_h sup_t M`. A validated root enclosure of `F(z)=0` is, by construction, an enclosure of the **true μ** — it does **not inherit the 0.380558 ceiling** because it never uses the cell-envelope relaxation at all. This is the legitimate reason the approach "sidesteps the relaxation ceiling."
- **It is two-sided.** Almost every other lever on file (3-point hierarchy, native Cohn–Elkies LP, dual-γ certificate) nibbles the **LB**. A Krawczyk enclosure delivers `[μ_lo, μ_hi]` simultaneously — the only route in the whole portfolio whose *limit object is μ itself*, two-sidedly. If it worked, it would close the gap, not just narrow it.
- **Cheap intermediate win is real and independently valuable.** A float-Newton KKT-tight `(h,γ,κ)` with residual `≪ 1e-9` would be **the first genuinely KKT-tight point** for μ. It (i) gives a primal `h` whose `max_t M` is a true UB that is `< 0.380871` (PRO-23 already proved Together's is slack), tightening the UB; (ii) **unblocks PRO-23 Step 4** (the analytical solve, blocked precisely for lack of a KKT-tight input); (iii) hands `ml_experimental_lens` Approach 1 a *clean learned dual support `S`* to place γ-atoms — closing the dual-side loop. This intermediate deliverable does **not** depend on the (obstructed) Krawczyk step and is the strongest reason to do *some* of this work.

So the "beat saturation" logic is sound *if* a finite, well-conditioned KKT reduction exists. The §3 data is the threat to that "if."

---

## 5. RISKS / why it might fail

Ordered by how lethal they are to the headline (two-sided enclosure) claim:

1. **[LETHAL to the certifier] Degenerate / large / near-continuum active set breaks the finite reduction.** Measured: 222 ±-symmetric magnitudes spanning |j|∈[22,256], 437 shifts agreeing to ~1e-9. The KKT system's `DF` is near-singular and several-hundred-dimensional ⇒ `K(Z)⊂int(Z)` will not hold; interval inflation explodes. The proposal's m=2–4 fallback is already falsified (§3). **This is the dominant risk and the §3 check shows it is *realized at Together's optimum*, not hypothetical.** Only escape: the *true* optimum (if Together is at a non-global degenerate stall) has a genuinely small, well-separated support — unproven, and PRO-23's evidence (residual barely improves as `S` grows) does not encourage it.
2. **[LETHAL to rigor if mishandled] Non-uniqueness / multiplicity of KKT points.** `sup_t M` is a sup of bilinear (indefinite) forms — non-concave in `h` (`min_overlap_report.md` flags this). The EL system can have **many** stationary points (the 437-fold degeneracy is a symptom). Krawczyk certifies a root is *unique in the box `Z`*, **not that it is the global optimizer**. Certifying a non-global KKT point gives a rigorous enclosure of the *wrong* number — an *invalid* μ bound presented as valid. This is the same epistemic trap family as the retracted Lasserre tail bound (`project_tail_bound_rigor_trap.md`): a locally-correct certificate of a globally-wrong object. Avoiding it requires an *independent* global-optimality argument (e.g. the dual-γ LB matching), which the approach does not supply on its own.
3. **[Blocks the intermediate win] Bang-bang gradient vanishing.** The objective gradient w.r.t. interior values carries a `(1−h)·h` factor (PRO-26 §4) → flat directions near the active sets; Newton/PINN can stall. Mitigation (log-barrier/entropy on `h`, solve dual `γ` first) is standard but adds tuning.
4. **[Tooling friction] No `arb`/INTLAB/`torch`.** Krawczyk must be hand-built on `mpmath.iv` (slow interval linear algebra at several-hundred dimension is painful); the "PINN" would need a `torch`/`jax` install or be replaced by numpy/`scipy` Newton (recommended anyway). Not fatal, but the "PINN" branding overstates what the tooling supports and what the problem needs.
5. **[Possibly fatal in principle] The optimizer may have infinitely many breakpoints / be aperiodic.** `OUT_OF_BOX_CROSS_DOMAIN.md` §3d–3e raises that the extremal `h` may not have finitely many breakpoints. If so, **no** finite `F(z)=0` reduction exists and the certifier is ill-posed in principle (not just in practice). The wide active band (|j| up to 256) is mildly consistent with rich structure.
6. **[Redundancy] Duplication with on-file proposals.** The proposer is ≈ `ml_experimental_lens` Approach 2; the rigorous-interval-over-finite-reduction idea is ≈ `computer_assisted_formal_global_optimization` Approach C. Spending effort here partly re-treads work already proposed; the *only* non-redundant new piece (Krawczyk-on-KKT) is the piece §3/§5.1 says is obstructed.

---

## 6. Recommendation (honest, scoped)

- **DO** the cheap *intermediate-win* solve **as part of / merged into `ml_experimental_lens` Approach 2 and Approach 1**, NOT as a standalone "PINN→Krawczyk" track: forward-solve a small parametric `(h, γ, κ, S)` by float-Newton/least-squares (numpy/`scipy`, no net) to drive the KKT residual `≪ 1e-9`. Decisive reads: (i) does a *self-consistent* optimum exist with support markedly smaller than 437? (ii) does its `max_t M` beat 0.380871 (tighten UB)? (iii) does it unblock PRO-23 Step 4? This is hours of work on existing infra and has independent value regardless of the certifier.
- **DO NOT** invest in the **interval-Krawczyk certifier** as the path to a two-sided μ enclosure **unless** the intermediate solve *first* discovers a small, well-separated, non-degenerate active set (which the current data argues against). Building several-hundred-dimensional `mpmath.iv` Krawczyk machinery against a near-singular `DF` is high-effort, low-probability.
- **If** a small clean support emerges, the natural rigorous endgame is **Approach C's interval B&B** (gap-closing by exhaustion, robust to multiplicity via the global search) rather than a single Krawczyk step (which certifies *a* root, not *the* optimum) — i.e. fold this into Approach C, supplying it the reduced support.

**EFFORT:** intermediate win ≈ low (1 focused session, existing tools). Full certifier ≈ very high (multi-session, new `mpmath.iv` Krawczyk + the unresolved degeneracy), with low success probability at the measured active-set structure.

---

## 7. One-line verdict

A scientifically *appealing* and genuinely two-sided idea whose **proposer half is a near-duplicate of an existing on-file approach** and whose **novel certifier half is structurally obstructed by the very degeneracy this repo already measured (222-atom near-continuum active set)**. Worth doing **only** for its cheap intermediate KKT-tight `h*` (which tightens the UB, unblocks PRO-23 Step 4, and seeds the dual-γ lane) — fold that into `ml_experimental_lens` Approaches 1–2; defer the Krawczyk certificate, and if a small clean support ever appears, route it into Approach C's interval B&B rather than a lone uniqueness step.
