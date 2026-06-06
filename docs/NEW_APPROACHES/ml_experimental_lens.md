# New Approaches — ML / Learning-Augmented & Experimental-Math Lens

**Date:** 2026-06-06
**Lens:** ml-experimental (learn a dual feasible certificate / symbolic-regress the optimal structure / PINN-solve the Euler–Lagrange)
**Scope:** NEW attacks on µ that could change the *value* of the LB (past the ~0.380558 SDP-framework ceiling) or give a fundamentally different proof. NOT a re-certification of the existing 0.380284. Proposals only — no claimed results.

## Why this lens is well-aimed

The UB side has already been won by *learning the primal* (AlphaEvolve, TTT-Discover evolve a good step-function `h`). The symmetric, *unexploited* move is to **learn the dual** — the lower-bound certificate — and then certify it rigorously. Three project facts make this the right target:

1. **A clean dual exists and is orthogonal to the mined SDP lane.** PRO-29 / PRO-26 establish the Fenchel–Sion dual `µ = sup_γ inf_h Σ_t γ_t M(h,t)`, where `γ` is a **probability measure over shifts t**. The current SDP never optimizes `γ` as a free variable — the cell-envelope (`white_full_convex.py:176–190`) is a *fixed-grid* surrogate for one piece of this, and PRO-22 showed it is load-bearing precisely because it enforces the `(a,b)↔(c,d)` Parseval link. The dual-`γ` measure-LP is a different object in a different variable.

2. **The "blocked" status is a misread.** PRO-23 tested whether *Together's primal h\** satisfies the KKT functional equation and found it does not (residual 1e-2 vs 1e-9) — proving `µ < 0.380871` strictly, but NOT that the dual route is dead. Nobody has ever *optimized γ directly as the bounding object*. The dual was named as the intended LB method in the very first kickoff (`min_overlap_report.md:211`: "maximize λ s.t. a nonnegative measure μ on shifts with λ ≤ ∫ J(t) dμ(t)") and in `OUT_OF_BOX_CROSS_DOMAIN.md:40`, then abandoned in favor of White's cell-envelope. **It was proposed and never built.**

3. **The residual gap is structural, not numerical (TOGETHER_DIAGNOSTIC).** The deviation between the SDP-optimal `f̃` and the true optimizer is **99.9% low-frequency** (smooth structural disagreement), not Gibbs ringing. This is exactly the regime where a *learned low-dimensional structural ansatz* on the optimizer (and its dual support) beats throwing more PSD constraints at a Fourier basis.

The retracted Lasserre attempt was a **primal-f-side** moment lift `(f²)̂(m)` with a bad tail bound. Everything below operates on the **dual t-side measure** or the **real-space free boundary** — disjoint objects.

---

## Approach 1 — Learn-and-certify the dual γ-measure on shifts ("learning to certify" the LB)

**Technique.** Treat the Fenchel–Sion dual `µ = sup_{γ ∈ Δ(shifts)} V(γ)` where, using `M(h,t) = 1 − ⟨h, T_t h⟩` (PRO-29) and `Σγ=1`,
```
V(γ) = inf_{ h:[0,2]→[0,1], ∫h=1 }  [ 1 − ∫ h(x) (Σ_t γ_t h(x+t)) dx ].
```
For **any** feasible `γ`, `V(γ)` is a rigorous lower bound on µ. The inner minimization is a concave-in-`γ`, box-constrained **indefinite QP** in `h` (kernel `K_γ(x,y) = Σ_t γ_t·[shift kernels]`); its global min can be **lower-bounded rigorously** by a Bochner/SOS certificate on the cosine-polynomial `1 − \hat K_γ` (positive-trig-polynomial cone — the *same* cone the project already trusts, but applied to the dual kernel, not the primal moment matrix). So the pipeline is: **(i) LEARN a good γ** by gradient ascent / evolutionary search on `V(γ)` over a parametric family (mixture of `atoms` at learnable shift locations + smooth background — motivated by PRO-23's finding that the active set is ~437 near-degenerate shifts), then **(ii) CERTIFY** the learned γ by solving the small fixed-γ SOS feasibility problem at high precision (SDPA-GMP), yielding a number-theoretically clean, *White-machinery-independent* LB. The learner only needs to find a *good* γ; rigor comes entirely from step (ii), so the network/optimizer is never in the trusted path (the "neural certificate" discipline: synthesize with ML, verify with a sound checker — Abate et al. 2024; de Klerk–Laurent GPM duality, arXiv:1811.05439).

**Citation.** de Klerk & Laurent, *A survey of SDP approaches to the Generalized Problem of Moments and their error analysis*, arXiv:1811.05439 (the dual-measure-LP + SOS-certificate framework with convergence/error guarantees). Lasserre, *A semidefinite programming approach to the generalized problem of moments*, Math. Prog. 2008. Neural-certificate "synthesize-then-verify" discipline: Abate et al. (neural certificates, 2024) and the soundness-of-certifier line (arXiv:2504.04542, 2025). Original-problem precedent that the dual is the intended LB object: this repo's own `min_overlap_report.md:211`.

**Why untested here (vs the ledger).** The ledger's "White §5 LP/SDP base" and "Bochner-PSD on f, 1−f" both live on the **primal/f-side**: they bound `Ω(f)` and enforce `f∈[0,1]` via PSD on `f` and `1−f`. They never instantiate `γ` (the measure over shifts) as a free optimization variable, never *learn* it, and never produce a standalone fixed-γ SOS certificate. The retracted Lasserre is a **primal** `(f²)̂(m)` lift — different variable, different tail issue. PRO-23's KKT work only *checked Together's primal* against the dual stationarity condition; it never *optimized γ*. This is the dual the kickoff proposed and the project never built.

**First experiment (cheap).** Single `git`-free script `_dual_gamma_poc.py`:
1. Discretize shifts on a coarse grid `t ∈ {k/R : k=0..2R}`, `R≈50`. Build the kernel `K_γ` as a function of weights `γ_k ≥ 0, Σγ=1`.
2. Compute `V(γ)` for two sanity γ's: (a) uniform γ, (b) γ = empirical histogram of PRO-23's 437 active shifts (already in `lp_research_state/data/together_active_set.json`). For each, lower-bound the inner indefinite QP two ways — a quick `cvxpy` SDP relaxation of the box-QP, and a direct CLARABEL solve of the convex *dual* — and read `V(γ)`.
3. **Decision metric:** does `V(γ_active) > 0.30`, and is `dV/dγ` informative (nonzero gradient pointing to a better γ)? If yes, scale to a 100–200-shift learnable mixture optimized by projected gradient ascent; final γ certified once at GMP precision. Total first-pass: a few hours, no new heavy infra (reuses the existing cvxpy/CLARABEL + SDPA-GMP stack).

**Why it could beat the SDP saturation.** The framework ceiling (`C_∞≈0.380558`, PRO-6) is a property of the *cell-envelope + f-side Bochner* relaxation specifically — PRO-22 showed the cell-envelope is necessary for *that* encoding's validity, i.e. the ceiling is intrinsic **to that encoding**. The dual-γ certificate does not contain the cell-envelope at all; its slack is governed by how well a measure-on-shifts + SOS certificate approximates `inf_h`, which is a *different* approximation with *no a-priori reason* to saturate at 0.380558. A measure with atoms placed by a learner exactly where the optimizer's autocorrelation binds (PRO-23 active set) can, in principle, realize the GPM value to the SOS-truncation error — and the GPM hierarchy is known to converge to the true optimum.

**Risk.** (a) The inner box-QP is indefinite; the SDP relaxation of it may itself have a relaxation gap, so `V(γ)` computed via relaxation could be *below* the true `V(γ)` — that's still a valid (if weaker) LB, but might not beat 0.380558. (b) The dual of an `inf-sup` over a non-convex inner set can have its *own* duality gap (the `sup_t M` is a sup of bilinear forms, non-concave in h — `min_overlap_report.md:222–227` flags exactly this), so `sup_γ V(γ)` might be `< µ` by a fixed amount. (c) Certifying at GMP precision for 100+ atoms may be heavy. Mitigation: the PoC's two sanity γ's tell us the gap size *before* any investment.

---

## Approach 2 — PINN / neural free-boundary solve of the KKT Euler–Lagrange (discover the true optimizer + its dual support)

**Technique.** PRO-23 derived the exact stationarity system the *true* optimizer must satisfy:
```
Σ_{t∈S} γ_t [h*(x+t) + h*(x−t)]   =  κ     where h*(x) ∈ (0,1)   (interior / Euler–Lagrange)
                                  ≥  κ     where h*(x) = 0       (lower-active)
                                  ≤  κ     where h*(x) = 1       (upper-active)
```
with `S = argmax_t M(t)`, `γ ∈ Δ(S)`, `κ` scalar. This is a **free-boundary (obstacle-type) variational problem**: the unknowns are the function `h*`, the *free boundaries* separating the bang-bang sets `{h=0}, {h=1}, interior`, the *active shift set* `S`, and the dual weights `γ`. Solve it with a **two-network PINN**: one net `h_θ(x)` for the density (with a smooth surrogate for the box `[0,1]` via a sigmoid head), one net/parametrization for the free boundaries and the shift-measure `γ_φ` (mixture of learnable atoms), trained on the residual of the coupled system (interior EL residual + complementarity on the active sets + the max-condition that `M(t)=κ_M` exactly on `S` and `≤` off it). PINNs for free-boundary / obstacle / bang-bang optimal-control problems are now standard (two-network "predict the free boundary + the state" architectures; Pontryagin-informed nets for bang-bang). The payoff is a **high-resolution, structurally-resolved candidate `h*`** that (a) tightens the UB toward the *true* µ (PRO-23 proved Together's is slack), and (b) hands Approach 1 a *learned γ-support* `S` — closing the loop.

**Citation.** Wang et al., two-network PINN for moving/free boundaries; "A physics-informed neural network framework for obstacle-related equations," arXiv:2304.03552; "Single-Loop Bilevel Deep Learning for Optimal Control of Obstacle Problems," arXiv:2601.04120; Physics-Informed Pontryagin Neural Networks for path-constrained / bang-bang control (AIAA J. GCD). Free-boundary EL structure for this exact problem: this repo's `LEVER_FUNCTIONAL_EQUATION.md` (PRO-23).

**Why untested here (vs the ledger).** The "spectral/translation-operator min-max" entry was a *Rayleigh-quotient* bound (4× loose, PRO-29) — a different object. PRO-23 *wrote down* the EL system but only **evaluated Together's h\* against it** (a verification, which failed); it never *solved it forward* for the unknown `(h*, free boundaries, S, γ)`. The Rechnitzer ansatz (ledger: "UB-side only") is a *fixed* `(1−4x²)^{j−1/2}` real-space basis for the *smooth* `ν₂²` problem and explicitly does **not** transfer to the bang-bang min-max (PRO-26 §3). No PINN / neural free-boundary solve has been attempted. Solving the EL system is categorically different from fitting a prescribed basis.

**First experiment (cheap).** `_kkt_pinn_poc.py`, 1-D, small:
1. Fix the active set to PRO-23's 437 shifts (or a 12-atom coarsening) and the box via a sigmoid; train `h_θ(x)` on a 2000-point grid to minimize the interior-EL residual `‖Σγ_t[h(x+t)+h(x−t)] − κ‖²` + complementarity penalties, with `γ, κ` co-optimized (a few thousand Adam steps — minutes on CPU).
2. **Decision metric:** does the trained `h_θ` reach KKT residual `< 1e-4` (PRO-23's best with Together's *fixed* h* was 7.6e-3)? Does `max_t M(h_θ;t) < 0.380871` (beat Together)? If residual stays `~1e-2`, that is itself a strong signal (supports the conjecture that the optimal active set is *not* 437 shifts — informs the structure).
3. If promising, let the free boundaries and `S` become learnable and re-solve.

**Why it could beat the SDP saturation.** The SDP framework's limit is a *relaxation* property; the EL system is the *exact* optimality condition. A converged PINN solution gives a primal `h*` whose `max_t M` is a true UB *and*, if it satisfies KKT tightly, pins the true µ from above to many digits — directly shrinking the open gap independent of the 0.380558 relaxation ceiling. Crucially it also *discovers the dual support* that Approach 1 needs, where blind shift-grids waste atoms.

**Risk.** (a) Bang-bang gradients vanish on the active sets (`(1−h)h` factor — PRO-26 §4 flags this), so training may stall on flat directions. Mitigation: log-barrier / entropy regularization on h, or train on the dual `γ` first. (b) PINNs give *approximate* solutions; to claim a rigorous UB you still must evaluate `max_t M(h_θ)` on a verified grid (cheap, FFT-based) — that part is sound, the PINN is only a *proposer*. (c) The optimizer may be genuinely aperiodic / infinitely-many-breakpoints (OUT_OF_BOX_CROSS_DOMAIN §3d–3e raises this), in which case a finite-net ansatz plateaus. (d) Non-uniqueness of the EL system's solutions (multiple KKT points) could mislead.

---

## Approach 3 — Symbolic regression on the optimizer's structure → exact low-parameter ansatz (FunSearch-style, on the dual/structural side)

**Technique.** Instead of hunting a closed form for the *scalar* µ (the ruled-out wide-basis PSLQ, NEGATIVE to 50 digits), **symbolic-regress the *functional structure* of the optimizer and its dual**: (i) the *breakpoint locations* and interior profile of the high-resolution `h*` (from Approach 2 or a refined Together construction), and (ii) the *atom locations / weights* of the dual measure `γ` (from Approach 1). Feed (breakpoint index → position) and (x → interior value) tables to a symbolic-regression engine (PySR / deep generative SR / LLM-guided SR) with a physics-flavored basis (algebraic, `√`, trig, rational) to propose a **closed-form ansatz with O(10) parameters**. A recognized structural law (e.g. breakpoints at `x_k = g(k)` for simple `g`, or interior `h*(x)=φ(x;a)`) collapses the infinite-dimensional program to a *finite* system that can be solved and **certified exactly** (rational SOS / interval arithmetic) — turning an empirical optimizer into a theorem. This is the FunSearch paradigm (LLM/evolution discovering a *construction* in extremal combinatorics — the cap-set program) redirected from "build a better step function" to "discover the *law* governing the extremal structure."

**Citation.** Romera-Paredes et al., *Mathematical discoveries from program search with LLMs* (FunSearch), Nature 2024 (cap-set extremal construction). Cranmer, PySR / symbolic regression for scientific laws; Deep Generative Symbolic Regression, arXiv:2401.00282; "Enhancing Symbolic Regression with Quality-Diversity and Physics-Inspired Constraints," arXiv:2503.19043; LLM+SR, arXiv:2505.07956. The repo's own structural conjectures to test: symmetry of optimizers and finite breakpoint count (`min_overlap_report.md:215–229`).

**Why untested here (vs the ledger).** The ruled-out entry is **"wide-basis PSLQ closed-form hunt (NEGATIVE to 50 digits)"** — that hunted an integer relation for the *single number* µ against a constant basis. Symbolic regression on the *function* `h*(·)` and the *measure* `γ` is a different target (a structural law, not a scalar identity) and a different tool (SR over expression trees with operators, not PSLQ over a fixed numeric basis). The "Rechnitzer real-space ansatz" entry is a *hand-derived, fixed* `(1−4x²)^{j−1/2}` family for the smooth L² problem (UB-side, doesn't transfer per PRO-26); SR *discovers* the ansatz from data rather than assuming it, and targets the bang-bang min-max optimizer's breakpoints. No automated structure-discovery on the optimizer has been run.

**First experiment (cheap).** `_symreg_structure_poc.py`:
1. Take Together's 600-cell `h*` (`lp_research_state/data/together_f_star.json`) and PRO-23's active-shift set. Extract: breakpoint positions (cell boundaries where h* crosses 0↔interior↔1), the interior sample profile, and the active-shift positions.
2. Run PySR (CPU, minutes) on three tiny datasets: `k ↦ breakpoint_k`, `x ↦ h*_interior(x)`, `j ↦ active_shift_j`. Restrict operators to `{+,−,×,÷,√,sin,cos,1/(·)}`.
3. **Decision metric:** does any expression fit with `R² > 0.999` at *low complexity* (≤ ~8 nodes) AND survive on a *held-out higher-resolution* optimizer (from Approach 2)? A law that is stable across resolutions is a real structural signal; one that drifts is overfitting (the F5 PSLQ failure mode — `LEVER_F5_EMPIRICAL_PHYSICS.md` — to be explicitly guarded against by the cross-resolution holdout).
4. If a law survives, instantiate the finite ansatz, solve for its O(10) parameters, and attempt an exact (rational-SOS / interval) certificate.

**Why it could beat the SDP saturation.** If the extremal structure has a simple closed form (the TOGETHER_DIAGNOSTIC's 99.9%-low-frequency finding *suggests* low structural complexity), a finite exact ansatz sidesteps *all* relaxation gaps — the SDP ceiling, the truncation tails, the cell-envelope — and yields µ to arbitrary precision from a closed system. This is the only route here that could in principle produce an *exact* µ (closing the gap, not just narrowing it).

**Risk.** (a) High prior probability the optimizer has *no* elementary closed form (the F5/PSLQ negatives, and the aperiodicity worry in OUT_OF_BOX §3d–3e). (b) SR's NP-hardness (arXiv:2207.01018) → may return only fitting noise; the cross-resolution holdout is essential but not foolproof. (c) Even a correct ansatz must be *certified* exactly, which can be hard for transcendental breakpoints. This is the highest-risk / highest-payoff of the three; best run *after* Approach 2 supplies a higher-resolution, less-locally-stuck optimizer than Together's.

---

## How the three compose (not three isolated bets)

They form a **closed learning loop on the dual/structural side**, mirroring the AlphaEvolve loop on the primal side:
- **A2 (PINN)** discovers a high-resolution optimizer `h*` and its active-shift support `S` → tightens the UB toward the true µ (PRO-23 proved headroom) and seeds the others.
- **A1 (dual γ)** uses `S` to place atoms, learns `γ`, and emits a *rigorous, White-independent LB* certificate (the actual goal: change the LB value).
- **A3 (SymReg)** reads `h*` and `γ` from A2/A1 and tries to compress them into an *exact* finite ansatz → potentially an exact µ.

The single highest-value, lowest-cost first move is **Approach 1's PoC** (evaluate `V(γ)` at the PRO-23 active-set γ): it is hours of work on existing infra, it directly tests whether the *dual lane* can clear 0.380558, and its outcome (the gap size) decides whether to invest in A2/A3.

## Honest cross-cutting caveat

All three are **proposers**, not provers: the LB rigor in A1 and the closed form in A3 must pass an independent sound checker (SDPA-GMP SOS / interval arithmetic / the project's `_independent` re-implementation policy) before any number is quoted — exactly the discipline the retracted Lasserre claim violated (tail-bound trap, `project_tail_bound_rigor_trap.md`). The shared failure mode to watch is **fitting noise masquerading as structure** (F5/PSLQ); the cross-resolution holdout (A3) and the fixed-γ GMP certificate (A1) are the designed guards.
