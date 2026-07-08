# Three best previously-untested attacks on µ (Erdős minimum overlap)

**Date:** 2026-06-06
**Workflow:** goal-trial-1 — final selection pass.
**Status:** PROPOSALS + VETTING ONLY. No bound is claimed. Each item has been adversarially
deep-dived (memos in this directory); decisive *cheap* checks were run, but no heavy solve.

**Current state.** LB µ ≥ 0.380284 (PRO-21, numerically certified augmented SDP); UB µ ≤ 0.380871
(Together 2026, shown NOT tight by PRO-23). Open gap 5.87 × 10⁻⁴. The framework ceiling
**C_∞ ≈ 0.380558** is *provably* a property of White's cell-envelope + Bochner cone (PRO-6:
`f_CB(N|n)=f_C(N)`; PRO-22: dropping the cell-envelope makes the bound INVALID — reconstructed f
exits [0,1], true Ω ≈ 5× reported). **The single most important ranking fact: in-cone moves are
capped at C_∞; a different cone, or solving the exact optimality condition, is not.**

**Selection principle.** The task asks for *diversity of technique*, not three flavors of one idea.
The three below sit in three genuinely different lanes:
1. **higher-arity convex relaxation** (a strictly larger moment cone) — the principled way to *pass* C_∞;
2. **exact optimality + validated numerics** (no relaxation at all) — the only *two-sided* route, a *fundamentally different proof*;
3. **representation theory / numerical *-algebra** (block-diagonalization) — the *enabler* that makes #1 computable at the scale where its gain appears.

All three are "promising"-rated in the deep dives. The four analytic certificate lanes
(NPA noncommutative-moment, Fisher-information/log-Sobolev, sign-uncertainty/Poisson, native
Cohn–Elkies LP) were each vetted to **weak** — they either collapse into the #1 cone, bound the
wrong functional, or lose the load-bearing realizability link; they are recorded as honorable
mentions with the reason each fails, so they are not silently re-proposed.

---

## #1 — Three-point / k-point (triple-correlation / bispectrum) SDP lift

**One-liner.** Lift White's purely 2-point Fourier program to a 3-point relaxation whose decision
object is the *bispectrum* `B(m₁,m₂)=f̂(m₁)f̂(m₂)conj(f̂(m₁+m₂))` of the single density f; impose
its PSD moment/Gram block **alongside** (never replacing) the existing cell-envelope, so the bound
lives in a strictly larger cone than the one proven to cap at C_∞.

**Technique & tools.** The Bachoc–Vallentin / de Laat–Vallentin t-point SDP hierarchy.
- Refs: Bachoc–Vallentin *JAMS* 2008 (arXiv:math/0608426) [3-point SDP]; **de Laat–Vallentin
  *Math. Program.* 2015 (arXiv:1311.3789)** [t-point hierarchy + *convergence theorem* — so the
  question is *at which level* the gain appears, not *whether* the hierarchy is sound];
  Cohn–de Laat–Salmon *Adv. Math.* 2024 (arXiv:2206.09876) [higher-arity proved Cohn–Elkies
  non-sharpness in d=12,16 — the exact precedent for a higher-arity SDP breaking a 2-point
  plateau]; three-point bounds arXiv:2206.15373; **Henrion/Korda–Rudi arXiv:2401.07734 Prop. 7**
  [periodic-Hausdorff truncation = the rigorous cubic-moment tail bound]; ClusteredLowRankSolver
  arXiv:2202.12077 / SDPB arXiv:1502.02033 [block-exploiting solvers for scale].
- Tools: cvxpy 1.8.2 + CLARABEL for the probe (all in-repo); the moment-lift cvxpy plumbing already
  exists in `lp_research_state/code/lasserre3_toeplitz.py` (a trilinear/quadrilinear bordered-LMI
  lift — known-good, re-usable). SDPA-GMP (`lp_research_state/bin/sdpa_gmp`) for the precision
  cross-check the project standardizes on (the kissing-number 3-point literature *required* GMP at
  scale — Mittelman–Vallentin, Machado–Oliveira). Symmetry reduction (#3) is the scaling enabler.

**Why novel vs the ledger.** The retracted **Lasserre level-2** lifted the *degree* of the *same
pairwise* object (it truncated `(f²)̂(m)` with a bad tail bound) — it never moved to a higher number
of points. The bispectrum is the genuinely-new axis: it constrains the **joint phases**
`arg f̂(m₁)+arg f̂(m₂)−arg f̂(m₁+m₂)` that the power spectrum `|f̂|²` — and therefore the *entire*
current program (PRO-22: the load-bearing cell-envelope is the degree-2 Parseval link
`M̂(m)=a_m f̂(m)−4|f̂(m)|²`) — leaves completely free. By Bartelt–Lohmann the bispectrum
magnitude+phase determines the signal up to translation, so B carries *strictly more* information
than `|f̂|²` (the necessary condition for any gain). Grep over `lp_research_state/code/` for
`bispectr|R3|B\(m1,m2|three-point|triple` returns only name-collisions — the genuine-B lift is
**unbuilt** (the existing `lasserre3_toeplitz.py` lifts the *inner aux vars* c_k,d_k, NOT the
bispectrum of f). OUT_OF_BOX_CROSS_DOMAIN §2a explicitly flagged the 3-point lift as never-built
but mis-routed it into a non-sharpness *diagnostic*; the refs were listed, never executed.

**First experiment (concrete next object).** A *cheap analytic precursor was already run this
session and is decisive: a "naive 3-point Gram" indexed by pairs (i,j) has entries
`E[χ_{i+j} conj χ_{k+l}] = f̂` at index *differences* ⇒ it is Hermitian-Toeplitz ⇒ it IS EXACTLY the
Bochner block already in `white_full_convex.py` (probe min-eig ≈ −1e-16, zero gain — the
Fejér–Riesz collapse the `lasserre3_toeplitz.py` docstring itself warns of).* **Consequence: the
approach CANNOT be tested by the easy matrix; you MUST lift the genuine cubic B.** The correct first
solve (≈ 0.5–1 day, a go/no-go): on the **torus** at T₀ = 6–8 (~36–64 index pairs → ~64-dim PSD
block, trivial for CLARABEL even joint with the N=300–3000 program), introduce complex variables
`b[m₁,m₂]`, impose (a) the **linking/localizing PSD tie** connecting b to the program's existing
`(a_m,b_m)`, (b) PSD of the Gram block built from b — **bolted onto the cell-envelope, not
replacing it**. Solve the row-4 center small-N. Three gates: (i) does the dual rise above the
2-point value at the same N,R? (ii) is the b-block **active** (nonzero dual)? (iii) **MANDATORY
PRO-22 validity gate** — reconstruct f, confirm f ∈ [0,1] and reported Ω ≤ true sup_t(f⋆f). The
torus probe is **tail-free**, so it detects any gain *before* re-entering the Lasserre tail trap.

**Payoff.** A *different and higher* ceiling than C_∞. M(k) is intrinsically a correlation count, so
triples encode *which inner shifts can simultaneously be low* — precisely the structure a min-over-
shifts / L^∞ object is expected to need, and exactly where the L² sibling (White autoconvolution
arXiv:2210.16437) is already almost-tight at degree 2 while the L^∞ object is not. If the b-block is
active and lifts the dual, this is the only in-portfolio lever with a principled mechanism to push
the LB past 0.380558 toward (and potentially past) the open gap — a genuine value change, not a
re-certification.

**Risk.** (1) **Transfer mechanism, under-weighted:** the sphere-packing 3-point gain is powered by
the *non-abelian* SO(n) symmetry (stabilizer-subgroup rep theory, SO(n−1) zonal kernels); the Erdős
symmetry is the *abelian* translation group on [0,2] + Z/2, where the extra content is "only" joint
phases and there is **zero positive precedent** for a bispectrum SDP improving a 1-D autocorrelation
*lower* bound (a targeted sweep of merit-factor/LABS/flat-Littlewood literature — Borwein–Choi–Jedwab
arXiv:1205.0626; BBMST *Ann. Math.* 2020 arXiv:1907.09464 — found none). (2) **Cubic-tail rigor
trap:** a rigorous µ-bound needs the truncated cubic moment sum to carry an analytic remainder —
structurally *identical* to the failure that retracted Lasserre-2; the probe must report gain *net
of* the Henrion–Rudi tail, not just at finite T₀. (3) **Modal outcome is bound-neutral:** most
likely the block is inactive at row 4 (no gain, ~1 day spent) or shows a finite-T₀ gain the rigorous
tail erases — a genuine net advance is the high-variance tail, not the mode. (4) O(T₀²)×O(T₀²) blocks
are heavy at useful T₀ without #3's symmetry reduction (itself unbuilt, and over an abelian group
buys less than in the SO(n) setting where it was developed).

---

## #2 — Solve the KKT / Euler–Lagrange system: float-Newton proposer → interval-Krawczyk certifier

**One-liner.** Forward-solve the *exact* free-boundary stationarity system of the original
`inf_h sup_t M` (PRO-23's Euler–Lagrange equations), then validate the unique root with a single
interval-Krawczyk step — a rigorous **two-sided** enclosure of the TRUE µ that never touches the
relaxation cone, hence never inherits the C_∞ ceiling.

**Technique & tools.** Verified-numerics for first-order-optimality systems.
- Refs: Krawczyk (1969); Moore, *Interval Analysis* (1966); **Rump, *Acta Numerica* 2010** [the
  uniqueness test, INTLAB `verifynlss`]; certified-root pipelines using the interval Krawczyk test —
  *J. Comput. Appl. Math.* 2019; "Certifying zeros of polynomial systems" arXiv:2011.05000;
  "Certified surface approximations using the interval Krawczyk test" arXiv:2602.07718 (2026);
  `NumericalCertification` in Macaulay2 arXiv:2208.01784. Sibling-constant feasibility evidence:
  **Rechnitzer, "The first 128 digits of an autoconvolution inequality," arXiv:2602.07292 (2026)**
  (real, confirmed) — rigorous certification of a *sibling L² constant* at extreme precision *is*
  possible. The EL system itself: repo PRO-23 / `docs/archive/LEVER_FUNCTIONAL_EQUATION.md`.
- Tools: numpy + scipy.optimize (Newton / least-squares) for the proposer — both in `.venv`;
  `mpmath 1.3.0` with `mpmath.iv` interval arithmetic (verified present) as the hand-rolled Krawczyk
  backend. (No torch/jax — the "PINN" framing is overkill; a plain parametric float-Newton on the
  finite reduction is the right and cheaper tool. No python-flint/arb, no INTLAB — Krawczyk is
  hand-built on mpmath.iv.) cvxpy/CLARABEL + SDPA-GMP for the LP feasibility re-test of any
  candidate (S, γ).

**Why novel vs the ledger.** This targets the *exact* KKT/EL condition of the original problem, not
a relaxation — it is categorically outside every ruled-out SDP/LP/moment lever, and is the only
route in the portfolio whose limit object *is* µ itself (everything else nibbles the LB from one
side). PRO-23 derived the EL system and ran an LP feasibility test but **only evaluated Together's
fixed h\*** against it (best residual 7.6 × 10⁻³ vs the 1 × 10⁻⁹ needed) — it never forward-solved
for the unknown (h\*, S, γ), and never attempted any interval certificate. The *forward solve* and
the *Krawczyk-on-KKT-root* step are both genuinely untested here. (Honest dedup: the proposer half
≈ `ml_experimental_lens.md` Approach 2; the rigorous-interval idea ≈
`computer_assisted_formal_global_optimization.md` Approach C as B&B. The *specific composition* —
Newton-seed → single Krawczyk uniqueness step on the finite KKT reduction for a two-sided enclosure —
is new and not on file.)

**First experiment (concrete next object).** A *cheap data inspection was already run this session
and sharpens the plan: `lp_research_state/data/together_active_set.json` shows the active set is a
near-continuum of 222 distinct shift magnitudes |j| ∈ [22,256] (444 with signs), with 437 shifts
agreeing to ~1e-9 — so the proposal's own "smallest m=2–4 atom reduction" fallback is already
falsified (PRO-23 got residual 7.6e-3 even with the full 437-set).* The genuine first experiment is
therefore a **SOLVE** (deferred per scope, but it is the gate): the *intermediate-win* float-Newton —
`min_{h,γ,κ} ‖interior-EL residual‖² + complementarity penalties` with h a free variable (PRO-23
held h fixed and solved only γ,κ; letting h move is a modest extension of existing Step-3 code),
seeking a self-consistent (h, γ, S) with residual ≪ 1e-9. The decisive question it answers: does a
*re-optimized* h\* (if Together sits at a degenerate non-global stall, as PRO-23's 437-fold tie
suggests) have a genuinely *small clean support*? If yes, the certifier becomes feasible; if no, the
two-sided certificate is blocked and the run still delivers the intermediate win.

**Payoff.** Two tiers. (a) **Intermediate win (no certifier needed):** a float-Newton KKT-tight
(h, γ, κ) with residual ≪ 1e-9 (vs PRO-23's 7.6e-3) is the *first genuinely KKT-tight point for µ*,
giving a primal h whose max_t M is a true UB < 0.380871 (PRO-23 already proved Together's UB is
slack) — it tightens µ *from above*, unblocks PRO-23 Step 4, and hands #1 / the dual-γ lane a clean
learned support S. (b) **Full certifier:** `K(Z) ⊂ int(Z)` ⇒ unique root ⇒ a rigorous
**[µ_lo, µ_hi] enclosure simultaneously** — the only lever that could *close the gap two-sidedly*
and the only *fundamentally different proof* of µ in the portfolio.

**Risk.** (1) **Lethal to the certifier (already realized in the data):** the measured 222-atom
near-continuum active set ⇒ a several-hundred-dim, near-singular Jacobian DF ⇒
`C ≈ DF(ẑ)⁻¹` is huge and the interval inflation `(I−C·DF(Z))(Z−ẑ)` explodes, so
`K(Z) ⊂ int(Z)` cannot hold; the m=2–4 mitigation is falsified by the data. (2) **Lethal to rigor
if mishandled:** Krawczyk certifies a root *unique-in-box*, **NOT global**. sup_t M is a sup of
indefinite bilinear forms (non-concave in h), so many KKT points exist (the 437-fold degeneracy is a
symptom); certifying a non-global KKT point gives a rigorous enclosure of the *wrong* number — the
same trap family as the retracted Lasserre tail bound. Avoiding it needs an independent
global-optimality argument the approach does not supply. (3) **Blocks even the intermediate win:**
the bang-bang objective gradient carries a (1−h)h factor (PRO-26 §4) ⇒ flat directions near active
sets; Newton can stall (mitigants: entropy/log-barrier on h, solve γ first). (4) **Possibly fatal in
principle:** if the optimizer has infinitely many breakpoints / is aperiodic
(OUT_OF_BOX_CROSS_DOMAIN §3d–3e), *no* finite F(z)=0 reduction exists and the certifier is ill-posed
— the wide active band (|j| up to 256) is mildly consistent with such rich structure.

---

## #3 — Representation-theoretic symmetry reduction of White's SDP (regular *-representation block-diagonalization)

**One-liner.** Block-diagonalize White's Fourier SDP unconditionally via the de
Klerk–Pasechnik–Schrijver regular-*-representation (Schur-decompose the group-averaging projector
`P_G=(1/|G|)Σ_g ρ(g)`), exploiting the reflection x↦2−x and the cosine/sine parity-sign symmetry of
the *constraints* — cutting PSD flops/memory so the same 4 GB hardware reaches the bochner_n ≥ 40–80
regime where #1's (and the existing levers') gains are projected to appear.

**Technique & tools.** Textbook symmetry-reduction for Fourier SDPs.
- Refs: **de Klerk–Pasechnik–Schrijver, "Reduction of symmetric SDPs using the regular
  *-representation," *Math. Program.* 109 (2007) 613–624** (matrix order drops to #orbits of the
  group action — the key number); Bachoc–Vallentin *JAMS* 2008 [the kissing-number precedent];
  Vallentin, "Symmetry in semidefinite programs," arXiv:0706.4233; Gatermann–Parrilo *JPAA* 2004;
  de Klerk, "Numerical block diagonalization of matrix *-algebras," *Math. Program.* 2011. The free
  complex→real U(1) factor: arXiv:2307.11599 and PICOS complex-SDP docs; invariant-SDP survey
  arXiv:1007.2905. It is THE standard tool for Cohn–Elkies-type Fourier SDPs.
- Tools: all in-repo — numpy/cvxpy 1.8.2 (which **natively supports `cp.Variable(hermitian=True)`,
  `H >> 0`** — the one-line lever below), mpmath/sympy, SDPA-GMP. Optional structure-aware export:
  ClusteredLowRankSolver.jl (a new dependency) for a solver that consumes block structure directly.

**Why novel vs the ledger.** The repo's *only* symmetry use is
`lp_research_state/code/symmetric_push.py`, which is a **different, weaker thing**: a *conditional*
even-f assumption (d=0, dlt=0, v=w) feasible only for rows 5,6 (h=0), explicitly NOT an
unconditional bound on µ — a manual variable drop on the optimizer, not an isotypic decomposition of
the *constraint* symmetry. No commutant/isotypic/irrep/regular-representation code exists anywhere
(grep confirms only name-collisions). This exploits symmetry of the *constraints* (the Bochner
moment matrices + cell-envelope are invariant kernels under x↦2−x), so it is **unconditional, all 7
rows** — categorically distinct from `symmetric_push.py`.

**First experiment (concrete next object).** *The decisive structural census was already run this
session (numpy, <1s):* (1) the Bochner real-form `[[Re,−Im],[Im,Re]]` of the (n+1)×(n+1) complex
Hermitian Toeplitz moment matrix has **every eigenvalue of multiplicity exactly 2** (n=6: 7 distinct
values each doubled) — the standard fact that the real embedding block-diagonalizes back into two
identical complex blocks under the SO(2)/U(1) action, a **free ~4× on PSD flops**
(2·(1/2)³ = 1/4) obtainable with **zero representation theory**; (2) the Hermitian Toeplitz block
carries a genuine constraint-level Z/2 centro/persymmetry — verified `S·RF·Sᵀ == RF` for
`S=[[J,0],[0,−J]]` (J = anti-diagonal reversal), another ~2× via centrosymmetric splitting (= the
DKPS #orbits count, ~(n+1)/2); (3) a cone census (get_problem_data/CLARABEL at N=300, T=200)
confirms the PSD blocks grow O(bn²) in svec and **dominate** as bn rises (bn=12 → svec 702;
bn=60 → svec 15006), while the nonneg cone stays ≈N (cell-envelope) and is invariant under the
Bochner group. **Net: realizable symmetry = U(1) × Z/2 ≈ 4–8× on the PSD blocks only; ~4× is free.**
The recommended *actual first step*: do the **one-line complex-Hermitian Bochner swap** in
`bochner.py` (replace the `[[Re,−Im],[Im,Re]]` real embedding with `cp.Variable(hermitian=True)`,
`H >> 0`), cross-check the dual to 10+ digits against `bochner_independent.py` + one SDPA-GMP spot
check, then **re-attempt binding row 4 at the bochner_n that previously OOM'd**. The Z/2
centrosymmetric split and re-expressing each isotypic sub-block as a *separate smaller* cvxpy PSD
constraint (so CLARABEL actually sees the smaller blocks) is the ~1–2-week deliverable.

**Payoff.** It does **not** change the relaxation, so it cannot by itself pass C_∞ (the bound still
asymptotes to ~0.380558 < µ_UB). What it changes is solver **cost**: a 4–8× cut on the PSD blocks
turns LEVER_F3's projected "+4–6 × 10⁻⁴ at bochner_n ≥ 40–80, intractable in 4 GB" from *research*
into an *engineering run* at binding row 4 — realizing a real ~+3 × 10⁻⁴ over the current 0.380284
headline (a landing near ~0.38056). More importantly, it is the **enabler** that makes #1 (the
genuine ceiling-breaker) computable at the T₀ where its gain appears — the highest-leverage
*multiplier* on every already-validated rigorous lever, and the cheapest decisive first step in the
whole portfolio.

**Risk.** (1) **Small group ⇒ modest win:** realized symmetry is only U(1) × Z/2 ≈ 4–8× *on the PSD
blocks*, not the order-of-magnitude the "large-T_max barrier" rhetoric implies — there is no hidden
dihedral/translation group of large order (Toeplitz gives persymmetry, the complex embedding gives
U(1), and that is the entire inventory). (2) **Half the win is free and needs none of the machinery:**
the ~4× U(1) factor is the one-line complex-Hermitian swap; if that is the only delivered gain, the
"regular *-representation/commutant" framing is overkill (a scoping issue, not a correctness one).
(3) **The cell-envelope cone is the floor and is barely touched:** nonneg ≈ N and the SOC block carry
only the v↔w reflection and do *not* shrink under the Bochner group, so at large N / modest bn they
dominate — the PSD reduction's payoff is non-uniform (helps only in the high-bochner_n regime, which
is, fortunately, the targeted regime). (4) **CLARABEL won't consume block structure for free:** it is
a general IPM, so the speedup must be plumbed by re-expressing each isotypic sub-block as a separate
PSD constraint or exporting to a structure-aware solver; without plumbing, forming the commutant
changes nothing the solver sees. (5) **Bookkeeping/rigor surface:** svec scaling, the
`[[Re,−Im],[Im,Re]]`↔complex map, and matching cvxpy canonicalization are classic off-by-√2/sign
traps that can silently corrupt `rigorous_dual_LB` — mitigate with the repo's `_independent`
re-derivation + 10-digit cross-check + one SDPA-GMP spot check. (6) **Ceiling not closed:** even a
perfect reduction at bochner_n=200 asymptotes to C_∞ ≈ 0.380558 < µ_UB = 0.380871 — a strong outcome
(~+3 × 10⁻⁴), not a closure.

---

## Honorable mentions (vetted to *weak*, or out-of-top-3 on tractability — recorded so they are not re-proposed)

- **NPA / Helton–McCullough noncommutative-moment certificate** for `µ = 1 − sup_h inf_t ⟨h, T_t h⟩`
  (Navascués–Pironio–Acín NJP 2008; Pironio–Navascués–Acín SIAM J. Optim. 2010). **Weak — collapses
  into #1's cone:** a session check verified the certified object is exactly the ≤L-point
  correlations of h (length-2k alternating words = k-point correlation), the SAME cone class as the
  bispectrum lift, dominated by it (lighter, repo-native validity gate). The "tail-free" claim is
  also *backwards*: the continuum→Z_m discretization makes Z_m *over*-estimate µ_dual hence
  *under*-estimate µ (an invalid-direction bound, PRO-22 failure shape); a Z_6 toy gave 0.90 vs the
  continuum 0.619. Memo: `NPA_Helton_McCullough_noncommutative_mom.md`.

- **Fisher-information / log-Sobolev convex surrogate** `J(A)=∫(A')²/A` for the L^∞ functional
  (Stam 1959; Blachman 1965; Carlen 1991; Gross 1975). **Weak — wrong functional + non-convex:** a
  session check on Together's h\* found ‖A‖_∞ = A(0) = ‖h‖₂² = 0.7747 ≈ 2× µ, because µ is the sup
  over shifts *bounded away from the origin* (sup_{|t|≥t₀≈2/3} A), which a *global* Fisher-info floor
  does not see; and the claimed DCP encoding fails — `quad_over_lin(A',A)` with A,A' both *quadratic*
  in f is not convex in f (verified UNKNOWN in CVXPY 1.8.2). Memo:
  `Fisher_information_log_Sobolev_convex_su.md`.

- **Sign-uncertainty / Poisson-summation certificate on the autocorrelation A** (absorbs the native
  Cohn–Elkies/Delsarte LP) (Bourgain–Clozel–Kahane 2010; Cohn–Gonçalves 2019;
  Gonçalves–Oliveira e Silva–Steinerberger arXiv:2003.10771). **Weak — no realizability link:** the
  certificate's only structural lever is `Â ≥ 0`, and a session LP showed that lever (plus mass 1,
  A(0) ≤ 1, support, A ≥ 0) admits autocorrelations no [0,1]-valued h realizes — giving invalid
  (µ ≥ 0.78) or vacuous (µ ≥ 0.0116) bounds (PRO-22's direct-sup_t failure in new clothing). The box
  h ≤ 1 is load-bearing and has no faithful A-space description; re-injecting it returns White's dual
  (capped at C_∞). The only sharp version needs a dim-1 Viazovska eigenfunction with no modular
  analogue for this transform. Memo: `Sign_uncertainty_Poisson_summation_certi.md`.

- **Flag-algebra / continuous-combinatorics certificate on the M(n)→µ limit.** Genuinely novel and
  certificate-grade, but the sup over a *growing* shift family (k ~ αn) needs a 1-parameter outer
  loop plus a flag-enumeration + SDP-assembly pipeline (~1–2 weeks to a first ℓ=5 number; flagmatic
  is graph/permutation-adapted), and flag bounds can stall below the truth. Just out of the top 3 on
  tractability vs the surviving new-cone SDPs — **promote if #1 stalls.**
  (Discussed in `RANKED_SHORTLIST_2026-06-06.md`.)

- **Magic-function / Viazovska transplant** (closed-form dual certificate as a ±1 Fourier
  eigenfunction). Lottery ticket — no modular-symmetry analogue in dim 1 for this transform, a PSLQ
  hit is not a proof; near-free probe (reuse `pslq_hunt.py` on an extracted dual φ\*), keep as a
  probe only.

---

## Selection notes

- **Why these three and not the filter's literal top-3** (which were #1 symmetry, #2 3-point, #3 KKT
  with symmetry ranked first): the task mandates *diversity of technique* and the goal stresses
  attacks that *change the value or give a fundamentally different proof*, beyond the
  convex-optimization lane already mined. Symmetry reduction, by its own honest verdict, **cannot
  change the value** — it is an *enabler*. So the two levers that can actually move the needle
  (3-point = the principled ceiling-breaker; KKT/Krawczyk = the two-sided different-proof) are
  ranked first and second, and symmetry reduction is kept as #3 in its true role: the cheapest,
  highest-confidence near-term move and the *multiplier* that makes #1 solvable at the scale its
  gain appears. The three lanes (higher-arity relaxation / exact-optimality validated numerics /
  representation theory) are maximally distinct — no two are flavors of one idea.
- **All decisive *cheap* checks have been run** (documented per-item above and in the per-approach
  memos); every *heavy solve* is deferred per task scope. No bound is claimed anywhere.
- **The two recurring rigor traps to honor on any follow-through:** (i) the PRO-22 validity gate
  (reconstruct f, confirm f ∈ [0,1] and Ω ≤ true sup_t(f⋆f)) — any new attack that drops the
  cell-envelope realizability link is invalid; (ii) the tail-bound trap (Lasserre/poly-moment) —
  any truncated higher-moment sum must carry an analytic remainder, and gain must be reported *net
  of* the tail, not at finite truncation.
