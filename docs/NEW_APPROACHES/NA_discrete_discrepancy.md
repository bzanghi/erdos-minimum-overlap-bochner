# NEW APPROACHES — Discrete / Additive-Combinatorics / Discrepancy lens

**Author lens:** discrete-discrepancy. **Date:** 2026-06-06.
**Status:** PROPOSALS ONLY — not results. Each vetted against the ruled-out ledger.

Setup recap (for self-containment):
- Discrete: partition `[2n] = A ⊔ B`, `|A|=|B|=n`; `M(k) = |A ∩ (B+k)|`; `M(n) = min_partition max_k M(k)`; `µ = lim M(n)/n`.
- Continuous (White): `µ = inf_{h:[0,2]→[0,1], ∫h=1} sup_t R_h(t)` where the overlap functional reduces (PRO-29 duality) to `µ = 1 − sup_h inf_t ⟨h, T_t h⟩`, the **two-point autocorrelation** `⟨h,T_t h⟩ = ∫ h(x)h(x+t)dx`.
- Current bracket: `µ ≥ 0.380284` (augmented two-point SDP) … `µ ≤ 0.380871` (Together). Framework ceiling conjectured ≈ `0.380558`.

The unifying observation across all three proposals: **every lever the project has pulled is a two-point object** (Bochner PSD on `f̂(m)`, M-side Bochner, Lasserre on `(f²)̂`, ellipse cover). The discrete problem `M(k)` is literally a two-point correlation count. The discrepancy lens says the *two-point* relaxation is provably non-sharp for problems of this type (Cohn–de Laat–Salmon proved exactly this for sphere packing), and the route past the ceiling is to add **three-point / triple-correlation** structure or to transfer a **discrete combinatorial floor** that the continuous two-point program cannot see. That is the through-line.

---

## Approach D1 — Three-point (triple-correlation) SDP lift via the Cohn–de Laat–Salmon discrete-reduction machinery

### Technique / object to compute
Replace White's two-point program with a **three-point relaxation**. The decision object is upgraded from the autocorrelation `R_h(t)=⟨h,T_t h⟩` (a function of one shift) to the **triple correlation**
`R_h(s,t) = ∫ h(x) h(x+s) h(x+t) dx`,
a function of two shifts. The three-point SDP imposes that the matrix of triple-correlation moments, symmetrized over the shift group, is PSD — exactly the Cohn–de Laat–Salmon construction that "symmetrizes a constraint on `H`-invariant pair distributions into a constraint on `G`-invariant triple distributions." Concretely:

1. Build the three-point moment tensor `T[(i,j)] = (\widehat{h·T_{s_i}h})(m_j)` over a finite grid of inner shifts `{s_i}` and Fourier modes `{m_j}`, using the **same cell-envelope cosine/sine representation** the repo already trusts for the two-point coefficients (so the load-bearing cell-envelope rigor — the thing PRO-22's "direct sup_t SDP" violated — is preserved).
2. Impose the symmetrized PSD blocks. This is a **clustered low-rank SDP**: the symmetry over the dihedral shift action factors the big PSD constraint into many small low-rank blocks.
3. Solve it at high precision with **ClusteredLowRankSolver.jl / SDPB** (de Laat–Leijenhorst, arXiv:2202.12077) — the same solver family that powers the three-point kissing-number bounds and is built for exactly this block structure. Read off the dual objective as the rigorous LB (same `rigorous_dual_LB` epistemics the repo already uses with CLARABEL).

### Citation
- H. Cohn, D. de Laat, A. Salmon, *Three-point bounds for sphere packing*, arXiv:2206.15373 (the discrete-reduction / symmetrization template).
- N. Leijenhorst, D. de Laat, *Solving clustered low-rank semidefinite programs arising from polynomial optimization*, arXiv:2202.12077; solver `ClusteredLowRankSolver.jl`; SDPB lineage (Simmons-Duffin arXiv:1502.02033).
- H. Cohn, D. de Laat, et al., *non-sharpness of the two-point LP/Cohn–Elkies bound in d=12,16* (arXiv:2206.09876) — the proof that two-point relaxations of autocorrelation problems are generically *not* tight, which is the motivation.
- Rational-certificate finishing for these solves: arXiv:2603.19901 (SDP bounds on quantum codes: rational certificates) — same de Laat-school certification pipeline.

### Why previously untested here (vs ledger)
The ledger rules out "White's §5 LP/SDP base; Bochner moment-matrix PSD on f,1−f; polynomial-moment Hausdorff cuts; Lasserre level-2; direct sup_t SDP." **Every one of those is a two-point object.** I verified this in `lp_research_state/code/white_full_convex.py:99–180`: the decision variables are `c,d` (cosine/sine Fourier coeffs of `f`) and the constraints act on `f̂(m)` and pairwise products — there is no triple-correlation variable, no symmetrized three-point block, anywhere in the program. The "direct sup_t SDP" that gave an INVALID bound (PRO-22) was an attempt to drop the cell-envelope, **not** a three-point lift — different failure mode. Cohn–de Laat–Salmon three-point machinery has *never* been applied to µ (PRO-25's lit mine surveyed it only as a non-sharpness *diagnostic* idea in OUT_OF_BOX 2a, and noted "a three-point analogue would see configurations the two-point LP cannot" — but no one built it). This is the single biggest unexplored structural lever.

### First experiment (cheap)
Before the full Julia/SDPB build: a **finite-n discrete sanity check** in Python that the three-point relaxation is strictly stronger than two-point. Take the exact `M(n)` optimizers already stored (`lp_research_state/data/Mn_optimizers_large.json`, n≤20) and compute, for the *continuous relaxation at small Fourier truncation T≈6 and a coarse inner-shift grid {s} of size ≤8*, both (a) the two-point dual bound and (b) the two-point + symmetrized-three-point dual bound, in cvxpy/CLARABEL. If (b) > (a) by more than solver noise at this toy scale, the lift has signal and justifies the SDPB port. This is a ~1-day cvxpy experiment reusing `build_problem`'s coefficient helpers; no new solver needed yet.

### Why it could beat the saturation
The framework ceiling (~0.380558) is a *two-point* ceiling: PRO-6's complementarity/KKT saturation proof is entirely about the two-point dual reaching its complementary-slackness fixed point. A three-point relaxation lives in a strictly larger constraint cone, so its ceiling is a *different* number — and Cohn–de Laat–Salmon *proved* the analogous jump is real and quantitative for sphere packing (two-point Cohn–Elkies is non-sharp in d=12,16; three-point sees the gap). The min-overlap M(k) is intrinsically a correlation count, so triple correlations encode genuinely new combinatorial information (which inner shifts can be simultaneously low) that the two-point program is blind to.

### Risk
- **Cost/conditioning:** three-point SDPs are expensive; the inner-shift grid blows up the variable count, and the cell-envelope residual analysis (Step E) must be re-derived for triple products — non-trivial new rigor work. The repo's whole rigor edifice (`rigorous_dual_LB`, cell envelopes) would need a triple-correlation analogue, and getting that wrong is exactly the historical overclaim trap (cf. Lasserre tail bound).
- **Possible null:** the three-point gain for *this* problem might be tiny even if nonzero (the toy experiment is the gate). If the first-experiment Δ is at noise level, abandon before the SDPB port.
- The symmetrization group for the line `[0,2]` is just reflection/translation (smaller than the orthogonal group in sphere packing), so the three-point gain could be smaller here than in the packing analogue.

---

## Approach D2 — Discrete-to-continuous transference floor via Kravitz's generalized-difference-set ↔ autocorrelation duality (manufacture a provable µ-floor from a finite combinatorial covering bound)

### Technique / object to compute
Kravitz (arXiv:2004.06611) proves an **exact duality** between a *discrete* covering quantity (minimum size of a "generalized difference set" covering `[N]`) and a *continuous* constant defined by "nonnegative-valued functions on ℝ with autocorrelation integral bounded below on `[0,1]`." That is structurally **the same shape as µ** (`µ = 1 − sup_h inf_t ⟨h,T_t h⟩` is exactly an autocorrelation-floor optimum over nonnegative `h`). The proposal is to **adapt his discrete↔continuous transference machinery to µ's own function class** (compact support `[0,2]`, `∫h=1`, **`h≤1` pointwise** — the constraint that PRO-32 showed kills the Barnard–Steinerberger transfer) to obtain:

> a **finite combinatorial covering/packing inequality** at level `n` whose value is a *rigorous lower bound on µ* — bypassing the SDP entirely.

The exact object: define the discrete dual of µ — a generalized-difference-set-style covering number `D(n)` for partitions of `[2n]` with the `h≤1` (indicator) constraint baked in — and prove, via Kravitz's step-function/limiting argument run **in reverse**, that `µ ≥ g(D(n))` for an explicit monotone `g`. Then compute `D(n)` exactly by ILP/SAT at the *small* `n` that is tractable, and let the inequality propagate to the limit.

### Citation
- N. Kravitz, *Generalized difference sets and autocorrelation integrals*, arXiv:2004.06611 — the discrete↔continuous *duality theorem* and its step-function transference proof (abstract: "the optimal constant in an analogous problem concerning nonnegative-valued functions on ℝ with autocorrelation integral bounded below on `[0,1]`").
- Cilleruelo–Ruzsa–Vinuesa and Martin–O'Bryant (the generalized-Sidon ↔ continuous-autoconvolution dictionary that Kravitz extends).
- Context: J. Cilleruelo, *Generalized Sidon sets* (the `B_2[g]` density side).

### Why previously untested here (vs ledger)
PRO-25 ("μ is not a renamed constant") surveyed Kravitz's τ and the generalized-difference-set world **only to check for a numerical collision** — found the constant ≈1.56 ≠ 0.38 and stopped. PRO-32 killed the *Barnard–Steinerberger* transfer because their function class lacks the `h≤1` ceiling. **Neither examined Kravitz's transference *method* as a tool to build µ's own discrete dual.** The ledger's "transference" entries are all about *importing an existing constant* (rejected) — this proposal instead **runs the transference construction natively for µ**, with the `h≤1` constraint included from the start (precisely the constraint whose omission broke B-S). The discrete dual `D(n)` of µ-with-ceiling has, to my reading of the archive, never been defined or computed. This is not the "ILP/SAT for exact M(n)" entry (ruled out as uninformative): that computes the *primal* `M(n)` whose ratio sits at 0.40 for small n; `D(n)` is a *dual covering* quantity that lower-bounds µ directly and need not have the 1/n rounding overhang that makes `M(n)/n` useless below n≈50.

### First experiment (cheap)
1. Read Kravitz §2–3 carefully and write down his duality theorem with the **exact** function class and the step-function map `(finite set) ↦ (indicator/step function)`.
2. On paper, specialize his construction to `h:[0,2]→[0,1], ∫h=1` and check whether the `h≤1` ceiling makes the discrete dual a *covering* problem (each integer hit ≤ once per shift) — if yes, that covering number is computable by the SAT harness already in `lp_research_state/code/_sat_Mn.py`.
3. Compute `D(n)` for n = 8…20 with that harness (reusing the existing pseudo-boolean encoding), and tabulate `g(D(n))`. If `g(D(n))` exceeds 0.380284 at any tractable n, that is an immediate rigorous improvement; if it converges to something `<0.3803` it bounds the *method's* reach (still a publishable structural result).

### Why it could beat the saturation
The SDP ceiling is a property of the *continuous Fourier relaxation*. A transference floor is a *combinatorial* lower bound that does not pass through the Fourier-truncation cone at all — so it is not subject to PRO-6's two-point complementarity ceiling. Kravitz's theorem is an *equality*, meaning the discrete and continuous optima genuinely coincide in the limit; if µ's ceiling-constrained dual `D(n)` is exactly computable, the limit *is* µ (or a clean lower bound on it), not a relaxation of it. Even a partial transference (a one-sided `µ ≥ g(D(n))` with explicit `g`) sidesteps the saturation because it is a different proof object.

### Risk
- **The duality may not survive the `h≤1` ceiling.** Kravitz's clean equality is for autocorrelation *bounded below on [0,1]* without an L∞ cap; adding `h≤1` may turn the equality into a one-sided inequality or destroy the discrete-dual interpretation. This is the same wall (in mirror image) that PRO-32 hit — honest risk it is fatal.
- **`D(n)` may inherit the same 1/n overhang** as `M(n)/n` (the reason small-n integer data is uninformative); if so, tractable n won't reach 0.3803.
- Substantial new theorem-writing, not a code probe — weeks of math, and the payoff might be only a *characterization* (µ = a covering limit) rather than a numeric improvement.

---

## Approach D3 — Eigenvalue/determinant discrepancy lower bound on the *finite-n* overlap matrix, certified at finite n and pushed to the limit

### Technique / object to compute
Treat `M(n)` as a genuine **two-coloring discrepancy** of the shift-incidence system and apply spectral/determinant discrepancy lower bounds (the "eigenvalue method" of combinatorial discrepancy, plus the linear-algebra/determinant bound of Lovász–Spencer–Vesztergombi and its modern strengthenings). Concretely, for fixed `n` form the `(2n−1) × 2n` **shift-incidence structure**: encode the coloring `χ ∈ {+1,−1}^{2n}` (the `A/B` partition, with `Σχ=0`) and the family of shifted overlap functionals `M(k)`. The overlap `M(k)` is a *bilinear* form `M(k) = ¼(n − χᵀ S_k χ) + (linear)`, where `S_k` is the (non-periodic) shift-by-k operator. Then:

1. Compute the **eigenvalue discrepancy bound**: `max_k |M(k) − E M(k)| ≥ c·σ_min`-type lower bounds from the spectrum of the Gram matrix of `{S_k}`. PRO-29 tried the *naive* Rayleigh bound on a single `A_j` and got 4× loose **because it used the unconstrained L² ball**; the fix is the **determinant/eigenvalue discrepancy bound applied to the `{±1}` cube with the `Σχ=0` linear constraint**, which is a fundamentally tighter inequality (it is sharp for Hadamard-like systems).
2. Compute the **Lovász–Spencer–Vesztergombi determinant lower bound** `disc ≥ |det(B)|^{1/m}/(2·max row norm)` for the best square submatrix `B` of the shift-incidence matrix — this gives a *certified* per-n floor on `max_k M(k)` over all colorings, hence on `M(n)`.
3. Track the floor as `n→∞` and compare its growth rate to `µn`.

### Citation
- L. Lovász, J. Spencer, K. Vesztergombi, *Discrepancy of set-systems and matrices*, Europ. J. Combin. 7 (1986) — the determinant lower bound `disc(A) ≥ det/...`.
- J. Spencer, *Six standard deviations suffice* (Trans. AMS 1985) and the eigenvalue / "linear algebra" discrepancy lower bounds (Matoušek, *Geometric Discrepancy*, ch. 4 — the eigenvalue method).
- Modern: the matrix-Chernoff / Banaszczyk lineage and Larsen's *Constructive discrepancy minimization* (arXiv:1711.02860) for the algorithmic eigenvalue bounds; the spectral lower-bound technique in additive shift systems (cf. the shift-operator polynomial method, arXiv:2311.08873, for the algebra of `S_k`).

### Why previously untested here (vs ledger)
The ledger rules out "spectral/translation-operator min-max reformulation (naive bound 4× loose)" — that is **exactly** the naive Rayleigh-quotient attempt (PRO-29), which bounded a *single* shift operator over the *unconstrained unit ball*. The discrepancy-theoretic lower bounds (determinant bound; eigenvalue method *with the `{±1}`-cube and zero-sum constraint*) are a **different inequality family**: they are lower bounds on the *coloring* discrepancy of the *whole shift system simultaneously*, sharp for structured systems, and they natively respect the integrality `χ∈{±1}` and the constraint `Σχ=0`. PRO-29 explicitly noted its looseness came from ignoring `‖h‖_∞≤1` and `‖h‖_1=1`; the determinant/eigenvalue discrepancy bounds are built precisely for the integral constrained problem. No archive memo computes a determinant lower bound or the constrained eigenvalue discrepancy bound. (It also is *not* the ruled-out "ILP/SAT for exact M(n)" — that computes the exact value by search; this computes a *certified floor* by linear algebra, cheaply, and is informative as a *rate*, not a single value.)

### First experiment (cheap)
For n = 10…30, build the shift-incidence matrix of `{S_k}` (dense, ≤ 60×60), and compute (a) the LSV determinant bound over a greedy-selected square submatrix, and (b) the constrained eigenvalue bound `min over zero-sum χ` via the second-smallest eigenvalue of the appropriate quadratic form. This is pure numpy/scipy, an afternoon. Plot the certified floor `/n` vs n. The decision criterion: does the floor/n trend *toward* something ≥ 0.3803 (win or framework-independent confirmation), or does it plateau below (then the determinant bound is too weak for this system and we stop)? Cross-check the floor against the known exact `M(n)` (must satisfy floor ≤ M(n)).

### Why it could beat the saturation
This is a **finite-n, integral** lower bound that never enters the continuous Fourier relaxation, so PRO-6's two-point asymptotic ceiling (a statement about the *continuous SDP dual*) does not bind it. Discrepancy lower bounds are *exact* for highly structured incidence systems (Hadamard, arithmetic progressions), and the shift system `{S_k}` of the min-overlap problem is extremely structured (Toeplitz/circulant-adjacent). If the determinant of a well-chosen shift submatrix is large, it forces *every* coloring to have a large overlap at some shift — a direct combinatorial floor on `M(n)` that is *constructive* and certifiable in exact arithmetic (SDPA-GMP not even needed; integer determinants). It attacks the *discrete* `µ = lim M(n)/n` directly rather than its continuous relaxation.

### Risk
- **The determinant bound may be quantitatively weak.** LSV/eigenvalue bounds are tight for Hadamard-type systems but can be far from sharp for general Toeplitz shift systems; the realized floor/n might sit well below 0.38 (e.g., near the old 0.25–0.29 classical bounds), giving no improvement. This is the most likely outcome and the first experiment is explicitly designed to detect it within an afternoon.
- The min-overlap functional is `max_k`, an `ℓ∞` discrepancy; converting the `ℓ2`/determinant bound into a sharp `ℓ∞` floor loses a `√(#shifts)` factor that could erase the gain.
- The limit `n→∞` of a finite-n determinant floor needs its own (possibly hard) asymptotic analysis; a good finite-n floor may not extrapolate to a good µ-floor.

---

## Cross-cutting note (lens summary)

The discrete-discrepancy lens converges on one diagnosis: **the project has exhausted the two-point relaxation and is now bumping its two-point ceiling.** The three escape routes are (D1) go to *three-point* correlations where the analogous sharpness gap is *proven* to exist (Cohn–de Laat–Salmon), (D2) *leave the relaxation entirely* via Kravitz-style discrete↔continuous transference to get a combinatorial floor, or (D3) bound the *discrete integral* problem directly with discrepancy-theoretic determinant/eigenvalue inequalities that the naive PRO-29 spectral attempt never reached. D1 is the highest-expected-value (real SDP machinery exists, gap proven elsewhere) but the most engineering; D3 is the cheapest to falsify (an afternoon) and the most likely to be a fast null; D2 is the most likely to yield a *fundamentally different proof* if the duality survives the `h≤1` ceiling. None re-treads the ledger: all prior spectral/transference/SAT entries were either two-point, unconstrained, primal-value, or numerical-collision checks.

A natural **flag-algebra framing** (Razborov; FlagAlgebraToolbox arXiv:2601.06590; Kiem–Pokutta–Spiegel 2024) sits underneath D1: the symmetrized three-point block *is* a flag-algebra SDP on the shift-density limit object. If D1's hand-built three-point lift shows signal, re-expressing it in flag-algebra form (densities of 3-point shift-configurations in the limit coloring) would give a principled, automatable hierarchy — but that is a follow-on, not a first experiment.
