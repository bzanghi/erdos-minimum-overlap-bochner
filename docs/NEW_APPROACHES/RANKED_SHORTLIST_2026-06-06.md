# Ranked shortlist of NEW attacks on µ (Erdős minimum overlap)

Date: 2026-06-06. Additive curation pass over the per-lens dossiers in this directory
(`NA_optimization_hierarchies.md`, `harmonic_analysis_lens.md`, `probabilistic_info_theory.md`,
`NA_discrete_discrepancy.md`, `computer_assisted_formal_global_optimization.md`,
`CROSS_DOMAIN_WILDCARD_2026-06-06.md`).

Goal bar: push the LOWER bound past the conjectured framework ceiling **C_∞ ≈ 0.380558**
(PRO-6), OR give a fundamentally different proof, OR close the gap from the UPPER side.
Current state: LB ≈ 0.380284 (PRO-21), UB ≤ 0.380871 (Together, shown NOT tight by PRO-23).

## Load-bearing facts established this pass (verified in-repo)
- **C_∞ ≈ 0.380558 is PROVABLY a property of White's cell-envelope + Bochner cone specifically.**
  PRO-6 shows `f_CB(N|n) = f_C(N)` (equality): the joint augmentation is *just* the
  cell-envelope at a higher Bochner level. PRO-22 shows dropping the cell-envelope makes
  the bound INVALID (reconstructed f exits [0,1], true Ω 5× the reported value). ⇒ the ceiling
  is intrinsic to *that* relaxation; a **different cone has a different ceiling**. This is the
  single most important lever for ranking: in-cone approaches are capped, new-cone ones are not.
- **All 19 candidate techniques are genuinely unbuilt.** Grep over `lp_research_state/code/`
  for de-laat/vallentin/three-point/triple/bispectrum/copositive/krawczyk/interval-newton/
  npa/noncommutative/flag-algebra/fisher/stam/entropic/sign-uncertainty/viazovska/second-moment/
  OGP/discrepancy/symmetry-reduction/isotypic/pinn/symbolic-regress returns only NAME COLLISIONS
  (`sweep_T5p.py`, `symmetric_push.py`, `pslq_hunt.py`, `lasserre3_toeplitz.py`). No real code.
- **PRO-23 KKT residual = 7.65e-3 at best** with Together's h*; its optimizer is degenerate
  (437 shifts tie within 1e-9, ~36% of cell-shifts) and is NOT a tight KKT point. So forward
  solving the Euler–Lagrange system attacks a genuinely-open object — but the degeneracy is the
  central risk for any finite active-set reduction.
- **`symmetric_push.py` is CONDITIONAL** (assumes f even ⇒ rows 5,6 only, NOT unconditional);
  a representation-theoretic isotypic reduction would be unconditional and is unbuilt.
- **Entropy was DCP-blocked on `cp.entr(f)` ONLY** (LEVER_B), never applied to the
  autocorrelation A, and Fisher-info via `quad_over_lin` was never attempted.
- **OUT_OF_BOX_INEQUALITIES ran only one-liners** (all collapse to 1/8); never re-ran the
  Martin–O'Bryant/Yu averaging with the L∞ box, never used the entropic UP, never used
  sign-uncertainty.
- **OUT_OF_BOX_CROSS_DOMAIN §2a explicitly flags the 3-point lift as never-built** ("a
  three-point analogue would see configurations the two-point LP cannot") but mis-routed it into
  a non-sharpness *diagnostic*; refs (Cohn–de Laat–Salmon arXiv:2206.09876; three-point bounds
  arXiv:2206.15373) are listed, unexecuted.
- **PRO-29 spectral is 4× loose because it dropped the polytope/L∞ constraints** (Rayleigh over
  the unconstrained L² ball) — precisely the gap an operator-localizer (NPA) would close.

## Dedup performed (near-duplicates merged across lenses)
- **3-point SDP lift × 3 framings → ONE.** (a) de Laat–Vallentin measure/GPM hierarchy +
  Henrion–Rudi tail [optimization-hierarchies], (b) "Three-point SDP via Cohn–de Laat–Salmon
  discrete-reduction" [discrete-discrepancy], (c) "Three-point/k-point hierarchy"
  [computer-assisted]. Same object (triple correlation R₃(s,t) / bispectrum B(m₁,m₂)), same cone,
  same precedent. Merged → **#2**, keeping the strongest validity gate (preserve cell-envelope +
  PRO-22 reconstruction check) and the cleanest rigor path (run finite small-N FIRST to detect
  any gain with NO tail bound; only then add the Henrion–Rudi periodic Hausdorff truncation).
- **KKT-system solve × 2 framings → ONE.** Interval-Krawczyk certifier [HA-3] + PINN/free-boundary
  proposer [ml]. Complementary, not redundant: PINN proposes the active set & optimizer profile,
  Krawczyk certifies the finite F(z)=0 reduction. Merged → **#3** as proposer→certifier pipeline.
- **Analytic dual certificate on A × 2 → ONE.** "Native Cohn–Elkies/Delsarte LP" folded into the
  **sign-uncertainty / Poisson-summation certificate** [#6]; the bare native-LP's own #1 risk is
  being equivalent-to-White, so the sign-uncertainty-geometry version (more clearly a different
  ceiling) is the representative.

---

## TOP 6 (ranked by novelty × promise-toward-bar × tractability)

### 1. Representation-theoretic symmetry reduction of White's SDP (regular *-representation block-diagonalization)
**Why #1:** Highest-confidence, highest-leverage, genuinely-unbuilt move. Pure linear algebra
(de Klerk–Pasechnik–Schrijver 2007 block-diagonalizer / Schur-decompose the averaging projector
P_G = (1/|G|)Σ ρ(g)) — no new theory, near-term. Attacks the *named* barrier (large-T_max
intractability) head-on: it doesn't change the relaxation, only the cost of solving it, so F3's
projected +4–6e-4 at bochner_n≥40–80 (flagged "research not engineering") becomes an engineering
computation. UNLIKE `symmetric_push.py` (conditional even-f, rows 5,6 only) this exploits symmetry
of the CONSTRAINTS ⇒ unconditional, all 7 rows. First object: form G=⟨x↦2−x⟩ action on the
cosine/sine + Bochner index sets at N=300, bn=12; verify reduced SDP reproduces the unreduced dual
to 10+ digits (repo cross-check standard), measure flop/memory ratio, then rerun binding row 4 at
bn=60–80. **Honest cap:** if the only symmetry is Z/2 the split is ~2×; the payoff hinges on a
larger effective (dihedral/translation) group in the cell-index set — must be checked. Even a
perfect 10× only buys ~0.38058 (approaches, can't pass C_∞) — it is the **enabler that makes #2
solvable at useful scale**, which is what transcends the ceiling.
Refs: de Klerk–Pasechnik–Schrijver Math.Program. 109 (2007); Bachoc–Vallentin JAMS (2008);
Vallentin arXiv:0706.4233; Gatermann–Parrilo JPAA (2004); de Klerk Math.Program. (2011).

### 2. Three-point / k-point (triple-correlation / bispectrum) SDP lift  [MERGED]
**Why #2:** The single most principled candidate to *pass* C_∞. The ceiling is *provably* a
2-point phenomenon (PRO-6/PRO-22 above); a 3-point relaxation lives in a strictly larger cone with
a different ceiling. The sphere-packing precedent is direct and quantitative: the
Cohn–Elkies-LP → Bachoc–Vallentin-SDP jump is exactly what broke the 2-point LP plateau and proved
non-sharpness in d=12,16. M(k) is intrinsically a correlation count, so triple correlations
R₃_h(s,t)=∫h(x)h(x+s)h(x+t)dx (Fourier: bispectrum B(m₁,m₂)) encode *which inner shifts can
simultaneously be low* — invisible to any degree of the 2-point program (this is the axis the
RETRACTED Lasserre-2 did NOT move: it lifted the *degree* of the same pairwise object, not the
number of points). First object: at N=300, R=5, coarse inner-shift grid |{s}|≤8, build the
symmetrized 3-point Gram/PSD block, BOLT IT ONTO the existing cell-envelope cone (do NOT replace it
— PRO-22), solve row-4 center in cvxpy. Decision gates: (i) does the dual rise above the 2-point
value at same N,R? (ii) is the 3-point block ACTIVE (nonzero dual)? MANDATORY validity gate:
reconstruct f, verify reported Ω ≤ true sup_t(f⋆f). **Run the FINITE small-N version first — it
needs NO tail bound, so it detects ANY gain before re-entering the retracted-Lasserre tail-bound
trap;** only if gain is confirmed, add the Henrion–Rudi periodic Hausdorff truncation
(arXiv:2401.07734 Prop. 7) for rigor, then scale with #1's symmetry reduction (the 3-point block
carries dihedral shift symmetry → block-diagonalizes). **Risk:** the triple cut may be IMPLIED by
2-point+cell-envelope for this pairwise objective (inactive block, zero gain — the "is it active?"
probe kills this fast); O(T²)×O(T²) blocks are heavy without #1.
Refs: Bachoc–Vallentin JAMS (2008)/arXiv:math/0608426; de Laat–Vallentin Math.Program. (2015)
arXiv:1311.3789; Cohn–de Laat–Salmon arXiv:2206.09876; three-point bounds arXiv:2206.15373;
Henrion/Korda–Rudi arXiv:2401.07734; ClusteredLowRankSolver arXiv:2202.12077 / SDPB arXiv:1502.02033.

### 3. Solve the KKT/Euler–Lagrange system: PINN proposer → interval-Krawczyk certifier  [MERGED]
**Why #3:** Only route that could *close the gap two-sidedly* rather than nibble the LB — it
sidesteps the relaxation framework entirely by solving the EXACT optimality condition (PRO-23's
free-boundary EL system: Σ_{t∈S} γ_t[h*(x+t)+h*(x−t)] = κ on the interior, with bang-bang plateaus
and S=argmax_t M). A validated enclosure of its unique root is a rigorous enclosure of the TRUE µ.
Rechnitzer's 128-digit ball-arithmetic solve of the sibling L² autoconvolution constant
(arXiv:2602.07292, 2026) is direct evidence of feasibility at extreme precision for this family.
PRO-23 only *evaluated* Together's h* against the system (residual 7.6e-3, degenerate) — it never
solved it forward. Pipeline: (a) two-network PINN (sigmoid head for the [0,1] box; learnable
free-boundaries + dual atoms γ) trained on the coupled residual to PROPOSE a structurally-resolved
h* and its active set S; max_t M(h_θ) gives a true UB (PRO-23 proved Together's is slack);
(b) seed a float-Newton then one interval-Krawczyk step on the finite reduction F(z)=0 in z=(break-
points, γ, {t_i}, κ, µ); K⊂int(Z) ⇒ unique root ⇒ rigorous two-sided µ enclosure. **Cheap
intermediate win:** a float-Newton KKT-tight point with residual ≪1e-9 (vs PRO-23's 7.6e-3) is the
first tight KKT point, immediately tightens µ from above and unblocks PRO-23 Step 4 — *and feeds S
to #2/the dual-γ route.* **Risk (shared by both halves):** the optimizer's active set may be
genuinely large/degenerate/Cantor-like (437 near-ties) ⇒ finite F(z)=0 reduction wrong, Krawczyk
has no isolated root; Jacobian ill-conditioned near the degenerate optimum inflates radii. Mitigate
by starting from the smallest symmetric m=2–4 reduction (cheap to falsify) and by training the
dual γ first (bang-bang gradients vanish on plateaus).
Refs: Moore (1966); Rump INTLAB / Acta Numerica (2010); Rechnitzer arXiv:2602.07292 (2026);
van den Berg–Lessard (rigorous numerics in dynamics); obstacle-PINN arXiv:2304.03552; repo PRO-23
(`docs/archive/LEVER_FUNCTIONAL_EQUATION.md`). Tooling note: no interval/Krawczyk lib in repo
(python-flint/arb not installed; only mpmath/sympy) — hand-implement one Krawczyk step in mpmath.

### 4. NPA / Helton–McCullough noncommutative-moment certificate for sup_t ⟨h, T_t h⟩
**Why #4:** Changes the certified object from "Fourier coefficients of continuous f" to
"noncommutative moments of the shift ACTION" — whose relaxation gap is unrelated to White's C_∞.
Two repo-specific reasons it could beat the saturation: (i) it inserts the L∞ bound LOSSLESSLY as
an OPERATOR inequality (localizing matrices [h]Γ⪰0, [1−h]Γ⪰0, ⟨h⟩=1) — exactly the polytope
restriction whose omission made PRO-29's spectral bound 4× loose and killed every B-S/M-R
autocorrelation transfer; (ii) it truncates by WORD LENGTH, which is finite-dimensional and
TAIL-FREE — it structurally dodges the recurring truncated-sum rigor trap that retracted Lasserre-2
and the poly-moment cuts. Genuinely noncommutative (T_t and mult-by-h don't commute), so it is NOT
a re-tread of any commutative f-side moment work nor of the naive Rayleigh spectral entry. First
object: discretize the shift group to cyclic Z_m (m=8–16), build level-1 then level-2 NPA moment
matrices over words in {T_1..T_{m−1}, h} (~80×80 at level-2, m=8 — trivial in cvxpy), impose
T-unitarity as moment equalities + the L∞ localizers, solve min⟨h T_k h⟩ jointly over the dual
shift-measure; check whether the LB clears 0.30 and RISES with m and level. **Risk:** the
continuous shift group t∈ℝ/2ℤ → discretization to Z_m needs a Bochner/Naimark extension argument
(a tail-flavored step, though here it controls a modulus-1 unitary character, much cleaner than an
f² Fourier tail); whether the L∞ localizers actually recover polytope tightness at level-2 is
unknown a priori. Needs an ncpol2sdpa/TSSOS-style encoding.
Refs: Navascués–Pironio–Acín NJP (2008); Pironio–Navascués–Acín SIAM J.Optim. (2010);
Helton–McCullough Positivstellensatz; arXiv:2402.02126 (2024); arXiv:2510.08427 (2025).

### 5. Fisher-information / log-Sobolev convex surrogate for the L∞ functional
**Why #5:** Cheapest genuinely-new analytic lane, and it resolves the EXACT blocker that killed
LEVER_B: ∫(A')²/A is the perspective of (u,v)↦u²/v, which IS jointly convex and IS DCP-encodable
via `cp.quad_over_lin` — UNLIKE the `cp.entr(f)` that LEVER_B correctly rejected (concave ≤ const).
It bounds *flatness* through a differential/gradient functional the Bochner moment hierarchy is
structurally blind to (the documented Lever-H blocker: sup A has no polynomial expansion in f̂),
and its convolution superadditivity (Stam/Blachman, A = f⋆f̃ a self-convolution) is a genuinely
nonlinear structural law with no finite-moment proxy. Every link (de Bruijn, Stam, LSI,
Gagliardo–Nirenberg) has a known sharp constant ⇒ best chance of a clean closed form. First object
(minutes): (i) on Together h* compute J(A)=∫(A')²/A and calibrate c in ‖A‖_∞ ≥ c·A(0)²·J(A)^(−1/2),
see how tight at the optimizer; (ii) solve the small convex program min_f c·A(0)²·J(A)^(−1/2) in
cvxpy via quad_over_lin at N=500, read the value. >0.293 is informative, →0.380 is a contender.
**Risk:** the sup-norm link ‖A‖_∞ ≥ c·‖A‖₂²/√J(A) is the least-standard step; the box-optimizer's
flat top → small A' on the plateau could soften the floor (mitigant: its steep shoulders make J(A)
large exactly where the bound bites). Experiment (i) settles slack before any theorem.
Refs: Stam Inf.&Control (1959); Blachman IEEE-IT (1965); Carlen JFA 101 (1991); Gross AJM (1975);
Costa IEEE-IT (1985); arXiv:1608.05431.

### 6. Sign-uncertainty / Poisson-summation certificate ON the autocorrelation A  [absorbs native Cohn–Elkies LP]
**Why #6:** Builds an ANALYTIC DUAL certificate directly on the continuous, globally
positive-definite A (Â_h = |f̂|² ≥ 0, even, mass 1, supp ⊆ [−2,2]) with NO cell-envelope — so no
validity hole (unlike the invalid direct sup_t SDP) and no finite-Fourier truncation ceiling (its
ceiling is certificate quality, a different limit than Bochner saturation). It treats min-overlap as
a +1-type sign-uncertainty instance and uses the sign-uncertainty refinements (modular-form-built g)
that have repeatedly hit SHARP constants in sibling problems and that beat the naive Logan first-zero
bound (which collapses to 1/8). Distinct from the ruled-out Beurling–Selberg majorants (those
majorize an INDICATOR; this bounds concentration of a sign-CONSTRAINED pair) and from PRO-29's naive
Rayleigh. First object (1–2 days): finite certificate LP — discretize an even g on [−2,2] in a
degree-~40 cosine basis, enforce ĝ ≤ 0 on the band and g(t) ≤ 0 for |t| ≤ t₀, maximize the ratio
that lower-bounds sup A; sweep t₀ (effective min-shift), read the floor; cross-check vs A_+(1) under
the GOSS normalization; use Poisson summation over the lattice 2Z (matching support) to upgrade to
an EXACT Cohn–Elkies-style inequality; evaluate at Together's h* to confirm no contradiction with
0.380871. **Risk:** the certificate LP could turn out to be a dual description of the SAME cone as
the SDP (no gain — but proving that equivalence is itself a worthwhile rigorous no-go); the
near-sharp modular-form g may need a Viazovska-style ansatz not available off-the-shelf in dim 1.
Medium risk, high informational value either way.
Refs: Bourgain–Clozel–Kahane (Ann.Inst.Fourier 2010); Cohn–Gonçalves (Invent.Math. 2019);
Gonçalves–Oliveira e Silva–Steinerberger arXiv:2003.10771; Carneiro–Quesada-Herrera arXiv:2006.00959;
Cohn–de Laat–Salmon arXiv:2206.09876; Cohn–Elkies Ann.Math. 157 (2003) arXiv:math/0110009.

---

## Dropped (one-line reasons)
- **Completely-positive / set-copositive reformulation of discrete M(n)** — author-rated lowest
  promise; CP hierarchies converge notoriously slowly, the n→∞ limit of a per-n CP bound has no
  uniform control, and max-over-shifts needs an extra epigraph lift. Logged for lens-completeness.
- **Discrepancy determinant/eigenvalue (Lovász–Spencer–Vesztergombi) on the shift-incidence matrix**
  — LSV/eigenvalue bounds are tight only for Hadamard-type, generically loose for Toeplitz shift
  systems; ℓ∞ (max_k) → ℓ²/det conversion loses √(#shifts); most-likely plateaus at 0.25–0.29.
- **Second-moment / OGP / Gamarnik–Kızıldağ on M(n)** — natively bounds the TYPICAL value, not the
  worst-case minimum; overlaps X_k share the same partition ⇒ effective # independent shifts likely
  O(1), collapsing the large-deviation gain to µ≥1/4 (1955-era). Highest-variance; the
  Cov(X_k,X_{k'}) probe is cheap but the prior is poor.
- **Symbolic regression on the optimizer structure (FunSearch/PySR)** — author-rated
  highest-risk/lottery; depends on a higher-resolution optimizer from #3 first; the value-side PSLQ
  is already NEGATIVE to 50 digits (weak evidence the function is ugly too); NP-hard SR may return
  fitting noise. A downstream tool, not a primary attack — revisit only if #3 yields a clean h*.
- **Native Cohn–Elkies/Delsarte LP (witness=proof)** — MERGED into #6; its own #1 risk is being
  equivalent to a weakening of White's program (Cohn–Elkies LPs are famously non-sharp off magic
  dimensions), so the sign-uncertainty-geometry version is the stronger representative.
- **Magic-function / Viazovska transplant (closed-form dual certificate as ±1 Fourier eigenfunction)**
  — author-rated lottery ticket; no modular-symmetry analogue in dim 1 for this transform, a PSLQ
  hit is not a proof, and closing it needs Romik-level modular-form-inequality work that may not
  exist. First experiment is near-free (reuses pslq_hunt.py on an extracted dual φ*) — keep as an
  honorable-mention probe, not a top-6 commitment.
- **Flag-algebra / continuous-combinatorics certificate on the M(n)→µ limit** — genuinely novel and
  certificate-grade, but the sup over a GROWING shift family (k~αn) needs a 1-parameter outer loop
  plus a flag-enumeration+SDP-assembly pipeline (~1–2 weeks to a first ℓ=5 number; flagmatic is
  graph/permutation-adapted), and flag bounds can stall below the truth. Just out of the top 6 on
  tractability vs the surviving new-cone SDPs; promote if #2 stalls.
