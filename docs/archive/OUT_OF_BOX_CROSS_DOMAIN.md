# Cross-Domain Analogy Hunt: Erdős Minimum-Overlap Problem

**Date:** 2026-05-10
**Author:** out-of-box research scan
**Scope:** Structural connections to other extremal mathematics whose machinery might apply to pushing the rigorous lower bound on `µ` past `0.379544`.

The problem in functional form:
```
µ = inf_{h: [0,2]→[0,1], ∫h=1} sup_k ∫ h(x)(1 − h(x+k)) dx
  = 1 − sup_h inf_k R_h(k),     R_h(k) = ∫ h(x) h(x+k) dx
```
We want a `[0,1]`-valued density on `[0,2]` with mass 1 whose **autocorrelation
is as flat as possible** (large infimum over shifts). Current rigorous interval:
`0.379544 ≤ µ ≤ 0.380871` (gap ~1.3×10⁻³). White's program is a finite-truncation
Fourier-analytic LP/SDP; the proven ceiling for the Bochner-PSD+ellipse-extension
family at currently-tractable SDP scale is exactly our `0.379544`.

---

## 1. TOP-1 RECOMMENDED BRIDGE — Autoconvolution / autocorrelation extremal theory (Cilleruelo–Ruzsa–Vinuesa / Martin–O'Bryant / Matolcsi–Vinuesa / Cloninger–Steinerberger), and in particular **Ethan P. White's own L² autoconvolution paper** (arXiv:2210.16437).

**Why this is the top-1.** This is not an analogy — it is the same author using
genuinely adjacent technique on a sibling extremal problem. The autoconvolution
problem `inf_f ‖f∗f‖_∞ / ‖f‖_1²` and the related `inf ‖f∗f‖_2` are extremal
problems for the **autocorrelation of a nonnegative density**, structurally
identical to the Erdős minimum-overlap problem (a function plus its reflection).
White (Canad. Math. Bull., arXiv:2210.16437) gave an **"almost-tight" L² result
determining the infimum to 0.0014%** and **proved uniqueness of the minimizer**
— a level of analytical control we do not yet have for `µ`. The follow-up arXiv
2508.02803 (2025) and 2506.16750 (2025) push the autoconvolution L∞ bound to
0.94136 of the conjectured value using **2,399-interval step functions plus a
4× upsampling**, conceptually the same Together-style sequential refinement that
gives our upper bound on `µ`, but in the analytical-bridge direction.

**Concrete proposed approach.**
1. Read arXiv:2210.16437 in detail and identify what White's L² proof exploits
   that our minimum-overlap SDP does *not* exploit. Two specific candidates:
   (a) The L² problem reduces to an explicit eigenvalue problem for a compact
   operator — the minimizer is the principal eigenfunction. If we can reformulate
   the inner `inf_k R_h(k)` via a *measure on shifts k* (Lagrangian dual) and
   take a related eigenvalue point of view, we may get an analytical handle
   beyond Bochner-PSD.
   (b) Uniqueness via convexity of an associated functional — if the minimum-
   overlap optimizer is unique, computing high-precision moments of the
   Together upper-bound step-function should converge to those of the analytical
   minimizer and let us close the gap from the construction side via certified
   smoothing.
2. Implement an SDP for the **L¹-autocorrelation extremal problem** in our
   framework as a sanity check: if our Bochner-augmented program reproduces
   White's L² constant to 10⁻⁵, then the encoding is sound and the technology
   should transfer.
3. Pose to White directly (`communications/email_to_ethan_white.md` already
   exists): can the L² techniques (perhaps via a smoothed `L²+ε‖·‖_∞` proxy)
   be applied to the min-overlap problem? White is the most likely person to
   know whether his autoconvolution machinery composes with his min-overlap
   machinery, and the email is half-written.

**Why this beats the alternatives.** Same author, same Fourier toolkit, same
genus of problem (extremal for autocorrelation of a nonnegative density on a
finite interval), recent progress in the same year, and an existing
communication channel. This is the only bridge where a concrete technique
transfer is plausible within a single research session rather than requiring
new theory.

Refs: White, *Canad. Math. Bull.* (arXiv:2210.16437); arXiv:2508.02803 (2025);
arXiv:2506.16750 (2025); Cilleruelo–Ruzsa–Vinuesa, *Proc. AMS* 145 (2017);
Matolcsi–Vinuesa (2010); Cloninger–Steinerberger (2017).

---

## 2. TOP-3 ALSO-RANS

### 2a. Delsarte / Cohn–Elkies LP duality and **non-sharpness theorems**

**What it is.** Delsarte (1973) LP bound for codes; Cohn–Elkies (Ann. Math. 2003,
arXiv:math/0110009) continuous LP bound for sphere packing; Viazovska (2017)
showed the bound is sharp in dimensions 8 and 24. Critically, Cohn–de Laat–Salmon
(arXiv:2206.09876, *Adv. Math.* 2024) developed **dual LP bounds via discrete
reduction** that proved Cohn–Elkies is *non-sharp* in d=12, 16.

**Specific connection.** Our SDP, like Cohn–Elkies, is a Fourier-side LP with
PSD constraints (the Bochner moment matrix is precisely the cone of positive
trigonometric polynomials — exactly what Delsarte uses). The **non-sharpness
machinery from Cohn–de Laat–Salmon** lets one prove an LP-bound technique
*cannot* close a gap by constructing a discrete-reduced dual witness that
exceeds the known optimum. This is exactly what is missing from our retracted
Lasserre analysis: we showed empirically the tail-bound kills the gain, but we
do not yet have a **non-sharpness theorem** for the Bochner+ellipse family.

**Applicability.** If we can adapt Cohn–de Laat–Salmon's discrete-reduction
argument to our setting, we could potentially **prove** that the Bochner-PSD
hierarchy has a strict positive duality gap to `µ` — making `0.379544` a *proved
ceiling* of the technique, not just a tractability ceiling. This converts a
diagnostic ("no obvious path forward") into a theorem ("this LP class cannot
close the gap"). Note: arXiv:2410.04800 (2024) proposes a "New LP method"
that breaks the Cohn–Elkies framework in some dimensions — there may be a
direct analogue using three-point bounds (arXiv:2206.15373) for our problem,
since the autocorrelation `R_h(k) = ⟨h, T_k h⟩` is intrinsically two-point and
a three-point analogue would see configurations the two-point LP cannot.

Refs: Cohn–Elkies, *Ann. Math.* 157 (2003), arXiv:math/0110009;
Cohn–de Laat–Salmon, arXiv:2206.09876 (2022/2024); arXiv:2211.09044 (d=6);
arXiv:2206.15373 (three-point bounds).

### 2b. Beurling–Selberg extremal majorants/minorants

**What it is.** Beurling (1930s) and Selberg (1970s) constructed entire
functions of exponential type that majorize/minorize indicator functions of
intervals while having Fourier transforms supported on a fixed band. The
"box minorant problem" (Carruth–Gonçalves, arXiv:1702.04579) is precisely the
question: what is the largest bandlimited function bounded below by `1_{[a,b]}`?

**Specific connection.** White's program parametrizes `f` by the first `T`
Fourier coefficients with tail corrections. A Beurling–Selberg-style argument
gives **certified bandlimited majorants/minorants of `1_{f>α}`** for any level
`α`, which converts an SDP solution into an analytical bound that survives
truncation rigorously. The natural use here is: rather than enforce `f ∈ [0,1]`
via Bochner-PSD on `f` and `1−f`, enforce it via certified Beurling–Selberg
majorants — potentially giving a stronger constraint per unit of SDP variable
budget. The Cohn–Elkies framework itself uses Beurling–Selberg type majorants.

**Applicability.** Moderate. The technique converts pointwise constraints to
Fourier-side constraints in a different (and provably tighter for some
problems) basis than Bochner-PSD, and could potentially give a finite-dimensional
exactness argument that our current encoding lacks. Survey: Vaaler/Graham
1985, Carneiro–Vaaler 2010, and the survey at mc.sbm.org.br on
Beurling-Selberg majorants.

Refs: Carruth–Gonçalves, arXiv:1702.04579; Beurling–Selberg survey at
*Matemática Contemporânea*; Cohn–Elkies (uses BS-style majorants).

### 2c. Low-autocorrelation binary sequences (LABS) / merit-factor problem

**What it is.** Find a ±1 sequence `s = (s_1, …, s_L)` minimizing the aperiodic
autocorrelation side-lobes; the **merit factor** is `F = L²/(2 Σ_{k>0} C_k²)`
where `C_k = Σ_i s_i s_{i+k}`. Golay (1982) conjectured `lim sup F = 12.32`;
Littlewood (1966) conjectured `∞`; Høholdt–Jensen (1988) conjectured 6;
Borwein–Choi–Jedwab (2012/2013, arXiv:1205.0626) achieve `F > 6.34` on a
specific family. Equivalent (up to scaling) to the Bernasconi spin-glass model.

**Specific connection.** This is **the discrete combinatorial analogue** of
our problem. The merit factor minimizes `‖s∗s̃‖_2² − L²` (sum of squared
autocorrelations), while we minimize `sup_k R_h(k)` (max autocorrelation).
L² vs L∞ extremals differ in general but are closely related: Littlewood's
flat-polynomial conjecture (proved by Balister–Bollobás–Morris–Sahasrabudhe–
Tiba, Ann. Math. 2020) gives uniform `δ√n ≤ |P(z)| ≤ Δ√n` for some ±1
polynomial `P`, which is exactly an L∞ flatness statement on the autocorrelation
generating function. The BBMST construction uses **Spencer-style discrepancy
theory + probabilistic methods** — completely orthogonal toolset to our SDP.

**Applicability.** Lower-promise. The connection is real but the BBMST proof
is an *existence* result (constants `δ, Δ` not explicit and likely far from
optimal), not a sharp constant — they prove flat polynomials exist but don't
extremize. For minimum overlap we already have sharp numerics; what we need
is an *analytical* technique that transfers. The most actionable transfer
would be Borwein–Choi–Jedwab style explicit family constructions to push
the *upper* bound (Together side) rather than the lower bound.

Refs: Balister–Bollobás–Morris–Sahasrabudhe–Tiba, *Ann. Math.* 192 (2020),
arXiv:1907.09464; Jedwab merit factor survey 2005; Borwein–Choi–Jedwab,
arXiv:1205.0626; Bernasconi 1987; Mertens "Ground States of the Bernasconi
Model".

---

## 3. RULED OUT

### 3a. Sidon sets / B_h sequences (Cilleruelo et al.) — **partial overlap, not a bridge**

The technique of bounding `|A + A|` for Sidon-like `A` uses additive
combinatorics (Plünnecke–Ruzsa); the autoconvolution-inequality work
(Cilleruelo–Ruzsa–Vinuesa) overlaps with our top-1 bridge above, but the
"Sidon/B_h" combinatorial side does not transfer — the minimum-overlap problem
is not a sumset cardinality question, it is a continuous extremal problem on a
density. We've already captured the relevant analytical-Cilleruelo content
via the autoconvolution bridge.

### 3b. Erdős–Ko–Rado / Sárközy / Ruzsa density theorems — **no concrete bridge found**

EKR is an intersection-system theorem (uniform set families), not an
autocorrelation extremal. Sárközy's results on difference sets containing
squares/primes use sieve+circle-method machinery on integers, not on continuous
densities. There is no extant paper connecting these to the minimum-overlap
program; we found none in 20 minutes of searching.

### 3c. Bombieri–Selberg sieve / explicit-formula methods — **structurally different**

Sieve methods bound the count of integers escaping a residue system; the
extremal functional involved is a *majorizer* of an indicator function on
shifted residues. Beurling–Selberg majorants (top-3 above) do come from this
lineage — but the sieve toolbox itself (large-sieve inequality, Selberg's
λ-weights) is set up for sums over primes/residues, not for `[0,1]`-valued
densities on a continuous interval. The relevant *piece* (BS majorants) is
already extracted; the rest does not transfer.

### 3d. Aperiodic order / quasiperiodic functions (Penrose / cut-and-project / Meyer sets) — **plausible but speculative**

It is conceivable the true minimum-overlap optimizer is aperiodic (not a
trigonometric polynomial), since the Fourier coefficients of any periodic
optimizer satisfy a finite-rank polynomial system that has no obvious
solution. Meyer-set / model-set theory describes diffraction patterns of
aperiodic structures via almost-periodic measures. However, **we found no
paper using cut-and-project schemes to solve an extremal problem of our
type**. This is a "possibly true, but no machinery to apply" item.

### 3e. Bohr almost-periodic functions — **probably the right function space, but no extremal theorems**

The Bohr compactification gives the natural setting in which `sup_k inf` over
arbitrary shifts is well-defined for non-periodic `h`. If the optimum is
strictly outside the trigonometric-polynomial / band-limited class, the
extremizer lives in some Bohr-almost-periodic class. But: we found no
**extremal theorem** in the Bohr-AP class that we could apply to compute or
bound `µ`. The space is suggestive; the techniques are not deployable.

### 3f. Quasi-Monte Carlo / discrepancy theory — **structurally different functional**

Low-discrepancy sequences minimize sup over hyperrectangles of |#points ∈ R −
N·vol(R)|; our problem maxes shift-autocorrelation. The Roth lower bound
`D_N ≥ c (log N)^{(d-1)/2}` uses orthogonal function method; the analogue here
would be a lower bound on `sup_k R_h(k) − (∫h)²`, which is a *fluctuation*
bound, not an extremal-density question. Different objective.

### 3g. Khintchine / Rademacher–Steinhaus concentration — **wrong direction**

These bound the probability that random `±1` sums deviate from typical
behavior. Our problem is deterministic — we want a *specific* `h` whose
autocorrelation is flat, not a probabilistic flatness statement. (BBMST does
use probabilistic methods to construct flat polynomials, but that goes the
other way — it's a way of *building* the optimizer, not bounding it.)

### 3h. Other Erdős cash-bounty problems — **no shared technique found**

Sum-free sets (Croot et al.), Erdős multiplication table (Ford), Erdős discrepancy
(Tao) — all are integer-combinatorial. The one *technique* that recurs is the
**polynomial method** (Croot–Lev–Pach, Ellenberg–Gijswijt for cap sets), but
the polynomial method is a finite-field / mod-p technique that does not apply
to a continuous functional inequality.

---

## 4. OPEN QUESTION (the specific theorem to look up)

**Does Cohn–de Laat–Salmon's discrete-reduction dual-bound machinery
(arXiv:2206.09876, *Adv. Math.* 2024) apply to White's program?**

Specifically: their method takes an *infinite-dimensional* LP (Cohn–Elkies for
sphere packing) and constructs, by discretizing the dual, a *finite* dual LP
whose value upper-bounds the primal — and they use this to prove **strict
duality gaps** in dimensions 12 and 16. White's program is already finite, but
it has a continuous limit (`N, T, R → ∞`) which is the actual `µ`-defining LP.
Translating Cohn–de Laat–Salmon to our continuous LP would either:

  (a) **Prove a strict duality gap** for the Bochner-PSD hierarchy → upgrade
      our diagnostic "no path forward at tractable scale" to the theorem "no
      path forward at any scale within the Bochner+ellipse family" — a
      publishable rigorous-impossibility result.

  (b) **Fail to find a gap** → empirically suggest the SDP hierarchy might be
      sharp in the limit, which would shift the strategic focus to scale
      (precision arithmetic, larger N) rather than new constraints.

Either outcome is more informative than the current state.

A second open question: **Can the L²-autoconvolution uniqueness theorem (White,
arXiv:2210.16437) be adapted to characterize the minimum-overlap optimizer?**
If `h*` is unique, then the Together upper-bound construction must converge
*structurally* (Fourier coefficients match) to `h*` as the step count grows.
A converse "smoothing" of the Together construction could provide a *certified*
upper bound matching our lower bound to higher precision.

---

## 5. HYPOTHESIS

If the goal is **a new lever on the lower bound**, the most likely productive
move is to send the existing email to White (top-1 bridge) and read
arXiv:2210.16437 carefully. If the goal is **a rigorous diagnostic** that
the current technique is at its ceiling, the most likely productive move is
to study Cohn–de Laat–Salmon (top-2a) and adapt their discrete-reduction
duality to White's program. The current research note already documents the
tractability ceiling (`+5.4e-4` over White) — converting that to a *proved*
duality-gap theorem would be the natural next paper, independent of any new
upper bound on `µ`.

---

## Sources (selected)

- White, *Acta Arith.* 2023, arXiv:2201.05704 — the program being augmented.
- White, *Canad. Math. Bull.*, arXiv:2210.16437 — L² autoconvolution, top-1
  bridge.
- arXiv:2508.02803 (2025), arXiv:2506.16750 (2025) — autoconvolution lower
  bound to 0.94136.
- Cilleruelo–Ruzsa–Vinuesa, autoconvolution inequalities literature.
- Cohn–Elkies, *Ann. Math.* 157 (2003), arXiv:math/0110009.
- Cohn–de Laat–Salmon, *Adv. Math.* 2024, arXiv:2206.09876 — top-2a bridge,
  non-sharpness via discrete-reduction duality.
- Carruth–Gonçalves, arXiv:1702.04579 — Beurling–Selberg box minorant.
- Balister–Bollobás–Morris–Sahasrabudhe–Tiba, *Ann. Math.* 192 (2020),
  arXiv:1907.09464 — flat Littlewood polynomials exist.
- Jedwab, merit factor survey, 2005; Borwein–Choi–Jedwab arXiv:1205.0626.
- Together Computer 2026 GitHub release; AlphaEvolve 2025; TTT-Discover 2026
  — current upper-bound constructions.
