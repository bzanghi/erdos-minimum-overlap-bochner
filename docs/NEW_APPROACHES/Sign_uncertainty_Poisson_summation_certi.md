# Deep-dive vetting: Sign-uncertainty / Poisson-summation certificate ON the autocorrelation A

**Date:** 2026-06-06. **Status:** ADVERSARIAL ASSESSMENT — not a result. Proposes/vets only.
**Approach under review:** RANKED_SHORTLIST item **#6** = HA-2 (`harmonic_analysis_lens.md`).
Absorbs the "native Cohn–Elkies/Delsarte LP" candidate.

**Verdict: WEAK.** The numerics below show the certificate, in the form proposed (a dual
certificate on the *continuous positive-definite A with no cell-envelope*), is either
**invalid** (gives `µ ≥ 0.76`, contradicting the known UB `µ ≤ 0.380871`) or, when run in the
correct direction, **vacuous** (`µ ≥ 0.0116`). The reason is structural and is the *same*
failure that killed PRO-22's direct sup_t SDP: the constraint set
`{Â ≥ 0, A ≥ 0, ∫A=1, A(0) ≤ 1, supp A ⊆ [−2,2]}` is a **strict superset** of the realizable
autocorrelations `{A = h⋆h̃ : h ∈ [0,1], ∫h=1, supp h ⊆ [0,2]}`. The box constraint `h ≤ 1`
has **no faithful description purely in A-space**, and it is *load-bearing*. Positive-definiteness
of A — the one piece of structure the sign-uncertainty machinery exploits — is far too weak to
reach the binding value. The approach is genuinely previously-untested here, and the vetting is
itself a worthwhile rigorous no-go (it converts the shortlist's "#1 risk: same cone" into a
sharper statement: *without* the envelope it is an INVALID larger cone, and *with* the envelope
re-added to restore validity it collapses to White's existing dual).

---

## 0. What the approach claims

Reformulate (PRO-29, verified below): `A_h(t) = (h⋆h̃)(t) = ∫h(x)h(x+t)dx`, so
`Â_h = |ĥ|² ≥ 0` (A is positive-definite), `A_h ≥ 0`, even, `∫A_h = 1`, `supp ⊆ [−2,2]`,
`A_h(0) = ‖h‖₂² ∈ [1/2,1]`, and
> `µ = inf_h sup_t [1 − A_h(t)] = 1 − sup_h inf_t A_h(t)`.

The proposal: treat min-overlap as a `+1`-type sign-uncertainty instance; build an even auxiliary
`g` with `ĝ` sign-controlled and `g ≤ A` outside a window `[−t₀,t₀]`, optimize `g` in a degree-~40
cosine basis to force a floor on the relevant extremum of A, sharpen with sign-uncertainty
refinements (modular-form `g`) that beat the Logan first-zero bound, and upgrade to an EXACT
inequality via Poisson summation over `2ℤ` (matching the support). Refs: Bourgain–Clozel–Kahane
(*Ann. Inst. Fourier* 2010); GOSS arXiv:2003.10771 / 2003.10765; Carneiro–Quesada-Herrera
arXiv:2006.00959; Cohn–Elkies *Ann. Math.* 157 (2003) arXiv:math/0110009; Cohn–de Laat–Salmon
arXiv:2206.09876.

---

## 1. FEASIBILITY with available tools

**High** — the first probes are cheap and were run here in minutes.
- The finite certificate LP (cosine basis on A, sign constraints, Parseval) is a small cvxpy/CLARABEL
  program (`D≈60`, grid 200–400) — solves instantly. Built and solved below.
- Sharp modular-form `g` (the only part that could be strong) would need a dim-1 Viazovska-style
  ansatz that is **not** off-the-shelf; that is a months-scale research artifact, not an
  engineering task. mpmath/sympy can manipulate candidate `g`, but constructing the eigenfunction
  is the hard open step (no dim-1 modular symmetry analogue exists for this transform — same wall
  the shortlist flagged for the Viazovska-transplant lottery ticket).
- Poisson summation over `2ℤ` is symbolic/numeric-easy; it is *not* the bottleneck.

So the cheap LP probe is fully feasible (done), and it is exactly the probe that settles the
question. The expensive modular-`g` route is feasible only as open research.

---

## 2. PRIOR ART

**The technique is real and active; none of it has touched this problem.**
- **Sign uncertainty:** Bourgain–Clozel–Kahane 2010 (origin); Cohn–Gonçalves *Invent. Math.* 2019
  (the sharp `−1`/+12-dim eigenfunction via modular forms); GOSS "New sign uncertainty principles"
  arXiv:2003.10771 and "regularity & mass concentration" arXiv:2003.10765 (best 1-D numerics:
  `A₊(1)` conjecturally ≈ 0.573, rigorous lower `1/√(2πe) = 0.2419`); Carneiro–Quesada-Herrera
  arXiv:2006.00959 (weighted/derivative variants). Fresh: "Sign uncertainty and low-degree
  polynomials" arXiv:2210.01684 (2022); "Fourier inequalities and sign uncertainty"
  arXiv:2505.15994 (2025). These are the right family and they are live.
- **LP-on-positive-definite-functions:** Cohn–Elkies 2003; Cohn–de Laat–Salmon arXiv:2206.09876
  (discrete reduction → finite dual bounds, the rigor template the proposal would borrow).
- **Applied to Erdős min-overlap / autocorrelation-sup: NOTHING.** White (arXiv:2201.05704) uses
  "elementary Fourier analysis → convex optimization," not sign uncertainty. The web/arXiv scan
  returns no paper joining sign-uncertainty/Poisson-summation certificates to min-overlap or to
  `sup_t (h⋆h̃)`. So the lever is genuinely unbuilt here.
- **In-repo:** the shortlist's grep audit confirms `sign-uncertainty/viazovska` returns only name
  collisions — no real code. The *closest* prior in-repo work that already bears on it:
  - **OUT_OF_BOX_INEQUALITIES item 10** (Logan/de Branges/Beurling): showed Logan-type first-zero
    bounds on a nonneg-Fourier-transform `A` give `sup A ≥ (∫A)/L_eff ≈ 1/8`. The sign-uncertainty
    refinement only improves the *constant* in front, not the structure.
  - **PRO-22** (`LEVER_SUPT_DIRECT.md` / `docs/.../PRO22`): dropping the cell-envelope and bounding
    M(t) directly gave an apparent `+4.4e-3` that was **invalid** — reconstructed `f` exited `[0,1]`,
    true Ω was 5× the reported value. The cell-envelope is "necessary for validity, not just a
    relaxation."
  - **BS_MR_NORMALIZATION / PRO-29 / PRO-32**: the autocorrelation reformulation `µ = 1 − sup_h inf_t A`
    is exactly here; the naive transfer to A-only bounds manufactured the spurious `µ ≥ 0.59`.

So: technique untested here ✔, but the in-repo PRO-22 / OUT_OF_BOX-10 results already predict the
failure mode, which the experiment below confirms quantitatively.

---

## 3. FIRST EXPERIMENT (run here — cheap analytic/LP probe, no heavy solve)

All on Together's near-optimal `h*` (`data/together_f_star.json`, 600-cell step on `[0,2]`,
`∫h=1`), plus a small certificate LP. Reproduced numbers:

**(a) Verify the reformulation & locate the binding object.**
```
∫h = 1.000000 ;  A(0) = ‖h‖₂² = 0.774675 ;  ∫A = 1.000000 ;  supp A = [−2,2]
Together UB  M_T = (2/n)·max corr(h,1−h) = 0.3808703   (matches claimed 0.380871)
Binding overlap shift  t* = 0.11 (cell 33),  A(t*) = 0.6191297
sup_{|t|≥0.15} A_{h*}(t) = 0.61913 = µ_dual = 1 − µ  ⇒  µ = 1 − 0.61913 = 0.38087  ✔
```
So the object a *correct* certificate must control is `sup_h inf_t A_h = µ_dual = 0.61913`,
attained at an **interior** shift `t* = 0.11` — NOT a first-zero / tail length scale.

**(b) The pure-A certificate, run in the CORRECT direction (LB on µ).**
To lower-bound `µ` we need an UPPER bound on `V = sup_h inf_t A_h`. Best single-atom certificate
`inf_t A_h ≤ A_h(t₀)`; bound `sup_h A_h(t₀)` using only the no-envelope structure
(`Â ≥ 0` ⇒ cosine coeffs `b_m ≥ 0`; `∫A=1` ⇒ `b₀=1/4`; `A(0)=b₀+2Σb_m ≤ 1`; `A ≥ 0` on grid):
```
min_{t₀≥0.05}  [ max_A A(t₀) ]  = 0.98839   at t₀ = 1.377
   ⇒  µ ≥ 1 − 0.98839 = 0.01161        (VACUOUS)
```
The relaxed A-cone contains autocorrelations that sit at ≈0.99 at essentially any single shift, so
no single-atom (and, by the same mechanism, no finite probe-measure) certificate built only on
positive-definiteness clears even `µ ≥ 0.30`.

**(c) The same LP run in the NAIVE direction exposes the invalidity.**
Minimizing `sup_{|t|≥t₀} A` over the relaxed A-cone (the "force A to spread its mass" framing the
memo literally proposes):
```
t₀=0.10:  min_A sup_{|t|≥t₀} A = 0.23072  ⇒  "µ ≥ 1 − 0.23072 = 0.76928"   ← FALSE
t₀=0.15:  min_A sup_{|t|≥t₀} A = 0.22068  ⇒  "µ ≥ 0.77932"                 ← FALSE
```
These contradict `µ ≤ 0.380871`. They are not a sign error in the experiment; they are the
**signature of an invalid relaxation** — the cone is too big.

**(d) Explicit witness of the realizability gap.**
`A(t) = ¼(1+cos πt)` on `[−2,2]` satisfies `Â ≥ 0` (atoms at `0, ±1`), `∫A=1`, `A(0)=½`,
`A ≥ 0`, and `inf_{|t|≥0.15} A ≈ 0` — so a pure-A certificate would "prove" `µ ≥ 1`. But its
`ĥ = √Â` is a sum of Diracs, so the only `h` realizing it is a sum of point masses — **not an
`L^∞[0,1]` density.** Inadmissible. Concrete proof that `{Â≥0, A≥0, ∫A=1, A(0)≤1}` ⊋
`{h⋆h̃ : h∈[0,1]}`.

**(e) Sign-uncertainty constant ceiling (closing the "modular g is much stronger" hatch).**
The sharpest 1-D sign-uncertainty constant bounds a *first-zero / concentration length* `r*`,
yielding a magnitude floor of order `1/L_eff`. Even a hypothetical `3×` tightening over Logan
gives `sup A ≳ 1/(4/3) = 0.75` of the support-scale — but that is a bound on `A(0)`-type
concentration, **not** on the interior value `A(t*) = 0.619` that sets `µ`. The binding number is
fixed by **box-constrained overlap geometry** (how much `h` and `1−h` can simultaneously avoid
each other under `h ≤ 1`), which no sign-uncertainty/first-zero constant sees. There is no path
from `A₊(1) ≈ 0.573` to the `0.61913` floor.

**Code:** ad-hoc, run inline (numpy `np.correlate` on `together_f_star.json`; cvxpy/CLARABEL
cosine-basis LP). Not persisted as a library (throwaway, per repo `_`-convention norms).

---

## 4. WHY IT COULD (in principle) BEAT THE ~0.380558 SATURATION — and why that path is blocked

The *hope* was legitimate and worth stating precisely:
- The conjectured framework ceiling `C_∞ ≈ 0.380558` is **provably a property of White's
  cell-envelope + Bochner cone specifically** (PRO-6: `f_CB(N|n)=f_C(N)`). A certificate that does
  **not** use that cone is not subject to *that* saturation; its ceiling would be "certificate
  quality," a different limit. So *a priori* a genuinely-different-cone analytic certificate could
  pierce 0.380558. This is the same logic that makes the 3-point/bispectrum lift (shortlist #2) the
  principled candidate.
- Sign-uncertainty certificates have repeatedly hit **sharp** constants in sibling problems
  (12-dim, sphere packing in 8/24), via modular-form `g`. If such a `g` existed for *this*
  constraint set, it could in principle give a clean, framework-independent (even closed-form)
  bound — the holy grail the PSLQ hunt sought from the wrong end.

**Why it is blocked (the core finding):** the certificate's *only* structural lever is
`Â_h ≥ 0` (positive-definiteness of A). The experiment shows that lever, plus `∫A=1`, `A(0)≤1`,
support, and `A ≥ 0`, is **nowhere near enough** — it admits autocorrelations that no `[0,1]`-valued
`h` can realize, so the resulting bound is either invalid (0.76) or vacuous (0.0116). The piece that
makes White's program *valid and tight* is the **realizability link** `A = h⋆h̃ with h ∈ [0,1]`,
encoded by the cell-envelope as a Parseval identity between the Fourier coefficients of `M` and of
`h` (White (W.1); see PRO-22). That link is **not** positive-definiteness, is **not** captured by
any sign-uncertainty refinement, and has **no representation purely in A**. To make the certificate
valid you must re-inject the realizability link — at which point you are back inside White's cone
(its dual), so you inherit `C_∞`, not a new ceiling. The approach therefore cannot beat the
saturation: it is invalid without the envelope and envelope-equivalent (capped at `C_∞`) with it.

This sharpens the shortlist's stated "#1 risk (could be a dual description of the same cone)":
the truth is worse than equivalence — **without** the envelope it is a strictly LARGER, INVALID
cone (PRO-22 redux); **with** the envelope it is White's dual. There is no intermediate regime
where it is both valid and stronger.

---

## 5. RISKS / why it might fail (consolidated)

1. **(Realized, fatal) Positive-definiteness ≠ realizability.** The box `h ≤ 1` is load-bearing
   and absent from A-space. Pure-A certificate → invalid (`µ ≥ 0.76`) or vacuous (`µ ≥ 0.0116`).
   This is PRO-22's exact failure mode in sign-uncertainty clothing. **This alone sinks the
   headline goal.**
2. **Wrong limiting object.** Sign-uncertainty bounds a first-zero/concentration *length*; the
   binding value `A(t*)=0.619` is an *interior magnitude* set by box-constrained overlap geometry.
   No sign-uncertainty constant (even sharp `A₊(1)≈0.573`) reaches it. The Logan one-liner already
   collapses to `1/8` (OUT_OF_BOX-10); the refinement only nudges the constant.
3. **Cone-equivalence on the valid sub-problem.** If one *does* restrict the certificate to valid
   autocorrelations (re-adding the realizability link), it becomes White's dual — capped at
   `C_∞ ≈ 0.380558`, no gain. (Still a *worthwhile rigorous no-go*: a clean proof "any
   positive-definite-only certificate on A is either invalid or ≤ White" closes a recurring
   temptation in this repo — PRO-22, PRO-29, PRO-32 all flirted with it.)
4. **Modular-`g` unavailability.** The only version with a shot at sharpness needs a dim-1
   Viazovska-style eigenfunction; no modular symmetry analogue exists for this transform
   (same wall as the dropped Viazovska-transplant ticket). High research cost, low odds.
5. **Poisson-summation step is a red herring for the gap.** Upgrading the finite LP to an exact
   inequality via `2ℤ`-summation is easy but only certifies whatever the (weak) finite certificate
   already gives — it does not add the missing box structure.

---

## 6. Bottom line / disposition

- **Verdict: WEAK.** Highest *informational* value is as a **no-go**: the experiment is a clean,
  cheap, rigorous demonstration that any certificate built on the positive-definiteness of `A`
  alone (the sign-uncertainty/Cohn–Elkies-on-A lever) cannot bound `µ` non-trivially without the
  realizability/cell-envelope link, and with that link it collapses to White's dual. This *kills
  the bare "native Cohn–Elkies/Delsarte LP" candidate too* (the shortlist already merged it here),
  with numbers.
- **Do NOT invest** in the modular-`g` / dim-1 Viazovska construction for min-overlap: even a
  sharp sign-uncertainty constant targets the wrong (length-scale) object and cannot reach the
  `0.619` interior floor.
- **Where the live value actually is** (redirect): the realizability link the certificate is
  missing is precisely what a *larger cone* supplies. The shortlist's #2 (3-point / bispectrum SDP
  lift, kept on the cell-envelope per PRO-22) and #4 (NPA operator-localizer that inserts `h ≤ 1`
  *losslessly* as an operator inequality) are the principled ways to add structure beyond
  positive-definiteness. This deep-dive's contribution is to confirm that the sign-uncertainty
  lane, attractive as it looks, is the PRO-22 trap and should not be promoted.

**Reproducibility note.** Numbers in §3 are from `data/together_f_star.json` via
`np.correlate` (autocorrelation and `corr(h,1−h)`) and a cvxpy/CLARABEL cosine-basis LP
(`D=60`, grid 200–400). The headline UB reproduces to `0.3808703` (vs claimed `0.380871`); the
binding interior shift is `t*=0.11`, `A(t*)=0.6191297`; the pure-A certificate ceiling is
`µ ≥ 0.0116` (correct direction) and the naive-direction value `µ ≥ 0.769` exhibits the
invalidity. No persisted code (throwaway probe).
