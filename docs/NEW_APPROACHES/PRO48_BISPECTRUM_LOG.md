# PRO-48 — Bispectrum / 3-point SDP lift: running log

Goal: break the conjectured 2-point ceiling **C_∞ ≈ 0.380558** by bolting a genuine
cubic bispectrum PSD block ONTO the cell-envelope (NOT replacing it, per PRO-22),
solving the row-4 center small-N, and reporting gain NET of the cubic-moment tail.

Object: triple correlation `R3(s,t)=∫f(x)f(x+s)f(x+t)dx ≥ 0` (f≥0), Fourier
coefficients = bispectrum `B(m,n)=y_m y_n conj(y_{m+n})`, `y_k=f̂(k)`.
Code: [`lp_research_state/code/_pro48_bispectrum.py`](../../lp_research_state/code/_pro48_bispectrum.py).

Gates (advance one per iteration): (i) dual rises above 2-point value; (ii) block
active (nonzero dual); (iii) PRO-22 validity — reconstruct f∈[0,1], Ω≤true sup_t(f⋆f).

---

## Iteration 1 (2026-06-06) — block built, validated, INACTIVE at small N

**Built.** `add_bispectrum_block(cons, c, d, L, Omega)` adds:
- B(m,n) as independent moment variables on the difference set of Λ={0..L}².
- **(PSD)** Bochner-2D moment matrix `G[a,b]=B(a−b) ⪰ 0` (R3≥0).
- **(symmetry)** B(−δ)=conj(B(δ)); 6-fold perm B(m,n)=B(n,m)=B(−m−n,n).
- **(marginal link to 2-point)** B(0,0)=1/8; B(m,0)=B(0,m)=B(m,−m)=½Pₘ with
  Schur slack Pₘ≥|yₘ|²=(cₘ²+dₘ²)/4.
- **(localizing link to OBJECTIVE)** R3≤Ω (since 0≤f≤1 ⇒ R3≤A(s)≤Ω), encoded as
  `G2[a,b]=Ω·δ−B(a−b) ⪰ 0`.

**Self-test PASS** (5 genuine nonneg f = |trig|² normalized to ∫f=½): both G and
G2 PSD (min eig ≈ 0.05–0.20), all symmetry/marginal/center identities to <1e-7.
⇒ block is a valid relaxation; the true f stays feasible ⇒ any rise is real, not
an over-constraint artifact.

**Gate (i)/(ii) at row-4 center, N=300,T=120,R=5,bochner_n=6, MOSEK:**

| config | Ω | gain vs 2-pt | dual L1 (PSD/LOC) |
|---|---|---|---|
| baseline 2-pt | 0.368121171 | — | — |
| +bispec L=3 (marginal only) | 0.368121168 | −3e-9 | 5e-10 / — |
| +bispec+LOC L=3 (default tol) | 0.368121261 | **+9.0e-8** | 1.3e-7 / 8.7e-8 |
| +bispec+LOC L=4 (default tol) | 0.368121257 | +8.5e-8 | 2.0e-7 / 1.4e-7 |
| +bispec+LOC L=3 (**tight tol** 1e-12) | 0.3681211679 | **−5.8e-10** | **2.1e-11 / 1.4e-11** |

**Verdict: the +9e-8 was MOSEK solver noise.** Under tight tolerances the gain
collapses to ~0 and both PSD blocks go inactive (dual ~1e-11). At small N the
**block is INACTIVE** — the 2-point + cell-envelope relaxation already admits a
PSD completion of the cubic moments.

**Why (diagnosis):** at N=300 the cell-envelope discretization slack ~π/(2N)≈5e-3
dominates Ω; the objective and every binding constraint are 2-point. The cubic
block couples only through (a) loose marginals and (b) the localizer R3≤Ω, which
is itself slack (R3 ≪ Ω typically). So the small-N detector cannot tell "idea is
implied" from "cell slack masks it."

**Next iteration (gate i retry):** strengthen the localizer from the scalar
`R3≤Ω` to the **localizing matrix** `R3(s,t)≤A(s)` (uses the autocorrelation
function, tighter), and/or reduce the base slack. If even that is inactive at
tight tolerance, the approach is decisively dead at small N.

---

## Iteration 2 (2026-06-06) — DECISIVE NEGATIVE: block is IMPLIED by 2-point + cell-envelope

**Added** the tighter **localizing matrix** `(A(s) − R3) ⪰ 0`
(`= ∫ f(x)f(x+s)(1−f(x+t))dx ≥ 0`, valid since 0≤f≤1), Fourier
`h^(m,n) = |y_m|² δ_{n,0} − B(m,n)` with `A_m` carried by the Schur slack `P_m`.
Fixed a self-test generator bug (random `|g|²` could exceed 1, violating f≤1 —
now `f = ½ + Σ a_k cos`, `Σ|a_k| ≤ ½`, guaranteeing f∈[0,1], ∫f=½). With valid
f∈[0,1], **self-test PASS**: G, G2[Ω−R3], G3[A−R3] all PSD.

**Two independent decisive tests, row-4 center, MOSEK tight tol:**

1. **Re-solve gain** (full block PSD + Ω-localizer + A-localizer): +6.5e-10 (L=3),
   +2.7e-11 (L=4) = ZERO; all dual norms ~1e-10.
2. **Separation probe** — fix (c*,d*,Ω*) to the 2-point optimum, ask if the
   bispectrum block is feasible (i.e. admits a valid completion). **Completion
   EXISTS at every N tested:**

   | N | bn | Ω(2-pt) | probe verdict |
   |---|---|---|---|
   | 300 | 6 | 0.368121 | COMPLETION EXISTS (inactive) |
   | 1500 | 12 | 0.376912 | COMPLETION EXISTS (inactive) |
   | 6000 | 16 | 0.379208 | COMPLETION EXISTS (inactive) |

   As N→ceiling (0.368→0.379, the 2-pt bound climbs most of the way to
   C_∞≈0.3806), the optimum STILL admits a bispectrum completion at T0=8.

**Airtight conclusion.** If the 2-point optimum (c*,d*,Ω*) admits a completion
(B*,P*), then (c*,d*,Ω*,B*,P*) is feasible for the augmented program at objective
Ω* ⇒ augmented_min = baseline_min. **The cubic block (T0≤8, marginal + R3≤A +
R3≤Ω coupling) CANNOT raise the bound.** Confirmed at the binding row across a
20× range of N spanning most of the gap to the ceiling.

**Why (structural).** The objective M = sup of the **2-point** autocorrelation;
the cubic block couples to it only via (a) marginals `B(m,0)=½P_m` (one index =0)
and (b) UPPER localizers `R3 ≤ A, Ω`. The interior moments `B(m,n)` (m,n,m+n≠0)
are FREE, and there is no natural **lower** bound on R3 in 2-point terms — so the
PSD-completion freedom always absorbs the constraints. Larger T0 only adds MORE
free interior entries ⇒ easier to complete, never harder. This is exactly the
shortlist's stated #1 risk ("the triple cut may be IMPLIED by 2-point+cell-
envelope for this pairwise objective — inactive block, zero gain").

**Tail bound (Henrion–Rudi): moot.** The finite, tail-free small-N test already
yields zero gain, so there is no gain for the cubic-moment tail to erase. The
protocol's purpose ("detect any gain before re-entering the tail trap") is served
— it detects none.

### VERDICT — gate (ii) fails decisively → loop STOPPED.
The bispectrum/3-point lift, as a PSD block bolted onto the 2-point cell-envelope
with marginal + upper-localizer coupling, does **not** break C_∞ at T0≤8. It is
implied by the existing relaxation at the binding row across N=300–6000. **No
overclaim made; µ headline and findings.md untouched.**

**What it would take to revive (research, not engineering):** a constraint that
genuinely makes the objective or feasible set depend on 3-point data the 2-point
program cannot complete — e.g. an objective reformulation sensitive to the
bispectrum, or a *lower* bound on R3 forcing the interior moments. PRO-49
symmetry-reduction (the enabler for scale) does NOT help: an inactive block stays
inactive at every scale.
