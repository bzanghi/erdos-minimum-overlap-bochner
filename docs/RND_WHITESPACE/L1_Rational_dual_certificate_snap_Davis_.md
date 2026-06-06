# L1 — Rational dual-certificate snap (Davis–Papp): ADVERSARIAL DEEP-DIVE

**Date:** 2026-06-03
**Workflow:** R&D whitespace scout → Phase-2 deep-dive on a single proposed direction.
**Scope:** Rigorously and adversarially assess the L1 direction (apply the Davis–Papp rational dual-certificate construction to White's binding-center dual, bypassing the SDPA-GMP serializer). **This memo SCOUTS and PLANS; it does NOT certify any bound.**

**Status of the direction before this memo:** *proposed only* — surfaced in the same-day external-literature scan (`docs/RND_WHITESPACE/2026-06-03_external_literature_scan.md`, lever **L1**) and in megamemory (`rd-whitespace-scout-2026-06-03-external-lit-track...`). It had **never been assessed against White's actual cone structure or Davis–Papp's stated preconditions.** That is what this memo does.

---

## TL;DR verdict

**MAYBE — but re-scope; the headline thesis rests on a false premise. Do not pursue as stated.**

The proposal's load-bearing claim is that *"White's dual is exactly [Davis–Papp's] featured class — univariate polynomial nonnegativity on a bounded interval after the cos-substitution."* **This is false.** White's program is an irreducibly **mixed LP + SOC + PSD** convex program whose *dominant* block is a ~28 000-variable LP encoding the discretized step-function `M(x)` — the load-bearing validity ingredient. Davis–Papp's algorithm certifies **one univariate polynomial nonnegative on one interval** (a single WSOS cone). There is no substitution reducing the whole program to that form.

Two further hard obstructions, both verified from the papers' text:

1. **Boundary obstruction.** Davis–Papp *requires the certified polynomial to be in the strict interior of the WSOS cone* and explicitly **cannot certify boundary points** (Scheiderer). White's binding bound is an **optimum with active constraints / complementary slackness** — i.e. on the boundary. You must back off to a strictly-weaker bound `c < c*`.
2. **It is not a float-snap.** The bit-size-bounded method (their Algorithm 1) **runs entirely in rational arithmetic** with a **barrier-Hessian oracle**; it does **not** consume the numerical dual. The "without rounding or projection" claim refers to not rounding a *primal* Gram matrix — not to a free numerical-dual → rational snap.

A **salvageable sub-target** exists (a Bochner-block-only exact certificate), but it collapses into directions **already scouted** (companion **D2** exact Fejér–Riesz dual; external **L4** Magron trigonometric exact-SOHS, which has a concrete published algorithm). The cheap, honest, *non-reductive* rigor upgrade is **L2 (Jansson/VSDP)**, which handles the mixed cone as-is. **Net: Davis–Papp is the wrong tool for the whole-program certificate; at best it is the symbolic finish on the Bochner sub-block, merged with D2/L4 — and it should be sequenced after L2.**

---

## 1. FEASIBILITY

### 1.1 What Davis–Papp actually requires (verified from full text)

I extracted both papers to text (`pdftotext`) and read the load-bearing statements verbatim:

**arXiv:2105.11369 (SIAM J. Opt. 2022) — existence.**
- *"rational polynomials on the boundary of the sums-of-squares cone may not have a rational sums-of-squares decomposition. On the other hand, polynomials in the interior of the sums-of-squares cone do have rational decompositions"* (lines 36–39).
- *"every rational polynomial in the interior of the WSOS cone has a rational dual certificate … (We cannot hope this to be true on the boundary; this follows from Scheiderer's seminal result …)"* (lines 711–714).
- **Theorem 3.5:** when `t − c*·1` is on the boundary of `Σ` (i.e. `c*` is the optimal WSOS lower bound), the gradient-certificate norm `‖1‖*_y` is related to the margin by `c* − c ≤ (C·‖1‖*_y)^{-1}` — so as `c → c*` the certificate norm **diverges**.
- A **dual certificate is a vector** `x ∈ (Σ*)°`, i.e. a point in the dual SOS cone (a moment/measure vector). *Not a scalar.*

**arXiv:2305.19039 (2023) — bit size + algorithm.**
- Bit-size bound (the explicit formula): `log(‖y‖∞) ≈ O(log(U) + log(‖t‖₂) + log(µ) + log(ν) + log(1/ε))` — **logarithmic** in the distance-to-boundary `ε`.
- *"… it is limited by the precision of the floating point arithmetic, and the bit sizes of the computed certificates cannot be bounded. The new algorithm proposed in this paper, Algorithm 1, runs entirely in [rational] arithmetic"* (lines 102–104). So the **bit-bounded** route is **not** a snap of the numerical dual; the snap-the-float variant (Algorithm 1 of [6]) is precision-limited and **has no bit-size guarantee.**
- **Algorithm 1 input:** *"A polynomial `t`; a tolerance `ε > 0`"*, plus **parameters: an oracle for computing the barrier Hessian `H` for `Σ`** and an initial certificate. It is an **exact-arithmetic interior-point-style iteration**, not a post-processing of CLARABEL's output.
- Featured special case: univariate polynomials nonnegative over the real line or a **bounded interval**, in several bases — but this is **one** WSOS cone with **known** Markov–Lukács weights.

### 1.2 White's program is NOT that class (verified by probe)

I built the binding-center program at a small size purely to inspect its cone structure (`build_problem(N=300, T=120, R=10, h=0.004, p=0.3875, q=±0.02, bochner_n=8)`, `white_full_convex.py`) and read the CLARABEL cone dimensions:

```
1 equality
2049 linear inequalities
SOC cones: forty 3-dim  + two 122-dim   (= 42 SOC blocks)
PSD cones: two 18×18                      (Bochner real-form, 2(n+1)=18 at n=8)
variables: 1163 (Ω, w[300], v[300], c[120], d[120], eps[10], dlt[10])
```

At **production** scale (`N≈10⁴–24 000`, `T=4000`, `bn≈40`) the **LP block alone** has ~2N+2T+… ≈ **28 000 variables** and tens of thousands of inequalities; the SOC ball constraints scale with `T`; the PSD blocks are ~`2(bn+1)≈82`.

Structurally (from `white_full_convex.py:143–262`):
- The LP block (`w,v ≥ 0`, `w,v ≤ Ω`, `L·Σ(w+v)=1`, the cell-consistency rows `lhs + 2·sq(a_m)+2·sq(b_m) − rhs ≤ 0`, the sine cell-consistency rows, the `eps/dlt` tail-bound rows, the box on `c,d`) **encodes the discretized step-function `M(x)`** via the cell averages `w_j, v_j`. This is exactly the ingredient the project's "direct sup_t SDP" experiment found is **load-bearing for validity** (dropping the cell-envelope produced an invalid bound — see ledger).
- The SOC cones come from `sum_squares(c)+sum_squares(d) ≤ 0.5` and the per-`m` quadratic terms.
- Only the **two PSD Bochner blocks** are honest trigonometric-WSOS objects.

**There is no `x = cos θ` substitution that turns this into "one univariate polynomial nonnegative on one interval."** The cos-substitution applies to a *single* trigonometric nonnegativity test; White's dual is a point in the dual of a **product cone** dominated by LP, not a univariate WSOS membership. The proposal conflates "White's *test family* `cos(πmx/2)` is trigonometric" with "White's *whole dual* is a univariate WSOS instance." It is not.

### 1.3 The dual is not even extracted as a vector

`dual_extractor.py` (read in full) parses CLARABEL's verbose **iteration log** and returns only the **dual objective scalar** `dual_obj` (the `dcost` column, printed to ~5 sig figs) at the last iteration with `dual_residual ≤ 1e-4`. It returns **no dual vector**. Davis–Papp needs the dual *vector* `x ∈ (Σ*)°`. So even for the Bochner sub-block, one would first have to (re)solve to obtain a numerical dual matrix — which is the same step D1/D2 require.

### 1.4 The boundary obstruction is real here

White's binding bound is an **optimum** — the whole PRO-6 "saturation / complementarity" line of work is precisely about which constraints are **active** at the optimum (`docs/archive/PRO6_COMPLEMENTARITY_PROOF.md`). Active constraints ⇒ the relevant nonnegativity is **tight** ⇒ the dual point is on the **boundary** of the cone, exactly the case Davis–Papp **cannot** certify rationally. One must certify a backed-off `c < c*`.

**Cheap probe — is the backoff affordable?** White's margin over White-2023 is `0.3802838 − 0.379005 = 1.279e-3`. To still beat White you need `ε < 1.279e-3`:

| `ε` (backoff) | certified LB | > White? | `log₁₀(1/ε)` (bit-scale) |
|---|---|---|---|
| 50% of margin (6.4e-4) | 0.379644 | yes | 3.2 |
| 10% of margin (1.3e-4) | 0.380156 | yes | 3.9 |
| 1% of margin (1.3e-5) | 0.380271 | yes | 4.9 |

So the **ε-backoff is benign in isolation** (logarithmic cost, and even 1% of the margin clears White). The blocker is **not** the bit size — it is (a) the program is the wrong cone class, and (b) the binding margin **over the project's own headline** (not over White) is only `+2e-4`, and Davis–Papp would certify something `ε` *below* that, eating into a margin the verification memos already flag as thin.

### 1.5 Tooling / effort reality

- **No released Davis–Papp code.** The only related software is `alfonso` (Papp–Yıldız, arXiv:2101.04274) — a generic **nonsymmetric conic solver**, not the rational-certificate Algorithm 1. The algorithm (incl. the **barrier-Hessian oracle** for the interval WSOS cone and the continued-fraction rational rounding) must be **re-implemented from the paper**.
- `mpmath`/`sympy` are present and adequate for the exact-arithmetic verification *step*, but not for the whole iteration.

**Feasibility verdict: the headline ("apply Davis–Papp to the binding dual, univariate-interval case") is NOT feasible** — wrong cone class, boundary optimum, no dual vector extracted, no code. A **Bochner-block-only** sub-application is feasible *in principle* but is dominated by existing scouted work (D2/L4) and still leaves the LP/SOC blocks to certify separately.

---

## 2. PRIOR ART

### 2.1 In-repo

- **Already proposed, not assessed:** `docs/RND_WHITESPACE/2026-06-03_external_literature_scan.md` lever **L1** (verbatim same arXiv refs, same thesis); megamemory `rd-whitespace-scout-2026-06-03-external-lit-track-cited-proof-grade-levers-l1-l4`. The proposal text under deep-dive is essentially a restatement of that scan. **This memo is the first adversarial technical check.**
- **Adjacent scouted directions it overlaps / collapses into:**
  - **D2** (`2026-06-03_whitespace_scout.md`): *exact Fejér–Riesz / trigonometric-SOS dual certificate for the binding center* — this is the honest home for the Bochner-block sub-certificate.
  - **L4** (external scan): **Magron et al.**, *Exact SOHS decompositions of trigonometric univariate polynomials* (arXiv:2202.06544) — a **concrete published algorithm with a bit-complexity bound** for exactly the trigonometric-circle WSOS object the Bochner block is. **This is a better-supported tool than Davis–Papp for the Bochner sub-block** (it is the literal Fejér–Riesz exact dual; it does not assume strict interior in the same fatal way because the perturb-and-compensate hybrid is designed for positive-on-the-circle polynomials).
  - **L2** (external scan): **Jansson/VSDP** rigorous a-posteriori SDP bound (Jansson–Chaykin–Keil, SIAM J. Numer. Anal. 2007) — handles **LP + SOC + SDP natively** by post-processing the **existing** solve. **This is the correct first rigor upgrade for the mixed cone**, with no reduction.
  - **D1/D5**: SDPA-GMP serializer + GMP re-solve + rational snap — the route Davis–Papp was pitched to *bypass*; but Davis–Papp does not bypass it (it still needs a (re)solve to get a dual vector).
- **Cohn–Elkies dual-certification culture** (external scan): the general "certify the dual witness in exact/interval arithmetic" practice. Davis–Papp and VSDP both live in this culture; **VSDP is the lighter member for a mixed cone.**

### 2.2 External (cited)

- Davis & Papp, *Dual certificates and efficient rational SOS decompositions for polynomial optimization over compact sets*, SIAM J. Opt. 2022 — arXiv:2105.11369. **Interior-only** existence (boundary excluded, Scheiderer).
- Davis & Papp, *Rational Dual Certificates for WSOS Polynomials with Boundable Bit Size*, 2023 — arXiv:2305.19039. Bit size `O(… + log(1/ε))`; **Algorithm 1 is exact-arithmetic with a barrier-Hessian oracle**, input `(t, ε)`, **not** the numerical dual; float-snap variant has **no** bit bound.
- Magron et al., *Exact SOHS decompositions of trigonometric univariate polynomials*, 2022 — arXiv:2202.06544 (the on-point tool for the Bochner block).
- Jansson, Chaykin & Keil, *Rigorous Error Bounds for the Optimal Value in SDP*, SIAM J. Numer. Anal. 46(1), 2007; VSDP-2012 (the mixed-cone a-posteriori bound).
- Powers; Scheiderer (boundary non-existence of rational SOS) — cited within 2105.11369.
- `alfonso` (Papp–Yıldız, arXiv:2101.04274) — nonsymmetric conic solver; **no** rational-certificate implementation.

**Prior-art verdict:** L1 is a *named, published* technique (true), but it is **already on the project's scouted list (un-assessed)**, and on inspection it is **subsumed** by D2/L4 (for the only block it can touch) and **dominated** by L2 (for the mixed cone as a whole). It is not the third independent route the proposal claims.

---

## 3. CONCRETE PLAN (if pursued anyway — as a Bochner-block finish, not the headline)

This plan honestly reflects that Davis–Papp can only touch the PSD/WSOS block. The whole-program certificate must come from L2 (or D1) for the LP/SOC parts.

1. **(Do L2 first.)** Port the Jansson/VSDP a-posteriori lower-bound formula and run it on the existing binding-center mixed-cone solve. This converts "numerically certified" → "verified" with **no reduction** and is the correct baseline. *(This is lever L2, not L1 — but it is the prerequisite that makes any L1 work meaningful.)*
2. **Obtain a numerical dual vector** for the Bochner block (re-solve the binding center, export the dual PSD matrix `Z`). *(Same step as D1/D2 — L1 does not avoid it.)*
3. **Choose the Bochner-block exact tool.** Prefer **Magron's trigonometric exact-SOHS (L4, arXiv:2202.06544)** over Davis–Papp here: it directly produces an exact weighted-sum-of-Hermitian-squares for a trig polynomial positive on the circle, which is what the Bochner block is, and has a concrete algorithm + complexity bound. Use Davis–Papp's **interval** algorithm only if recasting the block as univariate-on-`[−1,1]` is cleaner.
4. **Back off by `ε`.** Certify the Bochner block at a strictly-interior, ε-weaker point (Theorem 3.5 forces this at the optimum). Pick `ε` ≤ ~1% of the *White* margin so the block does not eat the `+2e-4` headline margin — and **first re-solve the binding center at `N ≥ 24 000`** (already recommended in `findings.md`) to widen the margin before snapping.
5. **Verify in exact arithmetic** (`sympy`/`mpmath`): check the rational certificate's PSD-ness (rational LDLᵀ / Fejér–Riesz factorization) and dual feasibility of the block.
6. **Combine with the LP/SOC verified bound from L2 (step 1)** to get a *single verified LB on the whole program*. The Bochner exact certificate alone does **not** certify the bound — it certifies one block of the dual.
7. **(Only then) consider Lean** for the finite rational PSD check (dovetails companion D5; `Matrix.PosSemidef` exists in Mathlib).

**Reality check on the plan:** steps 1–2 are L2/D1 work; step 3 is L4 work; steps 4–7 are the genuinely Davis–Papp-flavored part but apply to **one block**. So "pursuing L1" honestly means "do L2 + L4 + an ε-backoff exact factorization of the Bochner block." There is little that is *uniquely* Davis–Papp and *new* relative to the existing scout.

---

## 4. PAYOFF

- **Payoff class:** **1** (proof-grade / solver-independent certificate) *in principle*, but **only** if combined with L2/D1 for the LP/SOC blocks. Davis–Papp alone yields an exact certificate of **one PSD block of the dual**, which is **not** a bound on µ by itself.
- **"Truly meaningful" outcome it could produce (best case):** a single, solver-independent, rational/interval certificate that the *whole* mixed-cone dual is feasible at value `c`, with `c` strictly above White (and ideally above 0.380), making "µ ≥ c" a **theorem** with no floating-point solver in the trusted base. **But the Davis–Papp piece is a minor, replaceable component of that outcome; the heavy lifting is L2 (mixed-cone a-posteriori) and L4 (Bochner exact-SOHS).**
- **What it does NOT do:** it does not produce a stronger bound (class 3), nor an analytic structural theorem (class 2). The bound it certifies is **≤** the project's existing numerical bound by `ε`.

---

## 5. CHEAP FIRST-STEP PROBE — run, with finding

I ran three cheap probes (no heavy SDP solve):

1. **Cone-structure probe** (`build_problem(... bochner_n=8)`, `get_problem_data(CLARABEL)`): **42 SOC cones + 2049 linear inequalities + 2 PSD blocks** at tiny N=300. **Finding: the program is irreducibly mixed-cone, LP-dominated — not a univariate WSOS instance.** This *refutes the proposal's central premise.*
2. **Dual-extractor read** (`dual_extractor.py`): returns only the **dual objective scalar** (5-sig-fig log parse), **no dual vector**. **Finding: there is no extracted dual to "snap"; a re-solve to export the dual matrix is required regardless — Davis–Papp does not bypass that step.**
3. **ε-backoff arithmetic:** `log(1/ε)` bit cost; even `ε` = 1% of the White margin keeps the certified LB (0.380271) above White. **Finding: the bit-size / backoff is benign in isolation — so the obstruction is the cone class and the thin `+2e-4` margin over the project's own headline, NOT the certificate size.**

**Net first-step finding:** the premise is false at the cone level; the backoff is cheap; the unique-to-Davis–Papp content is small and replaceable.

---

## 6. FAILURE MODES (honest)

1. **(Primary, confirmed) Wrong cone class.** White's dual is mixed LP+SOC+PSD, LP-dominated; Davis–Papp certifies one univariate-interval WSOS cone. No reduction exists. The headline thesis is infeasible. *Mitigation:* restrict to the Bochner block — but then it's D2/L4, and the LP/SOC blocks still need L2/D1.
2. **(Confirmed) Boundary optimum.** The bound is an optimum (active constraints, complementary slackness); Davis–Papp cannot certify the boundary (Scheiderer). *Mitigation:* certify an ε-weaker `c < c*` — accepting a strictly weaker bound, and consuming part of the thin margin.
3. **(Confirmed) Not a float-snap.** The bit-bounded Algorithm 1 is exact-arithmetic with a barrier-Hessian oracle; the "from the numerical dual" hybrid has no bit bound. *Mitigation:* none that recovers the proposal's "lightweight snap" framing — the work is a re-implementation, not a post-processing.
4. **Thin headline margin.** Binding margin over the project's own 0.3802838 is `+2e-4`; an ε-backoff certificate sits below that. *Mitigation:* re-solve at `N ≥ 24 000` to widen the margin first (findings.md), reducing the usable `ε`.
5. **Re-implementation risk.** No released code; the interval-WSOS barrier Hessian and continued-fraction rounding are non-trivial to get exactly right. *Mitigation:* prefer Magron L4 (concrete algorithm) for the Bochner block.
6. **Dominated by alternatives.** L2 (VSDP) gives a *verified* mixed-cone bound now, with no reduction; L4 gives the Bochner exact-SOHS with a published algorithm. Davis–Papp adds little unique value. *Mitigation:* none — this is the reason to deprioritize.

---

## 7. Verdict

**MAYBE → re-scope; do not pursue as stated.** The headline ("apply Davis–Papp to White's binding dual, univariate-interval featured case") is **infeasible**: the dual is a mixed LP+SOC+PSD point (LP-dominated), not a univariate WSOS instance; the bound is a **boundary optimum** Davis–Papp cannot certify rationally; the bit-bounded method is **not** a float-snap; and the dual **vector** isn't even extracted today. The one block Davis–Papp could touch (Bochner/Toeplitz) is **already covered** by companion **D2** and external **L4 (Magron, with a concrete algorithm)**, and the mixed-cone bound as a whole is better served first by **L2 (Jansson/VSDP)**. **Recommendation: fold any Davis–Papp use into D2/L4 as a Bochner-block symbolic finish, sequenced after L2 — and reallocate the "L1 as a third independent route" framing, which does not hold up.**

---

## Sources (verified by full-text extraction)

- Davis & Papp, *Dual certificates and efficient rational SOS decompositions for polynomial optimization over compact sets*, SIAM J. Opt. 2022 — arXiv:2105.11369 (interior-only; Theorem 3.5; dual certificate = vector).
- Davis & Papp, *Rational Dual Certificates for WSOS Polynomials with Boundable Bit Size*, 2023 — arXiv:2305.19039 (bit size `O(…+log(1/ε))`; Algorithm 1 exact-arithmetic + barrier-Hessian oracle; float variant unbounded).
- Magron et al., *Exact SOHS decompositions of trigonometric univariate polynomials*, 2022 — arXiv:2202.06544.
- Jansson, Chaykin & Keil, *Rigorous Error Bounds for the Optimal Value in SDP*, SIAM J. Numer. Anal. 46(1), 2007; VSDP-2012.
- Papp & Yıldız, *alfonso*, arXiv:2101.04274 (generic conic solver; not the rational-cert algorithm).
- In-repo: `docs/RND_WHITESPACE/2026-06-03_external_literature_scan.md` (L1–L4), `2026-06-03_whitespace_scout.md` (D1–D5), `docs/archive/PRO6_COMPLEMENTARITY_PROOF.md`, `lp_research_state/code/white_full_convex.py`, `dual_extractor.py`, `bochner.py`.
