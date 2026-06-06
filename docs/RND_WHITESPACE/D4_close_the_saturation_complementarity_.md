# D4 — Close the saturation / complementarity result to an UNCONDITIONAL analytic theorem

**Author:** Claude (machine-assisted), whitespace deep-dive workflow
**Date:** 2026-06-03
**Scope:** Adversarial assessment of D4 (proposed in the 2026-06-03 scout, `docs/RND_WHITESPACE/2026-06-03_whitespace_scout.md`). SCOUTS and PLANS; certifies nothing.
**Verdict (one line):** **DROP the proposed framing (sub-target a), MAYBE on a reframed sub-target.** The named lemma D4 wants to prove (`|ξ| ≤ C·Ω`) is, even if proven exactly, **insufficient to make the ceiling theorem non-vacuous** — it is the wrong lever. The genuine open lemma (an N-uniform a-priori bound on `Σ m·λ_m`) is unattacked, hard, and the empirical data hints it may be **false**.

---

## 0. TL;DR for the orchestrator

D4's thesis: finish PRO-6 + LEVER-I-prime into an unconditional ceiling theorem ("the augmented SDP can't exceed C_∞ ≈ 0.380558 as N→∞"), blocked on (a) a-priori shadow-price bounds `|ξ|≤Ξ, τ≤T, ν₃≤V`, and (b) strict complementarity `r_CB ≤ max(r_C,r_B)` via the Bochner dual Z.

**Three findings, all adversarial:**

1. **The shadow-price route (sub-target a) does not reach the goal even if fully solved.** Theorem 2 of `LEVER_I_PRIME_THEOREM.md` bounds `Σ λ_m` (UNWEIGHTED). The ceiling needs `Σ m·λ_m` (m-WEIGHTED). The only a-priori bridge is `Σ m·λ_m ≤ 2R·Σ λ_m`, a **factor-20 (=2R) blowup**. Plugging the *exact empirical* `|ξ|=1.46Ω, τ=0.70Ω` into Theorem 2 gives `Σ m·λ_m ≤ 27.5`, i.e. a ceiling of **0.381209 at N=40,000 — ABOVE Together's UB 0.380871**. **Vacuous.** Proving `|ξ|≤C·Ω` rigorously buys *nothing* for the ceiling. (Probe 1 + 2 below; matches PRO-14 §"corrected ceiling formula": `Σ m·λ_m ≤ 72.8Ω` worst-case, 12× the empirical.)

2. **The scale-invariance argument D4 proposes for its first step is blocked** by exactly the obstruction PRO-14 already documented: under `f→s·f`, three constraints fail to scale — the anchor `c[0]∈[p₁,p₂]` (p fixed per row), `sum_squares(c)+sum_squares(d)≤0.5` (scales as s²), and the box `|c|,|d|≤2/π`. The rescaled problem is a *different, generally infeasible* row, so the envelope theorem doesn't transport `|ξ|/Ω`. (Probe 2.)

3. **The genuine open lemma is an N-UNIFORM bound on `Σ m·λ_m`, and the data suggests it may be FALSE.** `LEVER_I_PRIME_FINAL.md` §9 documents that the multipliers GROW with N at a near-linear rate (`Σ m·λ` sup-row: 5.97 at N=30K → 10.32 at N=40K). The residual `(π/2N)·Σ m·λ_m` therefore decays *much slower than 1/N*, and **`C_∞ ≈ 0.380558` is an extrapolation the authors themselves flag as unproven**. If `Σ m·λ_m` grows like `N^a` with `a≥1`, the residual does NOT vanish and there is **no finite ceiling at all** — the "ceiling theorem" target may be ill-posed. Neither sub-target (a) nor (b) touches this. (Probe 3.)

**Net:** D4 as scoped attacks a lemma that doesn't deliver the goal, via a method already known to be blocked, while the real blocker is elsewhere and possibly insurmountable. The honest current state — *"verifiable-per-solve, non-vacuous at measured N∈{30K,40K}"* (LEVER_I_PRIME_FINAL) — is already publishable and is NOT improved by the proposed work.

---

## 1. What is actually DONE (and solid)

- **Theorem 1 (KKT identity), `LEVER_I_PRIME_THEOREM.md` §2:** at an interior-`w₁` KKT point, `Σ λ_m^cos·α_m^-(1) = -2ξ + α₂^+(1)τ + 2Lν₃ + Δ_sin(1)`. Numerically verified to 0.4% at row4. **Exact algebra, stands.**
- **Tautological identity (‡), `PRO6_COMPLEMENTARITY_PROOF.md` §2:** `r_CB(N|n) = r_C(N) + r_B(N|n)` — this is an EQUALITY (set identity `K₀(n)∩F_B(N) = K₀(N)`), stronger than the conjectured inequality. **Done, rigorous, but tautological** (decomposition, not a bound).
- **Operative residual bound (Route B), `LEVER_I_PRIME_THEOREM.md` §3.3 + §3.1:** `ResidualGain ≤ (π/2N)·Σ m·λ_m^cos + O(1/N³)`, with the *corrected* per-m Case-A bound `πmL²/4` (PoC's Case-B `2L` was retracted, was 10⁸× too loose). **Rigorous given the measured `λ_m`.**
- **Non-vacuous saturation theorem, `LEVER_I_PRIME_FINAL.md`:** at N=30,000 and N=40,000, all 4 representative rows give `C_explicit < 0.380871` (sup `C_explicit = 0.380713` at N=40K, margin +1.58×10⁻⁴). This is a genuine, first-of-its-kind result: *the cell-envelope augmentation provably cannot match Together's UB at these N*. **CONDITIONAL on the measured `Σ m·λ_m` (verifiable-per-solve), and on the 4-row empirical sup.**

The gap between "done" and "D4's goal" is **purely**: turning the *verifiable-per-solve* `Σ m·λ_m` into an *a-priori, N-uniform* bound. That is the whole ballgame, and D4 mis-locates it.

---

## 2. Feasibility (adversarial)

### 2.1 Sub-target (a): a-priori shadow-price bounds — TRACTABLE BUT POINTLESS

**Is `|ξ|≤C·Ω` provable?** The empirical ratio is strikingly stable (Probe 2): `|ξ|/Ω = 1.456 ± 0.010`, 1.5% spread across 4 disparate centers. So D4's *premise* (it looks like a structural constant) is empirically sound. But:

- **The proof route (scale invariance) is blocked.** PRO-14 §A already tried `f→s·f` and found it breaks the `c[0]∈[p₁,p₂]` anchor. Probe 2 confirms two MORE non-scaling constraints (`sum_squares≤0.5`, box `2/π`). No clean group action survives.
- **Even granting `|ξ|≤1.46Ω` for free, the ceiling is vacuous** (Finding 1, Probe 1). Theorem 2 → `Σ λ_m ≤ 1.38`, then `Σ m·λ_m ≤ 2R·1.38 = 27.5`, giving ceiling **0.381209 @ N=40K ≫ UB**. The 2R blowup is the killer: the empirical mass sits at `m∈{3,4,5,6}` (PRO-14, LEVER_I_PRIME_THEOREM Obs 1), NOT at `m=2R=20`, but `2R·Σλ` assumes the worst. **Closing this requires bounding the m-DECAY, which is a spectral fact about the dual, not a shadow-price consequence.**

**Verdict on (a): tractable-but-useless.** Do not spend effort here.

### 2.2 Sub-target (b): strict complementarity `r_CB ≤ max(r_C,r_B)` via Bochner dual Z — HARD AND ORTHOGONAL

- PRO-6 §3 already shows the strict form holds *empirically within solver noise* (+2×10⁻⁵ apparent violation). Proving it would require a KKT-coupling argument extending Theorem 1 to the Bochner dual matrix Z (PRO-6 §7 step 2, "~1 week of math, deferred").
- **But (b) is about the cell-env↔Bochner residual SPLIT, not about the absolute size of `Σ m·λ_m`.** Even with strict complementarity proven, you still need to bound `r_C(∞) = (π/2N)·Σ m·λ_m` as N→∞ — same blocker as (a). Strict complementarity tells you the joint residual ≈ max of the two pieces; it does NOT tell you that piece converges.
- **Established theory (Nie, arXiv:1701.01549; matrix-poly version arXiv:2506.12579) gives a-priori multiplier bounds ONLY under strict complementarity + 2nd-order sufficiency + CQ at every global minimizer** — and PRO-23 (`LEVER_FUNCTIONAL_EQUATION.md`) showed the optimum is **massively degenerate** (437-shift active set, KKT residual 10⁻², not a strict/nondegenerate point). So the textbook machinery for a-priori multiplier bounds is itself blocked by the documented degeneracy.

**Verdict on (b): real math (~weeks), and even if done it does not deliver the N-uniform `Σ m·λ_m` bound.**

### 2.3 The REAL target: N-uniform bound on `Σ m·λ_m` — POSSIBLY FALSE

Probe 3 (using LEVER_I_PRIME_FINAL's own §9 data): `Σ m·λ` grows 5.97→10.32 from N=30K→40K. The precise exponent is noisy (row-instability: row4 is 10.74 at N=30K while row7 is 5.97), so the literal "N^1.90" is an over-read — but the *qualitative* signal (multipliers grow, decay is slower than 1/N) is **explicitly documented by the original authors** ("Σ m·λ and Σ m·σ continue to grow with N… the (π/2N) decay barely outpaces this"). Consequences:

- If `Σ m·λ_m → const`, the ceiling `C_∞` exists and the operative theorem (Route B) becomes unconditional once you prove that constant a-priori. **Nobody has.**
- If `Σ m·λ_m ~ N^a, a≥1`, `ResidualGain` does NOT vanish — **there is no finite framework ceiling**, and the "C_∞ ≈ 0.380558 theorem" is *ill-posed*, not just unproven.
- The data is in the ambiguous regime. **Resolving WHICH is the actual whitespace** — and it is genuinely hard (it asks for the asymptotic dual-multiplier growth rate of a degenerate SDP family).

---

## 3. Prior art (internal + external)

**Internal — heavily worked, all the pieces D4 cites are already in the repo:**
- `LEVER_I_PRIME_THEOREM.md` — Theorem 1 (KKT identity), Theorem 2 (verifiable shadow-price bound), Theorem 3 (conditional saturation), §2.2 explicitly flags `|ξ|,τ,ν₃` bounds "unproven, conjecture non-trivial."
- `PRO14_SHADOW_PRICE_AUDIT.md` — **disproved `|ξ|≤Ω`** empirically (`|ξ|≈1.46Ω`); built `_pro14_verifier.py`; §"corrected ceiling formula" already computed `Σ m·λ ≤ 72.8Ω` worst-case = vacuous. **D4 sub-target (a) is PRO-14's already-closed task, reopened.** PRO-14's own recommendation was "(B) verification-based, achievable now; (A) shadow-price proof remains research."
- `PRO6_COMPLEMENTARITY_PROOF.md` — tautological identity done; strict form deferred as "~1 week math."
- `LEVER_I_PRIME_FINAL.md` — non-vacuous theorem at N∈{30K,40K}; §9 documents the multiplier N-growth (the real obstruction).
- `LEVER_FUNCTIONAL_EQUATION.md` (PRO-23) — the optimum is degenerate (no strict complementarity), which **blocks the textbook a-priori-multiplier machinery**.

**External (WebSearch, 2026-06-03):** Nie, *Tight Relaxations for Polynomial Optimization and Lagrange Multiplier Expressions* (arXiv:1701.01549, Math. Prog. 2018) and the matrix-polynomial extension (arXiv:2506.12579, 2026) give Lagrange-multiplier-expression hierarchies with **a-priori bounds, but conditional on strict complementarity + 2nd-order sufficiency + CQ at every global minimizer** — hypotheses PRO-23 shows fail. No external work provides N-uniform multiplier bounds for a *degenerate* SDP family of this type. `arxiv_search.py` was HTTP-429 throttled (consistent with the same-day external scan note); WebSearch used instead.

---

## 4. Concrete plan (IF pursued anyway — reframed to the real target)

Do **not** execute sub-target (a) as written. If the orchestrator insists on a D4-flavored attack, the *only* version that could clear the bar is the reframed one:

1. **Settle the N-asymptotics of `Σ m·λ_m` empirically and cleanly first** (cheap, decisive). Re-extract `λ_m` and `σ_m` at a *fixed binding center* (cde_n30_iter3, the headline witness) across N ∈ {20K, 40K, 80K, 160K} with IDENTICAL config, to kill the row-instability noise in LEVER_I_PRIME_FINAL. Fit `Σ m·λ_m` vs N. **Decision gate:** if it plateaus → a finite ceiling exists, proceed; if it grows with exponent ≥1 → the ceiling theorem is ill-posed, **abandon D4 entirely and report the negative structural result** (itself meaningful: "the framework has no finite ceiling; the apparent saturation is a finite-N artifact").
2. **Only if plateau:** attempt an a-priori bound on `Σ m·λ_m` directly via the m-decay structure. The lever is NOT shadow prices — it is the spectral concentration at `m∈{3,4,5,6}`. Try to show the decay is a KKT consequence of the cosine-kernel curvature (PRO-14 §C2 "show the m-decay is itself a KKT consequence" — never attempted). This is the genuine open math.
3. **Sub-target (b) only as a sweetener:** extend Theorem 1 to include the Bochner dual Z (PRO-6 §7) to get the strict complementarity statement — publishable as a structural corollary, but it does NOT unblock the ceiling and should not be the headline.
4. Cross-check any analytic `Σ m·λ_m` bound against `_pro14_verifier.py` / `_lever_i_prime_lambda_m_all_rows.py` at ≥3 centers to 10+ digits.

**Effort:** Step 1 is ~1 day of solves (N up to 160K is memory-heavy, ~10–16 GB; feasible). Steps 2–3 are weeks of real math with **low success odds** (the textbook route is blocked by degeneracy; the bespoke m-decay-KKT argument is unattempted and speculative).

---

## 5. Payoff

- **Payoff class as proposed (sub-target a):** would-be class 2 (analytic theorem). **Actual delivered payoff: ~0** — proving `|ξ|≤C·Ω` leaves the ceiling vacuous; it does not produce a "truly meaningful" result.
- **Payoff class of the reframed target:** class 2 IF (and only if) `Σ m·λ_m` plateaus AND its limit is bounded a-priori → unconditional ceiling theorem `C_∞ < UB` = a clean negative structural result quantifying how much gap is beyond the framework. That clears the bar (class 2).
- **Consolation payoff (decision gate step 1):** if `Σ m·λ_m` grows, the finding "the framework has NO finite ceiling; apparent saturation is a finite-N artifact" is *also* a meaningful structural result (refutes the C_∞ ≈ 0.380558 folklore) — arguably class 2, and cheaply obtained.

**What "truly meaningful" would look like:** an unconditional statement "no augmentation in the cell-envelope + Bochner + poly-moment class can prove `µ ≥ C_∞` for an explicit `C_∞`, for ALL N" (or its refutation). NOT a proof of `|ξ|≤C·Ω`.

---

## 6. First-step probe — RUN, with findings

Three cheap probes run (no heavy SDP; reused cached `pro14_shadow_prices.json`, `lambda_m_scaled.json`):

**Probe 1 — two routes to the ceiling (`/tmp/_d4_probe.py`):**
| Route | Σ m·λ_m | ceiling @ N=40K | vs UB 0.380871 |
|---|---|---|---|
| A: a-priori shadow-price (Thm 2 + exact empirical ξ,τ) | 27.5 | **0.381209** | **+3.4×10⁻⁴ (VACUOUS)** |
| B: direct empirical Σ m·λ_m = 6.0 | 6.0 | 0.380364 | −5.1×10⁻⁴ (non-vacuous) |

Route A is **4.6× looser** than Route B. *Even with the exact empirical `|ξ|` plugged in*, the shadow-price route gives a ceiling above the UB. **Sub-target (a) cannot deliver a non-vacuous ceiling.**

**Probe 2 — `|ξ|/Ω` stability + scale-invariance obstruction (`/tmp/_d4_xi_probe.py`):**
- `|ξ|/Ω = 1.456 ± 0.010` (1.5% spread, 4 centers) — empirically a near-constant, as D4 claims.
- BUT `f→s·f` breaks `c[0]∈[p₁,p₂]` (p fixed), `sum_squares(c,d)≤0.5` (∝s²), box `|c|,|d|≤2/π`. **Scale-invariance argument blocked** — same obstruction PRO-14 documented.
- AND Theorem 2 bounds `Σ λ_m`, not `Σ m·λ_m`; the 2R-bridge loses everything. The m-concentration is spectral, not a shadow-price fact.

**Probe 3 — N-growth of `Σ m·λ_m` (`/tmp/_d4_smlambda.py`):**
- LEVER_I_PRIME_FINAL §9: `Σ m·λ` sup-row 5.97 (N=30K) → 10.32 (N=40K). Naive fit `~N^1.9` (noisy due to row instability, but the *growth* is real and author-documented).
- If growth exponent ≥1, `ResidualGain` does not vanish → **no finite ceiling**; `C_∞ ≈ 0.380558` may be an extrapolation artifact. The real open lemma (N-uniform `Σ m·λ_m`) may be **false**.

---

## 7. Failure modes (honest)

1. **The headline failure (confirmed, not hypothetical):** sub-target (a) is the wrong lever — proven `|ξ|≤C·Ω` ⇏ non-vacuous ceiling (2R blowup). D4 as written cannot reach its stated goal. *This is not a risk; it's a computed fact (Probe 1).*
2. **Scale-invariance blocked (confirmed):** the proposed first step doesn't get off the ground (Probe 2); PRO-14 already hit this wall.
3. **The real target may be ill-posed:** if `Σ m·λ_m` grows unboundedly with N (Probe 3 hint), there is no `C_∞` to prove a theorem about. Then the *only* deliverable is the negative "no finite ceiling" observation.
4. **Degeneracy blocks the textbook machinery:** Nie-style a-priori multiplier bounds need strict complementarity, which PRO-23 shows fails (437-shift degenerate optimum). The bespoke alternative (m-decay-as-KKT, PRO-14 §C2) is unattempted and speculative.
5. **Even full success is "only" class 2 and below the UB:** a ceiling theorem at `C_∞ ≈ 0.3806` does NOT close the gap to `µ` — it bounds the *framework*, not `µ`. It is meaningful (tells you to abandon D1/D2/D5's hope of reaching the UB from inside the framework) but does not itself produce a stronger bound on `µ`. The higher-leverage bar-clearing work remains D1 (rigor upgrade of the existing 0.380284 to proof-grade) per the same-day scout.

---

## 8. Recommendation

**DROP** sub-target (a) (proven useless here) and **DROP** the scale-invariance first step (blocked). **MAYBE** the reframed target, gated on a cheap decisive experiment: a single clean N-asymptotic extraction of `Σ m·λ_m` at a fixed center across N∈{20K…160K}. That ~1-day experiment either (i) shows a plateau → opens a real (hard, weeks-long, low-odds) class-2 attack on an N-uniform `Σ m·λ_m` bound, or (ii) shows growth → yields a cheap negative structural result ("the framework saturation ceiling is a finite-N artifact; no finite C_∞"), which is itself worth recording.

**Priority vs siblings:** D4 is the *weakest* of the scouted directions for clearing the bar. The same-day scout's D1 (cvxpy→SDPA-S serializer → GMP re-certification → rational dual) directly upgrades the *existing* 0.380284 to proof-grade (class 1) and is engineering, not blocked open math. Prefer D1. Use D4's decision-gate experiment only as a cheap side-probe to retire (or sharpen) the `C_∞ ≈ 0.380558` folklore.

---

## Appendix — artifacts

- Probes: `/tmp/_d4_probe.py`, `/tmp/_d4_xi_probe.py`, `/tmp/_d4_smlambda.py` (throwaway; reproduce from cached `lp_research_state/data/{pro14_shadow_prices,lambda_m_scaled}.json`).
- Key prior docs: `docs/archive/{LEVER_I_PRIME_THEOREM,LEVER_I_PRIME_FINAL,PRO6_COMPLEMENTARITY_PROOF,PRO14_SHADOW_PRICE_AUDIT,LEVER_FUNCTIONAL_EQUATION}.md`.
- External: Nie arXiv:1701.01549; matrix-poly extension arXiv:2506.12579.
