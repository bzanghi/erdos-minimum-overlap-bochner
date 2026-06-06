# D1 — cvxpy→SDPA-S serializer + arbitrary-precision (GMP) re-certification

**Date:** 2026-06-03
**Workflow:** R&D whitespace deep-dive (Phase 2). SCOUTS and PLANS; does **not** certify any bound.
**Author:** Claude (machine-assisted), adversarial assessment.
**Status of direction:** PURSUE — but **re-scoped**. The original "build serializer → GMP re-solve at production scale" is memory-blocked; the payoff-class-1 goal is reachable more cheaply via a route the original framing missed (Jansson verified post-processing). See verdict.

---

## 0. One-paragraph bottom line

The stated goal — turn the existing numerically-certified `µ ≥ 0.380284` into a **solver-independent** lower bound, removing reliance on CLARABEL's `reported − last_gap` log-parse trick — is genuinely worth pursuing and is **engineering, not open math**. But the specific mechanism in the proposal (serialize the production SDP to `.dat-s`, re-solve at GMP precision) **does not fit in memory**: at production `N=10000` the SDP scalarizes to ~36k variables / ~74k constraints and SDPA-GMP forms a **dense 52–73 GB GMP Schur matrix**. The correct, scalable mechanism is **Jansson–Chaykin–Keil verified-bound post-processing** (the method behind VSDP): a *single directed-rounding pass* over an approximate dual point yields a rigorous lower bound that accounts for all rounding errors, with no GMP re-solve and provably-finite cost because every primal block of White's program is box-bounded. GMP/SDPA then has a *narrow, valuable* role: crush the dual residual on the small **binding** instance so the Jansson correction term stays below the thin certification margin. Net: the keystone deliverable is a ~Jansson verifier in mpmath, not the full GMP serializer.

---

## 1. FEASIBILITY

### 1.1 SDPA-GMP is alive (after a trivial fix) — but the brief's "smoke-tested" status was stale
The binary `lp_research_state/bin/sdpa_gmp` is **dynamically linked** against `/tmp/sdpa_build/gmp-install/lib/libgmp.10.dylib`, which no longer exists (`/tmp` was wiped). `sdpa_gmp_wrapper.py` therefore returns all-`None` today (silent failure). It runs correctly with:

```
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib ./sdpa_gmp example1.dat-s out -p param.sdpa -pt 1
```

Homebrew has the exact GMP 6.3.0 build at `/opt/homebrew/lib/libgmp.10.dylib`. On `example1.dat-s` it reproduces the brief's claim: `p.feas.error = 2.3e-37`, relative gap `4.9e-31`, 30 digits. **Permanent fix:** `install_name_tool -change /tmp/...libgmp.10.dylib /opt/homebrew/lib/libgmp.10.dylib sdpa_gmp` (and the gmpxx dylib), or pass `env={"DYLD_FALLBACK_LIBRARY_PATH": "/opt/homebrew/lib"}` in the wrapper's `subprocess.run`. This is a 10-minute fix and a prerequisite for any GMP work. *(Flagged as a separate spin-off task.)*

### 1.2 The production SDP does NOT fit in GMP memory — the proposal's central mechanism is blocked
Scalarized cone dimensions (measured via `prob.get_problem_data(cp.CLARABEL)`):

| config | decision vars | eq rows | conic scalars | SOC cones | PSD blocks |
|---|---:|---:|---|---|---|
| N=1500,T=600,bn=10 | 5,483 | 11,560 | 9,729 | 40×`[3]` + 2×`[602]` | 2×`[22]` |
| N=5000,T=2000,bn=20 | 18,083 | 38,060 | 32,129 | 40×`[3]` + 2×`[2002]` | 2×`[42]` |
| **N=10000,T=4000,bn=20** | **36,083** | **74,060** | 64,129 | 40×`[3]` + 2×`[4002]` | 2×`[42]` |

SDPA-GMP forms a **dense** Schur-complement matrix of size `mDIM × mDIM` (`mDIM` = decision-var count) at GMP `precision 200` (≈ 40–56 effective bytes/entry):

- N=10000 → `36083²` dense GMP = **52–73 GB**. Infeasible on the 4 GB machine.
- N=5000 → **8–11 GB**. Also over budget.
- N≤~1500–2000 → ~**0.7–1.0 GB**. Tractable.

Two structural aggravators specific to SDPA's cone model:
1. **SDPA-S has no SOC cone — only LP (`R₊`) and PSD (`S₊`).** The two `sum_squares(c)+sum_squares(d)≤0.5`/box SOC cones (dimension `2T+1 ≈ 4002`) must be rewritten as PSD "arrow" blocks of size ~`(2T+1)`. That manufactures two ~4000×4000 *dense* PSD variables — extra cost the proposal didn't budget. (The 40 tiny `[3]` SOC cones from the per-`m` `cp.square` terms are cheap.)
2. The `w, v` variables (the `2N`-dimensional bulk) are "wide but shallow": they enter only ~25 linear moment functionals plus box bounds. A naive serializer pays the full `mDIM ≈ 2N` Schur cost for them; see §3.4 for the dual-reduction lever that avoids this.

**Conclusion:** "re-certify the production 0.380284 by GMP re-solve" is **not feasible** on available hardware. What is feasible at GMP precision is a **small instance** (N≤2000).

### 1.3 The small-instance bound at GMP-feasible scale
At `N=2000, T=800, bochner_n=20`, the binding center **row4** gives `primal_value_at_center = 0.37820` (from `parallel_results/path_b_N2000/row4.json`) — *below* White's 0.379005 at the bare center. The full-µ bound that clears White comes from the **ellipse-extension cover** (max over 7 row centers + Lipschitz grid term): `path_b.json` reports `coverage.grid_min_obj = 0.37957`, `headline.1e6_safety.min = 0.37955`. So **a GMP-certified small instance certifies White-beating µ only through the path_b cover machinery, not from a single solve.** This is the honest scope of a small-N GMP certificate: it would prove `µ ≥ ~0.3795` (beats White 0.379005), *not* `µ ≥ 0.3803` (which is intrinsically a large-N result).

### 1.4 The actually-correct, scalable tool: Jansson verified post-processing (VSDP)
This is the load-bearing prior-art the original D1 framing did not name. **Jansson, Chaykin & Keil, "Rigorous Error Bounds for the Optimal Value in Semidefinite Programming," SIAM J. Numer. Anal. 46(1):180–200 (2007)** — Theorem 3.2:

> Given any approximate dual point `ỹ`, set the **defect** (dual residual) matrices `Dⱼ := Cⱼ − Σᵢ ỹᵢ Aᵢⱼ`, compute rigorous lower eigenvalue bounds `dⱼ ≤ λ_min(Dⱼ)`, and let `xⱼ ≥ λ_max(Xⱼ)` be an a-priori upper bound on the primal block. Then
> **`p* ≥ inf{ bᵀỹ + Σⱼ sⱼ · dⱼ⁻ · xⱼ }`**, where `dⱼ⁻ = min{0, dⱼ}` and `sⱼ` is block dimension.

Properties that make this the right tool:
- **`ỹ` need not be dual-feasible.** The term `sⱼ·dⱼ⁻·xⱼ` exactly *absorbs* the infeasibility (this is what replaces the `reported − last_gap` trick with something rigorous).
- **Algorithm 3.1 terminates in ONE pass (step 3) when every `xⱼ` is finite.** Cost is "negligible compared to approximately solving the SDP." No GMP re-solve, no perturbation loop.
- **Every `xⱼ` IS finite and closed-form for White's program:** `0 ≤ wⱼ,vⱼ ≤ Ω ≤ 1` (so `xⱼ ≤ 1`); `‖(c,d)‖²≤0.5`; the Bochner real-form block is Toeplitz with diagonal `1/2`, hence `λ_max ≤ trace = nb+1 ≤ 21`. The bounded-primal hypothesis (the one nontrivial requirement) is **satisfied by construction**.
- Scales: VSDP "has rigorously solved problems up to 20 million variables"; median SDPLIB accuracy `2.2e-8`. It only post-processes a *float* solver's output (CLARABEL/MOSEK), so production `N` is fine.
- **Trusted base** shrinks to: directed-rounding arithmetic + a verified `λ_min` lower bound (Rump/Gershgorin-style, or a verified Cholesky) + the closed-form `xⱼ`. All implementable in `mpmath.iv` interval arithmetic; **no MATLAB/INTLAB strictly required** (VSDP uses INTLAB, but the algorithm is solver/language-agnostic).

**The catch — the quantitative tension (the real risk):** the correction term scales with the dual residual times the block size. For the two 42×42 Bochner blocks, the term is `42 · dⱼ⁻ · 21`. To keep the *total* correction under the binding full-space margin (`+1.1e-4` at R16 down to `+2e-5` at R6):
- residual `λ_min(D_Bochner) ≥ −1e-6` → correction `≈ −1.8e-3` per block → **blows the margin**.
- residual `≥ −1e-8` → correction `≈ −1.8e-5` → safe.

A typical CLARABEL solve has dual residual ~`1e-7..1e-6` — **insufficient for the binding row.** This is precisely where GMP earns its place: **SDPA-GMP on the small binding instance drives the residual to ~1e-30, making the Jansson correction ~1e-35 (utterly negligible).** The LP cone, by contrast, contributes ~`(#active-boundary vars)·residual` (the classical Neumaier–Shcherbina LP bound, where correct-sign reduced costs contribute 0), **not** `N·residual` — so the wide `N` is *not* the problem; the PSD-block residual quality is.

### 1.5 Feasibility verdict
- Jansson verifier (mpmath) on existing float duals: **tractable, days of work, scales to production N.**
- SDPA-GMP dylib fix: **trivial (10 min).**
- SDPA-GMP on the small binding instance (N≤2000): **tractable (~1 GB).**
- SDPA-GMP at production N: **infeasible (52+ GB).** Do not attempt.
- cvxpy→SDPA-S serializer with SOC→PSD arrow rewrite + dual bookkeeping: **medium (1–2 wks), but only needed for the GMP small-instance leg, not for the Jansson production leg.**

---

## 2. PRIOR ART

### 2.1 Inside this repo
- **`sdpa_gmp_wrapper.py`** explicitly flags the cvxpy→SDPA-S serializer as deferred-and-never-built ("the gating value is having the binary, not the interface"). Confirmed: no serializer exists; the binary handles only hand-written `.dat-s`.
- **`docs/RND_WHITESPACE/2026-06-03_whitespace_scout.md`** already lists this as direction **D1** (and D2 = exact Fejér-Riesz/SOS, D5 = small-N exact + Lean). This memo is the deep-dive on that D1; it adds the SDPA memory analysis and the Jansson reframe, which the scout did not have.
- **`dual_extractor.py`** is the current "central epistemic trick": it returns `dual_obj` at the last CLARABEL iteration with `dual_residual ≤ 1e-4` (`best_dual_lower_bound`, default `max_dual_residual=1e-4`). Crucially the threshold is **not zero** — even today's "rigorous_dual_LB" tolerates a 1e-4 per-constraint residual that is *not* rigorously absorbed. Jansson Thm 3.2 is exactly the principled replacement (it absorbs the residual via `sⱼ·dⱼ⁻·xⱼ` instead of hand-waving it).
- **`path_b_rigorous.py`** already documents a margin convention (`prob.value − 1e-6`) and admits the parsed dual is only "5 sig figs" exact. Same gap; same fix.
- **`docs/archive/FULLSPACE_RIGOR_MEMO.md`** (independent adversarial audit) confirms the core anchor `0.3802838` reproduces to `1e-11` across code paths — so the *number* is not in doubt; only its *certificate status* is. That is exactly what D1 targets.
- **`docs/archive/LEAN_LEMMA_INVENTORY.md`**: Mathlib lacks SDP duality / Bochner-Herglotz / SOS, so *general* formalization is 12+ mo. But `Matrix.PosSemidef` exists, so a *specific finite rational PSD certificate check* (the D2/D5 Lean target) is a much smaller, in-scope object.
- **MOSEK is installed** (`cp.installed_solvers()` → includes `MOSEK`). MOSEK's dual is far more reliable than CLARABEL's and is a cheap cross-solver feasibility check *before* any GMP work; it also reaches smaller residuals on the PSD blocks.

### 2.2 External literature (the methods D1/D2 are instances of)
- **Verified SDP bounds (D1's actual method):** Jansson–Chaykin–Keil (SIAM J. Numer. Anal. 2007); Jansson, *On verified numerical computations in convex programming* (Japan J. Indust. Appl. Math. 2009); **VSDP** toolbox (Jansson; vsdp.github.io). Companion LP method: **Neumaier & Shcherbina**, *Safe bounds in linear and mixed-integer programming* (Math. Prog. 2004). These are mature, peer-reviewed, and *exactly* our use case (rigorous LB from a float solver's output, all rounding accounted, bounded primal). **Decisive prior art that the bound-improvement is a known, validated technique.**
- **Exact rational SOS (D2's method):** **Peyrl & Parrilo**, *Computing sum of squares decompositions with rational coefficients* (SNC'07 / Theor. Comput. Sci. 2008) — numeric solve → rational rounding under a *strict-feasibility* assumption, with an explicit numerical-error-vs-rounding-tolerance relation. Extensions: Kaltofen–Li–Yang–Zhi (*Exact certification… via rationalizing SOS*, J. Symb. Comput. 2012); Magron–Safey El Din (*Dual certificates and efficient rational SOS*, SIAM J. Optim. 2021). The strict-feasibility/margin requirement is precisely the §1.4 tension: rational rounding needs the certified margin to exceed the rounding tolerance.
- **Erdős-µ-specific priority:** none of the above has been applied to the Erdős minimum-overlap SDP. The Together UB repo and White (2023) do not use verified-SDP post-processing. This is genuine whitespace *for this problem*, while resting on off-the-shelf, validated *method*.

---

## 3. CONCRETE PLAN (ordered)

**Track A — Jansson verifier (keystone; do first).**
1. **Fix the SDPA-GMP dylib** (`install_name_tool` or wrapper `env`); re-run `sdpa_gmp_wrapper.selftest()` and confirm `objValDual` parses. *(½ day; also spun off as a standalone task.)*
2. **Implement `jansson_lb(prob)` in mpmath interval arithmetic** operating on cvxpy's exported standard-form data (`get_problem_data`) + the solver's dual `ỹ`:
   - Build `Dⱼ = Cⱼ − Σᵢ ỹᵢ Aᵢⱼ` per block with directed rounding (`mpmath.iv`).
   - Rigorous `dⱼ ≤ λ_min(Dⱼ)` (verified Gershgorin first; upgrade to a verified Cholesky/Rump bound if Gershgorin is too loose on the Bochner block).
   - Closed-form `xⱼ`: `1` for `w,v` box; `√0.5`-derived for `(c,d)`; `nb+1` for each Bochner block; tail-bound values for `eps,dlt`.
   - Return `inf{ bᵀỹ + Σⱼ sⱼ·dⱼ⁻·xⱼ }` with outward rounding.
   - **Cross-check:** the LP-only part must reproduce the Neumaier–Shcherbina bound; sanity vs `path_b_rigorous`'s `value − 1e-6`.
3. **Validate at small N** where CLARABEL *and* SDPA-GMP both solve (N≤1500): run `jansson_lb` on the CLARABEL dual and on the GMP dual; confirm both are rigorous LBs and quantify how much the GMP residual tightens the correction term.

**Track B — small binding instance, GMP-grade certificate.**
4. **Serializer (cvxpy→SDPA-S)** for the small binding instance: emit LP + PSD blocks, **rewrite the SOC cones as PSD arrow blocks**, and the box/linear data. Validate the round-trip by checking SDPA-GMP's float-recomputed `objVal` matches CLARABEL on N≤1500 to printed digits.
5. **GMP-solve the binding center(s)** at N≤2000 (row4, and the thin gates R16/R17), residual ~1e-30; run `jansson_lb` on the GMP dual → a rigorous LB whose certificate is independent of any float solver's gap heuristic.
6. **Convert per-center GMP/Jansson LBs into a µ bound** via the existing `path_b_*` ellipse-extension cover (this is the step that lifts the sub-White center value to a White-beating `µ ≳ 0.3795`). Re-run `path_b_independent` with the verified `V_c` to keep the 10-digit cross-check policy.

**Track C — proof-grade hardening (optional, payoff toward "theorem").**
7. **D2/Peyrl-Parrilo:** snap the GMP dual on the (small) Bochner blocks to rationals, verify PSD by exact rational `LDLᵀ`; combine with rational LP duals (Neumaier–Shcherbina is already rational-friendly) → a fully rational, human/Lean-checkable certificate of the small-instance LB. Requires first **widening the binding margin** (N=24000 hardening, already recommended in `findings.md`) so the rational slack fits.
8. **D5/Lean (stretch):** formalize the *finite* rational PSD + feasibility check in Lean (`Matrix.PosSemidef`) for the smallest White-beating instance — "machine-checked `µ ≥ 0.3795`."

---

## 4. PAYOFF

- **Payoff class: 1** (proof-grade / solver-independent rigor).
- **"Truly meaningful" outcome, tiered by how far Track gets:**
  - *Track A alone* → the existing `µ ≥ 0.380284` (full-space) ceases to depend on CLARABEL's `reported − last_gap` log-parse; the LB is recomputed by a Jansson interval-arithmetic certificate that rigorously absorbs the dual residual. This **removes the project's single most-cited rigor caveat** ("central epistemic trick") and is *itself* the headline rigor upgrade — even if every other caveat (poly-moment, infeasibility exclusions) remains. **This is the high-value deliverable and it is achievable at production scale.**
  - *Track A+B* → a GMP-grade, residual-crushed certificate for the binding center, making the per-center value certified to ~30 digits; via path_b this yields a **solver-independent `µ ≳ 0.3795` (beats White 0.379005) with no float-gap heuristic anywhere in the trusted base.**
  - *Track C* → a **rational (and potentially Lean-machine-checked) certificate** — the literal word "theorem" on a White-beating bound.
- **What it does NOT do:** it does not push the *value* past the conjectured `C_∞ ≈ 0.380558` ceiling (that's a different, open-math axis), and Track A by itself does not upgrade the poly-moment-cut or infeasibility-exclusion caveats (those are D-other directions). It converts *rigor status*, not *strength*.

---

## 5. CHEAP FIRST-STEP PROBES (run; findings below)

1. **SDPA-GMP liveness.** Selftest returned all-`None`; root cause = stale `/tmp` rpath to `libgmp.10.dylib`. **Fixed in-session** with `DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib` → `p.feas.error 2.3e-37`, 30 digits, `pdOPT`. The "smoke-tested" claim is true *once the dylib path is repaired*; it is currently broken.
2. **Production memory model.** `get_problem_data` at N=10000 → 36,083 vars / 74,060 eq rows → dense GMP Schur **52–73 GB** ⇒ production GMP infeasible. N≤2000 fits (~1 GB).
3. **Small-instance scope.** `path_b_N2000/row4.json`: binding center at N=2000 = **0.37820** (sub-White at the bare center); `path_b.json` cover = **0.37955** (beats White only via the ellipse-extension machinery).
4. **Jansson `xⱼ` finiteness.** Verified every cone block has a finite closed-form `λ_max` bound (`w,v ≤ 1`; Bochner trace `≤ 21`; SOC slack `< 1.5`) ⇒ Algorithm 3.1 terminates in **one** pass — the cheap case.
5. **Quantitative margin tension.** Bochner-block correction `42·dⱼ⁻·21`: needs `λ_min(D_Bochner) ≥ −1e-8` to fit the `+1.1e-4..+2e-5` binding margin; CLARABEL's `~1e-7` residual is insufficient ⇒ GMP (or margin-widening) is *required* for the binding row, not optional.

---

## 6. FAILURE MODES (honest)

1. **Bochner-block residual vs margin (most likely failure).** If neither GMP-on-small-instance nor N=24000 margin-widening gets the Jansson correction below the binding margin, the *full-space* `0.380284` cannot be re-certified by this route — only a weaker constant (e.g. `0.380` or `0.3795`) clears. *Mitigation:* certify the weaker White-beating value (still a payoff-1 result); or widen the margin first.
2. **Verified `λ_min` looseness.** Gershgorin on the Bochner defect may be far looser than the true `λ_min`, inflating the correction. *Mitigation:* verified Cholesky / Rump's symmetric-eigenvalue enclosure (more code, standard).
3. **Production GMP is out (already established).** The proposal's literal "GMP re-solve the production SDP" fails on memory. *Mitigation:* this memo's reframe — Jansson at production, GMP only for the small binding instance.
4. **Serializer correctness (SOC→PSD).** A bug in the arrow-block rewrite or the dual sign/scaling bookkeeping would silently produce a wrong (possibly invalid) bound. *Mitigation:* the repo's standing policy — independent re-implementation agreeing to 10+ digits; round-trip-check SDPA's float `objVal` vs CLARABEL on N≤1500 before trusting GMP.
5. **Rational snap (D2) eats the margin.** Peyrl–Parrilo needs the certified margin > rounding tolerance; on a `+2e-5` margin the rational Gram may fail to stay PSD. *Mitigation:* margin-widening (N=24000) first; or keep the certificate at high-precision *interval* (Jansson) rather than exact-rational — still solver-independent, just not Lean-checkable.
6. **µ-conversion still needs path_b.** A per-center GMP certificate is not a µ bound until run through the ellipse-extension cover (which has its *own* Lipschitz-grid term, already rigorous). If that machinery has a latent issue it is inherited. *Mitigation:* path_b is independently re-implemented (`path_b_independent.py`) and audited (FULLSPACE_RIGOR_MEMO); reuse the cross-check.
7. **Diminishing returns on the *other* caveats.** Even a perfect Track A leaves the poly-moment-cut rigor and the solver-attested infeasibility exclusions (D3) untouched; a reader could still say the *full-space* bound isn't fully proof-grade. *Mitigation:* sequence D1→D3; D1 is necessary-not-sufficient but is the largest single rigor lever.

---

## 7. Verdict

**PURSUE — re-scoped.** The keystone is a **Jansson verified-bound post-processor in mpmath** (production-scale, removes the `reported − last_gap` trick), with **SDPA-GMP used only on the small binding instance** to crush the PSD-block residual below the certification margin. The literal "GMP re-solve at production scale" in the original framing is memory-blocked and should be dropped. Highest-leverage rigor upgrade in the project; rests on mature, peer-reviewed methods (Jansson/VSDP, Neumaier–Shcherbina, Peyrl–Parrilo) never applied to this problem.
