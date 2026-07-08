# Session synthesis — 2026-05-10

**Going-in state:** µ ≥ 0.3801279 (Phase 5 of CDE, saturated technique stack on Fourier basis at SDP scale N=10000 / T=4000 / bochner_n=30). Together's UB µ ≤ 0.380871. Open gap 7.43 × 10⁻⁴.

**Going-out state:** µ ≥ 0.3801279 unchanged. **What changed is our understanding of why pushing the LB is hard.** Three candidate levers were investigated and three were definitively ruled out as practicable next steps.

## The headline conclusion

The Erdős minimum overlap problem at the current SDP-relaxation regime has hit a wall that is not crossed by any of: more PSD constraint families, basis changes, switching to M-side encoding, or restricting to low-block-count step functions. Pushing the LB past 0.3801279 toward Together's 0.380871 requires *qualitatively* different mathematics — not more of the current style.

## What was investigated and what was learned

### 1. Together-as-primal SDP diagnostic ([TOGETHER_DIAGNOSTIC.md](TOGETHER_DIAGNOSTIC.md))

12-task implementation that loaded Together's piecewise-constant minimizer h\*, projected it into White's Fourier basis at T=4000, and read off every constraint slack in our Phase 5 SDP at row 4.

Key findings:
- **Ω_SDP(f_even, pinned at h\*) = 0.459311** — much larger than f_even's actual autocorrelation 0.387337 and much larger than Phase 5's 0.380128. The SDP's cell-kernel autocorrelation envelope is loose by 0.07 on this primal.
- **All three CDE augmentations (Bochner-PSD, poly_moment, Hankel-PSD) are healthily slack at Together's f.** None is binding. The binding constraint family is the *original* White cell-kernel envelope (`white_full_convex.py:176-190`), untouched by the CDE work.
- **The Phase 5 SDP-optimal f̃ is no longer Gibbs-dominated** (min −0.97, not the historical −3.78). The gap function `f̃ − f_even` is **99.9% low-frequency** — the deviation is smooth and structural, not high-frequency ringing.
- f_direct (asymmetric embedding) is **infeasible** in White's SDP, rejected by Bochner-PSD on (1−f) with λ_min = −0.88.

Original recommendation: pivot to **Lever D — structural restriction to O(1)-breakpoint step functions**.

### 2. Lever D pre-investigation ([LEVER_D_PRE_INVESTIGATION.md](LEVER_D_PRE_INVESTIGATION.md))

Tested the O(1)-breakpoint hypothesis directly. **REFUTED:**
- Together's h\* has 399 blocks at tolerance 10⁻³, 89 at 10⁻¹. Not low-complexity at any reasonable tolerance.
- 58% of h\*'s cells lie strictly in (0.05, 0.95) — h\* is *smoothly varying*, not bimodal toward {0, 1}.
- Discrete optimizer block counts grow ≈ 0.55·n + 2.3 (linear, not O(1)).

Pivot recommendation: reformulate Lever D as Lipschitz/BV restriction, or pivot to Lever E.

### 3. Lever decision — lifted discrete limits ([LEVER_DECISION.md](LEVER_DECISION.md))

Pushed exact M(n) brute force from n=12 to **n=18** (M(18) = 8) using symmetry-reduced branch-and-bound. Lifted each optimal discrete partition into a smoothed continuous density f_n and compared to a folded-to-[0,1] version of h\*.

**Result: DIVERGE.** L¹ distance plateaus at ≈ 0.25 across n=4..18, no decreasing trend; Pearson correlation stays at 0.5–0.65. The discrete optima (necessarily {0,1}-valued indicators) cannot represent h\*'s fractional structure (58% of cells in the strict interior). They live in different function classes.

This refutes both the Lipschitz reformulation of Lever D and motivates **Lever E — M-side SDP encoding**.

### 4. Lever E pretest ([LEVER_E_PRETEST.md](LEVER_E_PRETEST.md))

Three cheap experiments before committing weeks to building a new M-side SDP. **NOT PROMISING:**
- M(f̃) computed five different ways, all in [0.59, 0.71]. Our Phase 5 SDP-optimal f̃ is a *terrible* primal for Together's M-functional. Not a tighter UB.
- Existing `mside_bochner_n=8` at full Phase 5 scale: SDP value drops from 0.379653 to 0.379653 (Δ = −1.6 × 10⁻⁷, solver noise). The SOC-relaxed M-side constraint is *empirically vacuous*.
- Prior work (findings.md): `mside_bochner_n=5, 10` already known dead at Δ ≈ 10⁻⁹–10⁻⁸. The only exact M-side lift requires the retracted Lasserre level-2 approach.

## What is now ruled out as a near-term lever

| Lever | Status | Why |
|---|---|---|
| A — Lukács SOS / alt basis | ❌ Unlikely | Gibbs is already damped; gap is structural low-frequency |
| D — O(1)-breakpoint restriction | ❌ Refuted | h\* has 400+ blocks at fine tol, optimal block count grows linearly in n |
| D' — Lipschitz via discrete limit | ❌ Refuted | f_n and h_folded diverge; different function classes |
| E — M-side SDP encoding | ❌ Not promising | Convex M-side relaxations are vacuous; exact lift needs retracted Lasserre |

## What remains worth attempting

### Lever C — push exact M(n) to n=25-43 via SAT/IP

Scope: 1-2 weeks. Output: tighter upper bound on µ (some M(n)/n value < 0.380871, if it exists), plus better data on the limiting structure of discrete optima at scale.

Pros: clean infrastructure exists in [`lp_research_state/code/brute_force_Mn.py`](lp_research_state/code/brute_force_Mn.py) and [`_brute_force_Mn_extended.py`](lp_research_state/code/_brute_force_Mn_extended.py); strictly bounded scope; informative regardless of outcome.

Cons: doesn't directly improve the LB. Doesn't bridge to Together's 0.380871 unless their certificate is suboptimal.

### Speculative: structural theorem connecting the two formulations

What would actually close the gap is a *theorem* connecting White's Ω-functional and Together's M-functional at the level of µ. Specifically: a proof that

**inf_{f admissible} Ω(f) = inf_{h admissible} M(h) = µ**

is one thing (definitional, multiple proofs known). But for SDP purposes we want a *primal-side* connection: given an h with M(h) ≤ µ + ε, construct an f with Ω(f) ≤ µ + δ(ε). If we had that map, we could plug Together's h\* into our SDP via the bridge and get an LB approaching 0.380871.

The diagnostic suggests no such map exists with small δ — Ω(f_even) = 0.459 even though M(h\*) = 0.380871, so the natural bridges give very loose bounds.

This is a research mathematics question, not an engineering task.

### Stop and write up

The published-quality result µ ≥ 0.3801279 stands. The CDE Phase 1-5 result has not been written up as a paper yet — `communications/preprint_draft.tex` exists per CLAUDE.md. Time invested in writing up could be more valuable than further LB-pushing at the current technique frontier.

## Recommendation

I'd push you toward **two parallel tracks**, neither of which is a multi-week SDP-engineering investment:

1. **Track A (low-effort, 1-2 days):** Run Lever C — brute-force / SAT-based M(n) computation to n=25-30 to see if Haugland's ratio M(n)/n stays > 0.380871 (validating Together's certificate as a UB) or dips below (invalidating it). The infrastructure exists; this is mostly runtime.

2. **Track B (medium-effort, 1-2 weeks):** Write up the published-quality result. The diagnostic also surfaced a clean piece of new mathematical content: the **Ω-vs-M slack structure at the SDP scale**, with quantitative measurements of how far each functional is from the other on common primal points. This is publishable as a methodological note even without LB improvement.

If you want a stretch goal: investigate whether there's a *primal-side bridge* between the two formulations that's tighter than the trivial Cauchy-Schwarz estimates. That's the only path I can see to closing the [0.3801279, 0.380871] gap at this stack scale.

## Reproducing the session

Everything is committed on `main` (commits `cbe7978..9009a96`). Not pushed.

```bash
git log --oneline f5acb9e..HEAD     # 15 commits
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import aggregate_results; aggregate_results()"
```

Memos to read in order:
1. [TOGETHER_DIAGNOSTIC.md](TOGETHER_DIAGNOSTIC.md) — the main diagnostic
2. [LEVER_D_PRE_INVESTIGATION.md](LEVER_D_PRE_INVESTIGATION.md) — D refutation
3. [LEVER_DECISION.md](LEVER_DECISION.md) — D' refutation, pivot to E
4. [LEVER_E_PRETEST.md](LEVER_E_PRETEST.md) — E ruled out
