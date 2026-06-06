# PRO-25 Literature Mine — μ is not a renamed constant, but Rechnitzer technique-transfer is a strong lead

**Status:** Done. The Erdős minimum overlap constant μ is **not** known under a different name. Comprehensive search across the additive-combinatorics neighborhood produced no numerical collision. **Critical byproduct:** Rechnitzer (2026) computed 128 digits of an adjacent constant using a non-SDP method that may transfer to μ.

## Top-line conclusion

After surveying every plausible "dual" or "sibling" constant in additive combinatorics — autoconvolution L^p norms (Schinzel–Schmidt, Martin–O'Bryant, Matolcsi–Vinuesa, White, Rechnitzer), Cilleruelo–Ruzsa–Vinuesa's g-Sidon density σ, Kravitz's dual τ, the Barnard–Steinerberger / Fish–King–Miller / Madrid–Ramos autocorrelation problems, the distinct-distance constant, the Mian–Chowla reciprocal sum, and Martin–O'Bryant's continuous-Ramsey D(x) — **no constant in the literature lies in [0.379, 0.381] except μ itself**.

Terence Tao's `optimizationproblems` repository ([github.com/teorth/optimizationproblems](https://github.com/teorth/optimizationproblems)) is the authoritative curated index of such constants. μ is listed as **C_{1b}** (its own page); C_{1a} (Sidon-autocorrelation) is the *adjacent* constant at ≈ 1.28–1.50 — Tao curates them as siblings, not duplicates.

## The neighborhood (constants surveyed, all ≠ μ)

| Constant | Numerical bracket | Source | Match? |
|---|---|---|---|
| **C_{1b} = μ (Erdős min overlap)** | [0.379005, 0.380876] | Tao C_{1b} | identity |
| C_{1a} (Sidon-autocorrelation) | [1.2802, 1.502862] | Tao C_{1a}; arXiv:1403.7988 | none |
| ν_2² = inf ‖f∗f‖_2² | [0.574636, 0.574643] (128 digits known) | arXiv:2602.07292, arXiv:2210.16437 | none |
| ν_∞ = inf ‖f∗f‖_∞ | ≈ 0.64–0.69 | arXiv:0807.5121 | none |
| σ (g-Sidon density, Cilleruelo) | [1.147, 1.252] | arXiv:2004.06611 | none |
| τ (Kravitz dual difference) | (1.560, 1.643] | arXiv:2004.06611 | none |
| Barnard–Steinerberger autocorrelation | ≤ 0.411 ‖f‖₁²; cannot be ≤ 0.37 | arXiv:1903.08731, arXiv:2001.02326, arXiv:2003.06962 | **near-miss** (see below) |
| Distinct-distance constant | [2.16150003, 2.2473] | arXiv:2505.20851 | none |
| Mian–Chowla reciprocal sum | [2.158452685, 2.158532684] | Salvia 2014 | none |
| Martin–O'Bryant D(x) | varied, none near 0.38 | arXiv:math/0410004 | none |
| B_h[g] size constants σ_h(g) | all > 1 | arXiv:1601.00928, arXiv:1604.00661 | none |

## Strongest "near-miss": Barnard–Steinerberger

[arXiv:1903.08731](https://arxiv.org/abs/1903.08731): for f ∈ L¹(ℝ), `min_{0≤t≤1} ∫ f(x) f(x+t) dx ≤ 0.42 ‖f‖_{L¹}²`. Sharpened to 0.411 (Fish–King–Miller), then 0.4071 (Madrid–Ramos). The lower-bound side establishes the constant cannot be ≤ 0.37, so the true constant is in (0.37, 0.411).

**Structurally similar to μ** (min-over-translations of a correlation), and our μ ≈ 0.38 sits inside their bracket. But:
- Their f is in L¹(ℝ) — no compact support, no fixed [0, 2] domain, no `0 ≤ h ≤ 1` pointwise constraint
- Their bound is `min ≤ const · ‖f‖₁²`; μ is `min-max of overlap of (h, 1−h)` — different functional shape

Not the same constant, but the closest analog in the literature.

## **The Rechnitzer technique-transfer opportunity (critical)**

[arXiv:2602.07292](https://arxiv.org/abs/2602.07292) — Rechnitzer (2026), "128 digits of the autoconvolution L² constant".

What he did:
1. **Ansatz for the optimizer:** `f̂(k) = (-1)^k · a/√k · (1 + O(1/k))` at large k (PDE-aware form derived from the structure of the extremal function)
2. **Parameterize** the deviation from the ansatz by a small finite-dimensional vector
3. **Optimize** with BFGS in mpmath at 128-digit precision
4. **Levin acceleration** of the slow-convergent double sums for the objective evaluation

Result: 128 verified digits of ν₂² — *vastly* tighter than any SDP could produce.

**Why this might transfer to μ:**
- Our extremal h\* (Together's) has clean asymptotic structure near the boundary (h ≈ 0 region, h ≈ 1 region, and the smooth interior — see PRO-23 active-set analysis)
- The KKT functional equation `Σγ·[h(x+t)+h(x-t)] = κ` (LEVER_FUNCTIONAL_EQUATION) is exactly the kind of PDE-like structure Rechnitzer exploits
- If we can write h\*'s discrete optimum as `ansatz + small_deviation`, we can run mpmath BFGS at ~50+ digits on a *finite-dimensional* parameter space

**Possible payoff:** 30+ digit μ\_UB (or both bounds if the ansatz is two-sided) without needing SDPA-GMP at all. This is the contrarian path mentioned in the brainstorm — bypass the SDP infrastructure entirely.

## What the literature does NOT have

- A closed form for μ (confirmed by both PSLQ-negative and literature-negative)
- A non-SDP method tailored to μ specifically — Rechnitzer's technique is for ν₂², not μ, but the framework looks transferable
- Any tighter known bound than White (LB 0.379005) + Haugland/Together (UB 0.380871)

## Recommendations

1. **High-priority new task: investigate Rechnitzer technique transfer to μ.** Read his preprint carefully, identify the structural assumptions, see whether μ admits a similar ansatz. Spin as **PRO-26**.

2. **Don't drop PRO-11 (serializer).** Even with a working Rechnitzer-style optimizer, we still want SDPA-GMP for *verification* of the analytical result.

3. **Update preprint v2 (PRO-5)** to cite Tao's curated repo and acknowledge the Sidon/autoconvolution literature explicitly — currently the preprint doesn't situate μ in this neighborhood.

4. **Strategic ordering shift:**
   - 🥇 PRO-26: Rechnitzer technique transfer (new highest priority)
   - 🥈 PRO-11: cvxpy → SDPA-S serializer (still high; now for *verification* rather than *production*)
   - 🥉 PRO-5: preprint v2 update with new context
   - 4. PRO-8: White email v3

## Sources

All hyperlinked. See also: [Tao's `optimizationproblems` C_{1b} page](https://github.com/teorth/optimizationproblems/blob/main/constants/1b.md) — the canonical name for our constant μ.
