# Kickoff prompt — Erdős' minimum overlap problem

Paste the section between the `---` markers into a fresh session. Everything before that is context for Ben.

## Why this problem

Distinct flavor from the subset-sums work we just did: this is a **partition / convolution** problem with a **continuous reformulation** that turns it into step-function optimization. That makes it the rare open Erdős problem where AI-style numerical search has actually moved the SOTA in 2025–2026 (DeepMind AlphaEvolve, Stanford TTT-Discover, Together Computer). The lower-bound side is analytic (Fourier + convex programming, White 2023) and not yet AI-touched. So there are two genuine attack surfaces, and the gap is currently **about 0.0019**.

## What "AI on Erdős problems" actually does well — short version

From Tao's `teorth/erdosproblems` wiki and recent literature:

- **Strong:** numerical-construction problems. AlphaEvolve's main wins are extremal numerics and bound improvements (≈15 problems). FunSearch improved the cap-set lower bound in dim 8. Wagner-style cross-entropy RL refuted several spectral-graph conjectures by finding small explicit counterexamples.
- **Strong:** Lean formalization of known proofs (Aristotle dominates here; ~100+ problems formalized).
- **Weak:** producing rigorous new proofs of deep statements. Most "AI solved" entries on the wiki turn out to duplicate prior literature, attack a variant, or contain gaps. Tao explicitly notes that **most reasonably accessible Erdős problems have already been silently attempted by frontier models without success** — selection bias makes the success column look bigger than it is.
- **Failure mode to watch:** AI confusing one problem for another, presenting conditional results as unconditional, missing the paper that already solved it.

The minimum overlap problem fits the "strong" category: numerical, has well-defined certificates (a step function and an LP dual), every claimed bound can be verified rigorously by checking finitely many integrals.

---

# (Begin paste — kickoff prompt)

You are an expert mathematician working in the project "Erdős." Project rules apply: rigor, no bluffing, distinguish proved from heuristic, give the strongest partial result if a full one is unavailable.

## Problem (Erdős, 1955) — the minimum overlap problem

For each $n$, partition $[1,2n]=A\sqcup B$ with $|A|=|B|=n$. For each integer $k$, let $M_k(A,B)=|A\cap (B+k)|=|\{(a,b)\in A\times B:a-b=k\}|$. Set $M(A,B)=\max_k M_k(A,B)$, and
$$M(n)\;=\;\min_{A,B}M(A,B).$$
Erdős asked for the asymptotics of $M(n)/n$. Haugland (2016) proved the limit exists; call it $C_5\in(0,1)$.

**Equivalent continuous formulation (Haugland 2016).** Let $h:[0,2]\to[0,1]$ be a step function with $\int_0^2 h=1$. Then
$$C_5=\inf_h\;\sup_{k\in\mathbb R}\int_0^2 h(x)\bigl(1-h(x+k)\bigr)\,dx,$$
where $h$ is extended by $0$ outside $[0,2]$. Upper bounds on $C_5$ come from explicit step functions; lower bounds come from duality / Fourier-analytic obstructions to better step functions.

**Current state of the art (verify before quoting).**
- Lower bound: $C_5 \geq 0.379005$ — White 2023 (*Acta Arithmetica*), elementary Fourier reduction to a convex program.
- Upper bound: $C_5 \leq 0.380871$ — Together Computer, March 2026, sequential LP refinement of a 600-step function, building on TTT-Discover (Yuksekgonul et al., Jan 2026, $0.380876$) and AlphaEvolve (Georgiev–Gómez-Serrano–Tao–Wagner, May 2025, $0.380924$).

So the open gap is $C_5\in[0.379005,\,0.380871]$, about $1.9\times 10^{-3}$.

## What counts as progress (in priority order)

1. **A new lower bound** $C_5 \geq 0.379005 + \varepsilon$. The lower-bound side has *not* been AI-attacked publicly. White's bound comes from a finite-dimensional convex program over Fourier coefficients; tightening it likely means (a) using more frequencies, (b) better convexification, or (c) a structural lemma forcing the optimal $h$ to lie in a smaller class. Any rigorous improvement here is publishable.
2. **A new upper bound** $C_5 \leq 0.380871 - \varepsilon$ via an explicit step function. To be defensible: deliver $h$ as a JSON / array of (interval, value) pairs, verify $\int h=1$ and $\sup_k\int h(x)(1-h(x+k))\,dx \leq B$ to high precision, and show the verification reproduces with independent code.
3. **A meta-result** about the optimizer: e.g., "any optimal $h$ has at most $N$ steps" or "any optimal $h$ is symmetric about $x=1$." Structural lemmas of this kind would constrain the search and are valuable independent of bound improvements.
4. **Negative result / lower bound on the bound:** prove $C_5 \geq c$ for the strongest provable $c$ via a self-contained, short argument (this is partial credit; reproduces White-style bounds at lower constants).

A clean exposition of the current SOTA, with verified numerics, is also a valid output if no improvement materializes.

## Suggested attack plan

**Phase 1 — orient (low cost).**
- Reproduce the discrete $M(n)$ for small $n$ by brute force ($n\leq 20$ tractable). Compare to the continuous $C_5$. Tabulate $M(n)/n$ and the empirical rate of convergence; check it's consistent with $C_5\approx 0.3808$.
- Re-derive the continuous reformulation from scratch. Confirm signs and conventions match Haugland's; this is where 80% of follow-up errors enter.

**Phase 2 — pick a side.**
- *Lower-bound route.* Read White (Acta Arith. 2023) carefully. Set up the convex program in CVXPY/Mosek with $N=200,500,1000$ frequencies. Try to identify slack; consider auxiliary inequalities (e.g., monotonicity of the Fourier transform of $h$, sumset constraints). Report the exact LP/SDP and any provable improvement.
- *Upper-bound route.* Reimplement sequential LP refinement on step functions with 600–2000 steps, starting from the Together Computer 600-step solution if available, otherwise from the Haugland 51-step solution. Use double precision, then verify any candidate in interval arithmetic (e.g., `mpmath`) before claiming.

**Phase 3 — adversarial verification.** Whatever bound you produce, recompute $\int h(x)(1-h(x+k))\,dx$ over a *dense grid of $k$* and an analytic worst-case argument; never rely on the optimizer's reported value.

## Methods worth considering

- **Sequential / column-generation LP** for upper-bound step functions (proven productive: AlphaEvolve, TTT-Discover, Together used variants).
- **Fourier-analytic convex programming** for lower bounds (White's method; under-attacked computationally).
- **Symmetry exploitation:** if optimal $h$ is symmetric about $x=1$, the search dimension halves.
- **Interval arithmetic** for rigorous verification (`mpmath`, `flint`).
- **LLM-guided refinement** of step-function neighborhoods is the kind of thing AlphaEvolve does; a poor man's version is a local search + LLM critique loop.

## Anti-bluff guardrails

- Do not claim a numerical bound without an independent verification path.
- Do not claim a proof unless every step is rigorous and finite; if the proof is conditional, say so.
- Do not invent constants. If you cite a paper, quote the bound to the precision the paper states it.
- If you suspect the problem has been silently solved, search arXiv (`erdos minimum overlap`, `C_5 minimum overlap`, recent 2025–2026 preprints) before writing up.
- The Tao wiki note applies: assume frontier models have already poked at this problem; expect routine LLM ideas to fail; aim for an angle the SOTA constructions did not exploit (especially on the lower-bound side).

## Deliverables

A single markdown report with:
1. Restated problem + current SOTA (verified citations).
2. Whatever was attempted, with code and numerics in a reproducible block.
3. The strongest *rigorously defended* result obtained, clearly labeled "proved" vs "numerical" vs "heuristic."
4. A short "what I'd try next" section.

Begin by acknowledging that you have read this brief, then start with Phase 1.

# (End paste)

---

## References

- [Minimum overlap problem — Wikipedia](https://en.wikipedia.org/wiki/Minimum_overlap_problem)
- [Together Computer SOTA repo (Mar 2026)](https://github.com/togethercomputer/erdos-minimum-overlap)
- [White (2023), *Acta Arithmetica*, "A new bound for Erdős' minimum overlap problem"](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/online/115217/a-new-bound-for-erdos-minimum-overlap-problem)
- [Haugland (2016) "The minimum overlap problem revisited," arXiv:1609.08000](https://arxiv.org/abs/1609.08000)
- [Georgiev, Gómez-Serrano, Tao, Wagner (2025), AlphaEvolve, arXiv:2511.02864](https://arxiv.org/abs/2511.02864)
- [Yuksekgonul et al. (2026), TTT-Discover, arXiv:2601.16175](https://arxiv.org/abs/2601.16175)
- [Tao's `teorth/erdosproblems` wiki: "AI contributions to Erdős problems"](https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erd%C5%91s-problems)
- [Wagner (2021) "Constructions in combinatorics via neural networks," arXiv:2104.14516](https://arxiv.org/abs/2104.14516)
- [Romera-Paredes et al. (FunSearch, 2023), Nature](https://www.nature.com/articles/s41586-023-06924-6)
