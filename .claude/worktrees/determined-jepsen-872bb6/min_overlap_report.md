# Erdős minimum-overlap problem — Phase 1+2 working notes

**Date:** 2026-05-09  &nbsp;·&nbsp; **Author:** Claude (Erdős project)
**Status:** Reproduces methodology of state-of-the-art on both sides; **does not improve** SOTA.
Strongest rigorously verified results obtained in this session:
- **Upper bound:** `C_5 ≤ 0.3810916049…` (explicit step function on N=400 cells, three independent verification paths agree to 16 digits).
- **Lower bound:** `C_5 ≥ 1/4` (one-line pigeonhole; far from White's 0.379005).

Code, data, and reproducer scripts are in `outputs/erdos_minoverlap/`.

---

## 1 · Problem and current SOTA

### Discrete formulation
For each `n ≥ 2`, consider partitions `[1, 2n] = A ⊔ B` with `|A| = |B| = n`. For `k ∈ Z`,
`M_k(A,B) := |{(a,b) ∈ A × B : a − b = k}|`, `M(A,B) := max_k M_k`,
and `M(n) := min_{partitions} M(A,B)`.
Erdős (1955) asked for the asymptotic behavior of `M(n)/n`. Haugland (2016) proved the limit
exists; call it `C_5 ∈ (0,1)`.

### Continuous reformulation (corrected)
Let `h: [0,2] → [0,1]` be a step function with `∫h = 1`. Define
$$J(t) \;=\; \int_0^{\,2-t} h(x)\bigl(1-h(x+t)\bigr)\,dx, \quad t\in[0,2],$$
and similarly `J(-t) := ∫_0^{2-t} h(x+t)(1−h(x))dx`. Then
$$C_5 \;=\; \inf_h\;\sup_{t\in[-2,\,2]}\;J(t).$$
**Note on convention.** The brief states the formula with `h` "extended by 0 outside `[0,2]`" and
the integration on the full `[0,2]`. As a numerical check (script `derivation_check.py`), this
"extended" version gives `sup_t J = 1` for the rectangle `h = 1_{[1/2, 3/2]}`, while the discrete
problem and the corrected restricted-domain formulation both give `1/2`. The integration must be
restricted to `x` with both `x` and `x+t` in `[0,2]` — equivalently, `x ∈ [0, 2-t]` for `t ≥ 0`.
All code in this report uses the corrected version.

### State of the art (verified)
- **Lower bound:** White, Acta Arith. (2023), `C_5 ≥ 0.379005`, via finite-dim convex program
  on Fourier coefficients of `h`.
- **Upper bound:** Together Computer (March 2026), `C_5 ≤ 0.380871`, via sequential-LP
  refinement of a 600-step `h`, building on TTT-Discover (Yuksekgonul et al., Jan 2026,
  `0.380876`) and AlphaEvolve (Georgiev–Gómez-Serrano–Tao–Wagner, May 2025, `0.380924`).

Open gap: `C_5 ∈ [0.379005, 0.380871]`, width `≈ 1.9 × 10⁻³`.

---

## 2 · Phase 1 — orientation

### 2.1 Brute-force `M(n)` for small `n`
Code: `brute_force_Mn.py`. Approach: enumerate `A ⊂ [1,2n]` with `|A| = n` (fixing `1 ∈ A`),
compute `max_k M_k` via `np.correlate`, take the min over partitions.

| `n` | `M(n)` | `M(n)/n` | Optimal `A*` (one example) |
|---:|---:|---:|---|
| 2 | 1 | 0.500 | (1, 4) |
| 3 | 2 | 0.667 | (1, 2, 4) |
| 4 | 2 | 0.500 | (1, 2, 4, 8) |
| 5 | 3 | 0.600 | (1, 2, 3, 4, 7) |
| 6 | 3 | 0.500 | (1, 2, 3, 5, 8, 12) |
| 7 | 3 | 0.429 | (1, 2, 3, 6, 12, 13, 14) |
| 8 | 4 | 0.500 | (1, 2, 3, 4, 6, 10, 15, 16) |
| 9 | 4 | 0.444 | (1, 2, 3, 4, 7, 11, 16, 17, 18) |
| 10 | 5 | 0.500 | (1, 2, 3, 4, 5, 7, 11, 16, 18, 19) |
| 11 | 5 | 0.455 | (1, 2, 3, 4, 5, 8, 11, 16, 20, 21, 22) |
| 12 | 5 | 0.417 | (1, 2, 3, 4, 5, 9, 14, 20, 21, 22, 23, 24) |
| 13 | 6 | 0.462 | (1, 2, 3, 4, 5, 6, 9, 13, 18, 19, 24, 25, 26) |

This sequence matches the published Erdős minimum-overlap values (cross-check: `1, 2, 2, 3, 3,
3, 4, 4, 5, 5, 5, 6` for `n = 2, …, 13` — see Haugland's tabulations).

**Comment on convergence.** `M(n)/n` is **not** monotonically decreasing — it oscillates because
of the integer ceiling: `M(n)` can be no smaller than `⌈C_5 · n⌉`, and `⌈0.3808 n⌉ / n` itself
oscillates. For example, `M(12)/12 = 5/12 ≈ 0.417 ≈ ⌈0.3808 · 12⌉ / 12`, but
`M(13)/13 = 6/13 ≈ 0.462 > ⌈0.3808 · 13⌉ / 13`, so the asymptote is not yet "binding" at
`n = 13`.

### 2.2 Continuous reformulation — re-derivation and convention check
Done above; numerical confirmation in `derivation_check.py`. Two important practical consequences:
1. The autocorrelation `R̃(t) := ∫_0^{2-t} h(x) h(x+t) dx` is **piecewise-linear in t** when `h`
   is piecewise-constant on a uniform grid of width `Δ = 2/N`. Hence the supremum
   `sup_{t ∈ [0,2]} J(t)` is attained at some `t = jΔ` for `j ∈ {0, 1, …, N-1}`. Verifying the LP
   constraints at integer shifts is sufficient — no "between-grid" residual.
2. The discrete problem is symmetric under `(A, B) ↔ (B, A)`, so the continuous bound is
   `max(J^+_j, J^-_j)`. For `h` symmetric about `x=1` (`h_i = h_{N+1-i}`), `J^+_j = J^-_j` and
   only one direction is needed. **All numerical work below uses the symmetric ansatz.**

---

## 3 · Phase 2A — upper bound (verified construction)

### 3.1 Methodology
Discretize: `h ∈ [0,1]^N` with `Σ h_i = N/2` (so `∫h = 1`). Constraints:
$$J^+_j(h) = \tfrac{2}{N}\sum_{i=1}^{N-j} h_{i+j}(1-h_i), \qquad
  J^-_j(h) = \tfrac{2}{N}\sum_{i=1}^{N-j} h_i(1-h_{i+j}),$$
for `j = 1, …, N-1`. Objective: minimize `max_{j, σ} J^σ_j`.

This is bilinear in `h` with linear constraints. I attacked it with:
1. **Smoothed-max + L-BFGS-B**: replace `max` by `(1/β) log Σ exp(β J^σ_j)`, anneal `β` from
   10 to 100 000, projection onto `{Σh = N/2, 0 ≤ h ≤ 1}` between rounds.
2. **Multi-start** from constant `½`, triangle, rectangle, cosine bump, and 6 Gaussian random
   perturbations (`run_search.py`).
3. **Basin-hopping** with random kicks of std 0.05 followed by re-anneal (`basin_hop.py`).
4. **Sequential LP polishing** (`slp_step` in `basin_hop.py`): linearize each `J^σ_j` around
   the current `h`, solve the LP `min c subject to linearized J ≤ c, sum = N/2, |h - h_old|_∞ ≤ τ`
   with `scipy.optimize.linprog(method='highs')`. Trust radius `τ = 0.01`–`0.03`.
5. **Upsampling**: warm-start higher-`N` runs by linearly interpolating a refined lower-`N`
   solution.

### 3.2 Results
| Stage | `N` | best `bound(h)` | argmax j (sign) | runtime |
|---|---:|---:|---:|---:|
| L-BFGS multistart | 80 | 0.381533 | 18 (+) | <1s |
| L-BFGS multistart | 200 | 0.381623 | 59 (+) | 2.3s |
| Basin-hop + SLP | 200 | **0.381333** | 48 (+) | ~5s |
| Pure SLP polish from above | 200 | **0.381104** | 42 (+) | ~7s |
| Upsample to N=400 + L-BFGS + SLP | 400 | **0.3810916049** | 47 (+) | ~10s |

**Final verified upper bound** (this session): `C_5 ≤ 0.3810916049`. The optimizer's `h` is saved
in `refined_N400.json`; the optimal shift is `t* = 47 · (2/400) = 0.235`. The objective profile
`J(t)` over `t ∈ [0, 2]` is plotted in `h_and_J_curves.png`.

This is **above SOTA's `0.380871`** by `≈ 2.2 × 10⁻⁴`. With more compute I'd expect to close
this gap (see §6 for what would help), but I did not in this session.

### 3.3 Adversarial verification (Phase 3)
Three independent code paths in `verify.py`:

- **V1** Plain Python loop over `(i, j)`, no FFT, no vectorization: `0.3810916049…` ✓
- **V2** Continuous dense-`t` scan with `n_t = 20001` evaluations, walking each `J(t)` cell-by-
  cell with exact piecewise integration. Sup achieved exactly at `t = 0.23500` (= `47·Δ`) ✓.
  Confirms the piecewise-linearity argument: no between-grid residual.
- **V3** `mpmath` recompute at 50-digit precision: `0.3810916049136458917063676…`. Differs from
  V1 by `6.7 × 10⁻¹⁶` (float-rounding only) ✓.

The bound `0.3810916049` for the saved step function `h` is therefore rigorous to all printed
digits.

---

## 4 · Phase 2B — lower bound

### 4.1 Trivial pigeonhole bound (proved)
**Claim.** For all `n ≥ 2`, `M(n) ≥ ⌈ n² / (4n − 2) ⌉`. Hence `C_5 ≥ 1/4`.

*Proof.* `M_k(A,B) > 0` only for `k ∈ [-(2n-1), 2n-1] \ {0}`, so at most `4n − 2` values of `k`.
And `Σ_k M_k(A,B) = |A||B| = n²`. By pigeonhole, `max_k M_k ≥ n² / (4n − 2)`. Taking `n → ∞`
gives `C_5 ≥ 1/4`. ∎

Cross-check against brute force: every entry of the `M(n)/n` table satisfies the bound (column
"≥ trivial bound" in `lower_bound.py`).

### 4.2 Why simple LP relaxations don't help
Tried (in `lower_bound.py`): replace `J^+_j(h) = (2/N)(S^+_j − Σ_i h_i h_{i+j})` by its
McCormick relaxation `z_{i,j} ≤ min(h_i, h_{i+j})` and minimize `Σ_j λ_j J^+_j` with `λ_j` uniform.
For symmetric `h ≈ const`, `S^+_j ≈ Σ min(h_i, h_{i+j})`, so the relaxation's value collapses to
`0`. No improvement on `1/4`.

This is well known: the McCormick polytope is too loose for this problem. White's (2023)
`0.379005` uses a Fourier-analytic dual: replace the constraint `h ∈ [0,1]` by truncated Fourier
expansions and use specific positivity / Schur-test arguments on the resulting kernel. **I did
not reproduce White's bound in this session**; doing so requires reading the Acta Arith. paper
carefully and setting up the SDP/LP over Fourier coefficients. A faithful description of the
LP setup is in §6.

### 4.3 Range of validity
Combining: in this session I rigorously have `1/4 ≤ C_5 ≤ 0.3810916049`. Combining with cited
literature (which I have not independently verified) we have the open range
`C_5 ∈ [0.379005, 0.380871]` of width `≈ 1.9 × 10⁻³`.

---

## 5 · Strongest results, clearly labeled

| Claim | Status | Evidence |
|---|---|---|
| `C_5 ≥ 1/4` | **Proved** (rigorous) | One-line pigeonhole; `lower_bound.py`. |
| `C_5 ≤ 0.3810916049…` (explicit `h` on N=400 cells) | **Proved numerically** to 16 digits | `refined_N400.json`; three independent verification paths agree (`verify.py`). |
| Continuous reformulation as stated in brief is mis-formed at the boundary | **Proved by counterexample** | Rectangle `h` gives sup ≥ 1 with extension-by-0; `derivation_check.py`. The corrected restricted-domain version matches the discrete problem. |
| `M(n)` for `n = 2, …, 13` matches published values | Computational | `Mn_brute.json`. |
| White (2023): `C_5 ≥ 0.379005` | **Cited, not reproduced** | Methodology described in §6, not implemented. |
| Together Computer (2026): `C_5 ≤ 0.380871` | **Cited, not reproduced** | I would need their step-function data. |

I am **not** claiming any improvement on SOTA in either direction.

---

## 6 · What I'd try next

### 6.1 Closing the upper-bound gap (mine vs. SOTA)
My `0.3810916` reflects (a) too-low `N` (400 vs SOTA 600), (b) the smoothed-max + L-BFGS
optimizer getting stuck in basins that SLP-only methods avoid. I'd:

1. **Larger N + proper SLP loop, no smoothed max.** Implement only the trust-region SLP
   (`linprog` with `highs` is fast for `N ≤ 1000`; the LP has `2(N−1) ≈ 1200` rows). Run
   ~1000 iterations with adaptive trust radius (`τ` halves when no improvement). This is what
   AlphaEvolve / TTT-Discover / Together Computer use under the hood.
2. **Restart from published step functions.** If Together Computer's 600-step `h` is available
   (Github / preprint), warm-starting from it should let me match `0.380871` immediately and then
   try to refine.
3. **LLM-guided perturbations (the ollama suggestion).** Wire `ollama` (e.g. `qwen2.5-coder`,
   `llama3.1`) to propose structural moves: "merge two adjacent plateaus", "introduce a
   discontinuity at `x = 0.62`", etc. Score with the verified evaluator, accept on improvement.
   Trade-off: this works well only when the LP is already near-optimal and creative jumps are
   needed; before then, LP + SLP dominates.

### 6.2 Lower-bound side (more interesting)
The lower bound has been less attacked computationally (per the brief). I would:

1. **Reproduce White (2023) literally.** Set up the variables: Fourier coefficients
   `\hat{h}(ξ)` for `ξ ∈ {0, 1/2, 1, …, K/2}` (Nyquist-truncated), constraints
   `\hat{h}(0) = 1` and `h(x) ∈ [0,1]` enforced via positivity of certain trigonometric
   polynomials (LMI / SDP). Objective: maximize `λ` subject to existence of a non-negative
   measure `μ` on shifts such that `λ ≤ ∫ J(t) dμ(t)` for the dual. Solve with MOSEK / SCS.
2. **Tighten via more frequencies.** White uses `K ≈ 50–100`. Pushing to `K = 500–2000` with a
   structure-aware solver (Toeplitz / Hankel exploitation) could yield a small but rigorous
   improvement.
3. **Symmetry reduction.** Restrict to `h` symmetric about `x = 1`. Many SOTA upper-bound `h`'s
   are symmetric (worth checking); if true, a meta-result *"every optimal `h` is symmetric"*
   would halve the Fourier search dimension. **Open structural question**, see §6.3.
4. **Auxiliary inequalities.** White's program has slack; adding Cauchy-Schwarz on `(h ⋆ h)` or
   moment constraints (`∫ x h dx`, `∫ x² h dx`) can tighten the dual. Worth a sweep.

### 6.3 Structural / meta-results worth pursuing
- **Symmetry of optimizers.** Conjecture: every `h` achieving the inf in `inf_h sup_t J(t)` is
  symmetric about `x = 1`. **Attempted argument:** if `h` is optimal then so is `h(2−·)`, and
  by convexity in `h` of the smoothed objective, the average `(h + h(2−·))/2` is at least as
  good. But this requires convexity of `sup_t J(t)` in `h`, which fails (it's a sup of bilinear
  functions). So the convexity argument doesn't go through directly. **Status:** I think this
  is open; resolving it would simplify all numerical work by a factor of 2.
- **Finiteness of optimal step count.** Conjecture: the inf over step functions equals the inf
  over functions in `BV([0,2])` (finite total variation). True at any `N ≥ N_0`?
- **Equality at `t* = j*Δ`.** Empirically, my N=400 optimum has 1 active shift constraint
  (`j = 47`). At true optimum, we'd expect *multiple* active constraints (like a Chebyshev
  alternation). This suggests my optimizer hasn't reached the true min at N=400.

### 6.4 An "all-night ollama job" — concrete sketch
Given the working verifier, the productive overnight job is:

```
loop:
  proposal = ollama_propose(current_h, "Suggest a perturbation; one of:
                  - shift a plateau by Δ, - merge two cells,
                  - introduce a step at x = c, - reflect about x = 1")
  candidate = apply(current_h, proposal)
  candidate = lp_polish(candidate)   # SLP, ~50 iterations
  bound_new = verified_evaluate(candidate)
  if bound_new < bound_current - 1e-7:
      current_h = candidate
      log(proposal, bound_new)
```

Persist `current_h` to disk every iteration. The cost per iteration is ~5–30 s (LP-dominated).
Over 8 hours that's ~1000 iterations, plenty to explore.

---

## 7 · Honest assessment

**What this session delivered.** A correct, reproducible toolkit (`evaluator.py`, `verify.py`,
`refiner.py`, `basin_hop.py`, `lower_bound.py`) for the minimum-overlap problem. A verified
upper bound `C_5 ≤ 0.3810916049…` from an explicit step function. A clean trivial lower bound
`C_5 ≥ 1/4`. A documented convention bug in the brief's continuous reformulation. The
machinery is ready for the overnight ollama-guided job in §6.4.

**What this session did *not* deliver.** Any improvement on SOTA. My UB is `2.2 × 10⁻⁴` above
Together Computer's `0.380871`. My LB is the trivial `1/4`, not White's `0.379005`. The
techniques used here (smoothed max + L-BFGS, basin hopping, SLP polishing) are the same family
as the SOTA work — without their compute budget or their carefully-tuned LP loop, I cannot
match them in a single short session.

The realistic path to a *publishable* result, in order of cost-effectiveness:
1. Implement White's Fourier LP and look for slack (lower bound, under-attacked, high upside).
2. Prove a structural lemma about optimal `h` (e.g. symmetry).
3. Beat the upper bound by `1 × 10⁻⁵` (low priority — well-trodden, requires beating
   AlphaEvolve / TTT-Discover / Together Computer at their own game).

---

## 8 · Files

| File | Purpose |
|---|---|
| `brute_force_Mn.py` | Discrete brute force for `n ≤ 13`. |
| `derivation_check.py` | Cross-check of brief vs. corrected continuous reformulation. |
| `evaluator.py` | Fast & exact evaluators for `J^σ_j(h)`; sanity tests. |
| `fast_eval.py` | Vectorized FFT/correlate-based evaluator + analytic gradient. |
| `refiner.py` | Smoothed-max + L-BFGS-B refinement. |
| `run_search.py` | Multistart driver. |
| `basin_hop.py` | Basin-hopping + SLP step. |
| `verify.py` | Three-way independent verification (loop / dense-t scan / mpmath). |
| `lower_bound.py` | Trivial bound + McCormick experiment. |
| `refined_N400.json` | Best `h` (N=400) — `bound = 0.3810916049…`. |
| `Mn_brute.json` | Discrete `M(n)` table for `n = 2…12`. |
| `h_and_J_curves.png` | Plot of optimized `h(x)` and the `J^±(t)` profile. |

To reproduce the upper bound from scratch:
```
python3 run_search.py 200 6
python3 basin_hop.py 200
python3 -c "from basin_hop import slp_step; ..."   # SLP polish (see §3.2)
python3 verify.py refined_N400.json
```
