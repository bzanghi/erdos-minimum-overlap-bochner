# Lever D pre-investigation: does the optimal partition look "low-complexity"?

**Status:** complete. **Verdict:** HYPOTHESIS REFUTED.

This memo records an empirical test of the working assumption behind
Lever D, namely *"the optimal `f` for the continuous Erdős minimum
overlap problem lies near piecewise-constant densities with `O(1)`
breakpoints / low block-count structure."* Two independent data
sources are inspected: Together's claimed near-optimal certificate
`h*` (600 cells on `[0, 2]`, achieving `µ ≤ 0.380871`), and the
brute-force-optimal integer partitions `A*` for `n = 2..12`. Both
point the same direction. The hypothesis as stated is not supported.

Probe code: [`lp_research_state/code/_lever_d_structure_probe.py`](lp_research_state/code/_lever_d_structure_probe.py).
Data summary: [`lp_research_state/data/lever_d_structure_summary.json`](lp_research_state/data/lever_d_structure_summary.json).

## 1. Hypothesis under test

> The optimal density `f` for the continuous Erdős minimum overlap
> problem has `O(1)` breakpoints, i.e. it is well-approximated by a
> piecewise-constant function with a constant (or at most logarithmic)
> number of jumps as the discretization is refined.

This is the structural premise behind Lever D — that a structure
theorem restricting attention to densities with `O(1)` blocks should
be possible, and would cut the SDP feasible set drastically.

## 2. Evidence from Together's `h*`

Together's `h*` lives on `[0, 2]` with 600 equal-width cells, `Σ h_i =
300`, `h ∈ [0, 1]`. It is the best near-optimal step function in the
public literature (`µ ≤ 0.380871`). Plot:
[`lp_research_state/data/together_h_structure.png`](lp_research_state/data/together_h_structure.png).

| metric | value |
|---|---|
| cells (n) | 600 |
| nonzero cells (`h_i > 10⁻⁶`) | **436** (73% of cells) |
| "saturated to 1" cells (`h_i > 1 − 10⁻³`) | 54 |
| "in the middle" cells (`0.05 < h_i < 0.95`) | **347** (58%) |
| distinct values, rounded to 10⁻³ | 198 |
| distinct values, rounded to 10⁻² | 80 |
| distinct values, rounded to 5·10⁻² | 21 |
| distinct values, rounded to 10⁻¹ | 11 |
| total variation `Σ |h_i − h_{i-1}|` | **26.91** |

**Blocks** = maximal runs whose values agree to within tolerance:

| tolerance | block count |
|---|---|
| 10⁻³ | **399** |
| 10⁻² | 364 |
| 5·10⁻² | 206 |
| 10⁻¹ | 89 |
| 2·10⁻¹ | 13 |
| 5·10⁻¹ | 1 |

**Jumps** between adjacent cells (`|h_i − h_{i-1}| > θ`):

| θ | count |
|---|---|
| 0.01 | 363 |
| 0.05 | 205 |
| 0.10 | **88** |
| 0.20 | 12 |
| 0.50 | 0 |

Interpretation. `h*` is not bimodal toward `{0, 1}`: only 54 cells
saturate to ~1 and 164 cells sit below 10⁻³ (most of those are the
"tails" near `x = 0` and `x = 2`, where `h_i ≈ 10⁻¹¹` — almost
certainly numerical zeros from Together's solver). The bulk — 347 of
600 cells — sit strictly in the middle. Total variation is **27**,
roughly an order of magnitude larger than what a 5- or 10-block
indicator-like function could possibly produce on `[0, 1]`.

You can collapse `h*` down to roughly 13 blocks **only if you allow
tolerance 0.2** — i.e. only if you call values that differ by up to
±0.2 "the same block". By any tighter standard the function has at
least ~90 well-separated levels.

There is one charitable reading: `h*` is a near-optimum found by
stochastic optimization, not the true optimum, and its high-frequency
wiggle may be solver noise about a smoother underlying shape. The
top-panel plot does in fact look like a few smooth lobes overlaid
with a "bumpy" mid-section; if those bumps are noise, then the
effective complexity is closer to "a few smooth pieces" (which would
already not be `O(1)`-piecewise-constant, but would be `O(1)`-piecewise-smooth).
We treat that as inconclusive evidence at best — the Lever D
hypothesis is specifically about *step* functions with `O(1)`
*breakpoints*, and `h*` simply does not have that property at any
modest tolerance.

## 3. Evidence from small-n optimizers

Brute-force-optimal partitions `A* ⊂ {1, ..., 2n}` for `n = 2..12`
(taken from `min_overlap_session_2026-05-09/Mn_brute.json`). Converted
to 0/1 vectors of length `2n` and counted maximal runs ("blocks"):

| n | `M(n)` | `M(n)/n` | blocks | palindrome? | anti-palindrome? | A* as 0/1 |
|---|---|---|---|---|---|---|
| 2 | 1 | 0.500 | 3 | yes | no | `1001` |
| 3 | 2 | 0.667 | 4 | no  | yes | `110100` |
| 4 | 2 | 0.500 | 5 | no  | no  | `11010001` |
| 5 | 3 | 0.600 | 4 | no  | no  | `1111001000` |
| 6 | 3 | 0.500 | 7 | no  | no  | `111010010001` |
| 7 | 3 | 0.429 | 5 | no  | no  | `11100100000111` |
| 8 | 4 | 0.500 | 7 | no  | no  | `1111010001000011` |
| 9 | 4 | 0.444 | 7 | no  | no  | `111100100010000111` |
| 10| 5 | 0.500 | 10| no  | no  | `11111010001000010110` |
| 11| 5 | 0.455 | 9 | no  | no  | `1111100100100001000111` |
| 12| 5 | 0.417 | 7 | no  | no  | `111110001000010000011111` |

Plot: [`lp_research_state/data/Mn_optimizers_structure.png`](lp_research_state/data/Mn_optimizers_structure.png).

Growth fits (`blocks` against `n`):

- linear: `blocks ≈ 0.555 · n + 2.30`
- sqrt:   `blocks ≈ 2.78 · √n − 0.96`
- log:    `blocks ≈ 3.22 · log n + 0.33`

The linear fit is the most consistent visually (blocks = 3, 4, 5, 4,
7, 5, 7, 7, 10, 9, 7). There is real scatter — `n = 7, 12` are
abnormally clean (just 5 and 7 blocks); `n = 10` is unusually messy
(10 blocks) — but the trend is unambiguously *upward* and at least
proportional to `n` over this range. The block count is *not* `O(1)`
and is not even cleanly `O(log n)`.

Of the 11 optimizers exactly one is a palindrome (`n=2`) and one is
anti-palindromic (`n=3`). The other 9 are neither, ruling out a
clean "halved by reflection" structure theorem for the integer case.

Almost every optimizer has a prefix of consecutive `1`s and a suffix
of consecutive `1`s (e.g. `n=12`: `11111...11111`), with a structured
*middle* zone of length scaling with `n`. The middle zone is where
the block-count growth lives. That is consistent with the well-known
heuristic "A and B are dense at their respective ends" but it does
*not* shrink the middle to constant size.

## 4. Verdict — REFUTED

The hypothesis "optimal `f` has `O(1)` breakpoints / low block count"
is not supported by either evidence stream:

- Together's `h*` has on the order of 100 jumps that exceed 10%,
  total variation `≈ 27`, and 198 distinct values at 10⁻³ resolution.
  Even granting that some of the wiggle is solver noise, it is far
  from a step function with `O(1)` steps.
- Small-`n` optimizers have block counts growing at least linearly in
  `n` over `n = 2..12`. The slope (≈ 0.55) is *not* tending to 0.

A Lever D formulated as **"restrict the SDP feasible set to step
functions with `O(1)` breakpoints"** would therefore exclude the
optimum (or its best known approximation), and any bound it produces
would not be a bound on `µ`.

What the data *might* support instead, hedged by the small-`n` /
single-`h*` caveat:

- **`f` is dense almost everywhere** — `73%` of `h*` cells are
  nonzero. The optimum is *not* an indicator function. This already
  rules out the "Lever D" variant where one assumes `f` takes only
  values in `{0, 1}`.
- **`f` is bounded away from `{0, 1}` on a constant fraction of the
  support** — `58%` of cells sit in the strict interior `(0.05,
  0.95)`. So the right structural object, if any, is a *smooth* (or
  Lipschitz-bounded) function, not a step function.
- **The "structured middle, dense ends" pattern in the integer case**
  is the only observable that does have stable structure across `n`;
  but the middle zone scales with `n`, so this isn't an `O(1)`
  parameter family either.

## 5. Recommended next step

Drop Lever D as formulated. Two replacement levers are visible from
the same data:

1. **Lever C — push exact `M(n)` to larger `n`.** The block-count
   growth observed here is itself an open question: is it actually
   `Θ(n)`, or does it slow at `n = 15, 20, 25`? Extending the brute
   force / branch-and-bound from `n = 12` to `n ≈ 20` would tell us
   whether the structural picture stabilizes or keeps growing. Cost
   is bounded — `n = 20` is `C(39, 19) ≈ 6·10¹⁰` naively, but with
   symmetry-breaking and proper pruning is reachable in
   process-weeks. The output would also (independently of Lever D)
   refine `µ ≥ M(n)/n` lower bounds for n where the SDP isn't tight.

2. **Recast as a `BV` or Lipschitz constraint, not a step-function
   constraint.** If `f` is dense and varies smoothly (as Together's
   `h*` suggests modulo noise), then the right structural restriction
   is a bound on `||f'||₁` or `||f||_Lip`. That has the right
   character — it cuts feasible set, doesn't a priori exclude the
   optimum, and pairs naturally with Bochner-style moment constraints
   on `f̂`. This is a research question, not an engineering task; it
   needs a dedicated investigation.

Of the two, **Lever C** is the cleaner, lower-risk next step. The
infrastructure exists (`brute_force_Mn.py` is 40 lines; the
acceleration to `n = 16-20` is routine combinatorial optimization
work), and the result is informative regardless of which structural
direction Lever D-replacement eventually takes.
