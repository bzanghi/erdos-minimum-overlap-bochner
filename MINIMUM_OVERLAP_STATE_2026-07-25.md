# The Erdős minimum overlap problem: complete state, 2026-07-25

**Short answer to "solve it completely": the problem is open and this session did
not close it.** µ is not known exactly, no closed form is known, and no argument
in this repository — or in the 2026 literature — determines it. What follows is
the complete, verified state, what moved this session, and precisely what is left.

---

## 1 · The problem

Partition `{1, …, 2n} = A ⊔ B` with `|A| = |B| = n`. For `k ∈ ℤ` let
`M_k(A,B) = |{(a,b) ∈ A×B : a − b = k}|`, and set

```
M(n) = min over partitions of  max_k M_k(A, B),        µ = lim_{n→∞} M(n)/n.
```

The limit exists (Haugland). Equivalently (Moser–Murdeshwar; White §1), with
`h : [0,2] → [0,1]` measurable, `∫h = 1`, and `g := 1 − h` supported on `[0,2]`:

```
µ = inf_h  sup_{k ∈ ℝ}  ∫ h(x) g(x+k) dx.
```

Two facts organise everything:

- **Upper bounds are exhibitions.** Any admissible `h` gives `µ ≤ M(h)`. Step
  functions on `n` cells suffice, and `M_n := min` over them is non-increasing
  under cell doubling, with `µ = inf_n M_n`. So the UB side is a (nonconvex,
  bilinear) global optimisation — an AI-search benchmark target.
- **Lower bounds are exclusions.** They must rule out *every* `h` at once, via a
  certified relaxation. This is the hard side, and it is where this repository's
  contribution lies.

---

## 2 · State of the art — corrected

The bracket quoted throughout this repo (`CLAUDE.md`, `findings.md`,
`ub_core.ANCHOR`, the preprint drafts) is **stale at both ends**. Verified this
session by fetching primary sources, not summaries:

| | Value | Source | Date |
|---|---|---|---|
| Best known **UB** | **µ ≤ 0.380856** | SimpleTES ablation, [arXiv:2604.19341](https://arxiv.org/abs/2604.19341) §3.4.1 | 21 Apr 2026 |
| — leaderboard best | 0.3808591 | [Einstein Arena](https://einsteinarena.com/problems/erdos-min-overlap), lnzwz_AI4M_Agent | since Mar 2026 |
| — repo's anchor | 0.3808703105862199 | Together AI | Mar 2026 |
| Best **published** LB | **µ ≥ 0.37912** | Kim & Pilanci, [arXiv:2606.31182](https://arxiv.org/abs/2606.31182) | 30 Jun 2026 |
| — prior published LB | 0.379005 | White, *Acta Arith.* | 2023 |
| **This repo's LB** (unpublished) | **µ ≥ 0.3802838** | full-space cover, Bochner-PSD + ellipse extension | 2026 |

SimpleTES, verbatim: their headline is 0.380868 and *"an even better result was
found during our ablation study (0.380856), which is not included in the table
for fair comparison."* Four Einstein Arena entries beat Together's construction.

**Current true gap: `µ ∈ [0.3802838, 0.380856]`, width 5.72 × 10⁻⁴.**

Three consequences for this repo:

1. Its LB still leads the *published* LB by **1.16 × 10⁻³** — the advance is
   real and remains the strongest known lower bound. But Kim–Pilanci (June 2026)
   attack the same direction with the same rigor standard (interval arithmetic,
   outward rounding), so **priority is now time-sensitive**.
2. **PRO-34's central inference is falsified.** It argued that because Together's
   `h*` is a KKT point and a strict local minimum, "0.380871 is a serious
   candidate for µ itself". Since `µ ≤ 0.380856 < 0.3808703`, `h*` is *not*
   globally optimal. Its stationarity is a statement about its basin, nothing more.
3. The Lever I′ framework ceiling `C_explicit = 0.380713` is still **below** the
   new UB, so the saturation theorem is unaffected: the LB architecture remains
   provably unable to reach the true value.

---

## 3 · What moved this session

### 3.1 The upper bound is now certified in exact arithmetic — and the old one wasn't valid as stated

Both AI-search constructions I was able to obtain — Together's and SimpleTES's —
are **float64 evaluations of float64 iterates**, with neither the feasibility of
the iterate nor the value of the objective checked in exact arithmetic. (I have
not examined Haugland's 2016 construction, which may well have been rational;
this claim covers only the two witnesses checked here.)
`lp_research_state/code/evaluator.py` contains `bound_exact` /
`check_constraints_exact` for exactly this purpose — but nothing in the repo ever
calls them.

Doing it reveals that **Together's `h*` is not a feasible point**. Read as exact
dyadic rationals,

```
Σ h_i  =  185691005892807040627772388107 / 2^89  =  300 − 607096245493/2^89,
```

so `∫h = 1` fails by `9.81 × 10⁻¹⁶`. The float64 sum rounds to exactly 300.0,
which is why it was never caught. The defect is microscopic and harmless — a
mass repair moves `M` by at most `(2/n)·9.81e-16 ≈ 3 × 10⁻¹⁸` — but until it is
repaired the bound is not a theorem.

`lp_research_state/code/ub_certify.py` closes this. It snaps to an exactly
feasible rational point `h_i = a_i/2^60` (`Σa_i = nD/2` exactly, enforced by
integer repair) and evaluates all `2n−1` signed lags in integer arithmetic. A
priori rounding cost: `|∂M_j/∂h_k| ≤ 1` per cell, so snapping moves the scaled
objective by at most `1/D < 10⁻¹⁸`. Result:

```
µ ≤ 37969760969587463295692001413097862649
    ─────────────────────────────────────────  =  0.380870310586219904518562131404…
    99692099683868690467785529521025843200
```

with no floating-point step anywhere in the chain. Certified decimal, rounded
**up** as an upper bound must be: `µ ≤ 0.380870310586219904518562131405`,
checked `≥` the exact rational and tight to `9 × 10⁻³¹`.

**The same defect is systematic, not a Together quirk.** SimpleTES publish their
winning construction (`best_results/mathematics_extremal_analysis/erdos_minimum_overlap/`
in [wq-will/SimpleTES](https://github.com/wq-will/SimpleTES)) — a **2400-cell**
step function, which is what first suggested the finer-grid direction below. Scored
with the evaluator here (their data, not their code), it gives
`M = 0.3808676758273267`, reproducing their stored value bit-for-bit and confirming
their 0.380868 headline. It too is exactly infeasible — mass excess `+4.06 × 10⁻¹⁴`.

Certifying it the same way gives the **best exactly-certified upper bound
currently available**:

```
µ ≤ 607511972879059778836998329548027751357
    ─────────────────────────────────────────  =  0.380867675827326671773351376167…
    1595073594941899047484568472336413491200
```

certified decimal `µ ≤ 0.380867675827326671773351376168` (checked `≥` the exact
rational, slack `1.5 × 10⁻³¹`), improving the certified Together value by
`2.63 × 10⁻⁶`. The `0.380856` ablation remains the best *claimed* bound but has no
published witness, so it cannot yet be certified; `ub_certify.py` will do it in
~1 s once one appears.

### 3.2 The basin-uniqueness experiment

`PRO34_UB_REFINEMENT.md` ranks "basin-diversity search on the UB side" as
follow-on #2 and notes it had never been run: *"All known constructions may share
one basin … structured restarts would test whether 0.380871 is basin-unique."*

Run: `lp_research_state/code/ub_basin_sweep.py`, 1600 diverse starts at `n = 600`
across ten initialisation families (uniform, Bernoulli, block, structure-matched
fractional, smooth-Fourier, symmetric, antisymmetric, sparse, plus perturbations
and shuffles of `h*`), each by softmax-smoothed minimax descent, with the best 64
polished to first-order stationarity by trust-region SLP.

**Tier 1** (1600 descents, 160 per family) — best value reached by each family:

| family | best after descent |
|---|---|
| `perturb_hstar` | **0.381002933** |
| `smooth` | 0.381166876 |
| `blocks` | 0.381197519 |
| `antisym` | 0.381197768 |
| `uniform` | 0.381267643 |
| `symmetric` | 0.381312247 |
| `fractional` | 0.381377473 |
| `sparse` | 0.381415676 |
| `shuffle_hstar` | 0.381429654 |
| `bernoulli` | 0.381532769 |

The families span only `5.3 × 10⁻⁴` and every one of them is worse than `h*`.
Notably `shuffle_hstar` — a random permutation of `h*`'s own cell values, which
preserves the value *distribution* exactly and destroys only the arrangement —
lands at 0.381430, in the middle of the random pack. **What makes `h*` good is
its arrangement, not its value spectrum.**

**Tier 2** (top 64 SLP-polished to first-order stationarity): best
**0.380869107113270**, from `perturb_hstar`. This is `1.20 × 10⁻⁶` *below*
Together's own published value — a small but real improvement on their n=600
construction, obtained by perturb-and-repolish within their basin. It does not
approach SimpleTES's 0.3808676758 (`+1.43 × 10⁻⁶` above it). Exactly certified:

```
µ ≤ 0.380869107113271246758606524085
  = 7593928198568459773379995357147735909 / 19938419936773738093557105904205168640
```

**An important caveat about this run.** All 64 tier-1 points selected for deep
polish happened to come from `perturb_hstar`, because that family dominated the
descent ranking. So the run as executed shows that non-`h*` families are worse
*after descent* — it never polished one, and descent alone is far from converged
(residual certified gains ~6 × 10⁻⁶). `ub_basin_nonhstar.py` closes this by
regenerating the best tier-1 points of each non-`h*` family (the sweep is
seed-deterministic) and giving them the identical polish:

| family (best of 5, identically polished) | polished value |
|---|---|
| `uniform` | **0.380923941627155** |
| `symmetric` | 0.380941548201918 |
| `antisym` | 0.380942551946324 |
| `smooth` | 0.380948747087874 |
| `fractional` | 0.380951018397196 |
| `blocks` | 0.380968142045657 |
| `bernoulli` | 0.380974975146349 |
| `shuffle_hstar` | 0.380978578847031 |
| `sparse` | 0.380995007332619 |
| — | — |
| `perturb_hstar` (h\* basin) | **0.380869107113270** |

Under identical polish, **every non-`h*` family converges `5.48 × 10⁻⁵` worse than
the `h*` basin**, and they cluster tightly in `[0.380924, 0.380995]` — a broad
plateau of mediocre stationary points.

**Conclusion, stated at the strength the evidence supports.** 1600 diverse starts
with deep polish never produced anything competitive with `h*`. That is *not* a
proof of basin-uniqueness, and it cannot be: SimpleTES's n=2400 construction is
strictly better than `h*`, so better basins demonstrably exist. What the
experiment establishes is the weaker, still useful fact that **the good basin is
not reachable by unstructured multistart at this scale** — the landscape has a
wide, shallow plateau of near-`0.38095` optima that random restarts fall into.
This is consistent with every record construction (Haugland → AlphaEvolve →
TTT-Discover → Together → SimpleTES) descending from a shared lineage of
refinement rather than independent discovery, and it explains why UB progress has
come from search-procedure engineering rather than from re-rolling initialisations.

### 3.3 The current record construction is *not* first-order stationary

Applying the same trust-region SLP to SimpleTES's n=2400 construction shows it is
qualitatively less converged than Together's h*:

| point | radius | certified first-order gain |
|---|---|---|
| Together `h*` (n=600) | `10⁻⁴` | `1.94 × 10⁻¹⁰` |
| SimpleTES (n=2400) | `10⁻⁷` | `6.9 × 10⁻⁹` |

i.e. `h*` is ~35× tighter at a 1000× *larger* radius. Three SLP rounds moved
SimpleTES from `0.3808676758273267` to `0.3808676576575966`, a genuine
`1.8 × 10⁻⁸` improvement, so descent remains available in the published witness.

The reason it had never been polished out is structural: their construction
deliberately flattens the overlap profile across shifts (as their paper
describes), leaving **1580 of 4799 lags within `10⁻⁶` of the maximum**. That makes
the minimax LP massively degenerate — at a `10⁻³` trust radius it failed to
complete a single interior-point solve in 20 minutes. Only at radius `≤10⁻⁷`,
where the active set collapses to ~1200 lags, does it become tractable. This is
the same degeneracy PRO-34 flagged for `h*` at n=1200 ("HiGHS simplex effectively
hangs… always use IPM here"), one order worse.

Polishing it out (30 SLP rounds at radius `3 × 10⁻⁷`, 306 s) improves the published
witness by `1.30 × 10⁻⁷`, to `M = 0.3808675459609214`, and it is **still
descending** when the round budget runs out (residual certified gain
`2.7 × 10⁻⁹`/round, active set 1642 lags). Exactly certified:

```
µ ≤ 607511765732587328222592604544154703367
    ─────────────────────────────────────────  =  0.380867545960922320593700936552…
    1595073594941899047484568472336413491200
```

so the **certified UB ladder** at the end of this session is

| witness | certified `µ ≤` |
|---|---|
| Together n=600 | 0.380870310586219904518562131405 |
| this session's n=600 sweep best | 0.380869107113271246758606524085 |
| SimpleTES n=2400 (as published) | 0.380867675827326671773351376168 |
| **SimpleTES n=2400 + SLP polish (this session)** | **0.380867545960922320593700936552** |

Each row is an exactly-feasible rational witness with an integer-arithmetic
objective; each was re-checked as `≥` its own exact rational. None of them reaches
SimpleTES's *claimed* `0.380856`, whose witness is unpublished.


### 3.4 The lower-bound machinery reproduces

From a clean `.venv` (cvxpy 1.9.2 + CLARABEL, no Mosek license), the documented
single-row smoke test at `N=10000, T=4000, R=10, bochner_n=20` returns

```
rigorous_dual_LB = 0.379680000     (REPRODUCE.md expects >= 0.379653)   [86 s]
```

so the LB pipeline is intact and reproducible. Note this is the *smoke test*, one
row at reduced settings — not the headline. Reproducing `µ ≥ 0.3802838` end-to-end
requires the full 12-centre cover at `bochner_n=40, pm_k_max=20, hankel_n=6`
(~30 min, ~8 GB), which I did not run to completion here; that number remains
as documented, not as re-verified by me this session.

### 3.5 Independent replication of the stationarity certificate

Written from scratch against the mathematical definition (not by reading the
PRO-33 code), the trust-region SLP reproduces PRO-33's certified maximum
first-order gain at `h*` — **1.94 × 10⁻¹⁰** at radius `10⁻⁴` — to three
significant figures. PRO-33/34's stationarity finding replicates.

### 3.6 The lag-sign trap is live

The repo memory flags "sup evaluated over positive lags only" as a recurring bug
class that has previously poisoned a 50-digit anchor. My *first* draft of
`ub_core.py` hit it: `np.correlate(h, 1-h, 'full')` mirrors the lag axis — index
`m` is lag `(n−1) − m`, not `m − (n−1)`. The error was invisible in the objective
(a max over all lags is insensitive to mirroring) and only surfaced as an 8%
gradient/finite-difference mismatch. `ub_core` now exposes `lag_of_index` /
`index_of_lag` and a `selftest` that asserts the convention against a brute-force
double loop before any sweep is allowed to run. **The warning is worth keeping.**

---

## 4 · What would actually be required to solve it

Not achievable by more of what is here. Specifically:

**On the lower-bound side.** The Lever I′ saturation analysis says the current
architecture — cell-envelope LP + Bochner-PSD moment constraints + White's
ellipse extension — cannot certify past `C_explicit = 0.380713`. That is
`1.5 × 10⁻⁴` short of the current UB even if executed perfectly at infinite `N`.
Closing the problem needs an LB architecture that does *not* factor through
cell-envelope + Bochner-PSD duals. Nothing in the literature supplies one; ten
candidate transfers (Lasserre-2, White's L² autoconvolution machinery,
Li/Cohn–de Laat–Salmon discrete reduction, Beurling extremal functions,
noncommutative moment hierarchies, …) have been examined here and all failed.

**On the upper-bound side.** Progress is now driven by AI search over step
functions, moving roughly `10⁻⁵` per iteration of the benchmark. Nothing suggests
convergence to a recognisable closed form; PSLQ hunts on the high-precision
bracket found no relation.

**The structural obstruction.** `µ = inf_n M_n` where each `M_n` is a nonconvex
minimax over `n` variables. Finite `n` yields only upper bounds — exact `M(n)`
values (certified here through `n = 21`) cannot bound µ from below. So the two
sides are attacked by disjoint machinery with no known bridge, and the honest
statement is that the gap is not currently closable from either end.

### What *is* worth doing next, in priority order

1. **Publish the lower bound now.** `µ ≥ 0.3802838` leads the published record by
   `1.16 × 10⁻³`, and Kim–Pilanci (June 2026) demonstrates a second group working
   the certified-LB direction with the same rigor standard. This is the one place
   where this repo holds a large, defensible, time-limited lead. The blocking item
   is unchanged from `PROGRESS_AND_SIGNIFICANCE.md`: arbitrary-precision
   re-certification of the binding rows.
2. **Finish the polish of the record construction.** It is still descending at
   `2.7 × 10⁻⁹`/round after 30 rounds. Running it to convergence — and repeating on
   whatever construction currently tops Einstein Arena — is cheap and would put a
   certified value on the actual record. Use `slack ≤ 2 × 10⁻⁷`; larger trust radii
   make the LP intractably degenerate.
3. **Offer exact certification as a service to the UB community.** No AI-search
   entrant appears to check exact feasibility, and both witnesses examined here
   failed it. `ub_certify.py` settles any construction in ~1 s.
4. **Do not** run more Bochner/Lasserre level scans, more `N`-scaling, or more
   random multistart at n=600. All three are now measured dead ends.

---

## 5 · Reproducing this session's results

Environment: `python3 -m venv .venv && .venv/bin/pip install cvxpy clarabel numpy
scipy mpmath sympy highspy`. All commands from `lp_research_state/code`.

```bash
../../.venv/bin/python ub_core.py
```
```bash
../../.venv/bin/python ub_certify.py --out ../data/ub_certified_together.json
```
```bash
../../.venv/bin/python ub_certify.py --input ../data/simpletes_polished.json --out ../data/ub_certified_polished.json
```
```bash
../../.venv/bin/python ub_basin_sweep.py --n 600 --starts 1600 --deep 64 --rounds 150 --out ../data/ub_basin_sweep_n600.json
```
```bash
../../.venv/bin/python ub_basin_nonhstar.py --per-family 5 --rounds 150 --out ../data/ub_basin_nonhstar.json
```

New code: `ub_core.py` (objective, lag convention, `selftest` — run it before
trusting anything else here), `ub_local.py` (smoothed-minimax descent +
trust-region SLP), `ub_certify.py` (exact rational certification),
`ub_basin_sweep.py`, `ub_basin_nonhstar.py`, `ub_fine_search.py`.
New data: `ub_certified_{together,sweepbest,simpletes,polished}.json`,
`ub_basin_sweep_n600.json`, `ub_basin_nonhstar.json`,
`simpletes_construction.json`, `simpletes_polished.json`.

**Caveat on scope.** The headline `µ ≥ 0.3802838` was not re-derived end to end
this session — only the documented smoke test was run. Everything in §3.1–3.3 is
new computation reported here; §2 is literature verified against primary sources.
