# Lever C — M(n) reconnaissance + brute-force extension

**Date:** 2026-05-10
**Status:** DONE
**Verdict:** TOGETHER VALIDATED

## Question

Lever C of the diagnostic: can pushing M(n) computations further produce a tighter upper bound on µ than Together's `0.380871`, or confirm the published UB is the best available?

Recall: µ = lim inf_{n→∞} M(n)/n. Any specific n with `M(n)/n < 0.380871` would yield an immediate, integer-combinatorial upper bound tighter than the Together step-function construction.

## Known M(n) table (sorted by ratio, ascending)

All published exact values are for n ≤ 15. Wikipedia and HandWiki both stop there and attribute the values to Haugland (1996, 2016) and Guy *Unsolved Problems in Number Theory* (2004). Our session brute force (B&B + reflection symmetry, code at `lp_research_state/code/_brute_force_Mn_extended.py`) **independently reproduces every published value for n=2..15** and extends the table to n=18 (timed-out attempt at n=19 in progress; see below).

| n  | M(n) | M(n)/n   | Source                          |
|---:|-----:|---------:|---------------------------------|
| 15 |    6 | 0.400000 | Wikipedia / Haugland (published) |
| 12 |    5 | 0.416667 | Wikipedia / Haugland (published) |
| 17 |    7 | 0.411765 | session brute force              |
|  7 |    3 | 0.428571 | Wikipedia / Haugland (published) |
| 14 |    6 | 0.428571 | Wikipedia / Haugland (published) |
| 16 |    7 | 0.437500 | session brute force              |
|  9 |    4 | 0.444444 | Wikipedia / Haugland (published) |
| 18 |    8 | 0.444444 | session brute force              |
| 11 |    5 | 0.454545 | Wikipedia / Haugland (published) |
| 13 |    6 | 0.461538 | Wikipedia / Haugland (published) |
|  2 |    1 | 0.500000 | Wikipedia / Haugland (published) |
|  4 |    2 | 0.500000 | Wikipedia / Haugland (published) |
|  6 |    3 | 0.500000 | Wikipedia / Haugland (published) |
|  8 |    4 | 0.500000 | Wikipedia / Haugland (published) |
| 10 |    5 | 0.500000 | Wikipedia / Haugland (published) |
|  5 |    3 | 0.600000 | Wikipedia / Haugland (published) |
|  3 |    2 | 0.666667 | Wikipedia / Haugland (published) |
|  1 |    1 | 1.000000 | Wikipedia / Haugland (published) |

**Smallest known ratio:** M(15)/15 = **0.400000**.

Persisted: `lp_research_state/data/known_Mn_values.json` (published table + provenance), `lp_research_state/data/Mn_optimizers_large.json` (our brute-force data with optimizer sets and timings).

## Comparison to Together's UB 0.380871

| Metric                   | Value     |
|--------------------------|----------:|
| Smallest known M(n)/n    | 0.400000 (n=15) |
| Together's UB on µ       | 0.380871  |
| Gap                      | +0.019129 (the integer ratio is **looser**) |
| Does any known n beat Together? | **No** — by a margin of ~1.9 × 10⁻². |

**Zero values of n in the published table or in our extension to n=18 yield a ratio below 0.380871.** µ being a *lim inf* means individual M(n)/n can sit well above the limit; the published constant 0.380871 is a strictly better UB than any integer M(n)/n could give at currently-tractable n. The smallest known ratio at n=15 (0.40) is +5% above Together's bound.

## Our brute force extension

| n  | M(n) | M(n)/n | Wall time | Status |
|---:|-----:|-------:|----------:|--------|
| 16 |    7 | 0.4375 |  53.9 s   | certified |
| 17 |    7 | 0.4118 |  65.5 s   | certified |
| 18 |    8 | 0.4444 | 664.1 s   | certified |
| 19 |   ≤8 | ≤0.421 | 55 s — UB found, 1200 s — see Mn_optimizers_large.json | timed out (UB only) |

(Pre-existing data at `lp_research_state/data/Mn_optimizers_large.json` was retained; n=19 row was appended by this session's longer-budget attempt.)

Runtime scaling n=17 → n=18 was 10×. n=19 with the 20-minute budget completed without certifying optimality (an UB of M=8 was found but the search did not exhaust). Projecting n=20 from the geometric trend would require ~2 hours single-threaded, n=25 would be days-to-weeks, n=43 (Haugland-style reach via SAT) is overnight-cluster territory. **Even if achieved, the result would not move past the verdict above:** all known and observed M(n)/n ratios for n ≤ ~25 are dominated by floor effects (M(n) grows like ⌈µn⌉ which for µ ≈ 0.38 means the next M-jump happens roughly every n + 2..3 steps; ratios live in `[0.38, 0.44]` until n is large enough that the rounding error 1/n shrinks below `0.4 − 0.38 ≈ 0.02`, i.e. n ≥ 50).

## Verdict

**TOGETHER VALIDATED.**

- Sub-1: No published or computed M(n)/n is below 0.380871. The smallest known is 0.400 at n=15.
- Sub-2: Extended brute force certifies M(n) up to n=18 (matching all published values); n=19 budget exceeded.
- Sub-3: Combined with the prior diagnostic conclusions (Levers A, D, D′, E ruled out), and the lower-bound saturation at µ ≥ 0.3801279 (CDE Phase 5), Lever C confirms that **Together's `µ ≤ 0.380871` step-function construction stands as the best available UB**. The remaining open gap is `µ ∈ [0.3801279, 0.380871]`, width ≈ 7.4 × 10⁻⁴.

The framework — both on the LB side (SDP + Bochner + poly-moment + cover refinement) and on the UB side (continuous step-function constructions and integer brute force) — is at its current technical limit. Closing the remaining gap requires new mathematical levers, not more compute on the existing techniques.

## Honest caveats

- The integer route to µ is *asymptotically* informative, not bound-improving at small n. Even if M(n)/n at n=50 turned out to equal 19/50 = 0.380, that gives `µ ≤ 0.380` as a true UB — better than Together's 0.380871! Whether that ever happens is open. But: pushing the brute-force route from n=18 to n=50 is a research-grade open question (Haugland's 2016 paper reached n ≤ 43 only via SAT-style + symmetry tricks). Single-threaded B&B at our current implementation is not the path to n=50.
- A "needs more compute" verdict is *not* ruled out in principle. A SAT-encoded run, OEIS-style published computational tables (which we couldn't access through Cloudflare in this session), or a Haugland-2016 replication might extend the table. None changes the present session's conclusion — at every n we *can* check, Together's bound is tighter.
- We checked: `/tmp/together_repo/erdos-minimum-overlap/` contains **only step-function constructions** (`haugland_2016.py`, `alphaevolve_2025.py`, `together_ai_2026.py`, `ttt_discover_2026.py`) — *no* integer M(n) tables.

## Files touched

- `lp_research_state/data/known_Mn_values.json` (new) — published M(n) table + sources + comparison to Together UB
- `lp_research_state/data/Mn_optimizers_large.json` (updated with n=19 row)
- `LEVER_C_RESULT.md` (this file)

## Sources

- Wikipedia — Minimum overlap problem: https://en.wikipedia.org/wiki/Minimum_overlap_problem
- HandWiki — Minimum overlap problem: https://handwiki.org/wiki/Minimum_overlap_problem
- Haugland (1996), "Advances in the Minimum Overlap Problem," *J. Number Theory* 58:71-78
- Haugland (2016), "The minimum overlap problem revisited," arXiv:1609.08000
- Together (2026), GitHub: https://github.com/togethercomputer/erdos-minimum-overlap (UB 0.380871)
- White (2023), *Acta Arith.*; arXiv:2201.05704 (LB 0.379005, baseline)
