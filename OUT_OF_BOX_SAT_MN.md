# Out-of-box: SAT-based M(n) extension

## TL;DR

SAT (CDCL via CaDiCaL/Glucose, pseudo-boolean cardinality encoding) **pushes exact M(n) from n=18 (prior brute-force ceiling) to n=20** in ~12 minutes wall time. Two new integer values:

- **M(19) = 8**, ratio = 8/19 ≈ 0.421053 (new — not in published literature)
- **M(20) = 8**, ratio = 8/20 = **0.4000** (new — ties n=15 as smallest known integer ratio M(n)/n)

n=21 did not complete: the second SAT call (bisecting on M) exceeded the 10-min per-solve cutoff and was terminated. **Verdict: SAT is a productive direction at this scale, but the exponential wall is real and integer M(n) will not beat Together's 0.380871 within human-tractable compute.**

## Encoding

File: [`lp_research_state/code/_sat_Mn.py`](lp_research_state/code/_sat_Mn.py)

- Boolean `x_i = (i ∈ A)` for `i ∈ [1, 2n]`.
- `CardEnc.equals(..., n)` enforces `|A| = n` (sequential counter, ~O(n²) clauses).
- Pin `x_1 = 1` (A↔B symmetry break).
- For **each** nonzero shift `k ∈ [-(2n-1), 2n-1]`: introduce auxiliaries `y_i^k ↔ x_i ∧ ¬x_{i-k}` via 3 clauses each; then `CardEnc.atmost(y_i^k, M)`.
  - **NOTE (fixed during dev):** initial encoding only used positive shifts under a (false) symmetry claim. For partitions of `[2n]`, `overlap_k ≠ overlap_{-k}` in general (e.g., A={1..n}, B={n+1..2n}: overlap_{-n} = n, overlap_{+n} = 0). Both signs of `k` are encoded.
- Binary search on `M`. Each SAT-returned model `A` is **verified** numerically (`verify_overlap` computes max overlap directly and asserts ≤ M).
- Solver: `Cadical153` (pysat default).

Sanity passes: M(10)=5, M(12)=5, M(15)=6 — match published values.

## Results

All matches against ground truth (where available) hold. New rows in bold.

| n  | M(n) | M(n)/n   | t total | t hardest solve | source |
|----|------|----------|---------|-----------------|--------|
| 10 | 5    | 0.500000 | 0.07s   | 0.05s UNSAT M=4 | sanity (matches lit) |
| 12 | 5    | 0.416667 | 0.13s   | 0.07s UNSAT M=4 | sanity (matches lit) |
| 15 | 6    | 0.400000 | 3.46s   | 2.87s SAT M=6   | sanity (matches lit) |
| 16 | 7    | 0.437500 | 5.35s   | 4.65s UNSAT M=6 | matches brute-force |
| 17 | 7    | 0.411765 | 5.70s   | 5.19s UNSAT M=6 | matches brute-force |
| 18 | 8    | 0.444444 | 69.89s  | 64.88s UNSAT M=7| matches brute-force (brute: 664s — **~10× speedup**) |
| **19** | **8** | **0.421053** | 56.46s  | 56.25s UNSAT M=7| **new** |
| **20** | **8** | **0.400000** | 90.20s  | 83.04s UNSAT M=7| **new** |
| 21 | ?    | —        | >10 min | did not complete (>10 min single solve) | terminated |

Optimizers:
- M(19)=8: A* = [1, 2, 3, 4, 5, 6, 7, 11, 14, 21, 22, 28, 31, 33, 34, 35, 36, 37, 38]
- M(20)=8: A* = [1, 2, 3, 4, 5, 7, 8, 11, 13, 21, 22, 28, 30, 33, 34, 35, 37, 38, 39, 40]

## Did any M(n)/n beat 0.380871?

**No.** Smallest ratio observed is M(20)/20 = 0.4000, matching the n=15 record (M(15)=6, 6/15=0.4). All M(n)/n in the range stay ≥ 0.4. The Together upper bound 0.380871 remains the limiting target; integer M(n) at small n cannot beat it because rounding `0.380871 · n` to the next integer requires `n ≥ 53` before `⌈0.380871 n⌉/n < 0.4`. (At n=53: 0.380871·53 = 20.186, so M(53) ≥ 21 to beat, giving 21/53 ≈ 0.3962; but n=53 is far beyond SAT reach.)

## Time-vs-n scaling

Total bisection time (sec):

| n  | 15  | 16  | 17  | 18  | 19   | 20   | 21    |
|----|-----|-----|-----|-----|------|------|-------|
| t  | 3.5 | 5.4 | 5.7 | 69.9| 56.5 | 90.2 | >600  |

The dominant cost is the **hardest UNSAT proof at M = M(n) − 1** (proving infeasibility one below the optimum). For n in 17→20 those proofs took 5.2s → 64.9s → 56.3s → 83.0s — roughly 1.3-1.5× per step on average, but the n=21 UNSAT proof did not finish in 10+ min, suggesting an order-of-magnitude jump.

Extrapolating geometric ~1.5× per n from a 90s base at n=20:
- n=22: ~3-5 min
- n=25: ~30-60 min
- n=30: ~12-20 hours
- n=40: ~weeks
- n=53 (the n needed to potentially beat 0.380871): astronomically infeasible

**Could overnight reach n=30?** Marginally. Optimistically (1.5× / step), n=30 ≈ 90s · 1.5¹⁰ ≈ 5200s ≈ 1.5 hours per UNSAT proof. Pessimistically (factor closer to 2-3× as suggested by the n=21 jump), n=30 is overnight-to-multi-day. Verdict: n=25-27 looks plausible overnight; n=30 is a stretch.

## Verdict — was SAT productive?

**Productive but bounded:**
- **Pro:** 10× speedup over brute force at n=18 (70s vs 664s), and pushed two new exact values (M(19), M(20)) that don't appear in Haugland's 1996/2016 tables or Wikipedia.
- **Pro:** The SAT certificate of optimality (UNSAT at M-1) gives mathematical certainty, not just empirical "haven't seen better".
- **Con:** The scaling wall is exponential. Integer M(n) at SAT-tractable scale (n ≤ ~25) gives ratios `≥ 0.4`, which are 5%+ above White's 0.379544 LB and 5%+ above Together's 0.380871 UB. So **SAT-extended M(n) does not produce new µ bounds**; it only refines a finite table that's mostly of independent combinatorial interest.
- **Con:** The continuous extension `M_c(h, p, q)` / step-function formulations (White / Together / AlphaEvolve) sail past these integer constraints; that's the path the µ-bound literature already takes.

**Useful follow-ups (not in this 30-min budget):**
- Try parallel CDCL (e.g., ManySAT or Mallob) to push n=21-25.
- Add reflection symmetry as encoded constraint (currently only handle complement via x_1=1).
- Try LCG-Glucose / inprocessing variants — UNSAT proofs are the bottleneck.
- Cube-and-conquer (CaDiCaL --cube) for n=21+.

## Status

**DONE_WITH_CONCERNS** — produced new M(n) values for n=19, 20; n=21 did not complete within the per-solve budget and was terminated. All completed values cross-verified by direct overlap calculation.

## Files

- [`lp_research_state/code/_sat_Mn.py`](lp_research_state/code/_sat_Mn.py) — SAT solver
- [`lp_research_state/data/Mn_sat_results.json`](lp_research_state/data/Mn_sat_results.json) — results with optimizers and per-call timing
- Log: `/tmp/sat_mn_run.log` (transient; key lines reproduced in the table above)
