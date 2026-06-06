# PSLQ Tooling Notes — Approach A (closed-form hunt for μ)

**Status:** Tooling assessed and installed. **Verdict:** PSLQ approach is **NOT
productive at the current bracket width**. Need a factor-100 narrowing of the
[μ_LB, μ_UB] bracket before integer-relation finding can produce a definitive
hit.

---

## 1. Tooling state (post-install)

Project venv: `/Users/benzanghi/Documents/Claude/Projects/Erdos/.venv/bin/python`.

| Package | Pre-session | Post-install | Purpose |
|---|---|---|---|
| `numpy` | ✓ 2.4.4 | ✓ | numerical |
| `scipy` | ✓ 1.17.1 | ✓ | linear algebra |
| `cvxpy` | ✓ 1.8.2 | ✓ | SDP modeling |
| `clarabel` | ✓ 0.11.1 | ✓ | SDP solver |
| `scs` | ✓ 3.2.11 | ✓ | SDP solver (cross-check) |
| `mpmath` | ✗ | ✓ 1.3.0 | PSLQ + arbitrary-precision arithmetic |
| `sympy` | ✗ | ✓ 1.14.0 | symbolic algebra |
| `gmpy2` | ✗ | ✗ | (not installed; mpmath sufficient) |
| `sdpa-gmp` | ✗ | ✗ | (would need brew/source install; deferred) |

Tooling sanity-check passed: `mpmath.pslq([π, atan(1)])` recovers the
relation `[1, -4]` (π = 4·atan(1)) at 50 dps.

## 2. Achievable precision summary

| Source | Precision | Notes |
|---|---|---|
| μ_LB (Phase 5 rigorous, post-margin) | 7 digits | `0.3801279` |
| Row-4 SDP value (CLARABEL, N=3000, n=20) | 5 digits cross-solver | CLARABEL: 0.37879232684054, SCS: 0.37880187769190, diff ~10⁻⁵ |
| μ_UB (Together's certificate, mpmath) | ~16 digits | `0.3808703105862199036502394311513507222973` (limited by Together's float h* values) |
| Bracket width [μ_LB, μ_UB] | 7.4 × 10⁻⁴ | — |

**Bottleneck:** the SDP solvers (CLARABEL, SCS) at the highest tractable
scale (Phase 5: N=10000, T=4000, n=30) give the row-binding value to
~5-7 digits at best, with `optimal_inaccurate` status. The rigorous Phase 5
LB after margin (`0.3801279`) is conservative by ~10⁻⁶. Together's UB is
known to ~16 digits.

## 3. Why PSLQ fails at this width

PSLQ finds integer relations between real numbers given to *N digits*. To
constrain μ to a specific closed-form expression, we need the bracket
width to be **smaller than the typical separation between candidate
expressions**:

```
required-bracket-width  ≈  10^(-digits-of-precision-needed)
                       ≈  10^(-distinctness-of-closed-forms)
```

For closed-form candidates with coefficient magnitude `≤ K` in basis size
`B`, the number of distinct expressions per unit interval scales as
`K^B`. To distinguish them, the bracket must be `≤ 1/K^B`.

For modest `K=50, B=2` (e.g., `a·π + b·log(c)` with `a, b, c ≤ 50`), this
demands `width ≤ 10⁻⁵` at minimum. **Current bracket width `7.4 × 10⁻⁴` is
~100× too wide.**

## 4. Candidate enumeration (the actual finding)

Despite PSLQ not being decisive, we can list all "natural" candidate
closed forms within the current bracket. From a brute-force search over
small-coefficient expressions:

| Family | Candidates in [μ_LB, μ_UB] | Closest to midpoint |
|---|---|---|
| `arctan(p/q)`, p,q ≤ 100 | many (all equivalent to arctan(2/5)) | `arctan(2/5) = 0.3805063771...` |
| `log(p/q)`, p,q ≤ 200 | 8 distinct values | `log(79/54) = 0.3804638059...` |
| `n·π/d`, n ≤ 30, d ≤ 200 | 6 distinct | `19π/157 = 0.3801927415...` |
| `n·√k/d` | 3 distinct | `9·√3/41 = 0.3802062748...` |
| `n/(d·π)` | 4 distinct | `43/(36π) = 0.3802034752...` |

**Total natural closed forms in the bracket: ~25 distinct values.** This is
consistent with random density at the chosen complexity bound; no
candidate is statistically privileged.

The single most "elegant" candidate (smallest integers, simplest function):

> **`arctan(2/5) = 0.38050637711...`**

It sits comfortably inside the bracket (margin `+3.8 × 10⁻⁴` from LB,
`-3.6 × 10⁻⁴` from UB) and uses minimal coefficients. **There is no
mathematical reason to privilege this candidate** — it's only the
"shortest description" candidate in the list.

## 5. mpmath.identify() output is uninformative

Running `mpmath.identify()` on points in `[μ_LB, μ_UB]` with basis sets
`{π, e, log(2), sqrt(2), ...}` produces expressions with rational
coefficients (mostly denominators 7-333) for *every* point. These are
spurious matches — `identify()` is essentially fitting noise.

This confirms the bracket is too wide for tolerance-based identification.

## 6. Path forward (for future sessions, if attempted)

To make Approach A productive, one of:

1. **Push μ_LB up** by ≥ 5 × 10⁻⁴. Requires the saturation theorem from
   LEVER_I_PRIME_THEOREM.md to be proved at significantly larger `N`
   (likely > 30,000), AND the lever stack to extend with new constraint
   families. Discussed at length in SESSION_FINAL.md; ruled out as
   in-session-tractable.

2. **Push μ_UB down** by ≥ 5 × 10⁻⁴. Requires extending Together's
   step-function certificate from 600 to several thousand steps (more
   compute, no math). Not our work; would require
   replicating/contributing to Together's repo.

3. **Use structure-specific reasoning.** If μ has a closed form, it
   likely arises from a specific functional equation (e.g., the optimum
   of `inf_h max_t (h*h)(t)` over piecewise-constant `h` satisfies an
   Euler-Lagrange equation). Solving that equation analytically — if
   tractable — would give a candidate that PSLQ could verify even at
   current precision.

4. **Combine A with B (AI-driven constraint discovery).** If an
   ML-discovered constraint pushes `μ_LB` significantly, the bracket
   narrows and PSLQ becomes viable.

## 7. Honest summary

This approach has been correctly assessed as **non-productive at the
current bracket width**. The tooling investment (mpmath, sympy) is
*reusable* and was the cheapest possible way to determine this. The
candidate enumeration is a *useful baseline*: should future bracket
narrowing produce a value matching one of these candidates, the
match is *not unique* and a priori must be cross-checked with multiple
candidates.

**No closed-form identification is possible at this stage.**
