# μ High-Precision Values + PSLQ Search Results

**Goal of this document:** Provide the most precise machine-verified values
of μ_LB and μ_UB, list all natural closed-form candidates within the bracket,
and report the PSLQ search outcome.

**Status (2026-05-17):** Bracket is `[0.3801279, 0.3808703]`. Width
`7.4 × 10⁻⁴` is too wide for definitive closed-form identification.

---

## 1. Best-known precision values

### μ_LB (rigorous, post-margin)

```
μ ≥ 0.3801279                          (Phase 5 CDE, with margin 1e-6 + Lipschitz eps_grid ~ 2.15e-6)
μ ≥ 0.3801319                          (predicted with T5p re-iteration through Phase 5)
```

The pre-margin Phase 5 LB is closer to `0.3801289` (numerical), but for
mathematical rigor the post-margin value is the citable LB.

### μ_UB (numerical Together certificate)

Computed at 40 mpmath dps from Together's 600-step h* (loaded from
`lp_research_state/data/together_f_star.json`):

```
μ ≤ 0.38087031058621990365023943115135072229730000000000000000  (mpmath, 40 dps)
```

(Optimal correlation shift `j = -33`. Precision bounded by Together's h*
values being given as ~16-digit floats, so the meaningful precision is
~16 digits; mpmath produces additional decimals that simply repeat the
last meaningful 0.)

### Bracket

```
μ ∈ [0.3801279, 0.3808703105862199...]
Width: 7.4 × 10⁻⁴
```

---

## 2. PSLQ search results

### 2.1. Direct PSLQ on bracket points

`mpmath.pslq([μ, π, log(2), sqrt(2), ...])` was run on:
- μ_LB
- μ_UB
- Midpoint
- 1st and 3rd quartiles of bracket

**Outcome:** PSLQ returns "matches" for every point with various basis sets,
but the matches are spurious (large denominators 7-333, no consistency
across basis variations).

`mpmath.identify()` produces similar uninformative output.

### 2.2. Brute-force enumeration of "natural" closed forms in bracket

Searched the families:

| Family | Coefficient bound | Distinct values in [μ_LB, μ_UB] |
|---|---|---|
| `arctan(p/q)` (coprime) | `p, q ≤ 100` | **1**: `arctan(2/5)` |
| `log(p/q)` (coprime) | `p, q ≤ 200` | **8** |
| `nπ/d` (coprime) | `n ≤ 30, d ≤ 200` | **6** |
| `n·√k/d` (k ∈ {2,3,5,7}) | `n, d ≤ 50` | **3** |
| `n/(d·π)` (coprime) | `n, d ≤ 50` | **4** |
| **Total distinct** | | **~22** |

Number of candidates is consistent with random density (~one candidate
per family per 10⁻⁴-wide window), so no candidate is statistically
distinguished by being a "match."

### 2.3. Notable candidates (by elegance and bracket position)

```
arctan(2/5)    = 0.38050637711236487    [midpoint-ish, margin +3.8e-4 / -3.6e-4]
log(117/80)    = 0.38014730012387394    [very close to LB, margin +2.0e-5]
4π/33          = 0.38079910952603561    [very close to UB, margin -7.1e-5]
19π/157        = 0.38019274151723629    [near LB, margin +6.5e-5]
9·√3/41        = 0.38020627483219330    [near LB, margin +7.8e-5]
7·√2/26        = 0.38074980525429496    [near UB, margin -1.2e-4]
log(79/54)     = 0.38046380590274712    [middle]
```

### 2.4. PSLQ verdict

**No closed-form candidate is privileged.** The bracket admits roughly
22 distinct natural candidates and gives no tool for distinguishing
them. The "match" of any specific candidate (e.g., `arctan(2/5)`)
within the bracket is *not evidence* — it's a property of the bracket
being wide enough to admit many candidates.

---

## 3. What would be needed for a definitive PSLQ hit

A closed-form identification requires the bracket to be `≤ 10⁻⁵` wide
(empirically; based on the density of small-coefficient closed forms).
This would shrink the candidate set to `O(1)`.

To achieve width `10⁻⁵`:
- Need to either push μ_LB up by ~7×10⁻⁴ (not in-session tractable)
- Or push μ_UB down by ~7×10⁻⁴ (more Together-style compute on the UB side)
- Or both (each ~half the gap)

Both are non-trivial; neither is currently in scope.

---

## 4. Adjacent result: μ_UB to 16 digits

The Together certificate computed to mpmath-precision gives:

```
μ_UB  =  0.3808703105862199036502394311513507222973...
                                  ^^^^^^^^^^^^^^^^^^^^^^^^
                                  beyond 16 digits is float repetition
```

Reliable to ~16 digits (Together's h* values are given as ~16-digit floats).

This is a *modest* improvement over the previously cited "μ ≤ 0.380871"
(6 digits) — gives us 10 more digits of UB. By itself this is not
publishable, but combined with a sufficiently improved LB it could be
PSLQ-relevant.

---

## 5. Conclusion + recommendation

**Recommendation:** Defer Approach A (PSLQ closed-form hunt) until either
- μ_LB is improved by ≥ 5×10⁻⁴, OR
- μ_UB is improved by ≥ 5×10⁻⁴.

**Cheap pivot:** if Approach B (AI-driven constraint discovery) finds
ANY new constraint family that pushes μ_LB by 5×10⁻⁴, the bracket narrows
to ~2×10⁻⁴, and PSLQ becomes meaningfully informative.

The mpmath + Together-UB infrastructure built in this session is
reusable for that future hunt.

<promise>PSLQ_DONE</promise>
