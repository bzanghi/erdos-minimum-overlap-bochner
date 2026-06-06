# PRO-22: Direct sup_t SDP — NEGATIVE RESULT

**Status:** FAILED. The contrarian attack surface from the round-3 analysis turned out to be invalid as a pure replacement for the cell-envelope.

---

## What was tried

Replaced the cell-envelope cosine + sine constraints (lines 176-190 of `white_full_convex.py`) with direct constraints `M(t) ≤ Ω` at a grid of 201 shifts `t ∈ [0, 2]`, where M(t) is computed from the (a_m, b_m) Fourier variables for m=1..2R.

The hypothesis: cell-envelope is the bottleneck of the framework saturation theorem (PRO-6's `C_∞ ≈ 0.380558`), so replacing it might break the ceiling.

## What happened

Build_problem_supt at small scale (N=200, T=200, R=5, bn=8):

| SDP | Ω | Δ vs original |
|---|---|---|
| Original (cell-envelope cell-min) | 0.366133 | — |
| **PRO-22 (direct sup_t)** | **0.370557** | **+4.42 × 10⁻³** |

The `+4.42e-3` apparent improvement looked promising — substantially larger than the open gap.

**T-stability check** (does the SDP exploit the missing high-m tail?):

| T | Ω_supt |
|---|---|
| 50 | 0.3705570 |
| 100 | 0.3705570 |
| 200 | 0.3705570 |
| 400 | 0.3705570 |
| 800 | 0.3705570 |

Stable across 16× variation in T — suggesting the existing constraints (Parseval, |c|/|d| ≤ 2/π, Bochner-PSD) already implicitly control the tail.

## The killing test: validity check

Solved the supt-SDP, extracted optimal (c, d), reconstructed f via Fourier, computed actual `M(t) = (f * f)(t)` numerically, and compared to the SDP's claimed Ω.

```
supt-SDP reported Ω        =  0.370557
True Ω from reconstructed f =  0.929000    ← 5× LARGER
diff                        =  +0.558
```

The reconstructed f has range `[-1.52, +1.69]` — completely outside `[0, 1]`. The Bochner-PSD constraint at `bochner_n=8` was insufficient to enforce `f ∈ [0, 1]`.

**Conclusion:** the supt-SDP is **NOT** a valid lower bound. The cell-envelope was doing essential WORK that the direct M(t) representation does not replicate: it constrains (a, b, c, d) to come from a valid f ∈ [0, 1], not just bounds M's cell averages.

## Why the cell-envelope is multi-purpose

Going back to White's derivation: the cell-envelope constraint

```
(L/2) · α_m^-(j) · (w + v) + 2(a_m² + b_m²) - (4/πm) sin(πm/2) · a_m ≤ 0
```

is the SDP encoding of (W.1):

```
∫_0^2 M(x) cos(πmx/2) dx = (4 sin(πm/2)/πm) · a_m - 2(a_m² + b_m²)
```

which is a Parseval-like identity relating `M's m-th Fourier coefficient`
to `f's m-th Fourier coefficient via |f̂(m)|²`. **Removing this constraint
breaks the link** between (a, b) [encoded as Fourier of M] and (c, d) [Fourier
of f]. Without it, (a, b) are unconstrained by (c, d), and the SDP can choose
(a, b) values that make M(t) small (satisfying the direct sup_t constraints)
while (c, d) implies a different M with much larger sup.

## Why T-stability is a red herring

T-stability LOOKED like evidence of validity, but it actually means: the SDP
is exploiting the FREEDOM in (a, b) that's decoupled from (c, d). The
"high-m tail" isn't the problem — the LOW-m a, b values themselves are wrong
(don't match the actual Fourier coefficients of M for the SDP's chosen f).

## What this teaches us about the framework ceiling

PRO-6's saturation theorem said: at fixed Bochner level n, the cell-envelope
augmentation has residual `r_C(n) → 0` as N→∞. With the cell-envelope
DROPPED entirely, the SDP becomes invalid — so the cell-envelope is the
bottleneck NOT because of its relaxation gap, but because it's NECESSARY for
validity.

The framework ceiling `C_∞ ≈ 0.380558` is therefore a more fundamental limit
than I initially thought. It reflects the tightest LB achievable by:
- Keeping the cell-envelope (necessary for validity)
- Making it as tight as possible (exact integral vs cell-min)
- With Bochner-PSD truncated at any n

The ~3.1×10⁻⁴ "beyond-framework" portion of the open gap is genuinely beyond
the reach of this constraint structure.

## Implications for further work

1. **The cell-envelope cannot simply be replaced.** Any new attack surface
   must either keep the cell-envelope or find an equivalent (a, b) ↔ (c, d)
   link that establishes validity.

2. **A HYBRID approach is the right idea:** keep the cell-envelope as-is,
   ADD direct M(t) ≤ Ω at a t-grid AS AUXILIARY CONSTRAINTS. The new
   constraints would only tighten the SDP (the cell-envelope is unchanged,
   so validity is preserved). But the existing cell-envelope ALREADY
   implies M(t) ≤ Ω + O(L) at all t in any cell, so the auxiliary
   constraints would only tighten by O(L) — small improvement at our N.

3. **The cell-envelope's role is essentially Parseval enforcement.** White's
   (W.1) identity is THE link between Fourier coefficients of M and of f. We
   can't avoid this link without violating validity.

4. **The framework ceiling is genuine.** Closing the beyond-framework 55% of
   the open gap requires methodology FUNDAMENTALLY OUTSIDE the SDP+PSD+Fourier
   relaxation paradigm.

## Honest summary

- **Goal:** test whether the cell-envelope can be replaced by direct M(t) ≤ Ω.
- **Result:** No. The supt-SDP gives invalid LBs (off by 5× of the open gap).
- **Lesson:** the cell-envelope is necessary for validity, not just a relaxation.
- **Status:** PRO-22 → Done with negative result. The framework ceiling
  `C_∞ ≈ 0.380558` is more robust than expected.

Code: `lp_research_state/code/white_full_convex_supt.py` retained as a
*demonstration of the failure mode*, with explicit validity-check script
provided for future reference. Should NOT be used as a production LB tool.
