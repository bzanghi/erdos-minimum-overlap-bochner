# Six attempts at a new lower-bound architecture: five dead, one terminal

2026-07-25. Six independent lenses were pointed at the question "find a lower-bound
architecture that does not factor through the saturated cell-envelope + Bochner
duals", then each proposal was handed to an independent adversary instructed to
kill it. Five died. The survivor is valid, new, and worth about `4 × 10⁻⁵`.

The useful output is not the survivor. It is that **the f-side of White's
architecture is now closed by measurement**, and that the repo's own published
ceiling was not sound as derived.

---

## 1 · The headline correction: `C_explicit = 0.380713` is not a theorem

`docs/archive/LEVER_I_PRIME_FINAL.md` is the source of the claim, repeated
throughout this repo and in my own write-ups today, that the architecture
"provably cannot certify past 0.380713". Under audit it does not stand up:

- its own §2 table gives row7 `Ω = 0.381586` plus residual `5.85e-4` = `0.382171`,
  which is **above the certified upper bound** — so the per-row statement cannot
  be a valid ceiling as written;
- the published number substitutes the min-over-rows Phase-5 anchor `0.3801279`
  for the theorem's per-row right-hand side, which is not the quantity the
  derivation bounds.

**Replacement, which needs no sensitivity linearisation and no per-row sup.** By
feasible-set containment, the value of the discretised program at any `N` is at
most the value of the continuum program `P_∞` at fixed `(T, R, bochner_n, pm)`.
That limit was pinned two independent ways at the record witness's own parameter
point: Richardson on an N-ladder gives `0.380653`; the complementary-slackness
table with the entire cosine-envelope slack removed gives `0.380658`. Call it

> **`≈ 0.38065 ± 2e-5` at the witness's parameters**, versus the certified UB
> `0.3808591` — so the architecture is short by `≈ 2.1 × 10⁻⁴` even at the single
> parameter point where the extremal function lives, before any cover loss.

That is a weaker and more honest statement than the one it replaces, and it is
still enough to conclude the architecture cannot close the gap.

## 2 · The f-side is closed

Two independent routes converged on the same number, which is why I believe it.

**Route A — an exact f-side slack meter.** Replace `c, d = Variable(T)` by the
exact moment map of a density `y ∈ [0,1]` on `[-1,1]` with unit mass. Production
has *no* M-side augmentation, and every f-side constraint it does have (`|c| ≤ 2/π`,
`Σ(c²+d²) ≤ 1/2`, T5', Bochner-40, poly-moment) is implied by `y` being a genuine
density. So this program's value is a rigorous **upper bound on production plus
any f-side cut family whatsoever** — present or future. Measured at N=20000:

| center | certified `p_lo` | exact-f ceiling | headroom |
|---|---:|---:|---:|
| `cde_n30_iter3` (binding) | 0.38031004 | 0.38034774 | **+3.77e-5** |
| `cde_n30_iter1` | 0.38035443 | 0.38039237 | +3.79e-5 |
| `row4` | 0.38037981 | 0.38045974 | +7.99e-5 |
| `cde_n30_iter2` | 0.38040896 | 0.38044698 | +3.80e-5 |
| `row6` | 0.38044526 | 0.38048169 | +3.64e-5 |

The headroom also **shrinks** with the only lever that raises the bound:
`+1.2e-4` at N=6000 → `+3.8e-5` at N=20000. So Bochner-40 + poly-moment has
already harvested roughly 90% of all available f-side slack.

**Route B — the surviving proposal, costed.** A half-integer (period-4) Bochner
constraint: Toeplitz PSD on the interleaved *odd* Fourier modes, together with
the exact indicator complement. Validity verified at 50 digits on five admissible
indicator sets; verifiably absent from the repo (the even-sublattice submatrix
identity is exact to `4e-17`); a strict addition of a valid constraint, hence
monotone. Its author measured `+1.73e-4`; the adversary found that number wrong
by `4.6×` and put the true value at `≈ 3.8e-5` — **the same figure Route A
predicts as the cap on any f-side addition**.

Conclusion, and it is a real one: **no further f-side proposal should be funded
without first running the meter.** Recommend committing it as
`lp_research_state/code/fside_ceiling.py`; a single 19 s solve at N=20000 bounds
the entire family, versus minutes for a bn=40 SDP, because it drops both 82×82
PSD blocks. Note it is *not* itself a bound on µ — piecewise-constant `y` is an
inner approximation of the moment body — it is a diagnostic.

## 3 · The five deaths, and why each one is instructive

**Discretization-error LB** (`µ ≥ M_n − ε_n`), the direction the brief ranked most
promising. Dead on both halves.
*The ε side* is settled by an exact identity: `C_{P_n h}(md) = C_h(md) + A_e(md)`
for every lattice lag, and `m ↦ A_e(md)` is a **positive-definite sequence**, so
its maximum is at `m = 0` and equals the full detail energy `‖h − P_n h‖²`. This
answers the question I posed to the lens ("can positivity of `|ĥ|²` *sign* the
error?") in the negative — positive-definiteness *forces* the error to be
maximally adverse. The sharpest lattice-only bound is `G_n ≤ 0` for every `n`,
with equality at `c ≡ 1/2`. Needed `ε_n < 6.7e-4` at n=256; provable `ε_n ≈ 0.36`.
*The `M_n` side*: moment relaxations of the n-cell QCQP saturate near `0.310` and
**decrease** in `n`, with Lasserre level 2 returning bit-identical values to
level 1 — a min-max randomization gap, which does not close with degree.

**Extreme-point / positive-definite-weight reduction** — the structural fact I
derived and handed to the lens. It is *correct*: `∫w C_h` is concave in `h` when
`ŵ ≥ 0`, so the infimum sits at an extreme point, and the extreme points are
indicators. It is also **vacuous**, with an exact rational certificate over 51
explicit measure-1 adversaries: the ceiling is `0.3294738` for *any* single
probability weight (signed ones included, being dominated by their normalised
positive part), and exactly `1/4` for positive-definite ones. The reason is
crisp — collapsing `sup_k` into a single `∫w` costs at least `0.05`
unconditionally, while the extreme-point structure you gain is free information
(indicators are weak-* dense and the functional is weak-* continuous, so the
infimum over the whole set already equals the infimum over indicators). You pay
`0.05` and receive nothing.

**Argmax-disjunctive branching.** Dead to a one-line vacuity theorem: if
`P ⊆ ⋃ D_j` then `min_j inf{f : P ∩ D_j} = inf{f : P}` exactly. The branches were
indexed by the argmax, which always exists, so the union is everything and the
gain is identically zero at every `N`, level, solver and tolerance.

**Bathtub (perspective) cell-envelope.** Mathematically correct — the sharp
mass-and-cap bound `Ω·ψ(w_j/Ω)` strictly dominates White's cell-min line, with
`t = 0` reproducing it exactly, and no truncation anywhere. Dead anyway: it is
contained in the same continuum program `P_∞`, so it is **a Richardson
extrapolator for the N-discretisation dressed as a new constraint**. Plain
N-scaling reaches the same place for about a tenth of the compute and without
adding ~40N epigraph variables that every one of which would need an interval
bound in the Jansson verifier.

**Bathtub cuts (Markov/Krein moment-body separation).** Also mathematically
clean — valid, no tail-bound trap, genuinely outside the Bochner cone at every
finite level. Dead on magnitude: capped by Route A's meter at `+3.8e-5`, and the
author's measured gains were taken against deliberately weakened stacks
(`bochner_n = 8` and `16`), where a zero-cut `bn=40` run already beats them by
`+2.1e-4`.

## 4 · What is actually left

The f-side is closed and the cover/ellipse machinery is understood. That leaves
exactly one direction with measured headroom: **`N`**.

The repo's graveyard entry "N sweeps don't help" was measured in the 10k–40k band.
Outside it, the 1/N tail is still worth roughly `1e-4` and nobody has banked it —
solve time is only about linear in `N`, and memory is dominated by `T` and
`bochner_n`, not `N`. Measured today at the certified cover's own binding point
`(h, p, q) = (0.002812, 0.392188, ±0.02)`:

| N | T | bochner_n | value | secs |
|---:|---:|---:|---:|---:|
| 20000 | 4000 | 40 | 0.380343028 | 175 |
| 20000 | 600 | 12 | 0.380105962 | 16 |
| 48000 | 600 | 12 | 0.380217878 | 46 |

**The single next experiment**, and it needs none of the dead proposals' code:
solve the base program at `N = 96000` with production `(T=4000, bochner_n=40,
pm_k_max=20)` at the binding point, then run `_jansson_verify.jansson_lower_bound`
on that exact problem. Two things to read off:

- **(a)** does the `~1e-4` N-tail survive at production `(T, bn)` and at the
  box-min parameters rather than the witness's?
- **(b)** does the Jansson interval pass still certify at 5× the cell count, or
  does the interval widening eat the gain?

If both hold, re-run `_jansson_reanchor.py` over all 12 centres and the headline
moves about `1e-4` with **zero new mathematics**. If (b) fails, the correct
conclusion is that the certified LB is **Jansson-limited, not envelope-limited** —
a more useful diagnosis than either the bathtub or Lever I′, and one nobody has
made.

## 5 · Byproducts worth keeping

- `exact_fside.py` — the f-side slack meter (§2). Commit as `fside_ceiling.py`.
- Certified global optima of the *continuous* n-cell minimax via a corrected
  epigraph-Lasserre level 2, tight to `1e-8` at n = 4, 6, 8, 10 and near-tight at
  n = 12: `M_4 = 0.4000000`, `M_6 = 0.3888889`, `M_8 = 0.3850717`,
  `M_10 = 0.3824271`, `M_12 = 0.3822141`. These are certified **upper** bounds on
  µ, and they certify global optimality of the local descent optimiser at those
  sizes — the first such certificates in this repo.
- A graveyard entry, stated at full strength: *single-weight lag-averaging bounds
  — including the concavity/extreme-point reduction, every Fejér and every
  positive-definite weight, and signed weights — are capped at `0.3294738`, with
  `1/4` attained exactly by the triangle weight.*
