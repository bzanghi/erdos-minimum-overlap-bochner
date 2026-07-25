# The upper-bound side is closed to local methods — measured, 2026-07-25

Everything here was verified by running. The conclusion is negative and I think
it is worth having in writing, because it retires a whole class of attempts.

**Claim.** At laptop scale, no local method reaches below the current record
`µ ≤ 0.3808590568145606` (Einstein Arena, `lnzwz_AI4M_Agent`, n = 512) by more
than `2 × 10⁻¹⁰`. Beating it needs a genuinely different construction, which in
practice means the large-scale evolutionary search the AI-search groups run.

---

## 1 · What was tried, and what it returned

| method | scale | result | verdict |
|---|---|---|---|
| Cell doubling + SLP, SimpleTES | n = 2400 → 4800 | 1.2e-9 / round, 60 s LPs | stationary |
| Cell doubling + SLP, arena | n = 512 → 1024 | 8.6e-12 / round | stationary |
| Cell doubling + SLP, arena | n = 1024 → 2048 | no gain, 94 s LPs | stationary |
| Non-commensurate resample + polish | n = 512 → 600…1200 | best 0.3808818 | **worse than record** |
| Symmetry projection (3 subgroups) | n = 512 | best 0.3814097 | **catastrophic** |
| Structured basin-hop, 6000 proposals | n = 512 | best 0.3808590566148070 | same basin |
| Random multistart, 1600 starts (prior session) | n = 600 | all ≥ 0.380924 | far worse |

The basin-hop is the strongest of these and deserves its numbers spelled out.
6000 proposals, 7 move families (`gauss`, `block`, `swap`, `reflect`, `mirror`,
`lowfreq`, `spike`), perturbation scale swept log-uniformly over `[1e-3, 2e-1]`,
each screened with 12 SLP rounds and the best 20 deep-polished with 150. 245
proposals beat the starting point — every one of them by between `1e-11` and
`2e-10`, i.e. all of them are polish artifacts *inside* the record basin. The
best, `0.3808590566148070`, arrives via a `reflect` move at scale `1.6e-3` and
improves the arena witness by `1.998 × 10⁻¹⁰`.

Nothing escaped the basin. Not once in 6000 tries.

## 2 · Why cell doubling cannot work, in hindsight

`M_n` is non-increasing under cell doubling and `ub_core.cell_double` preserves
the objective to 0 ulp, so refining a witness is free and strictly enlarges the
feasible set. That makes it *sound* like a lever. It is not, because both record
witnesses are already first-order stationary in the enlarged space: the doubled
point inherits stationarity, and the new degrees of freedom (splitting each cell
into two) are not descent directions. Measured across three doublings and two
independent witnesses.

The corollary is the useful part: **record progress on this problem has never
come from resolution.** The arena's 512-cell construction beats SimpleTES's
2400-cell one. Neither is anywhere near its own `M_n`. The binding resource is
search quality, not grid size.

## 3 · The symmetry group, and that the optimum breaks it

Verified numerically on the record witness — all three act on the feasible set
and preserve `M` exactly (`Δ = 0`, every digit):

| map | effect on the objective |
|---|---|
| `R : h(x) → h(2−x)` | reverses the lag axis, `sup` unchanged |
| `K : h → 1 − h` | reverses the lag axis, `sup` unchanged |
| `RK : h → 1 − h(2−x)` | fixes each lag pointwise |

`RK` is checked by substituting `u = 2 − x − k` in the overlap integral. So a
Klein four-group acts on the optimizer set.

The optimum is not a fixed point of any of it. `K`'s fixed subspace is the single
function `h ≡ 1/2`, which scores `0.5`. And although the record witness sits at
correlation `0.975` with its own reflection, projecting it onto the `R`-symmetric
subspace costs `5.5 × 10⁻⁴` even after 120 polish rounds. The small antisymmetric
component carries disproportionate weight — which is why symmetry-restricted
search (halving the dimension, an attractive idea) is a trap here.

## 4 · What would actually be needed

The record constructions descend from a shared lineage of refinement
(Haugland → AlphaEvolve → TTT-Discover → Together → SimpleTES → the arena
entries), not from independent discovery. Every local method tried here returns
to whichever basin it started in. So the open question on this side is not "can
this point be improved" — it cannot, materially — but "how is a better basin
found", and the only demonstrated answer is large-scale evolutionary program
search with a lot of compute.

Two things that *are* worth doing from here, neither of which is search:

1. **Exact certification as a service.** Three of three witnesses examined fail
   exact feasibility (mass off by 1e-16 to 4e-14), and one published record was
   a normalization artifact worth `9.3e-5`. `ub_certify.py` settles a 2400-cell
   construction in about a second.
2. **Watch the leaderboard.** It moves weekly, the witnesses are downloadable
   through `/api/solutions/best?problem_id=1`, and certifying whatever tops it
   is a one-second job.

## 5 · Code

`ub_search.py` (structured basin-hopping, multiprocessing, checkpointing),
`ub_ladder.py` (cell-doubling ladder), `ub_refine.py` (trust-region SLP with an
exact line search along the LP step — the stock `ub_local.slp_polish` discards
the entire LP when the full step fails, which near a degenerate optimum is most
of them).

```bash
cd lp_research_state/code && ../../.venv/bin/python ub_search.py \
    --input ../data/arena_lnzwz_n512.json --iters 6000 --workers 6 \
    --out ../data/ub_search_best.json
```
