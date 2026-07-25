# Erdős minimum overlap: state after the certification round, 2026-07-25

Supersedes `MINIMUM_OVERLAP_STATE_2026-07-25.md` (same day, earlier session).
Everything below is labelled **[RAN]** (verified by executing code or fetching a
primary source in this session), **[DOC]** (read in a repo document, not
re-verified), or **[INF]** (inferred).

---

## 1 · Headline

| | Value | Status |
|---|---|---|
| **UB, certified here** | **µ ≤ 0.380859056614806899090596051448** | [RAN] exact rational, arena witness + polish |
| UB, widely cited | ~~0.380856~~ | [RAN] **NOT A BOUND** — normalization artifact |
| **LB, this repo** | **µ ≥ 0.3802946** | [RAN] 12/12 core anchors Jansson-certified |
| LB, previous repo headline | µ ≥ 0.3802838 | [RAN] reproduced end-to-end |
| LB, best published | µ ≥ 0.37912 | [RAN] Kim & Pilanci, ICML 2026 |

**Bracket: `µ ∈ [0.3802946, 0.3808590567]`.**

The UB improves the previous best certified value (0.380867545960922320593700936552,
SimpleTES n=2400 + SLP) by **8.49 × 10⁻⁶**.

---

## 2 · The cited upper-bound record is not a bound

`0.380856` (SimpleTES / *Evaluation-driven Scaling for Scientific Discovery*,
[arXiv:2604.19341](https://arxiv.org/abs/2604.19341) §3.4.1, 21 Apr 2026) is
quoted in this repo's own docs and across the literature as the best known upper
bound. It is not one. **[RAN]**

The witness exists — as the *initial-commit* version of the SimpleTES repo JSON,
later overwritten. Downloaded from commit `406fc651` and re-evaluated:

```
stored score      0.38085596768904106      (n = 4096, mass exactly 2048, h in [0,1])
honest objective  0.3809489501030183       (dx = 2/4096)
discrepancy       -9.298241397726059e-05
```

The stored value is reproduced **bit-for-bit** by dividing by `4096.9999999999`
rather than `4096`, and the evolved program in the same commit contains

```python
epsilon   = 0.9999999999      # as close to 1 as possible, still < 1
n_points  = float(n) + epsilon
```

with `dx = 2.0 / n_points`, while the feasibility check uses `int(n_points)`.
Verified numerically: `max(corr)*2/4096` gives 0.3809489501030183,
`*2/4097` gives 0.38085596768903174, `*2/4096.9999999999` gives exactly the
stored 0.38085596768904106.

The authors caught it. Commit `6eb2ca0a` (2026-05-23), verbatim: *"fix a
potential hack possibility with n_points not being integer"*, replacing the
artifact with an honest `0.3808676758` construction. **The arXiv paper has never
been revised (still v1 as of 2026-07-25)**, so the number keeps propagating.

The true value of that construction, 0.3809490, is worse than every entry in the
Einstein Arena top 14.

### The real record, certified

Best UB with a public, downloadable, verified witness: **Einstein Arena** entry
`lnzwz_AI4M_Agent`, submission 2407, 512 cells. **[RAN]**

```bash
curl -s "https://einsteinarena.com/api/solutions/best?problem_id=1&agent_name=lnzwz_AI4M_Agent&limit=1"
```

Re-evaluated here: `M = 0.38085905681456067` vs reported `0.3808590568145606`
(diff 5.6e-17); mass exactly 256.0 = n/2; `h ∈ [0,1]`. Exactly certified by
`ub_certify.py`:

```
µ ≤ 32399905329033719063044589842123849041 / 85070591730234615865843651857942052864
  = 0.380859056814560651295303328196          (decimal rounded UP)
argmax lag = 100
```

Polishing that witness improves it further. A cell-doubling ladder to n=1024
gains `1.63 × 10⁻¹⁰`; a 6000-proposal structured basin-hop at n=512 does slightly
better at `2.00 × 10⁻¹⁰`, giving the **best certified upper bound available**:

```
µ ≤ 129599621248162196650621937272674201309 / 340282366920938463463374607431768211456
  = 0.380859056614806899090596051448          (decimal rounded UP)
n = 512, argmax lag = 69
```

Both polishes stay inside the record basin — see
[UB_SEARCH_NEGATIVE_2026-07-25.md](UB_SEARCH_NEGATIVE_2026-07-25.md).

This improves the previous best certified UB (0.380867545960922320593700936552,
SimpleTES n=2400 + SLP polish) by **8.49 × 10⁻⁶**. Files:
`data/ub_certified_arena_n512.json`, `data/ub_certified_best.json`.

**Third witness in a row that fails exact feasibility.** Together's `h*` was
short of mass by 9.8e-16, SimpleTES's over by 4.1e-14, the arena's short by
6.4e-16. All three were only ever checked in float64. `ub_certify.py` repairs to
an exactly feasible rational and evaluates all `2n−1` signed lags in integer
arithmetic; it settles any construction in about a second.

---

## 3 · The lower bound: what was wrong with the certificate

### 3.1 Two defects in `dual_extractor.py` **[RAN]**

`rigorous_dual_LB` was the project's central epistemic claim, quoted in
`CLAUDE.md`, `REPRODUCE.md` and `findings.md`. Two independent defects:

1. **Wrong column.** CLARABEL's iteration table is
   `iter pcost dcost gap pres dres k/t μ step` — captured verbatim from a live
   solve this session. The old regex took five groups and named the fifth
   `dual_residual`. The fifth column is **`pres`, the primal residual**. The
   eligibility gate that decided which iterations were trustworthy was reading
   the wrong quantity.
2. **No feasibility margin.** The docstring conceded that a strict bound needs
   zero dual residual and that the residual "can be absorbed into a margin".
   None ever was. The returned number is the raw dual objective of an
   approximately-feasible point — not a lower bound on anything.

Both fixed: the regex now captures `pres` and `dres` separately, the gate uses
`dres`, and the return dict carries `is_certificate: False` with a pointer to
the real tool. The name `rigorous_dual_LB` is kept for caller compatibility but
is documented as a misnomer.

### 3.2 The anchor convention was not a theorem **[RAN]**

The core (5.16) floor is `min` over the (h,p) box of `max_c Φ_c`, with
`Φ_c(θ) = anchor_c + shift_c(θ)` and `anchor_c := primal_c − 1e-5`. That haircut
on a solver-reported value was the load-bearing quantity. Before this round,
**2 of the 12** core anchors had an actual certificate (`docs/RND_WHITESPACE/L2_PROD.json`).

### 3.3 The fix, which also tightened the bound **[RAN]**

`_jansson_reanchor.py` re-solves each core center and returns the
Jansson-Chaykin-Keil a-posteriori bound `p_lo` in directed-rounding interval
arithmetic **together with the model duals read from the same solve**, so the
pair `(p_lo, λ)` is internally consistent:

```
p_lo = -b(θ_c)ᵀz + pen_zs − pen_Dx        (pen_zs ≤ 0, pen_Dx ≥ 0)
Φ_c(θ) := p_lo + shift_c(θ)  ≤  dual_obj(z; θ)  ≤  µ | θ
```

Validity across the box needs both penalties to be θ-independent. `pen_Dx`
involves `D = c + Aᵀz` and fixed model box bounds — no `b`, so it is. `pen_zs`
is θ-independent exactly when it **vanishes**, i.e. `z` lies in the dual cone;
the driver asserts `penalty_zs == 0` at every center and refuses to emit one
otherwise. It held at all 12.

All twelve anchors certified, and **every certified value came out above the
convention it replaced**:

| center | old anchor (`primal − 1e-5`) | Jansson `p_lo` | gain |
|---|---:|---:|---:|
| row1 | 0.380698115 | 0.381394387947 | +6.96e-4 |
| row2 | 0.380677211 | 0.381464078052 | +7.87e-4 |
| row3 | 0.381042228 | 0.382132253805 | +1.09e-3 |
| row4 | 0.380338114 | 0.380379812943 | +4.17e-5 |
| row5 | 0.380611839 | 0.380618330942 | +6.49e-6 |
| row6 | 0.380436448 | 0.380445257236 | +8.81e-6 |
| row7 | 0.382003546 | 0.384594385921 | +2.59e-3 |
| cde_n30_iter1 | 0.380346079 | 0.380354430229 | +8.35e-6 |
| cde_n30_iter2 | 0.380377646 | 0.380408959276 | +3.13e-5 |
| cde_n30_iter3 | 0.380302283 | 0.380310043153 | +7.76e-6 |
| cde_n30_iter4 | 0.380606798 | 0.381207476917 | +6.01e-4 |
| cde_n30_iter5 | 0.380393675 | 0.380663716920 | +2.70e-4 |

`row4` and `cde_n30_iter3` reproduce the June run in `L2_PROD.json` to 7 and 2
significant figures respectively, from an independent solve.

So the honest anchor is the stronger anchor. Core floor:
`0.3802837846529683` → **`0.3802946016`**.

Two conventions, both rigorous, differing only in how the Lipschitz cell error is
paid. A single 4001×4001 grid over the whole core box gives
`0.3802897673` (`grid_min` 0.3802946921 − `eps_grid` 4.92e-6). The certified duals
are steeper than the old ones (`L_max` 0.1491 → 0.3378), so a single grid pays
more than twice the eps it used to, eating half the gain. Adaptive subdivision
recomputes `L_max` per sub-box and drives eps to 4.1e-8, recovering it:

| method | floor | eps | binding point |
|---|---:|---:|---|
| single 4001² grid | 0.3802897673 | 4.92e-6 | (0.00285, 0.39225) |
| adaptive, depth 12 | 0.3802942806 | 3.68e-7 | (0.00281–0.00375, 0.39219–0.39375) |
| adaptive, depth 16 | 0.3802945525 | 9.20e-8 | — |
| **adaptive, depth 20** | **0.3802946016** | 4.12e-8 | (0.002812, 0.392188) |

Binding witness `cde_n30_iter3` (it was `row4` before re-anchoring).

### 3.4 The two-tier statement

**Tier 1 — shortest chain (what I would defend hardest): µ ≥ 0.380000.**
White's Table 2 column is headed *"Optimum lower bound"* and he uses the printed
values as rigorous bounds: *"The data of the first line of Table 2 shows that
either µ ≥ 0.38 or E(M\*) ≤ 0.75"* [RAN, read from arXiv:2201.05704 v1]. Take
his published table for regions 1–17; our 12 certified anchors cover the
residual region (5.16) **and** lift his weakest strip (region 18, his 0.37925)
to 0.380894 with witness `row6` [RAN]. Nothing else is load-bearing.
This is **+8.8 × 10⁻⁴** over the best published LB.

**Tier 2 — working frontier: µ ≥ 0.3802946.** Replaces White's Table 2
with our own 121-center full-space cover. The ~109 non-core centers are still
anchored at the uncertified `primal − 1e-5` convention, so this is a working
result, not a theorem. Re-anchoring them is mechanical (~109 × 4 min).

### 3.5 Reproduced end-to-end, twice **[RAN]**

The previous headline was re-derived from scratch this session, which the
earlier session explicitly had *not* done — and the driver regenerated
`fullspace_promote_final.json` **byte-identically** (`git status` reports it
unmodified):

```
_fs_recompute.py  ->  INDEPENDENTLY-CERTIFIED floor 0.3802838, binding region = core
                      regions still white-reliant: []   (all 18 clear 0.380000 on our own)
_eval_r6_box.py   ->  R6 floor 0.380309 (ceiling 0.380780), reproduced exactly
```

Re-run against the certified centers (`LP_DUALEXT=...dualext_reanchored.json`):

```
_fs_recompute.py  ->  INDEPENDENTLY-CERTIFIED floor 0.3802898, binding region = core
                      binding corner (h, p, q) = (0.00285, 0.39225, 0.0)
                      regions still white-reliant: []
```

so the core stays binding, and the full-space floor is the core floor. Tightest
outside regions, all comfortably clear of it: R6 0.3803090, R17 0.3803351,
R16 0.3803547, R9 0.3803667, R7 0.3805539. The headline `0.3802946` is the
adaptive core value; `0.3802898` is the same bound under the simpler
single-grid convention `_fs_recompute` uses.

**Gate-region floors are target-limited, not true infima.** The adaptive
subdivision in `_eval_r*_box.py` stops as soon as a sub-box clears its target,
so a reported region floor means "at least this", not "this is the infimum". R6
reports 0.380309 against a ceiling of 0.380780. Raising `$LP_TARGET` makes them
subdivide further and report higher, still-rigorous floors — which is why a
raised core floor is not capped by the stored region numbers.

---

## 4 · Negative results, measured

### 4.1 Grid refinement does not improve the record constructions **[RAN]**

`M_n` is non-increasing under cell doubling, so refining a witness is free and
strictly enlarges the search space. It does not help — both current records are
already first-order stationary:

| witness | doubling | descent/round | LP/round | verdict |
|---|---|---:|---:|---|
| SimpleTES n=2400 | → 4800 | 1.2e-9 | 60 s | stationary |
| Arena n=512 | → 1024 | 8.6e-12 | 0.1 s | stationary |
| Arena n=1024 | → 2048 | 8.8e-12 | 94 s | stationary |

To move the sixth decimal at those rates takes thousands of LP solves. **Local
refinement is not the lever**; record progress comes from better search, not
finer grids. New code: `ub_refine.py` (trust-region SLP with an exact line
search along the LP step — the stock `ub_local.slp_polish` discards the whole LP
when the full step fails), `ub_ladder.py` (the doubling ladder).

### 4.2 The architecture ceiling still binds

The cover + Bochner-PSD approach cannot certify past `C_explicit = 0.380713`,
which is `1.5 × 10⁻⁴` below the certified UB. **[DOC]** — not re-derived this
session. Closing the gap needs an LB architecture that does not factor through
these duals.

---

## 5 · New / changed code

| file | what |
|---|---|
| `_jansson_reanchor.py` | **new** — certify + re-anchor core centers; `--all`, `--new-center`, `--emit-dualext`, `--evaluate` |
| `_jansson_all12.py` | **new** — certify the 12 anchors against the old convention |
| `ub_refine.py` | **new** — line-searching trust-region SLP + grid refinement |
| `ub_ladder.py` | **new** — cell-doubling ladder |
| `_finalize_reanchor.sh` | **new** — one-command pipeline from certificates to both tiers |
| `dual_extractor.py` | **fixed** — column mis-parse; honest naming; `is_certificate: False` |
| `_fullspace_eval.py` | `load_centers(path)` / `$LP_DUALEXT`; `anchor_value(., 'p_lo')`; `$LP_TARGET` |

Data: `parallel_results/jansson_core12_reanchored.json`,
`phase5_N20K_bn40_dualext_reanchored.json`,
`data/arena_lnzwz_n512.json`, `data/ub_certified_arena_n512.json`.

Reproduce:

```bash
cd lp_research_state/code && ../../.venv/bin/python _jansson_reanchor.py --all
```

```bash
cd lp_research_state/code && bash _finalize_reanchor.sh
```

---

## 6 · What the next rung takes

1. **Publish.** The LB leads the published record by 1.16e-3 and a second group
   (Kim & Pilanci) is working the same direction to the same standard.
2. **Re-anchor the ~109 outside-region centers** the same way. Mechanical, ~7 h
   of compute, converts Tier 2 from a working frontier into a theorem.
3. **Correct the record publicly.** `0.380856` is in circulation and is not a
   bound; the arXiv paper is unrevised.
4. **Offer exact certification to the UB community.** Three of three witnesses
   examined here fail exact feasibility; no search system checks its own
   arithmetic outside float64.
5. **Do not** run more Bochner/Lasserre level scans, more N-scaling, more random
   multistart at n=600, or more grid refinement of record witnesses. All
   measured dead ends.
