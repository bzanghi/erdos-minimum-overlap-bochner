# Linear sync queue — 2026-07-25 certification round

**Not synced.** The `plugin:productivity:linear` MCP server needs OAuth, which
cannot run in a non-interactive session, and there is no authenticated Linear
CLI or token on this machine. Everything below is pending manual entry.

Commit: `c218b15` (Erdős repo) · `d810d24` (Portfolio, branch
`ben/erdos-minimum-overlap-july-2026`, unpushed).

---

## Team PRO (Erdős)

### Update existing

**PRO-38** — full-space promotion, µ ≥ 0.3802838. *Keep Done; add comment:*

> Re-derived end-to-end for the first time this session — `_fs_recompute.py`
> regenerates `fullspace_promote_final.json` byte-identically, confirming
> 0.3802838 with binding region = core and no region White-reliant. Now
> superseded by the Jansson-re-anchored floor 0.3802946 (binding still core,
> tightest outside region R6 at 0.3803090).

**PRO-5** — preprint. *Add comment:*

> Headline moved twice. LB 0.3802838 → **0.3802946** (all 12 core anchors now
> carry a Jansson-Chaykin-Keil interval-arithmetic certificate, with duals read
> from the same solve). Best certified UB is
> **0.380859056651254094841565810818**, not 0.380871. The preprint also needs a
> two-tier split: Tier 1 `µ ≥ 0.380000` resting only on our 12 certificates plus
> White's published Table 2 (his column is headed "Optimum lower bound" and he
> writes "either µ ≥ 0.38 or E(M\*) ≤ 0.75"), and Tier 2 `µ ≥ 0.3802946` using
> our own 121-center cover, whose ~109 non-core anchors are still the
> uncertified `primal − 1e-5` convention.

**PRO-34** — UB refinement / "0.380871 a serious candidate for µ". *Keep closed;
add comment:*

> Further falsified — certified UB is now 0.3808590567. Separately, a measured
> negative result: **grid refinement of record constructions is dead.** Cell
> doubling is free (M non-increasing), but both current records are already
> first-order stationary — SimpleTES 2400→4800 descends 1.2e-9/round with 60 s
> LPs; the Einstein Arena witness 512→1024 gains 8.6e-12/round and 1024→2048
> nothing. Against a 1e-5 target that is thousands of solves. New code:
> `ub_refine.py` (SLP with an exact line search along the LP step),
> `ub_ladder.py`.

### Create

**New — Todo, High.** *Re-anchor the ~109 outside-region centers with Jansson
certificates.*

> The core (5.16) anchors are now certified, but the stage-2 / halo / promote
> centers that certify White's Table-2 regions still use the uncertified
> `primal − 1e-5` convention. Re-running them through `_jansson_reanchor.py`
> (~109 × 4 min ≈ 7 h, serialized, ~8 GB each) converts the Tier-2 working
> frontier into a theorem. On this round's evidence it will also raise it: every
> certified core anchor came out above the convention it replaced, by +6.5e-6 to
> +2.6e-3. Highest-value remaining item.

**New — Todo.** *Correct the public record: 0.380856 is not a bound.*

> SimpleTES arXiv:2604.19341 §3.4.1's ablation value is a normalization
> artifact — the evolved program reports its bin count as
> `float(n) + 0.9999999999` and divides by that. The witness's honest value is
> 0.3809489501030183, worse than Together's. The authors fixed it themselves in
> commit `6eb2ca0a`; the arXiv paper is still v1, so the number keeps
> propagating. Draft note sits unsent at
> `communications/email_to_simpletes_authors.md` — **needs Ben's approval before
> sending.**

**New — Done.** *`dual_extractor.rigorous_dual_LB` was never a certificate.*

> Two independent defects, both fixed: the eligibility gate parsed CLARABEL's
> 5th iteration column and called it `dual_residual`, but the table is
> `iter pcost dcost gap pres dres k/t μ step` — the 5th column is `pres`, the
> primal residual. And no dual-infeasibility margin was ever subtracted despite
> the docstring saying one must be. The function now returns
> `is_certificate: False` and points at `_jansson_verify.jansson_lower_bound`.

---

## Team ZAN (benzanghi.com)

**SHIPPED — create the issue retroactively and set it straight to Done.**

Ben approved the deploy explicitly. [PR #26](https://github.com/bzanghi/Portfolio/pull/26)
squash-merged to `main` as `4ae8112`; CI green (Typecheck · Lint · Build, 1m43s);
Vercel deployed. Live and content-verified at
<https://www.benzanghi.com/blog/erdos-minimum-overlap>.

> Rewrote `src/content/blog/erdos-minimum-overlap.mdx` for July 2026 — every
> number in the May version was superseded. Headline LB 0.3802838 → 0.3802946
> (all 12 core anchors certified in interval arithmetic), best certified UB
> 0.380859056651254094841565810818, and the widely-cited 0.380856 refuted as a
> normalization artifact. A follow-up commit credits Ethan White's framework
> (Acta Arithmetica 208, 2023) properly and makes clear that none of the
> certificate failures in the post are his — his Appendix II already contains
> the a-priori floating-point feasibility argument the repo's own code skipped.

**Protocol gap to flag, not paper over.** That repo's CLAUDE.md says *"Claim
before code — run `/claim-ticket` before touching files."* No ticket was claimed,
because Linear was unreachable from the authoring session. This is noted in the
PR body as well. The rest of the protocol was followed: branched off `main`, no
stacked PRs, CI green, squash-merge, branch deleted.

Ben's unrelated `ben/blog-personal-tech-radar` branch was preserved intact
(`19761de`) — that work landed separately as PR #25.

---

# Addendum — architecture-audit round (same day, later)

Commits `42f3373`, `50cdcb7`, `97fe8ad`. Still unsynced, same reason.

### Correct, do not merely close

**PRO-6 / Lever I′** — *`C_explicit = 0.380713` is not sound as derived.*

> `docs/archive/LEVER_I_PRIME_FINAL.md` §2's own table gives row7
> `Ω = 0.381586` plus residual `5.85e-4` = `0.382171`, which is **above the
> certified upper bound** — so the per-row statement cannot be a valid ceiling as
> written. The published number also substitutes the min-over-rows Phase-5 anchor
> `0.3801279` for the theorem's per-row RHS. Replace with a feasible-set
> containment argument: the discretized program at any `N` is at most the
> continuum program's value, pinned two independent ways at **≈0.38065 ± 2e-5**
> at the witness's parameters. Weaker, but true, and still enough to conclude the
> architecture cannot close the gap. **This number has been quoted repeatedly in
> repo docs and in my own summaries; it needs correcting wherever it appears.**

### Create

**New — Todo, High.** *N-scaling is the only lever with measured headroom.*

> The graveyard entry "N sweeps don't help" was measured in the 10k–40k band
> only. At the certified cover's binding point `(0.002812, 0.392188, ±0.02)` with
> production `(T=4000, bochner_n=40, pm_k_max=20)`: `N=20000 → 0.380343028`,
> **`N=48000 → 0.380458890`** — `+1.16e-4` in 378 s at 5.7 GB peak. That is ~4×
> the entire remaining f-side headroom, for zero new mathematics. A 12-centre
> re-certification at N=48000 was launched at session end →
> `parallel_results/jansson_core12_N48000.json`. **Open question that decides
> this:** does the Jansson interval pass still certify at 2.4× the cell count, or
> does interval widening eat the gain? If it eats it, the correct diagnosis is
> that the certified LB is **Jansson-limited, not envelope-limited** — which
> nobody has established either way and is more useful than either ceiling story.

**New — Done.** *The f-side of White's architecture is closed.*

> Two independent measurements agree. An exact moment-body program (replacing
> `c,d = Variable(T)` by the exact moment map of a density) rigorously
> upper-bounds production **plus any f-side cut family, present or future**,
> because production has no M-side augmentation and all its f-side constraints
> are implied by a genuine density. Headroom at the binding anchor
> `cde_n30_iter3`: **+3.77e-5**, and it *shrinks* with N (+1.2e-4 at N=6000).
> Bochner-40 + poly_moment has already taken ~90% of all f-side slack. Action:
> commit the meter as `lp_research_state/code/fside_ceiling.py` and gate every
> future f-side proposal on running it first — 19 s versus minutes for a bn=40
> SDP, because it drops both 82×82 PSD blocks. Note it is a *diagnostic*, not a
> bound on µ (piecewise-constant `y` inner-approximates the moment body).

**New — Done.** *Five LB architectures killed, one terminal survivor.*

> Writeup: `LB_ARCHITECTURE_AUDIT_2026-07-25.md`. Dead: discretization-error
> `µ ≥ M_n − ε_n` (the aliased detail spectrum is a positive-definite sequence,
> which *forces* the projection error to be maximally adverse at lag 0 — so
> positivity cannot sign the error; sharpest lattice-only bound `≤ 0` for every
> n, needed `ε_n < 6.7e-4`, provable `≈0.36`); the extreme-point / PD-weight
> reduction (correct but capped at **0.3294738** by an exact rational certificate
> over 51 measure-1 adversaries, exactly 1/4 for positive-definite weights);
> argmax-disjunctive branching (one-line vacuity theorem: `P ⊆ ⋃D_j` implies the
> branched min equals the unbranched inf); bathtub cell-envelope (correct, but a
> Richardson extrapolator for N in disguise); bathtub Markov/Krein cuts (correct,
> capped at ≤3.8e-5). Survivor: half-integer period-4 Bochner — valid, novel,
> absent from the repo, worth ~3.8e-5. Terminal, not a lever.

**New — Done.** *UB side closed to local methods.*

> 6000 structured basin-hop proposals at n=512 (7 move families, scale over
> `[1e-3, 2e-1]`) never left the record basin; 245 "beats" were all 1e-11 to
> 2e-10 polish artifacts inside it. Non-commensurate resampling loses more than
> polish recovers. Symmetry projection costs 5.5e-4 — and the symmetry group is
> now verified: `R: h(2−x)`, `K: 1−h`, `RK: 1−h(2−x)` each preserve M exactly,
> and the optimum is a fixed point of none of them. Best certified UB is now
> **`µ ≤ 0.380859056614806899090596051448`** (n=512, argmax lag 69). Writeup:
> `UB_SEARCH_NEGATIVE_2026-07-25.md`.

**New — Todo, Low.** *Bank the certified continuous n-cell optima.*

> Corrected epigraph-Lasserre level 2 certifies global optimality of the
> *continuous* n-cell minimax, tight to 1e-8: `M_4 = 0.4000000`,
> `M_6 = 0.3888889`, `M_8 = 0.3850717`, `M_10 = 0.3824271`, `M_12 = 0.3822141`.
> First such certificates in the repo; certified **upper** bounds on µ. Worth
> committing as a table with the generating script.

### Site

The live post quotes the UB as `0.3808590567`, which is still correct to the
digits shown, but the exact certified rational improved after publication to
`0.380859056614806899090596051448`. Not worth a redeploy on its own; fold into
the next site update.

---

# ADDENDUM — everything after `c218b15`

Commits: `42f3373` (UB local search closed) · `50cdcb7` (certified UB + README)
· `97fe8ad` (LB architecture audit) · `a23d7c9` (N=48000 re-certification).
Portfolio: merged and deployed as `4ae8112` — https://www.benzanghi.com/blog/erdos-minimum-overlap

## Numbers that changed again

| | was (earlier today) | now |
|---|---|---|
| certified LB, core | 0.3802946 | **0.3803954** (N=48000) |
| certified LB, full-space | 0.3802946 | **0.3803351** (binding = R17, stale anchors) |
| certified UB | 0.380859056651254094841565810818 | **0.380859056614806899090596051448** |

## New issues to create

**PRO-49 — Re-run gate regions R17/R16/R9/R7 at the N=48000 anchors.** *Todo.*
The only remaining purely-mechanical LB work. Use
`LP_DUALEXT=parallel_results/dualext_reanchored_N48000.json` and a raised
`LP_TARGET`, with `max_depth` above 8 — R6 improved 0.3803090 → 0.380344 and
stopped on depth, not mathematics (its region ceiling is 0.380814). Pure
evaluation, no SDP, ~1–2 h. Expected headline ≈ 0.38040.

**PRO-50 — Bank the rest of the N lever.** *Todo, needs hardware.*
1/N fit at the binding point (production T, bn): `value(N) = 0.380541649 −
3.9724/N`, fitted on N=48000/96000 and reproducing N=20000 to 2e-6. +8.28e-5
remains beyond N=48000; ceiling ≈0.380542 at the binding point. Memory is the
constraint, not time: 5.7 GB at N=48000, roughly linear ⇒ ~23 GB at N=192000,
~60 GB at N=500000. One large-memory cloud instance exhausts this in a day.

**PRO-51 — Commit `fside_ceiling.py` (the f-side slack meter).** *Todo.*
An exact moment-body program that rigorously upper-bounds production **plus any
f-side cut family whatsoever**, in one ~19 s solve at N=20000. Verdict:
Bochner-40 + poly-moment has already harvested ~90% of all f-side slack;
headroom at the binding anchor is +3.77e-5 and *shrinks* with N. No f-side
proposal should be funded again without running this first.

## Issues to close as dead, each with a proof

**Discretization-error LB (µ ≥ M_n − ε_n).** *Close.* The aliased detail
spectrum is a positive-definite sequence, which forces the projection error to
be maximally adverse at lag 0. Sharpest lattice-only bound is `G_n ≤ 0` for
every n, equality at `c ≡ 1/2`. Needed ε_n < 6.7e-4 at n=256; provable ε_n ≈
0.36. The M_n side independently saturates near 0.310 and *decreases* in n, with
Lasserre level 2 bit-identical to level 1.

**Single-weight lag-averaging / extreme-point reduction.** *Close.* Correct and
vacuous: certified ceiling **0.3294738** for any single probability weight
(signed included), exactly **1/4** for positive-definite ones, by an exact
rational certificate over 51 measure-1 adversaries. Collapsing `sup_k` into one
`∫w` costs ≥0.05 unconditionally; the extreme-point structure gained is free
information.

**Argmax-disjunctive branching.** *Close.* One-line vacuity theorem: if
`P ⊆ ⋃ D_j` then `min_j inf{f : P ∩ D_j} = inf{f : P}` exactly.

**Bathtub cell-envelope / bathtub cuts.** *Close.* Both mathematically correct.
The envelope is a Richardson extrapolator for N in disguise (plain N-scaling
gets the same for ~1/10 the compute); the cuts are capped at +3.8e-5 by the
f-side meter.

**UB local search.** *Close.* 6000 structured basin-hops at n=512 across 7 move
families never left the record basin; best gain 2.0e-10. Cell doubling
(3 levels, 2 witnesses), non-commensurate resampling, and symmetry projection
all measured dead. Record progress needs large-scale evolutionary search.

## Correction to log against Lever I′

**`C_explicit = 0.380713` is not sound as derived.**
`docs/archive/LEVER_I_PRIME_FINAL.md` §2's own table gives row7 0.381586 +
5.85e-4 = 0.382171, above the certified UB, and the published figure substitutes
the min-over-rows Phase-5 anchor for the theorem's per-row RHS. Replace with the
feasible-set containment statement: the discretised program at any N is at most
the continuum program's value, pinned two independent ways at **≈0.38065** at
the witness's parameters. Weaker, true, and still enough to conclude the
architecture cannot close the gap.

## Byproduct worth an issue

Certified **global** optima of the continuous n-cell minimax via corrected
epigraph-Lasserre level 2, tight to 1e-8: M_4 = 0.4000000, M_6 = 0.3888889,
M_8 = 0.3850717, M_10 = 0.3824271, M_12 = 0.3822141. First global-optimality
certificates in this repo. Pushing to n ≈ 30–40 would measure how far the record
constructions sit from optimal at their own n — the most valuable open unknown
on the UB side.
