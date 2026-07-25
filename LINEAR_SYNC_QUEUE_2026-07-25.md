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

**Protocol gap to flag, not paper over.** That repo's CLAUDE.md says *"Claim
before code — run `/claim-ticket` before touching files."* The post rewrite was
branched and committed **without a ticket**, because Linear was unreachable.

Needs: a ZAN issue created, moved to In Progress retroactively, then In Review
when the PR opens.

> Rewrote `src/content/blog/erdos-minimum-overlap.mdx` for July 2026. The May
> version's numbers are all superseded. Branch
> `ben/erdos-minimum-overlap-july-2026` off `main`, `npm run build` passes 35/35,
> not pushed (production deploys only when Ben asks).

The rest of the protocol was followed: branched off `main`, no stacked PRs,
build green.
