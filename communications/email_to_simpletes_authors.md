# Draft — note to the SimpleTES authors (arXiv:2604.19341)

**STATUS: DRAFT. NOT SENT. Do not send without Ben's explicit approval.**

**To:** corresponding author(s) of *Evaluation-driven Scaling for Scientific Discovery*
(arXiv:2604.19341); repo owner of github.com/wq-will/SimpleTES
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** Erdős minimum overlap: the 0.380856 ablation figure

---

Hello,

I work on lower bounds for the Erdős minimum overlap constant, and I maintain
exact-arithmetic certification code for upper-bound constructions. While
re-checking the current record I ran into something in §3.4.1 that I think you
would want flagged, though I suspect you already know about it.

The paper reports 0.380856 from a Best-solution Restart run, excluded from
Table 10 for fairness but quoted in the abstract and Figure 1 as the best
result. I fetched the construction from the initial commit of the repository
(`406fc651`, `best_results/mathematics_extremal_analysis/erdos_minimum_overlap/`)
and re-evaluated it directly:

```
stored score        0.38085596768904106
recomputed value    0.3809489501030183     (dx = 2/4096)
difference          9.298e-05
```

The array is genuinely feasible — 4096 cells, mass exactly 2048, values in
[0, 1] — so this is purely a normalization question. The stored value is
reproduced bit-for-bit by dividing by 4096.9999999999 rather than 4096, which
matches the evolved program in that same commit:

```python
epsilon  = 0.9999999999    # as close to 1 as possible, still < 1
n_points = float(n) + epsilon
```

with `dx = 2.0 / n_points` in the objective and `int(n_points)` in the shape
check. So the reported figure is about 2.4e-4 relative below the construction's
true value, and 0.380856 is not an upper bound on µ. The honest value of that
construction, 0.3809490, is above Together's 0.3808703.

You caught this yourselves — commit `6eb2ca0a` (2026-05-23) says "fix a
potential hack possibility with n_points not being integer" and replaces the
artifact with a 2400-cell construction scoring 0.3808676758. I re-evaluated
that one too and it reproduces bit-for-bit; it is a real result and it was the
best available witness until recently.

The reason I am writing is that the arXiv paper is still at v1, so 0.380856 is
propagating. I found it quoted as the record in my own project's notes, and
it appears in downstream summaries. A short v2 correcting §3.4.1, the abstract
and Figure 1 would stop that.

Two things that may be useful to you either way:

1. I have exact-arithmetic certification code for this problem — it snaps a
   float construction to an exactly feasible rational point and evaluates all
   2n−1 signed lags in integer arithmetic, returning a certified decimal. It
   runs in about a second on a 2400-cell witness. Every construction I have
   examined (Together's, yours, and the current Einstein Arena leader) is
   exactly infeasible by between 1e-16 and 4e-14 in mass — harmless, but it
   means none of the published values were ever bounds as literally stated. I
   am happy to share the code or to certify constructions on request.

2. The current best witness I am aware of is the Einstein Arena entry from
   `lnzwz_AI4M_Agent` (512 cells, 0.3808590568). I have certified it, and a
   short polish improves it to 0.3808590566512541 exactly.

For context on the other side of the gap: I have a lower bound of about 0.3803
built on White's convex program with Bochner moment constraints, certified in
interval arithmetic, which I am preparing for publication.

Happy to be wrong about any of the above — please tell me if I have
misread the normalization.

Best,
Ben Zanghi
ben@benzanghi.com
https://www.linkedin.com/in/bzanghi

---

## Verification trail for this note (all run locally, 2026-07-25)

```bash
curl -s "https://raw.githubusercontent.com/wq-will/SimpleTES/406fc65137a4d56cd53b39a0509c4af302060cba/best_results/mathematics_extremal_analysis/erdos_minimum_overlap/erdos_minimum_overlap_best_construction.json" -o simpletes_v0.json
```

```
stored 0.38085596768904106 | honest 0.3809489501030183
max(corr)*2/4096              = 0.3809489501030183
max(corr)*2/4097              = 0.38085596768903174
max(corr)*2/4096.9999999999   = 0.38085596768904106   <- exact match
```

Tone note: state the mechanism, not a motive. The authors fixed it themselves,
which is the strongest evidence it was a search artifact rather than intent.
