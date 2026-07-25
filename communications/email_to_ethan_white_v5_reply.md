# Reply to Ethan White — v5 (response to his 2026-05-31 message)

> **SUPERSEDED — do not send.** Use `email_to_ethan_white_v6.md`. This draft also carries the
> wrong recipient: Ethan's address is **[personal address - see .local-contacts.md]**, confirmed 2026-07-25, not the
> ubc.ca address transcribed here.


**Status:** DRAFT — for Ben to review and send. Do not send programmatically.
**Context:** White replied to our v3 outreach (PRO-8): he validated the Bochner-PSD constraint, flagged two corrections (5.6/5.7 coeff 8→4; 5.8/5.9 should use 2m−1), and said neither is material. This is the reply to *that* message. (The earlier `email_to_ethan_white_v4.md` was a refined re-pitch written before his reply arrived — now superseded by this reply.)
**Verify** his current email/affiliation before sending.

**To:** (verify) ethan.white@ubc.ca
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** Re: Erdős minimum-overlap — your corrections applied, and a full-space result

---

Dear Ethan,

Thank you — this was generous and genuinely useful, and faster than I had any right to expect.

**On the Bochner constraint and your instinct about interval-averages of f.** You're exactly right that one *could* add variables for the average value of `f` on small intervals (mirroring the `v_j, w_i` for `M`) and derive many implied constraints — and I think your hunch is the heart of it: the moment-matrix PSD condition is the compact way to encode `f, 1−f ∈ [0,1]` *without* that variable/constraint blow-up. Bochner gives `f ≥ 0` (and `1−f ≥ 0`) as a single PSD block on the existing Fourier coefficients `f̂(k)`, rather than O(N) interval-average variables plus their linking inequalities. So it adds the same information more cheaply, which is why it composes so easily with your §5 program.

**Your two corrections — applied, and I can confirm your "not material."**
- *5.6/5.7 (8 → 4 in the RHS numerator):* applied. I had inherited the `8` verbatim; interestingly its real-part sibling already carried the correct `4`, so the asymmetry was a giveaway in hindsight. I made it a parameter and measured both: at the binding center the change is ~`−3×10⁻⁸` (pure solver noise), so your bound and mine are both unaffected — exactly as you said. (One small curiosity: at an *off*-binding center the corrected `4` raises the local objective by ~`+1.6×10⁻³`, i.e. the old `8` was a valid-but-looser constraint there, never an over-tight one — reassuring for rigor.)
- *5.8/5.9 (m → 2m−1 on the RHS):* it turned out my encoding already used `2m−1` there, matching the LHS indexing, so no change was needed.

**The main development since I wrote.** I've been able to promote the bound from your core residual region (5.16) to your **entire** `(E(M), c₁, d₁)` parameter space — i.e. across all of the Table-2 outside regions, not just Table 3. Concretely, each augmented dual point gives a globally-valid lower bound `Φ_c(h,p,q)`, and the cover `max_c Φ_c` clears the core value over every Table-2 region (the wide regions need adaptive subdivision to control the Lipschitz grid-error, and a couple of deep-`d₁` corners are simply infeasible). The upshot is that

  **µ ≥ 0.380284 now holds unconditionally over the whole space**, with your core region (5.16) remaining the binding one —

so the improvement no longer rests on the published per-region "0.38" floors for the outside regions; it's certified internally. (The core headline under my tighter tail-bound convention is `0.3802973`; the full-space certificate uses a slightly more conservative anchor, hence `0.380284`.)

**Two honest caveats I'd state in any writeup**, and where you could help most: the full-space bound is *load-bearing* on the polynomial-moment cuts and on a handful of fresh dual centers in two narrow `d₁`-strips (the core anchors alone fall ~`3×10⁻⁵` short there), and the binding margin is thin (~`10⁻⁴`). I'm currently re-solving those centers at higher `N` to widen the margin.

**One concrete ask, if it's easy for you.** My only external dependency in the full-space argument was reading your Table-2 entries' "0.38" as literally `≥ 0.380000`. If you still have the **exact per-region optimum values** behind that table (especially the regions adjacent to the core — your `E(M) ≤ 0.06`, `c₁ ∈ [0.33,0.45]`, `|d₁| ∈ [0.02,0.025]` strips), they'd let me cross-check my internal certificate against your originals and remove even that reading as an assumption. No worries at all if they're long gone.

I'd still love your read on whether the cell-envelope `sup_t (h⋆h)(t)` representation is the "right" one or just convenient — a direct sup-`t` SDP I tried produced an *invalid* lower bound, which suggests the cell-envelope is load-bearing for validity, not just convenience. And I remain very interested in writing this up with you involved in the framing if you'd like; otherwise I'll cite your work as the foundation it is.

Thank you again — this genuinely moved the work forward.

Best wishes,
Ben

---

*Send notes:* verify email/affiliation; this is a thread reply (keep "Re:"); if he offers the Table-2 values, fold them into `FULLSPACE_VERIFICATION.md` as an independent cross-check and it removes the last literal-floor assumption. Keep the caveats in — White is exactly the reader who would (rightly) catch an overclaim, and stating them is a credibility asset.
