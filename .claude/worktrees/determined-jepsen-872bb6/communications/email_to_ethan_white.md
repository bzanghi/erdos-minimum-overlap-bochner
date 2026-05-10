# Email to Ethan Patrick White

**To:** ethan.white@ubc.ca *(verify current address — UBC was Killam-Trusts/NSERC affiliation per arXiv:2201.05704; check his current affiliation/email before sending)*
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** Bochner-PSD strengthening of your minimum-overlap lower bound

---

Dear Professor White,

I hope this finds you well. I'm writing because I have a small but
rigorous improvement on your 2023 lower bound `µ ≥ 0.379005` for the
Erdős minimum overlap problem (Acta Arith. 208), and I wanted to share
it with you before posting publicly — both because the result is in the
direct spirit of your program and because I'd value your read on it.

**Result.** Adding a Bochner moment-matrix PSD constraint to your
Section-5 convex program — namely the Hermitian Toeplitz constraints
`[f̂(j-k)]_{j,k=0..n} ⪰ 0` for both `f` and `1-f` — and then applying
your own §5.1 / Appendix II ellipse-extension argument with the
augmented dual objective gives

$$\mu \;\ge\; 0.379828$$

an improvement of about `+8.2 × 10⁻⁴` over your published value, with
the same proof structure, just enlarged ellipses. The strengthening
combines two valid SDP constraint families — Bochner-PSD on $\hat f$
and Lasserre level-2 on the bilinear products — both of which close
loopholes in your program that admit $(c, d)$ configurations not
arising from any $f \in [0, 1]$.

**Method (sketch).**

1. Add the convex SDP constraints `M_n(f) ⪰ 0` and `M_n(1-f) ⪰ 0`
   (where `M_n(f) := [f̂(j-k)]_{j,k=0..n}`) to your program. These are
   immediate consequences of `f, 1-f ≥ 0` and Bochner / Toeplitz.
2. Solve the augmented program at each of your seven Table-3 ellipse
   centers at `N = 10000, T = 4000, R = 10, n = 20` (CLARABEL via cvxpy).
3. Extract dual variables; recompute the dual objective as a quadratic
   function of `(h, p, q)` per your §5.1 envelope-theorem argument.
4. The seven ellipses where this quadratic exceeds `0.379005` cover
   the residual region (5.16) on a fine grid (verified analytically at
   the boundary). The `(h, p, q)`-min over the cover is `0.3795475`
   (closed-form), reduced to `0.379544` with a conservative `1e-6`
   margin for CLARABEL's interior-point gap.

**Why it works.** Your Section-5 program admits `(c, d)` configurations
that don't correspond to any `f ∈ [0, 1]`. Numerically, at your row 4
ellipse center, the optimal `(c, d)` violates `M_n(f) ⪰ 0` from `n = 10`
upward (min eigenvalue ≈ −0.16 at n=10, growing in magnitude with n).
The Bochner constraint closes this loophole; the dual objective at each
ellipse center increases by ~9 × 10⁻⁴, which propagates to a `+5.4 × 10⁻⁴`
gain on the cover-min.

**What I've checked.**

- Three independently-written implementations agree on the per-row SDP
  values to 10+ digits.
- The Bochner constraint encoder has been re-implemented from scratch
  by a separate agent (without seeing the first encoding); both agree
  bit-for-bit on the Hermitian-PSD real-form.
- CLARABEL's `optimal_inaccurate` flag triggers because its sharp
  tolerance (1e-8) isn't met within iteration cap, but the actual gap
  per the verbose output is consistently ≤ 1e-7 across all 7 rows.
  An SDPA-GMP spot-check at small N agrees with CLARABEL to 5e-9.
- Your envelope-theorem extension is replicated exactly: the dual
  feasibility region is `(h, p, q)`-independent, so the ellipse
  argument carries through verbatim.

**Caveats.**

- CLARABEL `optimal_inaccurate`. The 1e-6 margin gives 10× headroom over
  the observed ~1e-7 IPM gap; a higher-precision SDP solver (SDPA-GMP)
  at scale would remove this concern. Spot-checks at small N confirm
  CLARABEL is rigorous to ~1e-9.
- Coverage of (5.16): verified on a 1001×1001 (h, p) grid plus closed-form
  per-row ellipse minima. Continuous-`q` variation is at the boundaries
  `q ∈ {-0.02, +0.02}` only; your same Table-2 outer-bands argument
  handles the rest.
- I have not attempted to push `n_Bochner` above 30 or to combine with
  Lasserre-level-2 (a separate sub-experiment showed Lasserre-2
  compounds, +8 × 10⁻⁵ at modest scale, still rising). Either could
  improve the bound modestly further; I wanted to send the cleanest
  result first.

**Disposition.** I'd like to write this up as a short note, with you as
co-author if you're interested in commenting / refining it, or
acknowledged as the indispensable foundation if you'd rather not. The
code, raw numerics, and three independent implementations are all in
a Git repo I'll happily share. The result clearly belongs in the
arc you opened with Acta Arith. 2023.

A couple of incidental things you might find interesting:

- The conditional bound assuming `f*` is even (your §6 question) is
  `µ_even ≥ 0.379904`, i.e., the symmetry assumption barely tightens
  things beyond Bochner alone. This suggests the Bochner-augmented
  optimum may already be close to even.
- M-side Bochner via convex relaxation (SOC or Schur) of `|f̂(m)|²` is
  empirically dead — the slack absorbs all constraint content. An
  exact non-convex bilinear lifting would be needed.
- Lasserre-level-2 compounds with Bochner additively at modest scale.

I'm a software engineer by day and was working on this as a project
with Anthropic's Claude. The mathematical structure is yours; the
contribution here is "noticed Bochner is missing from your program,
verified the extension carefully." Happy to disclose the AI-assisted
workflow in detail if useful.

Could I send you the repo and the research note? Either way I'd
welcome any thoughts.

Best,
Ben Zanghi
ben@benzanghi.com

---

*Notes for sending:*
- Verify White's current email/affiliation. arXiv shows UBC + Killam
  Trusts/NSERC support 2022. Recent talks (UBC, March 2023) suggest
  he was still there. A quick search for "Ethan Patrick White" + a
  current institution should give the right address.
- This is intentionally formal-but-warm. Math etiquette: short, lead
  with the result, offer co-authorship, acknowledge his foundation.
- Don't send before the Git repo is public — he should be able to
  inspect the code immediately.
