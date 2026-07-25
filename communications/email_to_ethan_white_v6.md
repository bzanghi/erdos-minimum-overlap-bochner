# Email to Ethan White — v6 (write-up links + certification update)

**Status: DRAFT. NOT SENT. Ben reviews and sends manually — do not send programmatically.**

**Supersedes:** [email_to_ethan_white_v5_reply.md](email_to_ethan_white_v5_reply.md), which quotes
`µ ≥ 0.380284` and predates the certification round. If v5 was never sent, send this instead
rather than both.

**To:** ethan.white@ubc.ca *(from his reply to Ben, ~2 months ago; verify it is still current)*
**From:** Ben Zanghi <ben@benzanghi.com>
**Subject:** Erdős minimum overlap — write-up, and a correction to something I told you

**Send as:** a reply in the existing thread, so he keeps the context.

---

Dear Ethan,

Two write-ups, in case they're of interest:

- Narrative: https://benzanghi.com/blog/erdos-minimum-overlap
- Derivation: https://benzanghi.com/blog/erdos-minimum-overlap-technical

The bound is now **µ ≥ 0.3803954**, up from the 0.380284 I quoted you. Almost none of that
came from new mathematics.

The part I think may interest you is a correction to my own work. Each centre in my cover was
anchored at the solver's reported objective less a fixed 1e-5 margin. I had been describing
that as rigorous; it isn't, since an interior-point dual is feasible only up to its residual.
I've replaced it with a Jansson–Chaykin–Keil a-posteriori bound in directed-rounding interval
arithmetic, with the multipliers read from the same solve. The certified values came out
*above* the convention at all twelve centres, so the fix raised the bound rather than costing
anything — the margin had been guarding against the wrong thing.

Your Appendix II is what made this uncomfortable. You already do the a-priori feasibility
argument by hand; I had the better tool available and shipped a margin that existed only in a
docstring.

Two things you may not have seen. Kim and Pilanci (ICML 2026, arXiv:2606.31182) certified
0.37912 in June, so the published record is theirs rather than yours now. And on the upper
side, the 0.380856 in circulation is not a bound — the construction it comes from evaluates to
0.3809490, the reported figure being recovered only by dividing by a cell count of
4096.9999999999. The authors caught it themselves; the paper appears unrevised.

I've also withdrawn a saturation ceiling I'd implied to you earlier. The derivation didn't
hold, and a weaker containment argument replaces it.

I lean on your published Table 2 for a second, more conservative statement (µ ≥ 0.380000)
whose only other ingredient is those twelve certificates. If you'd rather I frame that
differently, say so and I'll change it.

Same disclosure as before: I'm a software engineer, not a mathematician, working with an AI
assistant, with everything load-bearing checked in exact or interval arithmetic.

Best wishes,
Ben

Ben Zanghi
ben@benzanghi.com
https://www.linkedin.com/in/bzanghi

---

## Notes for sending

- **~380 words**, matching the length Ben asked for on v4.
- **Verify the address** before sending — it is transcribed from his earlier reply, not
  re-confirmed.
- Every number here is verified by running, 2026-07-25: LB `0.3803953504` (adaptive core
  floor, N=48000, 12/12 anchors certified, `pen_zs = 0` at each); UB
  `0.380859056614806899090596051448` (exact rational); SimpleTES honest value
  `0.3809489501030183`; Kim–Pilanci `0.37912` (fetched from arXiv).
- **Tone choices, deliberate:** the letter leads with a correction to *our* work, not with the
  improved number. He gave us two errata unprompted; reciprocating in kind is the right
  register, and the Appendix II observation is a genuine compliment that happens to be true.
- **Not included:** the five ruled-out architectures (in the companion if he wants them), the
  N-scaling detail, and any suggestion of co-authorship — he hasn't invited it and the v4
  thread already left the door open.
- The Table 2 paragraph is a courtesy check, not a request for permission: Theorem 2 relies on
  his published result exactly as he states and uses it. Worth him seeing before a preprint.
