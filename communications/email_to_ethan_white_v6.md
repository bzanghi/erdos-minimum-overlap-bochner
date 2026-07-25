# Email to Ethan White — v6 (write-up links + certification correction)

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

The bound is now µ ≥ 0.3803954, up from the 0.380284 I quoted you — though the interesting
part is a correction to my own work rather than the number.

Each centre in my cover was anchored at the solver's reported objective less a fixed 1e-5
margin, and I had been calling that rigorous. It isn't: an interior-point dual is feasible only
up to its residual. I've replaced it with a Jansson–Chaykin–Keil bound in interval arithmetic.
The certified values came out *above* the old convention everywhere, so the fix raised the
bound rather than costing anything.

Your Appendix II is what made that uncomfortable to discover. You do the a-priori
floating-point argument by hand, in a paper from 2022; I had better tools and shipped a margin
that existed only in a docstring. The more time I spend inside your Section 5, the more the
care in it shows — the errata you sent me were both things I'd have had no way to find from
outside, and the covering argument is doing far more work than it appears to on first reading.
Everything here is a bolt-on to your framework.

I should say plainly: I'm an AI engineer, not a mathematician, and this is a personal curiosity
project. I'm conscious that's a different thing from your actual research programme, and that
correspondence with an enthusiastic amateur may not be a good use of your time. Genuinely no
offence taken if you'd rather I stop — just say the word.

With thanks,
Ben

Ben Zanghi
ben@benzanghi.com
https://www.linkedin.com/in/bzanghi

---

## Notes for sending

- **~290 words** (was ~380 in the previous version).
- **Verify the address** before sending — transcribed from his earlier reply, not re-confirmed.
- Numbers verified by running, 2026-07-25: LB `0.3803953504` (adaptive core floor, N=48000,
  12/12 anchors certified, `pen_zs = 0` at each); prior figure quoted to him was `0.380284`.
- **Tone choices, deliberate:**
  - Leads with a correction to *our* work, not the improved number. He gave two errata
    unprompted; reciprocating in that register is the right response.
  - The Appendix II line and the "care in it shows" paragraph are genuine and specific, not
    flattery — his hand-derived feasibility margin is exactly the discipline our code lacked.
  - The closing gives him a graceful exit. He is a working mathematician and owes an amateur
    correspondent nothing; making that explicit costs nothing and removes any awkwardness in
    him not replying.
- **Cut from the previous draft** to make room: Kim & Pilanci holding the published record, the
  `0.380856` upper-bound correction, the withdrawn saturation ceiling, and the Table 2 courtesy
  note. All four are in the technical companion if he follows the link. If Ben wants any of
  them back, the `0.380856` correction is the most useful to him.
- **Not included:** co-authorship. He hasn't invited it, and the v4 thread already left that
  door open.
