"""PRO-35: A/B test of White's (5.6)/(5.7) erratum (author email, 2026-07).

E. White (private communication, July 2026): in the published program,
constraints (5.6)/(5.7) [labels fullsinep/fullsinem, the sine cell-average
sandwich] have an 8 in the RHS numerator that should be a 4 (matching
B_m = -(4/(m pi)) sin(m pi/2) b_m, eq. (Bm) in the paper); and in
(5.8)/(5.9) the m's on the RHS should be 2m-1. Our transcription
(white_full_convex.py) already used 2m-1 in the tail bounds (valid, looser
direction) but inherited the 8 at the sine sandwich (line ~188).

Measured at row4 (binding center), N=3000, T=1200, R=10, bochner_n=20,
CLARABEL:
    factor-8 (paper/our code): 0.378792326841  (reproduces recorded baseline)
    factor-4 (corrected):      0.378794084451
    delta = +1.758e-06
The corrected program is TIGHTER (higher optimum), so previously reported
lower bounds remain valid a fortiori; a full corrected re-run should
microscopically improve the headline. Full 12-center Phase-5 re-run with the
4 is the required follow-up before adopting any new headline number.

Method: two copies of white_full_convex.py differing only in the sine-sandwich
factor (sed 8.0 -> 4.0), solved back-to-back with identical params/solver.
"""
print(__doc__)
