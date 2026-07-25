"""
ub_certify — turn a float64 step function into an EXACTLY certified upper bound
on the Erdos minimum overlap constant mu.

Why this is needed
------------------
The published upper bound mu <= 0.380871 (Together AI, March 2026) and the
16-digit anchor 0.3808703105862199 used throughout this repo are both float64
evaluations of a float64 iterate.  Neither the feasibility of that iterate nor
the value of its objective has been checked in exact arithmetic.  They do not
survive it as stated: Together's h*, read as exact dyadic rationals, has

    sum_i h_i  =  300 - 607096245493 / 2^89   !=  300 = n/2,

so int h = 1 fails by about 9.8e-16 and h* is not, strictly, a feasible point
of the problem it is quoted as bounding.  The defect is microscopic and easily
repaired, but until it is repaired the bound is not a theorem.

What this module produces
-------------------------
Given any float64 h, it emits an exactly-feasible rational point

    h_i = a_i / D,   a_i integers in [0, D],   sum_i a_i = n*D/2  (exactly),

and the exact rational value

    M = (2/n) * max_j sum_i a_i (D - a_{i+j}) / D^2,

computed in integer arithmetic over ALL signed lags.  Since such an h is an
admissible competitor, mu <= M is then a rigorous inequality with no
floating-point step anywhere in the chain.

Rounding cost is controlled a priori: |dM_j/dh_k| <= 1 for every cell, so
snapping each cell to the 1/D grid moves the scaled objective by at most 1/D.
With the default D = 2^60 that is below 1e-18 and cannot affect the 16th digit.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction

import numpy as np

import ub_core as U

__all__ = ["snap_to_rational", "exact_overlap", "certify"]


def snap_to_rational(h: np.ndarray, D: int) -> list[int]:
    """Round h onto the 1/D grid, then repair the mass to sum(a) == n*D/2 exactly.

    Returns integer numerators a with 0 <= a_i <= D and sum(a_i) == n*D//2.
    """
    n = h.size
    if (n * D) % 2:
        raise ValueError("n*D must be even for an exact mass of n/2")
    target = n * D // 2

    a = [int(round(float(x) * D)) for x in h]
    a = [min(D, max(0, ai)) for ai in a]

    # Repair the mass to hit `target` exactly.  The input's mass defect can be
    # anywhere from ~1e-16 (a converged float64 iterate) to ~1e-12 (a projected
    # optimiser output), i.e. up to ~1e6 units at D = 2^60, so the repair is
    # distributed proportionally to available headroom rather than unit by unit.
    for _ in range(64):
        diff = target - sum(a)
        if diff == 0:
            break
        room = [(D - ai) if diff > 0 else ai for ai in a]
        total = sum(room)
        if total == 0:
            raise RuntimeError("no headroom to repair mass")
        sgn = 1 if diff > 0 else -1
        need = abs(diff)
        moved = 0
        for i in range(n):
            if room[i] == 0:
                continue
            give = min(room[i], need * room[i] // total)
            a[i] += sgn * give
            moved += give
        need -= moved
        # sweep the small remainder one unit at a time
        i = 0
        while need > 0 and i < 64 * n:
            k = i % n
            if 0 <= a[k] + sgn <= D:
                a[k] += sgn
                need -= 1
            i += 1
    if sum(a) != target:
        raise RuntimeError("could not repair mass")
    assert all(0 <= ai <= D for ai in a)
    return a


def exact_overlap(a: list[int], D: int) -> tuple[Fraction, int]:
    """Exact (M, argmax_lag) for h = a/D, over ALL signed lags.

    M = (2/n) * max_j sum_i a_i (D - a_{i+j}) / D^2, integer arithmetic only.
    """
    n = len(a)
    b = [D - ai for ai in a]  # numerators of 1 - h
    best, best_j = None, None
    for j in range(-(n - 1), n):
        lo, hi = max(0, -j), min(n, n - j)
        s = 0
        for i in range(lo, hi):
            s += a[i] * b[i + j]
        if best is None or s > best:
            best, best_j = s, j
    return Fraction(2 * best, n * D * D), best_j


def decimal_ceil(q: Fraction, digits: int) -> str:
    """Exact decimal string >= q, with `digits` places.  No float anywhere.

    Rounding UP is what preserves the direction of a *upper* bound.
    """
    scaled = -((-q.numerator * 10 ** digits) // q.denominator)  # ceil division
    s = str(scaled).rjust(digits + 1, "0")
    return f"{s[:-digits]}.{s[-digits:]}"


def certify(h: np.ndarray, D: int = 2 ** 60, label: str = "") -> dict:
    n = h.size
    a = snap_to_rational(h, D)
    M, j = exact_overlap(a, D)
    return {
        "label": label,
        "n": n,
        "D_log2": D.bit_length() - 1,
        "M_exact_num": M.numerator,
        "M_exact_den": M.denominator,
        "M_float": float(M),
        "M_certified_decimal": decimal_ceil(M, 30),
        "argmax_lag": j,
        "float_M_of_input": U.overlap_value(h),
        "snap_shift": float(M) - U.overlap_value(h),
        "a": a,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=str, default="",
                    help="JSON with a top-level 'best_h' list; default = Together h*")
    ap.add_argument("--Dlog2", type=int, default=60)
    ap.add_argument("--out", type=str, default="")
    args = ap.parse_args()

    U.selftest(verbose=False)

    if args.input:
        h = np.asarray(json.loads(open(args.input).read())["best_h"], dtype=np.float64)
        label = args.input
    else:
        h = U.load_together()
        label = "together_h_star"

    # exact feasibility of the raw float64 input, for the record
    fr = [Fraction(float(x)) for x in h]
    raw_sum = sum(fr)
    print(f"input: {label}   n={h.size}")
    print(f"  exact sum(h)  = {raw_sum}")
    print(f"  == n/2 ?        {raw_sum == Fraction(h.size, 2)}"
          f"   (defect {float(raw_sum - Fraction(h.size, 2)):+.3e})")

    res = certify(h, D=2 ** args.Dlog2, label=label)
    print(f"  float64 M     = {res['float_M_of_input']:.17f}")
    print(f"  EXACT   M     = {res['M_float']:.17f}   (snap shift {res['snap_shift']:+.3e})")
    print(f"  certified     mu <= {res['M_certified_decimal']}")
    print(f"  exact rational = {res['M_exact_num']} / {res['M_exact_den']}")
    print(f"  argmax lag     = {res['argmax_lag']}")
    if args.out:
        with open(args.out, "w") as f:
            json.dump(res, f)
        print("wrote", args.out)


if __name__ == "__main__":
    main()
