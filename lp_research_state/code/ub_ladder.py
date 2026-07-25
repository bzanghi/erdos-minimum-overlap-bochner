"""
ub_ladder — climb the cell-doubling ladder from a coarse witness.

`M_n` is non-increasing under cell doubling (`ub_core.cell_double` preserves the
objective exactly), so 512 -> 1024 -> 2048 -> ... is free, and each level has
strictly more freedom than the last.  Polish at every rung with the
line-searching trust-region SLP and keep the best.

This is worth doing on a COARSE witness and not on a fine one: measured here,
SimpleTES's n=2400 point is already first-order stationary to ~2e-9 and doubling
it to n=4800 unlocked only ~1e-9/round.  A 512-cell witness has far more slack.

    ../../.venv/bin/python ub_ladder.py --input ../data/arena_lnzwz_n512.json \\
        --levels 3 --rounds 40 --out ../data/ub_ladder.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ub_core import cell_double, overlap_value, overlap_profile
from ub_refine import slp_polish_ls, load_witness


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--levels", type=int, default=3)
    ap.add_argument("--rounds", type=int, default=40)
    ap.add_argument("--r0", type=float, default=1e-5)
    ap.add_argument("--slack", type=float, default=1e-5)
    ap.add_argument("--level-budget", type=float, default=900.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    h, _meta = load_witness(args.input)
    M = overlap_value(h)
    print(f"start   n={h.size:6d}  M={M:.16f}", flush=True)
    rungs = [{"n": int(h.size), "M": M, "stage": "input"}]
    best_h, best_M = h.copy(), M

    # polish at the input resolution first, then climb
    for level in range(args.levels + 1):
        if level > 0:
            h = cell_double(h)
            Md = overlap_value(h)
            assert abs(Md - M) < 1e-14, "cell doubling changed the objective"
            print(f"\ndouble  n={h.size:6d}  M={Md:.16f}  (exact, drift "
                  f"{Md - M:+.1e})", flush=True)

        prof = overlap_profile(h) * (2.0 / h.size)
        deg = int((prof >= prof.max() - 1e-6).sum())
        print(f"polish  n={h.size:6d}  degenerate lags within 1e-6: {deg}/{prof.size}",
              flush=True)
        t0 = time.time()
        h, M, gain, hist = slp_polish_ls(
            h, rounds=args.rounds, r0=args.r0, slack=args.slack,
            verbose=True, time_budget_s=args.level_budget)
        print(f"        n={h.size:6d}  M={M:.16f}  gain={gain:.3e}  "
              f"({time.time() - t0:.0f}s, {len(hist)} rounds)", flush=True)
        rungs.append({"n": int(h.size), "M": M, "stage": f"level{level}",
                      "gain": float(gain), "rounds": len(hist),
                      "degenerate_lags": deg})
        if M < best_M:
            best_M, best_h = M, h.copy()

    print(f"\nBEST  n={best_h.size}  M={best_M:.16f}  "
          f"vs input {rungs[0]['M']:.16f}  ({best_M - rungs[0]['M']:+.3e})", flush=True)

    if args.out:
        Path(args.out).write_text(json.dumps({
            "best_h": best_h.tolist(), "M": best_M, "n": int(best_h.size),
            "source": f"ub_ladder from {Path(args.input).name}",
            "rungs": rungs,
        }, default=float))
        print(f"saved -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
