"""
ub_fine_search — try to BEAT the best known upper bound on mu by re-optimising
on grids that are not commensurate with Together's n = 600.

Motivation.  PRO-34 showed h* is first-order stationary under cell DOUBLING
(600 -> 1200), which preserves M exactly.  But a 600-cell step function is not
representable on an 800- or 1000-cell grid, so re-projecting h* there and
re-optimising explores genuinely different function space.  The Einstein Arena
leaderboard (which superseded Together's 0.3808703 with 0.3808591) has public
discussion referencing n = 800 constructions, which is what suggested this.

Seeds per grid: the cell-average re-projection of h*, perturbations of it at
several amplitudes, and the diverse random families from ub_local.

Usage:
  python ub_fine_search.py --grids 800,1000,1200 --per-grid 40 --out fine.json
"""
from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import ub_core as U
import ub_local as L

BEST_KNOWN = 0.3808591          # Einstein Arena leaderboard, fetched 2026-07-25
TOGETHER = U.ANCHOR             # 0.3808703105862199


def _job(spec):
    n, kind, seed, amp, iters, lr, rounds = spec
    rng = np.random.default_rng(seed)
    hstar = U.load_together()
    if kind == "resample_hstar":
        h0 = U.resample(hstar, n)
    elif kind == "resample_perturb":
        h0 = U.project_feasible(U.resample(hstar, n) + rng.normal(0, amp, n))
    elif kind == "resample_smooth":
        r = U.resample(hstar, n)
        k = int(rng.integers(1, 6))
        h0 = U.project_feasible(np.convolve(r, np.ones(k) / k, mode="same"))
    else:
        h0 = L.random_start(n, kind, rng)
    h, M1 = L.smooth_descent(h0, iters=iters, lr=lr)
    h, M, gain = L.slp_polish(h, rounds=rounds, r0=2e-3, r_min=1e-9)
    return {"n": n, "kind": kind, "seed": int(seed), "amp": float(amp),
            "M_descent": float(M1), "M": float(M), "gain": float(gain),
            "h": h.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grids", type=str, default="800,1000,1200")
    ap.add_argument("--per-grid", type=int, default=40)
    ap.add_argument("--iters", type=int, default=8000)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--rounds", type=int, default=200)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    U.selftest(verbose=False)

    kinds = ["resample_hstar", "resample_perturb", "resample_perturb",
             "resample_perturb", "resample_smooth", "fractional", "symmetric",
             "blocks", "smooth"]
    amps = [0.0, 0.01, 0.03, 0.08, 0.0, 0, 0, 0, 0]

    specs = []
    for n in [int(g) for g in args.grids.split(",")]:
        for i in range(args.per_grid):
            k = kinds[i % len(kinds)]
            a = amps[i % len(kinds)] or 0.02
            specs.append((n, k, 50_000 + 977 * n + i, a, args.iters, args.lr, args.rounds))

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(_job, specs, chunksize=1))
    res.sort(key=lambda r: r["M"])

    print(f"{len(res)} runs in {time.time()-t0:.0f}s")
    by_n = {}
    for r in res:
        by_n.setdefault(r["n"], []).append(r)
    for n in sorted(by_n):
        b = min(by_n[n], key=lambda r: r["M"])
        print(f"  n={n:5d}  best={b['M']:.13f}  ({b['kind']})  gain={b['gain']:.2e}")
    best = res[0]
    print(f"\nBEST OVERALL   {best['M']:.15f}   n={best['n']} kind={best['kind']}")
    print(f"Together       {TOGETHER:.15f}   delta {best['M']-TOGETHER:+.3e}")
    print(f"Leaderboard    {BEST_KNOWN:.15f}   delta {best['M']-BEST_KNOWN:+.3e}")
    print(f"beats Together: {best['M'] < TOGETHER}   beats leaderboard: {best['M'] < BEST_KNOWN}")

    with open(args.out, "w") as f:
        json.dump({
            "best_M": best["M"], "best_n": best["n"], "best_kind": best["kind"],
            "best_h": best["h"], "together": TOGETHER, "leaderboard": BEST_KNOWN,
            "beats_together": bool(best["M"] < TOGETHER),
            "beats_leaderboard": bool(best["M"] < BEST_KNOWN),
            "all": [{k: r[k] for k in ("n", "kind", "seed", "M_descent", "M", "gain")}
                    for r in res],
        }, f)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
