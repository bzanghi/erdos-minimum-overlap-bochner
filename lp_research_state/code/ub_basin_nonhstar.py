"""
ub_basin_nonhstar — close the gap left by ub_basin_sweep.

In the n=600 sweep, the top-64 tier-1 points selected for deep SLP polish were
*all* from the `perturb_hstar` family, because that family dominated the descent
ranking. So the sweep showed that non-h* basins are worse AFTER DESCENT, but
never polished one — leaving open whether a random basin would catch up under
the full polish.

This script regenerates the best tier-1 points from each NON-h*-derived family
(the sweep is seed-deterministic, so they reproduce exactly) and gives each the
same deep SLP polish, making the comparison apples-to-apples.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import ub_core as U
import ub_local as L
from ub_basin_sweep import _tier1

EXCLUDE = {"perturb_hstar"}


def _job(spec):
    n, kind, seed, iters, lr, rounds = spec
    r = _tier1((n, kind, seed, iters, lr, 0.02))
    h = np.asarray(r["h"], dtype=np.float64)
    h2, M, gain = L.slp_polish(h, rounds=rounds, r0=3e-3, r_min=1e-8)
    return {"kind": kind, "seed": seed, "M_descent": r["M"], "M": float(M),
            "gain": float(gain), "h": h2.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", type=str, default="../data/ub_basin_sweep_n600.json")
    ap.add_argument("--per-family", type=int, default=5)
    ap.add_argument("--rounds", type=int, default=150)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    U.selftest(verbose=False)
    d = json.load(open(args.sweep))
    n = d["n"]

    by_kind = collections.defaultdict(list)
    for r in d["tier1_summary"]:
        if r["kind"] not in EXCLUDE:
            by_kind[r["kind"]].append(r)
    specs = []
    for kind, rows in by_kind.items():
        rows.sort(key=lambda r: r["M"])
        for r in rows[: args.per_family]:
            specs.append((n, kind, r["seed"], 6000, 1e-2, args.rounds))
    print(f"polishing {len(specs)} non-h* points from {len(by_kind)} families", flush=True)

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(_job, specs, chunksize=1))
    res.sort(key=lambda r: r["M"])

    best_by = {}
    for r in res:
        if r["kind"] not in best_by or r["M"] < best_by[r["kind"]]["M"]:
            best_by[r["kind"]] = r
    print("\nbest polished value per non-h* family:")
    for k in sorted(best_by, key=lambda k: best_by[k]["M"]):
        print(f"   {k:14s} {best_by[k]['M']:.15f}   (descent {best_by[k]['M_descent']:.9f})")
    print(f"\nbest non-h*      {res[0]['M']:.15f}  ({res[0]['kind']})")
    print(f"best perturb_h*  {d['best_M']:.15f}")
    print(f"Together         {U.ANCHOR:.15f}")
    print(f"gap non-h* vs h* basin: {res[0]['M'] - d['best_M']:+.3e}")

    json.dump({"best_nonhstar_M": res[0]["M"], "best_kind": res[0]["kind"],
               "best_h": res[0]["h"],
               "all": [{k: r[k] for k in ("kind", "seed", "M_descent", "M", "gain")}
                       for r in res]}, open(args.out, "w"))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
