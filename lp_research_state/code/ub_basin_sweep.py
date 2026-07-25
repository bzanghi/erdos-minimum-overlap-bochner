"""
ub_basin_sweep — is Together's mu <= 0.3808703105862199 basin-unique?

This is the experiment ranked #2 in PRO34_UB_REFINEMENT.md's follow-on list
("Basin-diversity search on the UB side"), which had never been run.  PRO-33/34
established that h* is a numerically exact KKT point and a strict second-order
local minimum under 2x grid refinement, so any UB progress must come from a
DIFFERENT basin, not from local refinement.  This sweep tests for one.

Two tiers:
  tier 1  wide   — many diverse starts, smoothed-minimax descent only (cheap)
  tier 2  deep   — the best tier-1 survivors, polished to first-order
                   stationarity with the trust-region SLP

Every reported M is the exact float64 objective from ub_core.overlap_value.

Usage:
  python ub_basin_sweep.py --n 600 --starts 800 --deep 48 --out sweep.json
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


def _tier1(job):
    n, kind, seed, iters, lr, anchor_perturb = job
    rng = np.random.default_rng(seed)
    if kind == "perturb_hstar":
        h0 = U.load_together()
        if n != 600:
            h0 = U.cell_double(h0, n // 600) if n % 600 == 0 else h0
        h0 = U.project_feasible(h0 + rng.normal(0, anchor_perturb, h0.size))
    elif kind == "shuffle_hstar":
        h0 = U.load_together()
        h0 = U.project_feasible(rng.permutation(h0))
    else:
        h0 = L.random_start(n, kind, rng)
    h, M = L.smooth_descent(h0, iters=iters, lr=lr)
    return {"kind": kind, "seed": int(seed), "M": float(M), "h": h.tolist()}


def _tier2(job):
    h, rounds, r0, r_min = job
    h = np.asarray(h, dtype=np.float64)
    h2, M2, gain = L.slp_polish(h, rounds=rounds, r0=r0, r_min=r_min)
    return {"M": float(M2), "gain": float(gain), "h": h2.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--starts", type=int, default=800)
    ap.add_argument("--deep", type=int, default=48)
    ap.add_argument("--iters", type=int, default=6000)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--rounds", type=int, default=120)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    U.selftest(verbose=False)  # never sweep on an unverified objective

    kinds = list(L.START_KINDS) + ["perturb_hstar", "shuffle_hstar"]
    jobs = []
    for i in range(args.starts):
        kind = kinds[i % len(kinds)]
        amp = [0.02, 0.05, 0.15, 0.4][(i // len(kinds)) % 4]
        jobs.append((args.n, kind, 10_000 + i, args.iters, args.lr, amp))

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        tier1 = list(ex.map(_tier1, jobs, chunksize=1))
    t1 = time.time()
    tier1.sort(key=lambda r: r["M"])
    print(f"tier1: {len(tier1)} starts in {t1-t0:.0f}s   best={tier1[0]['M']:.12f} "
          f"({tier1[0]['kind']})   median={tier1[len(tier1)//2]['M']:.6f}", flush=True)

    deep_jobs = [(r["h"], args.rounds, 3e-3, 1e-8) for r in tier1[: args.deep]]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        tier2 = list(ex.map(_tier2, deep_jobs, chunksize=1))
    t2 = time.time()
    for r, d in zip(tier1[: args.deep], tier2):
        d["kind"] = r["kind"]
        d["seed"] = r["seed"]
        d["M_tier1"] = r["M"]
    tier2.sort(key=lambda r: r["M"])
    print(f"tier2: {len(tier2)} polished in {t2-t1:.0f}s   best={tier2[0]['M']:.15f} "
          f"({tier2[0]['kind']})", flush=True)
    print(f"anchor (Together)                    {U.ANCHOR:.15f}")
    print(f"best - anchor = {tier2[0]['M'] - U.ANCHOR:+.3e}")

    out = {
        "n": args.n,
        "starts": args.starts,
        "anchor": U.ANCHOR,
        "best_M": tier2[0]["M"],
        "best_kind": tier2[0]["kind"],
        "beats_anchor": bool(tier2[0]["M"] < U.ANCHOR),
        "tier1_summary": [
            {k: r[k] for k in ("kind", "seed", "M")} for r in tier1
        ],
        "tier2": [
            {k: r[k] for k in ("kind", "seed", "M", "M_tier1", "gain")} for r in tier2
        ],
        "best_h": tier2[0]["h"],
        "seconds": t2 - t0,
    }
    with open(args.out, "w") as f:
        json.dump(out, f)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
