"""
ub_search — structured basin-hopping around a GOOD witness.

WHY THIS AND NOT `ub_basin_sweep.py`
------------------------------------
That sweep ran 1600 random multistarts at n=600 and every family landed in
[0.380924, 0.380995] — far worse than the record. Conclusion recorded at the
time: the good basin is not reachable by unstructured multistart. Fine. But two
things were never tried, and both are cheap:

 1. Perturb-and-repolish around the *arena* witness. The previous sweep's
    `perturb_hstar` family perturbed Together's n=600 h* (0.3808703). The arena
    n=512 witness is a strictly better basin (0.3808591) and was not available
    then.
 2. Non-commensurate resampling. A 512-cell step function is not representable
    on a 600- or 800-cell grid, so `ub_core.resample` lands in a genuinely
    different function space -- unlike cell doubling, which is exactly
    value-preserving and (measured, this repo) unlocks nothing.

The n=512 scale is what makes this affordable: one trust-region LP is ~0.4 s
there, versus ~60 s at n=4800. Thousands of attempts, not dozens.

MOVES
-----
Each proposal is `perturb -> project -> polish`, keeping the best exact
objective ever seen. Moves are deliberately structural, not just noise, because
noise at the scale that escapes a basin also destroys the arrangement that makes
the witness good (the sweep's `shuffle_hstar` family established that h*'s value
multiset is worthless without its arrangement):

  gauss     iid noise, scale swept over decades
  block     resample one contiguous block uniformly
  swap      exchange k random cell pairs
  reflect   reverse a contiguous block in place
  mirror    apply h -> 1 - h(2-x) on a block. This map is an exact symmetry of
            the objective (substitute u = 2-x-k in the overlap integral), so it
            moves within the level set globally, and applying it to a *block*
            perturbs while respecting the symmetry the optimum seems to have.
  lowfreq   add a low-order sinusoid
  spike     push a few cells to 0 or 1 (the optimum is partly bang-bang)

Every reported M is the exact `ub_core.overlap_value`, never a model value.

    ../../.venv/bin/python ub_search.py --input ../data/arena_lnzwz_n512.json \\
        --iters 400 --workers 8 --out ../data/ub_search_best.json
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np

from ub_core import overlap_value, project_feasible, resample
from ub_refine import slp_polish_ls

MOVES = ("gauss", "block", "swap", "reflect", "mirror", "lowfreq", "spike")


def perturb(h, kind, scale, rng):
    n = h.size
    x = h.copy()
    if kind == "gauss":
        x = x + rng.normal(0.0, scale, n)
    elif kind == "block":
        w = max(2, int(n * scale))
        i = int(rng.integers(0, max(1, n - w)))
        x[i:i + w] = rng.random(w)
    elif kind == "swap":
        k = max(1, int(n * scale))
        a = rng.integers(0, n, k); b = rng.integers(0, n, k)
        x[a], x[b] = x[b], x[a]
    elif kind == "reflect":
        w = max(2, int(n * scale))
        i = int(rng.integers(0, max(1, n - w)))
        x[i:i + w] = x[i:i + w][::-1]
    elif kind == "mirror":
        # h -> 1 - h(2-x) is an exact symmetry of M; apply it on a block only
        w = max(2, int(n * scale))
        i = int(rng.integers(0, max(1, n - w)))
        x[i:i + w] = 1.0 - x[i:i + w][::-1]
    elif kind == "lowfreq":
        k = int(rng.integers(1, 12))
        t = np.arange(n) / n
        x = x + scale * np.sin(2 * np.pi * k * t + rng.random() * 2 * np.pi)
    elif kind == "spike":
        k = max(1, int(n * scale))
        idx = rng.integers(0, n, k)
        x[idx] = (rng.random(k) < 0.5).astype(float)
    else:
        raise ValueError(kind)
    return project_feasible(np.clip(x, 0.0, 1.0))


def _one(args):
    h0, seed, screen_rounds, r0, slack = args
    rng = np.random.default_rng(seed)
    kind = MOVES[int(rng.integers(len(MOVES)))]
    scale = float(10.0 ** rng.uniform(-3.0, -0.7))
    hp = perturb(h0, kind, scale, rng)
    h, M, gain, _ = slp_polish_ls(hp, rounds=screen_rounds, r0=r0, slack=slack)
    return {"M": float(M), "kind": kind, "scale": scale, "seed": int(seed),
            "h": h.tolist()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--iters", type=int, default=400)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    ap.add_argument("--screen-rounds", type=int, default=12)
    ap.add_argument("--deep-rounds", type=int, default=120)
    ap.add_argument("--deep", type=int, default=12)
    ap.add_argument("--r0", type=float, default=1e-5)
    ap.add_argument("--slack", type=float, default=1e-5)
    ap.add_argument("--resample-to", type=int, nargs="*", default=[],
                    help="also try non-commensurate grids, e.g. 600 640 768 800")
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    d = json.loads(Path(args.input).read_text())
    h0 = np.asarray(d.get("best_h") or d["values"], dtype=np.float64)
    M0 = overlap_value(h0)
    print(f"start n={h0.size}  M={M0:.16f}", flush=True)

    results = []
    t0 = time.time()

    # ---- non-commensurate resampling (cheap, qualitatively different space) --
    for m in args.resample_to:
        hm = project_feasible(resample(h0, m))
        pre = overlap_value(hm)
        h, M, gain, _ = slp_polish_ls(hm, rounds=args.deep_rounds,
                                      r0=args.r0, slack=args.slack)
        print(f"resample {h0.size} -> {m:5d}: {pre:.16f} -> {M:.16f}", flush=True)
        results.append({"M": float(M), "kind": f"resample{m}", "scale": 0.0,
                        "seed": -1, "h": h.tolist()})

    # ---- screening pass ------------------------------------------------------
    jobs = [(h0, args.seed0 + i, args.screen_rounds, args.r0, args.slack)
            for i in range(args.iters)]
    if args.workers > 1:
        import multiprocessing as mp
        with mp.Pool(args.workers) as pool:
            for i, r in enumerate(pool.imap_unordered(_one, jobs, chunksize=1)):
                results.append(r)
                if r["M"] < M0:
                    print(f"  [{i:4d}] BEAT START  M={r['M']:.16f}  "
                          f"{r['kind']}@{r['scale']:.2e}  ({r['M']-M0:+.2e})", flush=True)
                elif i % 25 == 0:
                    best = min(x["M"] for x in results)
                    print(f"  [{i:4d}] best so far {best:.16f}  "
                          f"({time.time()-t0:.0f}s)", flush=True)
    else:
        for i, j in enumerate(jobs):
            results.append(_one(j))

    results.sort(key=lambda r: r["M"])
    print(f"\nscreen best {results[0]['M']:.16f} via {results[0]['kind']} "
          f"@ {results[0]['scale']:.2e}   ({time.time()-t0:.0f}s)", flush=True)

    # ---- deep polish of the survivors ---------------------------------------
    best_h, best_M, best_tag = h0, M0, "input"
    for r in results[:args.deep]:
        h = np.asarray(r["h"])
        h, M, gain, _ = slp_polish_ls(h, rounds=args.deep_rounds,
                                      r0=args.r0, slack=args.slack)
        tag = f"{r['kind']}@{r['scale']:.1e}"
        flag = ""
        if M < best_M:
            best_h, best_M, best_tag = h, M, tag
            flag = "  <== NEW BEST"
            # checkpoint immediately: a 13-minute screening pass is expensive to
            # lose, and this process has been killed mid-run before
            if args.out:
                Path(args.out).write_text(json.dumps({
                    "best_h": best_h.tolist(), "M": best_M, "n": int(best_h.size),
                    "source": f"ub_search from {Path(args.input).name} via {tag}",
                    "input_M": M0, "partial": True}, default=float))
        print(f"  deep {tag:24s} -> {M:.16f}{flag}", flush=True)

    print(f"\nBEST  M={best_M:.16f}  via {best_tag}  vs input {M0:.16f}  "
          f"({best_M - M0:+.3e})", flush=True)

    if args.out:
        Path(args.out).write_text(json.dumps({
            "best_h": best_h.tolist(), "M": best_M, "n": int(best_h.size),
            "source": f"ub_search from {Path(args.input).name} via {best_tag}",
            "input_M": M0, "iters": args.iters,
            "screen_top": [{k: r[k] for k in ("M", "kind", "scale", "seed")}
                           for r in results[:25]],
        }, default=float))
        print(f"saved -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
