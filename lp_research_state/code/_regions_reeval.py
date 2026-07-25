"""
_regions_reeval.py — re-certify the GATE regions against a new anchor set.

WHY
---
Each gate region has its own `_eval_r*_box.py` with the subdivision depth and the
target hardcoded, and each stores a floor produced against whatever core anchors
existed at the time.  Two consequences bit us after the N=48000 re-certification:

 1. The stored floors are STALE.  They remain *valid* — the new anchors are
    strictly higher, so every cover only improved — but they are the numbers
    `_fs_recompute.load_certified_region_floors()` reads, so the full-space
    minimum is pinned to an old, lower value.  After N=48000 the core floor is
    0.3803954 while the stored region floors sit at 0.3803090 (R6) through
    0.3805539 (R7), so a REGION binds rather than the core.

 2. The evaluators STOP as soon as a sub-box clears `target`.  A reported region
    floor therefore means "at least this", never "this is the infimum".  Raising
    the target makes them subdivide further and report higher, still-rigorous
    floors.  R6 re-run at the new anchors went 0.3803090 -> 0.380344 and stopped
    on `max_depth=8`, not on mathematics — its own ceiling (the eps-free grid
    min) is 0.380814.

So the fix is not new mathematics, it is re-running the same rigorous routine
with the new anchors, a raised target, and more depth.

WHAT THIS DOES
--------------
For each gate region: core anchors (from $LP_DUALEXT) + that region's fresh
promotion centers, fed to `_eval_r6_box.adaptive_boxmin` — the same rigorous
grid+Lipschitz box-min with per-sub-box `L_max` that produced the stored numbers.
Nothing here is a new bound mechanism; only the depth and the anchors change.

    LP_DUALEXT=../parallel_results/dualext_reanchored_N48000.json \\
      ../../.venv/bin/python _regions_reeval.py --target 0.3803954 --max-depth 14
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

from _fullspace_eval import load_centers, WHITE_TABLE2          # noqa: E402
from _eval_r6_box import adaptive_boxmin                        # noqa: E402

PR = CODE.parent / "parallel_results"
POOL = []

# region -> (promote file, centers to drop).  The drops are the ones the original
# evaluators excluded to keep L_max controlled: spiky deep-q centers whose
# gradient inflates eps_grid over the whole box while covering <2% of the grid.
GATE = {
    6:  ("fullspace_promote_R6.json",  ("R6_c2_p05_q05", "R6_c3_p05_q20")),
    7:  ("fullspace_promote_R7.json",  ()),
    9:  ("fullspace_promote_R9.json",  ()),
    16: ("fullspace_promote_R16.json", ()),
    17: ("fullspace_promote_R17.json", ()),
}


def fresh_centers(region):
    """That region's own promotion centers.

    NOTE the key varies by file: R7 stores `dual_lb_raw` where the others store
    `dual_lb`.  Filtering on `dual_lb` alone silently dropped all 7 of R7's
    centers and evaluated the region on the core anchors only, which is exactly
    the documented failure case (core-alone puts these corners at 0.3802561).
    What actually matters is that the center carries `duals` and a `primal`, so
    that is what we require.
    """
    fname, drop = GATE[region]
    p = PR / fname
    if not p.exists():
        return []
    d = json.loads(p.read_text())
    out = []
    for c in d.get("centers", []):
        if c.get("label") in drop:
            continue
        if not c.get("duals") or c.get("primal") is None:
            continue
        out.append(c)
    return out


def smooth_pool():
    """Core + halo + every promotion center, excluding the stage2 box-LP leaves.

    `_fs_recompute` calls this the `smooth` subset: the stage2 leaves are spiky
    (their con_513 duals inflate L_max ~45x), which blows up eps_grid over a
    whole box and LOWERS the rigorous floor on narrow regions even though their
    pointwise cover is identical.  Evaluating on a wider, smoother pool is what
    lets the wide regions (R6/R7/R9) clear at all -- core-alone does not.
    """
    import _fs_recompute as R
    centers, _cfg, _src = R.harvest_centers()
    centers = R.dedupe(centers)
    return [c for c in centers if c.get("_src") != "stage2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", type=float, required=True,
                    help="subdivide until each sub-box clears this")
    ap.add_argument("--max-depth", type=int, default=14)
    ap.add_argument("--base", type=int, default=21,
                    help="per-axis grid points inside each sub-box")
    ap.add_argument("--regions", type=int, nargs="*", default=sorted(GATE))
    ap.add_argument("--pool", action="store_true", default=True,
                    help="also try the smooth pool (core+halo+all promote centers)")
    ap.add_argument("--no-pool", dest="pool", action="store_false")
    ap.add_argument("--out", default=str(PR / "gate_regions_reeval.json"))
    args = ap.parse_args()

    core, cfg = load_centers()
    global POOL
    POOL = smooth_pool() if args.pool else []
    print(f"smooth pool: {len(POOL)} centers", flush=True)
    print(f"core anchors: {len(core)}  config={cfg}", flush=True)
    print(f"target={args.target}  max_depth={args.max_depth}  base={args.base}\n",
          flush=True)

    results = {}
    for reg in args.regions:
        hr, pr_, qr, _white = WHITE_TABLE2[reg - 1]
        fresh = fresh_centers(reg)
        # Evaluate on nested subsets and keep the BEST floor.  Max of valid
        # lower bounds is a valid lower bound, so a bigger pool can only help --
        # but only through the max, because a spikier pool can raise L_max and
        # thereby LOWER a given box's rigorous floor.
        subsets = [("core+fresh", core + fresh)]
        if args.pool:
            subsets.append(("smooth_pool", POOL + fresh))
        best = None
        t0 = time.time()
        for name, sub in subsets:
            r = adaptive_boxmin(sub, hr, pr_, qr, args.target,
                                max_depth=args.max_depth, base=(args.base,) * 3)
            if best is None or r[0] > best[1][0]:
                best = (name, r)
            if r[0] >= args.target:
                break
        mech, (lb, gm, pt, wit, leaves) = best
        dt = time.time() - t0
        clears = lb >= args.target
        print(f"[R{reg:2d}] h{hr} p{pr_} q{qr}", flush=True)
        print(f"      floor={lb:.7f}  ceiling(grid_min)={gm:.7f}  "
              f"leaves={leaves}  {dt:.0f}s  "
              f"{'CLEARS' if clears else 'BELOW TARGET'}", flush=True)
        print(f"      worst @ (h={pt[0]:.4f}, p={pt[1]:.4f}, q={pt[2]:.4f}) "
              f"wit={wit}  via {mech}  ({len(core)} core + {len(fresh)} fresh)\n",
              flush=True)
        results[reg] = {"region": reg, "floor": lb, "ceiling_grid_min": gm,
                        "worst_point": {"h": pt[0], "p": pt[1], "q": pt[2]},
                        "witness": str(wit), "leaves": leaves,
                        "clears_target": bool(clears), "secs": dt,
                        "mechanism": mech,
                        "n_core": len(core), "n_fresh": len(fresh)}
        Path(args.out).write_text(json.dumps(
            {"target": args.target, "max_depth": args.max_depth,
             "base": args.base, "core_config": cfg, "regions": results},
            indent=2, default=float))

    worst = min(results.values(), key=lambda r: r["floor"])
    print("=" * 70)
    print(f"gate-region minimum: {worst['floor']:.7f}  at R{worst['region']}")
    print(f"target was          : {args.target:.7f}")
    print(f"all gate regions clear the target: "
          f"{all(r['clears_target'] for r in results.values())}")
    print(f"saved -> {args.out}")


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()
