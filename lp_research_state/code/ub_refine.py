"""
ub_refine — push the certified UPPER bound down by grid refinement + a
line-searching trust-region SLP.

WHY
---
`mu = inf_n M_n` and `M_n` is non-increasing under cell doubling
(`ub_core.cell_double` preserves `M` exactly).  So refining a witness from `n`
to `2n` cells is FREE — it cannot make the bound worse — and it strictly
enlarges the search space.  The best published witness (SimpleTES, n=2400,
polished here to M = 0.3808675459609214) is therefore a valid starting point at
n = 4800, 9600, ... and any descent found there is pure gain.

WHAT THIS ADDS OVER `ub_local.slp_polish`
-----------------------------------------
`slp_polish` computes an LP step `delta`, tries `t = 1`, and on failure shrinks
the trust radius and re-solves the LP from scratch.  Near a degenerate minimax
optimum (SimpleTES has ~1580 of 4799 lags within 1e-6 of the max) most full
steps fail, so almost every expensive LP is discarded.

`slp_polish_ls` keeps the LP but adds an exact backtracking line search along
`delta`.  The objective along the ray,

    M(h + t*delta) = (2/n) max_j [ M_j(h) + t*(g_j . delta) - t^2 * q_j(delta) ],

is a max of one-dimensional quadratics, so evaluating a whole grid of `t` costs
one FFT-batched `batch_overlap_values` call.  A rejected full step becomes a
short step instead of a discarded LP solve.

The trust-region validity check from `ub_local` is preserved verbatim: a lag
excluded from the LP cannot overtake the modelled maximum, since a step of
radius `r` moves any scaled `M_j` by at most `4r` (|grad|_1 <= 4).  That check
is what makes the reported `M` an honest upper bound rather than a model value.
Every `M` reported here is the exact `ub_core.overlap_value`, never the model.

USAGE
    ../../.venv/bin/python ub_refine.py --input ../data/simpletes_polished.json \\
        --factor 2 --rounds 40 --r0 3e-7 --out ../data/ub_refined_n4800.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ub_core import (
    cell_double,
    grad_lag,
    index_of_lag,
    lag_of_index,
    overlap_profile,
    overlap_value,
    project_feasible,
)
from ub_local import batch_overlap_values

DATA = Path(__file__).resolve().parent.parent / "data"


# ------------------------------------------------------------------ line search

def _line_search(h, delta, M0, n_t=24):
    """Best exact M along h + t*delta for t in (0, 1], via one batched eval.

    Returns (t_best, M_best).  `M_best` is the exact objective, not a model.
    """
    ts = np.linspace(1.0, 1.0 / (1 << (n_t - 1)), n_t) if n_t > 1 else np.array([1.0])
    # geometric backtracking is the useful regime near a degenerate optimum
    ts = 2.0 ** (-np.arange(n_t, dtype=np.float64))
    H = np.clip(h[None, :] + ts[:, None] * delta[None, :], 0.0, 1.0)
    # renormalise each row exactly; project_feasible is cheap enough per row
    H = np.stack([project_feasible(row) for row in H])
    vals = batch_overlap_values(H)
    k = int(np.argmin(vals))
    if vals[k] >= M0:
        return None, M0, None
    return float(ts[k]), float(vals[k]), H[k]


# --------------------------------------------------------------------- polish

def slp_polish_ls(h0, rounds=40, r0=1e-6, r_min=1e-11, slack=2e-7,
                  verbose=False, time_budget_s=None, log=None):
    """Trust-region SLP with an exact line search along the LP step.

    Returns (h, M, certified_first_order_gain, history).
    """
    from scipy.optimize import linprog

    h = project_feasible(np.asarray(h0, dtype=np.float64))
    n = h.size
    M = overlap_value(h)
    r = r0
    gain = np.nan
    hist = []
    t_start = time.time()

    for it in range(rounds):
        if r < r_min:
            break
        if time_budget_s is not None and time.time() - t_start > time_budget_s:
            if verbose:
                print(f"  [time budget {time_budget_s}s reached at round {it}]", flush=True)
            break

        prof = overlap_profile(h) * (2.0 / n)
        Mmax = prof.max()
        # a step of radius r moves any scaled M_j by at most 4r (|grad|_1 <= 4)
        keep = np.nonzero(prof >= Mmax - (4.0 * r + slack))[0]
        A = np.empty((keep.size, n))
        for row, mi in enumerate(keep):
            A[row] = grad_lag(h, lag_of_index(int(mi), n))

        A_ub = np.hstack([A, -np.ones((keep.size, 1))])
        b_ub = Mmax - prof[keep]
        A_eq = np.hstack([np.ones((1, n)), np.zeros((1, 1))])
        b_eq = np.zeros(1)
        lo = np.maximum(-r, -h)
        hi = np.minimum(r, 1.0 - h)
        bounds = list(zip(lo, hi)) + [(-1.0, 1.0)]
        c = np.zeros(n + 1)
        c[-1] = 1.0

        t0 = time.time()
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=bounds, method="highs-ipm")
        lp_s = time.time() - t0
        if not res.success:
            r *= 0.5
            continue

        delta, u = res.x[:n], res.x[-1]
        gain = -u

        # validity: a lag dropped from the LP cannot overtake the modelled max
        dropped = np.delete(prof, keep)
        if dropped.size and dropped.max() + 4.0 * r >= Mmax + u:
            r *= 0.5
            continue

        t_best, M_new, h_new = _line_search(h, delta, M)
        if t_best is not None and M_new < M - 1e-16:
            h, M = h_new, M_new
            # grow the radius only when the full step was taken
            r = min(r0, r * 1.6) if t_best > 0.5 else max(r_min, r * 0.8)
            accepted = True
        else:
            r *= 0.4
            accepted = False

        rec = {"round": it, "M": M, "gain": float(gain), "r": float(r),
               "lags": int(keep.size), "t": t_best, "accepted": accepted,
               "lp_s": lp_s}
        hist.append(rec)
        if verbose:
            print(f"  [{it:3d}] M={M:.16f} gain={gain:.3e} r={r:.2e} "
                  f"lags={keep.size} t={t_best} lp={lp_s:.1f}s", flush=True)
        if log is not None:
            Path(log).write_text(json.dumps(
                {"n": n, "M": M, "history": hist}, indent=1, default=float))

    return h, M, gain, hist


# ----------------------------------------------------------------------- main

def load_witness(path):
    d = json.loads(Path(path).read_text())
    for k in ("best_h", "h", "values", "construction"):
        if k in d:
            return np.asarray(d[k], dtype=np.float64), d
    raise SystemExit(f"no recognised h key in {path}; keys={list(d)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=str(DATA / "simpletes_polished.json"))
    ap.add_argument("--factor", type=int, default=2)
    ap.add_argument("--rounds", type=int, default=40)
    ap.add_argument("--r0", type=float, default=3e-7)
    ap.add_argument("--slack", type=float, default=2e-7)
    ap.add_argument("--time-budget", type=float, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    h, meta = load_witness(args.input)
    n0 = h.size
    M0 = overlap_value(h)
    print(f"loaded n={n0}  M={M0:.16f}  (from {args.input})", flush=True)

    if args.factor > 1:
        h = cell_double(h, args.factor)
        M1 = overlap_value(h)
        print(f"refined n={h.size}  M={M1:.16f}  "
              f"(doubling drift {M1 - M0:+.2e}, must be ~1 ulp)", flush=True)
        assert abs(M1 - M0) < 1e-14, "cell doubling changed the objective"

    log = (args.out + ".partial") if args.out else None
    t0 = time.time()
    h, M, gain, hist = slp_polish_ls(
        h, rounds=args.rounds, r0=args.r0, slack=args.slack,
        verbose=True, time_budget_s=args.time_budget, log=log)
    wall = time.time() - t0

    print(f"\nFINAL n={h.size}  M={M:.16f}  gain={gain:.3e}  "
          f"improvement vs input {M - M0:+.3e}  ({wall:.0f}s)", flush=True)

    if args.out:
        Path(args.out).write_text(json.dumps({
            "best_h": h.tolist(), "M": M, "n": int(h.size),
            "source": f"{Path(args.input).name} x{args.factor} + slp_polish_ls",
            "input_M": M0, "certified_first_order_gain": float(gain),
            "rounds": len(hist), "wall_s": wall,
        }, default=float))
        print(f"saved -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
