"""
_jansson_all12.py — certify ALL 12 core (5.16) anchors with the Jansson
interval-arithmetic a-posteriori lower bound at the production config.

WHY THIS EXISTS
---------------
The headline `mu >= 0.3802838` is the min over the core box of the max-of-12
quadratic envelope built from the 12 anchors in
`parallel_results/phase5_N20K_bn40_dualext.json`.  Each anchor is used at the
CONSERVATIVE value `primal - 1e-5`.  That convention is only justified if, for
every anchor, the true SDP optimum at that center is >= `primal - 1e-5`.

Before this script, exactly TWO of the twelve had a rigorous certificate
(`docs/RND_WHITESPACE/L2_PROD.json`: row4 and cde_n30_iter3).  The other ten
rested on CLARABEL's own reported value and on `dual_extractor.py`, whose
"rigorous_dual_LB" is the solver's dual objective with NO correction for dual
infeasibility (its own docstring concedes a margin "can be absorbed" but none
ever is).  That is not a proof.

This driver closes the gap: it runs `_jansson_verify.jansson_lower_bound` at
each of the 12 centers and records whether `p_lo >= primal_dualext - 1e-5`.
Every anchor that clears is rigorously justified; any that fails invalidates
the envelope and must be re-anchored at its own `p_lo`.

DISK-FIRST and MEMORY-AWARE: one solve per invocation of `run_one`, appended to
the output JSON immediately, so an OOM kill loses at most one center.

    ../../.venv/bin/python _jansson_all12.py --all
    ../../.venv/bin/python _jansson_all12.py --center row1
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

from _jansson_verify import jansson_lower_bound  # noqa: E402

REPO = CODE.parent.parent
DUALEXT = CODE.parent / "parallel_results" / "phase5_N20K_bn40_dualext.json"
OUT_JSON = CODE.parent / "parallel_results" / "jansson_core12.json"

ANCHOR_MARGIN = 1e-5  # the envelope's anchor convention: primal - 1e-5


def load_anchors():
    """The 12 core anchors, verbatim from the file the envelope actually reads."""
    d = json.loads(DUALEXT.read_text())
    cfg = d["config"]
    out = []
    for c in d["centers"]:
        out.append({
            "label": c["label"],
            "h_c": c["h_c"], "p_c": c["p_c"],
            "q1": c["q1"], "q2": c["q2"],
            "primal_dualext": c["primal"],
            "anchor": c["primal"] - ANCHOR_MARGIN,
        })
    return out, cfg


def load_db():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"meta": {
        "charter": "Jansson interval-arithmetic p_lo for all 12 core (5.16) anchors",
        "anchor_convention": "primal_dualext - 1e-5",
        "headline": 0.3802838,
        "dualext_source": str(DUALEXT.name),
    }, "runs": []}


def save_db(db):
    OUT_JSON.write_text(json.dumps(db, indent=2, default=float))


def run_one(label, N, T, R, bn, pm):
    anchors, cfg = load_anchors()
    a = next((x for x in anchors if x["label"] == label), None)
    if a is None:
        raise SystemExit(f"unknown center {label!r}; have {[x['label'] for x in anchors]}")

    db = load_db()
    key = dict(center=label, N=N, T=T, R=R, bochner_n=bn, pm_k_max=pm)
    if any(r.get("key") == key and r.get("ok") for r in db["runs"]):
        print(f"[skip] already certified {label} at {key}", flush=True)
        return

    print(f"\n### {label}  h_c={a['h_c']} p_c={a['p_c']} q=({a['q1']},{a['q2']})"
          f"  N={N} T={T} bn={bn} pm={pm} ###", flush=True)
    t0 = time.time()
    try:
        r = jansson_lower_bound(N, T, R, a["h_c"], a["p_c"], a["q1"], a["q2"],
                                bn, pm, use_T5p=False, slack_infl=1.0, verbose=True)
        p_lo = r["p_lo"]
        clears = bool(p_lo >= a["anchor"])
        rec = {
            "key": key, "ok": True,
            "h_c": a["h_c"], "p_c": a["p_c"], "q1": a["q1"], "q2": a["q2"],
            "primal_dualext": a["primal_dualext"],
            "anchor": a["anchor"],
            "prob_value": r["prob_value"],
            "p_lo": p_lo,
            "p_lo_minus_anchor": p_lo - a["anchor"],
            "anchor_justified": clears,
            "penalty_total": r["penalty_total"],
            "defect_inf_rigorous": r.get("defect_inf_rigorous"),
            "self_checks": {k: r["self_checks"][k] for k in
                            ("|c@x - obj_val|", "|-b@z - obj_val_dual|",
                             "r_prim", "r_dual", "status")},
            "solve_time_s": r["solve_time_s"],
            "wall_s": time.time() - t0,
        }
        print(f"[{label}] p_lo={p_lo:.12f}  anchor={a['anchor']:.12f}  "
              f"margin={p_lo - a['anchor']:+.3e}  "
              f"{'JUSTIFIED' if clears else '*** ANCHOR NOT JUSTIFIED ***'}", flush=True)
    except Exception as e:
        rec = {"key": key, "ok": False, "error": f"{type(e).__name__}: {e}",
               "tb": traceback.format_exc(), "wall_s": time.time() - t0}
        print(f"[{label}] ERROR: {e}", flush=True)

    db = load_db()
    db["runs"].append(rec)
    save_db(db)
    print(f"[saved] {OUT_JSON}", flush=True)


def summary():
    anchors, _ = load_anchors()
    db = load_db()
    by = {}
    for r in db["runs"]:
        if r.get("ok"):
            by[r["key"]["center"]] = r
    print(f"\n{'center':<16} {'anchor':>14} {'p_lo':>16} {'margin':>12}  verdict")
    worst = None
    for a in anchors:
        r = by.get(a["label"])
        if r is None:
            print(f"{a['label']:<16} {a['anchor']:>14.9f} {'--':>16} {'--':>12}  NOT CERTIFIED")
            continue
        m = r["p_lo_minus_anchor"]
        worst = m if worst is None else min(worst, m)
        print(f"{a['label']:<16} {a['anchor']:>14.9f} {r['p_lo']:>16.12f} "
              f"{m:>+12.3e}  {'ok' if r['anchor_justified'] else 'FAIL'}")
    n_ok = sum(1 for a in anchors if by.get(a["label"], {}).get("anchor_justified"))
    print(f"\n{n_ok}/{len(anchors)} anchors rigorously justified"
          + (f"; worst margin {worst:+.3e}" if worst is not None else ""))
    return n_ok == len(anchors)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--center", type=str, default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--summary", action="store_true")
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=40)
    ap.add_argument("--pm_k_max", type=int, default=20)
    args = ap.parse_args()

    if args.summary:
        summary()
    elif args.all:
        anchors, _ = load_anchors()
        for a in anchors:
            run_one(a["label"], args.N, args.T, args.R, args.bochner_n, args.pm_k_max)
        summary()
    elif args.center:
        run_one(args.center, args.N, args.T, args.R, args.bochner_n, args.pm_k_max)
    else:
        ap.error("pass --center LABEL, --all, or --summary")
