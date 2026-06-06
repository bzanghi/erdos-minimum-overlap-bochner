"""
_jansson_prod.py — PRODUCTION-N Jansson rigorous-LB driver (one point at a time).

Extends the N<=3000 ladder (_jansson_ladder.py) toward production scale
(N=10000, then 20000) at the binding centers, using the SAME augmentation config
as the verified full-space cover:  T~0.4N, R=10, bochner_n=40, pm_k_max=20.

DISK-FIRST: appends each result to docs/RND_WHITESPACE/L2_PROD.json immediately.
MEMORY-AWARE: runs ONE solve per invocation (pass --N and --center); the caller
serializes invocations so only one ~8GB solve is resident at a time.

Each point reports p_lo (rigorous Jansson LB on SDP_opt(center)), prob.value,
penalty breakdown, and margins vs White 0.379005 and vs headline 0.380284.
"""
from __future__ import annotations
import sys, json, time, warnings, argparse
from pathlib import Path
warnings.filterwarnings("ignore")

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _jansson_verify import jansson_lower_bound

REPO = CODE.parent.parent
OUT_JSON = REPO / "docs" / "RND_WHITESPACE" / "L2_PROD.json"

CENTERS = {
    "row4": dict(h_c=0.004, p_c=0.3875, q1=-0.02, q2=0.02),
    "cde_n30_iter3": dict(h_c=0.000045, p_c=0.39015, q1=-0.02, q2=0.02),
}


def load():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"meta": {"charter": "L2 Jansson rigorous a-posteriori LB; PRODUCTION N",
                     "WHITE": 0.379005, "PRIOR_PUB": 0.379544, "HEADLINE": 0.380284},
            "runs": []}


def save(db):
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(db, indent=2, default=float))


def run_point(cname, N, T, R, bn, pm, use_T5p=False, slack_infl=1.0):
    db = load()
    key = dict(center=cname, N=N, T=T, R=R, bochner_n=bn, pm_k_max=pm,
               use_T5p=use_T5p, slack_infl=slack_infl)
    if any(r.get("key") == key for r in db["runs"]):
        print(f"[skip] already have {cname} N={N} T={T} bn={bn} pm={pm}")
        return
    cc = CENTERS[cname]
    print(f"\n### {cname}  N={N} T={T} R={R} bn={bn} pm={pm} ###", flush=True)
    t0 = time.time()
    try:
        r = jansson_lower_bound(N, T, R, cc["h_c"], cc["p_c"], cc["q1"], cc["q2"],
                                bn, pm, use_T5p=use_T5p, slack_infl=slack_infl,
                                verbose=True)
        rec = {"key": key, "ok": True,
               "prob_value": r["prob_value"], "p_lo": r["p_lo"],
               "penalty_total": r["penalty_total"],
               "penalty_Dx_upper": r["penalty_Dx_upper"],
               "penalty_zs_lower": r["penalty_zs_lower"],
               "defect_inf": r["defect_inf"],
               "defect_inf_rigorous": r.get("defect_inf_rigorous"),
               "defect_1": r["defect_1"],
               "margin_vs_white": r["margin_vs_white"],
               "margin_vs_prior_pub": r["margin_vs_prior_pub"],
               "margin_vs_headline": r["margin_vs_headline"],
               "psd_lambda_min_lowers": [bb["d_lower"] for bb in r["block_report"]
                                         if bb["kind"] == "psd"],
               "xbar_flags": r["xbar_flags"],
               "self_checks": {k: r["self_checks"][k] for k in
                               ("|c@x - obj_val|", "|-b@z - obj_val_dual|",
                                "r_prim", "r_dual", "status")},
               "solve_time_s": r["solve_time_s"],
               "wall_s": time.time() - t0}
    except Exception as e:
        import traceback
        rec = {"key": key, "ok": False, "error": f"{type(e).__name__}: {e}",
               "tb": traceback.format_exc(), "wall_s": time.time() - t0}
        print("ERROR:", e, flush=True)
    db = load()            # reload in case of concurrent writes
    db["runs"].append(rec)
    save(db)
    print(f"[saved] {OUT_JSON}", flush=True)
    return rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--center", type=str, default="row4", choices=list(CENTERS))
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--T", type=int, default=None)   # default 0.4*N
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=40)
    ap.add_argument("--pm_k_max", type=int, default=20)
    ap.add_argument("--use_T5p", action="store_true")
    ap.add_argument("--slack_infl", type=float, default=1.0)
    args = ap.parse_args()
    T = args.T if args.T is not None else int(round(0.4 * args.N))
    run_point(args.center, args.N, T, args.R, args.bochner_n, args.pm_k_max,
              use_T5p=args.use_T5p, slack_infl=args.slack_infl)
