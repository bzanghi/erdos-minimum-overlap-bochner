"""
Run the Jansson rigorous-LB ladder over escalating N at binding center(s) and
write incremental results to docs/RND_WHITESPACE/L2_RESULT.json AS WE GO
(disk-first: nothing lost on crash).

Caps N<=3000 per the L2 charter (production-scale solves are a separate step).
"""
from __future__ import annotations
import sys, json, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
from _jansson_verify import jansson_lower_bound

REPO = CODE.parent.parent
OUT_JSON = REPO / "docs" / "RND_WHITESPACE" / "L2_RESULT.json"

CENTERS = {
    "row4": dict(h_c=0.004, p_c=0.3875, q1=-0.02, q2=0.02),
    "cde_n30_iter3": dict(h_c=0.000045, p_c=0.39015, q1=-0.02, q2=0.02),
}

# (N, T, R, bochner_n, pm_k_max) ladder.  T grows with N (T ~ 0.4 N), R fixed-ish.
LADDER = [
    (300, 120, 6, 6, 14),
    (1000, 400, 8, 10, 14),
    (2000, 800, 10, 12, 14),
    (3000, 1200, 10, 16, 14),
]


def load():
    if OUT_JSON.exists():
        return json.loads(OUT_JSON.read_text())
    return {"meta": {"charter": "L2 Jansson rigorous a-posteriori LB; N<=3000",
                     "WHITE": 0.379005, "PRIOR_PUB": 0.379544, "HEADLINE": 0.380284},
            "runs": []}


def save(db):
    OUT_JSON.write_text(json.dumps(db, indent=2, default=float))


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--centers", nargs="+", default=["row4", "cde_n30_iter3"])
    ap.add_argument("--use_T5p", action="store_true")
    ap.add_argument("--slack_infl", type=float, default=1.0)
    ap.add_argument("--max_N", type=int, default=3000)
    args = ap.parse_args()

    db = load()
    for cname in args.centers:
        cc = CENTERS[cname]
        for (N, T, R, bn, pm) in LADDER:
            if N > args.max_N:
                continue
            key = dict(center=cname, N=N, T=T, R=R, bochner_n=bn, pm_k_max=pm,
                       use_T5p=args.use_T5p, slack_infl=args.slack_infl)
            # dedupe
            if any(r.get("key") == key for r in db["runs"]):
                print(f"[skip] already have {cname} N={N}")
                continue
            print(f"\n### {cname}  N={N} T={T} R={R} bn={bn} pm={pm} ###", flush=True)
            t0 = time.time()
            try:
                r = jansson_lower_bound(N, T, R, cc["h_c"], cc["p_c"],
                                        cc["q1"], cc["q2"], bn, pm,
                                        use_T5p=args.use_T5p,
                                        slack_infl=args.slack_infl, verbose=True)
                rec = {"key": key, "ok": True,
                       "prob_value": r["prob_value"],
                       "p_lo": r["p_lo"],
                       "penalty_total": r["penalty_total"],
                       "penalty_Dx_upper": r["penalty_Dx_upper"],
                       "penalty_zs_lower": r["penalty_zs_lower"],
                       "defect_inf": r["defect_inf"],
                       "defect_1": r["defect_1"],
                       "margin_vs_white": r["margin_vs_white"],
                       "margin_vs_prior_pub": r["margin_vs_prior_pub"],
                       "margin_vs_headline": r["margin_vs_headline"],
                       "psd_lambda_min_lowers": [
                           bb["d_lower"] for bb in r["block_report"]
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
            db["runs"].append(rec)
            save(db)   # <-- disk-first, after every solve
            print(f"[saved] {OUT_JSON}", flush=True)
    print("\nDONE.")


if __name__ == "__main__":
    main()
