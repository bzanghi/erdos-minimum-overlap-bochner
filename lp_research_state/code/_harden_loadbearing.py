"""PRO-38 hardening: re-solve the LOAD-BEARING R16/R17 poly-moment promotion centers
at N=24000 (up from the N=20000 they were certified at) to widen the thin 6th-decimal
margin and confirm the centers hold at higher N (SDP min rises monotonically with N, so
N=24000 >= N=20000 value is expected; a DROP would be a red flag). Disk-first: each
result is written to parallel_results/harden_loadbearing_N24K.json as it completes, so a
crash/OOM mid-run loses nothing. Sequential (one heavy solve at a time) to bound memory."""
import sys; sys.path.insert(0, '.')
import json, time, traceback
from path_b_with_polymoment import solve_with_pm

# (label, h_c, p_c, q1, q2, claimed primal at N=20000)
CENTERS = [
    ("R16_h0.0_p0.3877_qm0.02",      0.0,     0.3877, -0.02,  -0.02,  0.38044787),  # binding winner (+1.2e-4)
    ("R16_h0.0_p0.39_qm0.025",       0.0,     0.39,   -0.025, -0.025, 0.38069201),  # deep-q corner
    ("R16_h0.00375_p0.3915_qm0.025", 0.00375, 0.3915, -0.025, -0.025, 0.38073020),
    ("R17_h0.0_p0.39_q0.02",         0.0,     0.39,    0.02,   0.02,  0.38063000),  # R17 mirror
]
N, T, R, BN, PMK = 24000, 4000, 10, 40, 20
TARGET = 0.3802838
OUT = "../parallel_results/harden_loadbearing_N24K.json"
results = {}

def save():
    json.dump({"config": {"N": N, "T": T, "R": R, "bochner_n": BN, "pm_k_max": PMK,
                          "mside_sin_coeff": 4.0}, "target": TARGET, "results": results},
              open(OUT, "w"), indent=2)

print(f"=== PRO-38 hardening: re-solve {len(CENTERS)} load-bearing centers at "
      f"N={N}, bochner_n={BN}, pm_k_max={PMK} (corrected coeff 4.0) ===", flush=True)
for label, h, p, q1, q2, claimed in CENTERS:
    t0 = time.time()
    try:
        res = solve_with_pm(N, T, R, h, p, q1, q2, BN, PMK)
        val = float(res["value"]); st = res["status"]
        results[label] = {"h": h, "p": p, "q1": q1, "q2": q2, "claimed_N20K": claimed,
                          "N24K_value": val, "status": st, "delta_vs_N20K": val - claimed,
                          "margin_vs_target": val - TARGET, "sec": round(time.time() - t0)}
        print(f"{label}: N24K={val:.7f}  (N20K {claimed:.7f}, delta {val-claimed:+.2e})  "
              f"margin_vs_target={val-TARGET:+.2e}  status={st}  {time.time()-t0:.0f}s", flush=True)
    except Exception as e:
        results[label] = {"error": str(e)[:300], "trace": traceback.format_exc()[-600:]}
        print(f"{label}: ERROR {e}", flush=True)
    save()
print("DONE", flush=True)
