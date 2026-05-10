"""Test M-side Bochner via Lasserre-lifted bilinears at White Table-3 row 4.

Compares (per task spec):
  (B0) f-side Bochner_n=20 only                                   ~0.378187
  (B1) + Lasserre T_max=20                                        ~0.378267
  (B2) + Lasserre T_max=20 + M-side Bochner_n=10  (NEW)
  (B3) + Lasserre T_max=20 + M-side Bochner_n=20  (NEW)

Row 4: h=0.004, p=0.3875, q∈[-0.02, 0.02], N=2000, T=200, R=10.

Reports the Δ over (B1) — the f-side Bochner + Lasserre baseline — for the
two M-side configurations.
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")
import json
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from white_full_convex import solve_full_program

# Row 4 specs.
H = 0.004
P = 0.3875
QM, QP = -0.02, 0.02

# Program dimensions.
N = 2000
T = 200
R = 10
BOCH_N = 20            # f-side Bochner level
LAS_T = 20             # Lasserre lift cutoff
LAS_TLOC = 10          # Localizing-matrix order (smaller than T_max for solve speed)


def run(label, **kwargs):
    t0 = time.time()
    res = solve_full_program(
        N, T, R, H, H, P, P, QM, QP,
        cell_mode="exact", solver="CLARABEL", verbose=False,
        **kwargs,
    )
    dt = time.time() - t0
    print(f"  [{label:55s}] Ω* = {res['value']:.7f}   "
          f"({res['status']}, {dt:.1f}s)")
    return res, dt


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    print(f"=== Row 4 M-side Bochner via Lasserre test ===")
    print(f"  N={N}, T={T}, R={R}")
    print(f"  bochner_n={BOCH_N}, lasserre_T_max={LAS_T}, "
          f"lasserre_T_loc={LAS_TLOC}")
    print(f"  h={H}, p={P}, q∈[{QM},{QP}]")
    print()

    out = {
        "config": {
            "N": N, "T": T, "R": R, "h": H, "p": P, "q": [QM, QP],
            "bochner_n": BOCH_N, "lasserre_T_max": LAS_T,
            "lasserre_T_loc": LAS_TLOC,
        },
        "runs": {},
    }

    # (B0) f-side Bochner_n=20 only.
    print("(B0) f-side Bochner only (n_b=20):")
    resB0, tB0 = run(f"bochner_n={BOCH_N}", bochner_n=BOCH_N)
    out["runs"]["B0_bochner_only"] = {
        "value": resB0["value"], "status": resB0["status"], "time_s": tB0,
    }

    # (B1) + Lasserre T_max=20 (the task says T_max=20; T_loc smaller for cost).
    print(f"\n(B1) + Lasserre level-2 (T_max={LAS_T}, T_loc={LAS_TLOC}):")
    resB1, tB1 = run(
        f"bochner_n={BOCH_N} + lasserre_T_max={LAS_T}",
        bochner_n=BOCH_N, lasserre_T_max=LAS_T, lasserre_T_loc=LAS_TLOC,
    )
    out["runs"]["B1_lasserre"] = {
        "value": resB1["value"], "status": resB1["status"], "time_s": tB1,
    }

    # (B2) + M-side Bochner via Lasserre at n_M=10.
    print(f"\n(B2) + Lasserre + M-side Bochner via Lasserre (n_M=10):")
    resB2, tB2 = run(
        f"bochner_n={BOCH_N} + lasserre + mside_lasserre_n=10",
        bochner_n=BOCH_N, lasserre_T_max=LAS_T, lasserre_T_loc=LAS_TLOC,
        mside_bochner_lasserre_n=10,
    )
    out["runs"]["B2_mside_n10"] = {
        "value": resB2["value"], "status": resB2["status"], "time_s": tB2,
    }

    # (B3) + M-side Bochner via Lasserre at n_M=20 (full T_max).
    print(f"\n(B3) + Lasserre + M-side Bochner via Lasserre (n_M=20):")
    resB3, tB3 = run(
        f"bochner_n={BOCH_N} + lasserre + mside_lasserre_n=20",
        bochner_n=BOCH_N, lasserre_T_max=LAS_T, lasserre_T_loc=LAS_TLOC,
        mside_bochner_lasserre_n=20,
    )
    out["runs"]["B3_mside_n20"] = {
        "value": resB3["value"], "status": resB3["status"], "time_s": tB3,
    }

    print()
    print("=== Summary ===")
    print(f"  (B0) Bochner only           :  {resB0['value']:.7f}    ({resB0['status']})")
    print(f"  (B1) + Lasserre T_max={LAS_T}    :  {resB1['value']:.7f}    "
          f"Δ_B0 = {resB1['value'] - resB0['value']:+.3e}    ({resB1['status']})")
    print(f"  (B2) + M-side n_M=10        :  {resB2['value']:.7f}    "
          f"Δ_B1 = {resB2['value'] - resB1['value']:+.3e}    ({resB2['status']})")
    print(f"  (B3) + M-side n_M=20        :  {resB3['value']:.7f}    "
          f"Δ_B1 = {resB3['value'] - resB1['value']:+.3e}    ({resB3['status']})")
    print()

    out["summary"] = {
        "Delta_B2_over_B1": resB2["value"] - resB1["value"],
        "Delta_B3_over_B1": resB3["value"] - resB1["value"],
        "Delta_B3_over_B0": resB3["value"] - resB0["value"],
    }

    if args.out:
        with open(args.out, "w") as f:
            json.dump(out, f, indent=2, default=float)
        print(f"  Wrote {args.out}")
