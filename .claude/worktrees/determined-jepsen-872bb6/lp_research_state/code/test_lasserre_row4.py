"""Test Lasserre level-2 augmentation at White Table-3 row 4 binding cell.

Row 4: h = 0.004, p = 0.3875, q ∈ [-0.02, 0.02].

Compares four configurations at the binding row 4 center:
  (A) Baseline White program (no extra tightenings).
  (B) Bochner-on-f at level n_b=20 (existing).
  (C) Lasserre level-2 alone (T_max=10).
  (D) Bochner-on-f at level n_b=20 + Lasserre level-2 (T_max=10).

Reports the deltas vs baseline and vs Bochner-only.
"""
from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")
import time
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from white_full_convex import solve_full_program

# Row 4 specs.
H = 0.004
P = 0.3875
QM, QP = -0.02, 0.02

# SDP dimensions.
N = 2000
T = 200
R = 10
LAS_T = 10        # Lasserre lift cutoff
LAS_TLOC = 10     # Localizing-matrix order
BOCH_N = 20       # f-side Bochner level


def run(label, **kwargs):
    t0 = time.time()
    res = solve_full_program(
        N, T, R, H, H, P, P, QM, QP,
        cell_mode="exact", solver="CLARABEL", verbose=False,
        **kwargs,
    )
    dt = time.time() - t0
    print(f"  [{label:30s}] Ω* = {res['value']:.7f}   "
          f"({res['status']}, {dt:.1f}s)")
    return res


if __name__ == "__main__":
    print(f"=== Row 4 Lasserre level-2 test ===")
    print(f"  N={N}, T={T}, R={R}, lasserre_T_max={LAS_T}, "
          f"lasserre_T_loc={LAS_TLOC}, bochner_n={BOCH_N}")
    print(f"  h={H}, p={P}, q∈[{QM},{QP}]")
    print()

    # (A) Baseline.
    print("(A) Baseline White program:")
    res_A = run("baseline")

    # (B) Bochner-on-f only.
    print("\n(B) Bochner-on-f (n_b=20) only:")
    res_B = run(f"bochner_n={BOCH_N}", bochner_n=BOCH_N)

    # (C) Lasserre level-2 alone.
    print(f"\n(C) Lasserre level-2 (T_max={LAS_T}, T_loc={LAS_TLOC}) only:")
    res_C = run(f"lasserre T_max={LAS_T}", lasserre_T_max=LAS_T, lasserre_T_loc=LAS_TLOC)

    # (D) Both.
    print(f"\n(D) Bochner-on-f (n_b=20) + Lasserre level-2 (T_max={LAS_T}):")
    res_D = run(
        f"bochner+lasserre",
        bochner_n=BOCH_N, lasserre_T_max=LAS_T, lasserre_T_loc=LAS_TLOC,
    )

    print()
    print("=== Summary ===")
    print(f"  (A) baseline             :  {res_A['value']:.7f}")
    print(f"  (B) Bochner only         :  {res_B['value']:.7f}     "
          f"Δ_A = {res_B['value'] - res_A['value']:+.7f}")
    print(f"  (C) Lasserre only        :  {res_C['value']:.7f}     "
          f"Δ_A = {res_C['value'] - res_A['value']:+.7f}")
    print(f"  (D) Bochner + Lasserre   :  {res_D['value']:.7f}     "
          f"Δ_A = {res_D['value'] - res_A['value']:+.7f}     "
          f"Δ_B = {res_D['value'] - res_B['value']:+.7f}")
    print()
    print(f"KEY REPORT:  Δ over Bochner-only at row 4 center = "
          f"{res_D['value'] - res_B['value']:+.7e}")
    print(f"             (positive Δ means Lasserre-2 tightens beyond Bochner)")
