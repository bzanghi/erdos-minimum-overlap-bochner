"""Exact M(n) for the Erdős minimum overlap problem via HiGHS ILP.

Convention (Haugland; matches lp_research_state/data/known_Mn_values.json):
    Full set [1..2n], |A| = |B| = n, B = [2n] \\ A.
    overlap_k(A) = |{i : i ∈ A, i - k ∈ B, i ∈ [1, 2n], i - k ∈ [1, 2n]}|
    M(n) = min_A max_{k ∈ Z \\ {0}} overlap_k(A).
    µ = lim M(n)/n ∈ [0.380128, 0.380871].

ILP formulation
---------------
Binary x_i ∈ {0,1}, i = 1..N (N = 2n); x_i = 1 ⇔ i ∈ A.
    Σ x_i = n.
For each shift k ∈ Z \\ {0} with overlap-domain nonempty:
    overlap_k = Σ_{i in I_k} y_{i,k},   y_{i,k} = x_i · (1 - x_{i-k})
                                     = x_i - x_i · x_{i-k}.
Linearize with z_{i,k} := x_i · x_{i-k} (the AND), then y_{i,k} = x_i - z_{i,k}.
Each z_{i,k} needs only: z ≥ x_i + x_{i-k} - 1, z ≥ 0  (z ≤ x_i, z ≤ x_{i-k}
are unnecessary in a minimization with maximization adversary since increasing
z DECREASES y, which DECREASES t — the solver will saturate z to its lower
bound). [Equivalently we keep all four, harmless.]

Equivalently and more directly: y_{i,k} ∈ [0, 1] with
    y_{i,k} ≥ x_i - x_{i-k}        (y ≥ 1 if x_i=1, x_{i-k}=0; else ≥ 0 or negative)
    y_{i,k} ≥ 0
Then minimizing t with t ≥ Σ y_{i,k} forces y to its lower envelope, which is
exactly max(0, x_i - x_{i-k}) = x_i(1 - x_{i-k}) for integer x. So y can be
*continuous*, no need to bound from above.

We use the latter: half as many auxiliary continuous vars, no integrality on y.

Objective: minimize t  (integer).
Symmetry break: fix x_1 = 1 (breaks A ↔ B complement symmetry).

Returns M(n), optimal A (1-indexed), MIP gap at termination, wall time.
"""
from __future__ import annotations

import time
from typing import Optional

import highspy
import numpy as np


def _max_overlap(A: tuple[int, ...], N: int) -> int:
    """Direct verification: max over nonzero k of overlap_k(A) in [1, N]."""
    Aset = set(A)
    Bset = set(range(1, N + 1)) - Aset
    best = 0
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        count = sum(1 for i in Aset if (i - k) in Bset)
        if count > best:
            best = count
    return best


def solve_Mn_ilp(
    n: int,
    time_limit: float = 1800.0,
    threads: int = 0,
    verbose: bool = False,
    fix_x1: bool = True,
) -> dict:
    """Solve M(n) exactly with HiGHS.

    Parameters
    ----------
    n : int
        Half the set size (set = [1, 2n], |A| = n).
    time_limit : float
        Wall-clock seconds. If hit, returns the best feasible solution found
        and the proven LB; the returned `optimal` flag indicates whether the
        gap was closed.
    threads : int
        0 = HiGHS default.
    fix_x1 : bool
        Break A↔B complement symmetry by fixing x_1 = 1. Safe because if (A, t)
        is feasible, so is ([2n]\\A, t) (overlap_k(A) = overlap_{-k}([2n]\\A)
        plus boundary corrections — actually exactly equal because shifting the
        full set is symmetric; this remains a valid symmetry break).
    """
    N = 2 * n
    h = highspy.Highs()
    if not verbose:
        h.silent()
    # h.setOptionValue("output_flag", verbose)
    if threads:
        h.setOptionValue("threads", threads)
    h.setOptionValue("time_limit", time_limit)
    # Tighter MIP tolerances since the problem is purely integer/0-1.
    h.setOptionValue("mip_feasibility_tolerance", 1e-9)
    h.setOptionValue("primal_feasibility_tolerance", 1e-9)
    h.setOptionValue("dual_feasibility_tolerance", 1e-9)

    lp = highspy.HighsLp()

    # Variable layout:
    # [0 .. N-1]: x_i (binary)
    # [N]:        t (integer, [0, n])
    # [N+1 ..]:   y_{i,k} (continuous, [0, 1]) — one per (k, i) with i ∈ I_k.
    # We track y indices per k for the t ≥ Σ y_{i,k} constraint.

    # Enumerate all (k, i_lo, i_hi).
    shift_blocks = []  # list of (k, i_lo, i_hi, y_start_idx)
    y_count = 0
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        i_lo = max(1, 1 + k)
        i_hi = min(N, N + k)
        if i_lo > i_hi:
            continue
        shift_blocks.append((k, i_lo, i_hi, N + 1 + y_count))
        y_count += (i_hi - i_lo + 1)

    n_vars = N + 1 + y_count
    lp.num_col_ = n_vars
    col_cost = np.zeros(n_vars)
    col_cost[N] = 1.0  # minimize t
    col_lower = np.zeros(n_vars)
    col_upper = np.ones(n_vars)
    col_upper[N] = float(n)
    # x integer, t integer, y continuous
    integrality = [highspy.HighsVarType.kContinuous] * n_vars
    for i in range(N + 1):  # x_1..x_N and t
        integrality[i] = highspy.HighsVarType.kInteger
    lp.col_cost_ = col_cost
    lp.col_lower_ = col_lower
    lp.col_upper_ = col_upper
    lp.integrality_ = integrality
    lp.sense_ = highspy.ObjSense.kMinimize

    # Constraints (build CSR / row-wise sparse):
    # 1) Σ x_i = n.
    # 2) Optionally x_1 = 1 (handled via bounds).
    # 3) For each (k, i ∈ I_k): y_{i,k} ≥ x_i - x_{i-k}, i.e.
    #       -x_i + x_{i-k} + y_{i,k} ≥ 0.
    # 4) For each k: t - Σ_i y_{i,k} ≥ 0  i.e.  Σ y - t ≤ 0.

    if fix_x1:
        col_lower[0] = 1.0
        col_upper[0] = 1.0
        lp.col_lower_ = col_lower
        lp.col_upper_ = col_upper

    rows_start = [0]
    rows_index = []
    rows_value = []
    row_lower = []
    row_upper = []

    # (1) Σ x = n
    for i in range(N):
        rows_index.append(i)
        rows_value.append(1.0)
    row_lower.append(float(n))
    row_upper.append(float(n))
    rows_start.append(len(rows_index))

    # (3) y_{i,k} + x_{i-k} - x_i ≥ 0
    for (k, i_lo, i_hi, y_start) in shift_blocks:
        for off, i in enumerate(range(i_lo, i_hi + 1)):
            j = i - k
            y_idx = y_start + off
            # -x_i + x_j + y ≥ 0
            rows_index.extend([i - 1, j - 1, y_idx])
            rows_value.extend([-1.0, 1.0, 1.0])
            row_lower.append(0.0)
            row_upper.append(highspy.kHighsInf)
            rows_start.append(len(rows_index))

    # (4) For each k: t - Σ y_{i,k} ≥ 0
    for (k, i_lo, i_hi, y_start) in shift_blocks:
        idxs = list(range(y_start, y_start + (i_hi - i_lo + 1)))
        coefs = [-1.0] * len(idxs)
        idxs.append(N)  # t column
        coefs.append(1.0)
        rows_index.extend(idxs)
        rows_value.extend(coefs)
        row_lower.append(0.0)
        row_upper.append(highspy.kHighsInf)
        rows_start.append(len(rows_index))

    lp.num_row_ = len(row_lower)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
    lp.a_matrix_.start_ = np.array(rows_start, dtype=np.int32)
    lp.a_matrix_.index_ = np.array(rows_index, dtype=np.int32)
    lp.a_matrix_.value_ = np.array(rows_value, dtype=np.float64)
    lp.row_lower_ = np.array(row_lower)
    lp.row_upper_ = np.array(row_upper)

    t0 = time.time()
    h.passModel(lp)
    h.run()
    wall = time.time() - t0

    model_status = h.getModelStatus()
    info = h.getInfo()
    sol = h.getSolution()
    obj_val = info.objective_function_value
    # MIP gap: HiGHS exposes mip_gap (relative) and mip_node_count, mip_dual_bound.
    mip_dual_bound = info.mip_dual_bound
    mip_gap = info.mip_gap
    optimal = (model_status == highspy.HighsModelStatus.kOptimal)

    M_value = int(round(obj_val)) if obj_val is not None else None
    A_star = None
    if sol.col_value and len(sol.col_value) >= N:
        xv = sol.col_value[:N]
        A_star = tuple(i + 1 for i, v in enumerate(xv) if v > 0.5)

    # Verification
    verified_overlap = None
    if A_star is not None:
        verified_overlap = _max_overlap(A_star, N)

    h.clear()

    return {
        "n": n,
        "N": N,
        "M": M_value,
        "A_star": list(A_star) if A_star is not None else None,
        "verified_overlap": verified_overlap,
        "optimal": optimal,
        "mip_gap": mip_gap,
        "mip_dual_bound": mip_dual_bound,
        "wall_time_sec": wall,
        "ratio": M_value / n if M_value is not None else None,
        "solver": "HiGHS",
        "status": str(model_status),
    }


if __name__ == "__main__":
    # Quick self-test on small known cases.
    for n, expected in [(5, 3), (8, 4), (10, 5), (12, 5), (15, 6)]:
        r = solve_Mn_ilp(n, time_limit=120, verbose=False)
        ok = r["M"] == expected and r["verified_overlap"] == expected
        print(
            f"n={n}: M={r['M']} (expect {expected}) verified={r['verified_overlap']} "
            f"t={r['wall_time_sec']:.2f}s opt={r['optimal']} gap={r['mip_gap']} {'OK' if ok else 'FAIL'}"
        )
