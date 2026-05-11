"""
Phase 1.2 — Brute-force M(n) for n in {13..16} via branch-and-bound +
symmetry reduction.  Companion to brute_force_Mn.py (n ≤ 12 exhaustive).

Problem
-------
A ⊔ B = {1, ..., 2n}, |A| = |B| = n.  For k ∈ Z,
   M_k(A, B) = |{(a, b) ∈ A × B : a - b = k}|.
M(n) = min over A of max over k of M_k(A, B).

Symmetries used
---------------
1. A ↔ B (complement): both achieve the same M.  Fix 1 ∈ A.
2. Reflection r(x) = 2n + 1 - x.  M(r(A), r(B)) = M(A, B); reduces space.
   We canonicalize by requiring A ≤ r(A) lexicographically.
3. Branch-and-bound pruning: build A element-by-element (after fixing 1).
   At each partial state, compute a LOWER bound on the eventual max-corr
   that ANY completion can achieve, and prune if ≥ current best.

Lower bound used in B&B
-----------------------
Fix the partial A ⊂ {1..2n}.  Complete to a full A by adding elements.
The eventual A and B = {1..2n} \ A induce M_k = (1_A ⋆ 1_B)[k].
For any k ≠ 0, we have
   M_k = Σ_{i ∈ A} 1_{i + k ∈ B}   (treating A, B as subsets of {1..2n})
       = |{i ∈ A : i + k ∈ B}|.

Lower bound at partial state: for each candidate shift k, the
*minimum possible* contribution from already-decided elements is the
count of fixed (a, b) pairs where a ∈ A_fixed, b = a + k ∈ B_fixed.
Future additions can only increase M_k.  If this minimum ≥ current_best,
prune.

Empirically: at n=15, naive enumeration is C(29, 14) ≈ 78M; with
reflection + B&B pruning, runtime is a few minutes single-threaded.

Output
------
Appends to lp_research_state/data/Mn_optimizers_large.json:
   list of {n, M, A_star, num_blocks (in 2n-grid), time_seconds}.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
import numpy as np


def num_blocks(A_vec: np.ndarray) -> int:
    """Number of constant runs in 0/1 vector (treated as cells, no wrap)."""
    if len(A_vec) == 0:
        return 0
    diffs = np.diff(A_vec.astype(np.int8))
    return int(1 + np.count_nonzero(diffs))


def M_full(A_vec: np.ndarray) -> int:
    B_vec = 1 - A_vec
    return int(np.correlate(A_vec, B_vec, mode='full').max())


def reflect(A_set, N):
    return frozenset(N + 1 - x for x in A_set)


def is_lex_le(A_set, A_ref_set):
    """A_set ≤ A_ref_set lexicographically as sorted tuples?"""
    a = sorted(A_set)
    b = sorted(A_ref_set)
    return a <= b


def correlation_full(A_vec: np.ndarray):
    """Return full correlation vector."""
    B_vec = 1 - A_vec
    return np.correlate(A_vec, B_vec, mode='full')


def M_n_bb(n: int, time_budget_s: float = 600.0):
    """Branch-and-bound search for M(n) with 1 ∈ A and reflection
    canonicalization."""
    N = 2 * n
    # Domain: positions 1..N.  1 ∈ A.  We choose remaining n-1 from {2..N}.
    best = N + 1
    best_A = None

    # Stronger upper bound from known sequence (Haugland-style).  Optional.
    # We just initialise best = n (since trivially M ≤ n always; in fact
    # known M(n) grows roughly as 0.4n).
    best = n  # initial UB

    start = time.time()
    timed_out = False

    # State: chosen A elements (as sorted list), candidate positions.
    # We do DFS: choose remaining elements in increasing order.
    # That naturally enumerates each A set once.
    A_init = (1,)

    # Precompute correlation contributions of pos as we add/remove.
    # We track A_vec and B_partial.  But B is "everything not in A yet decided".
    # During DFS, "undecided" positions are those > current cursor.
    # B_decided are positions < cursor not in A.

    # M_k lower bound: count fixed A-elements a such that a+k is fixed in B.
    # Future adds to A reduce B; the contribution of an (a, a+k) pair where
    # a is in A_partial but a+k is currently undecided could go either way.
    # We instead lower-bound: only count CERTAIN B (fixed B).
    #
    # We'll proceed via DFS scanning positions 2..N, deciding A or B for each.
    # State: pos = next position; A_count, B_count; arrays.

    # Implement as recursion with explicit stack for speed.
    # For clarity start with recursion.

    A_arr = np.zeros(N, dtype=np.int8)
    A_arr[0] = 1
    B_arr = np.zeros(N, dtype=np.int8)  # 1 if position is committed to B

    # For pruning we need a lower bound on max-correlation given current
    # commitments. Use the partial correlation between A_arr and B_arr:
    # for each shift k, M_k >= Σ_i A_arr[i] * B_arr[i - k] (already fixed).
    # We track running partial corr as we add positions.

    partial_corr = np.zeros(2 * N - 1, dtype=np.int32)
    # np.correlate(A, B, full) has output[k] = Σ_i A[i] B[i - lag]
    # where lag = k - (n_B - 1).  Let's just maintain it via incremental.

    # When we set A[p] = 1, partial_corr[k] += B[p - lag_from_k]; equivalently
    # for each j with B[j]=1, partial_corr[k = (p) - j + (N-1)] += 1.
    # Similarly when we set B[p] = 1.

    nA_target = n
    nA_current = 1  # position 1 in A
    nB_current = 0

    nonlocal_best = [best]
    nonlocal_bestA = [None]

    # For reflection canonicalization: we want A ≤ reflect(A) lexicographically.
    # Reflection sends position p -> N+1-p.  So sorted(reflect(A)) is
    # [N+1-p for p in sorted(A, reverse=True)].
    # We'll check at completion (or eagerly prune when we know A > reflect(A)).

    def recurse(pos):
        nonlocal nA_current, nB_current
        if time.time() - start > time_budget_s:
            return False  # time out signal

        remaining = N - pos + 1  # positions pos..N still to decide
        needed_A = nA_target - nA_current
        needed_B = nA_target - nB_current  # |B| = n too
        if needed_A < 0 or needed_B < 0:
            return True
        if needed_A > remaining or needed_B > remaining:
            return True
        if pos > N:
            # Completed assignment.
            # Reflection canonicalization
            A_set = frozenset(np.where(A_arr == 1)[0] + 1)  # 1-indexed positions
            R = frozenset(N + 1 - x for x in A_set)
            if not is_lex_le(A_set, R):
                return True
            # Compute actual M
            m = M_full(A_arr)
            if m < nonlocal_best[0]:
                nonlocal_best[0] = m
                nonlocal_bestA[0] = tuple(sorted(A_set))
            return True

        # Pruning: check partial correlation lower bound
        # partial_corr already reflects committed (A, B) overlap contributions.
        # If max(partial_corr) >= best (strict), we cannot improve.
        if partial_corr.max() >= nonlocal_best[0]:
            return True

        # Branch order: try adding to A first if needed (since structure
        # of optimal solutions tends to be A-heavy at small positions).
        # Doesn't affect correctness.

        # Option 1: put pos in A
        if needed_A > 0:
            # Add pos to A.  Update partial_corr: for each j with B_arr[j]=1,
            # increment partial_corr at shift index pos-1 - j + (N - 1).
            # Simpler: vectorized.
            A_arr[pos - 1] = 1
            nA_current += 1
            # Find B positions and update partial_corr
            B_idx = np.where(B_arr == 1)[0]
            if len(B_idx) > 0:
                # shift index = (pos-1) - B_idx + (N - 1)
                shift_idx = (pos - 1) - B_idx + (N - 1)
                np.add.at(partial_corr, shift_idx, 1)
            cont = recurse(pos + 1)
            # Undo
            if len(B_idx) > 0:
                shift_idx = (pos - 1) - B_idx + (N - 1)
                np.add.at(partial_corr, shift_idx, -1)
            A_arr[pos - 1] = 0
            nA_current -= 1
            if not cont:
                return False

        # Option 2: put pos in B
        if needed_B > 0:
            B_arr[pos - 1] = 1
            nB_current += 1
            A_idx = np.where(A_arr == 1)[0]
            if len(A_idx) > 0:
                shift_idx = A_idx - (pos - 1) + (N - 1)
                np.add.at(partial_corr, shift_idx, 1)
            cont = recurse(pos + 1)
            if len(A_idx) > 0:
                shift_idx = A_idx - (pos - 1) + (N - 1)
                np.add.at(partial_corr, shift_idx, -1)
            B_arr[pos - 1] = 0
            nB_current -= 1
            if not cont:
                return False

        return True

    completed = recurse(2)
    elapsed = time.time() - start

    return {
        "n": n,
        "M": nonlocal_best[0],
        "A_star": list(nonlocal_bestA[0]) if nonlocal_bestA[0] else None,
        "time_seconds": elapsed,
        "timed_out": not completed,
    }


def main():
    out_path = Path(__file__).parent.parent / "data" / "Mn_optimizers_large.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    results = []

    # Sanity check: replicate n=10..12 to validate B&B against existing data
    sanity_ns = [10, 11, 12]
    expected_M = {10: 5, 11: 5, 12: 5}
    print("=== Sanity checks (replicate known M(n)) ===")
    for n in sanity_ns:
        t0 = time.time()
        r = M_n_bb(n, time_budget_s=300.0)
        dt = time.time() - t0
        N = 2 * n
        A_vec = np.zeros(N, dtype=np.int8)
        for p in r["A_star"]:
            A_vec[p - 1] = 1
        nb = num_blocks(A_vec)
        r["num_blocks"] = nb
        ok = (r["M"] == expected_M[n])
        print(f"n={n:2d}  M={r['M']}  expected={expected_M[n]}  ok={ok}  "
              f"A_star={r['A_star']}  blocks={nb}  t={dt:.2f}s  timed_out={r['timed_out']}")
        if not ok:
            print("ABORT: sanity check failed.")
            return
        results.append(r)

    # Extended runs
    print("\n=== Extended runs (n=13..16) ===")
    budgets = {13: 600.0, 14: 1200.0, 15: 1500.0, 16: 1800.0}
    cumulative_used = 0.0
    HARD_TOTAL_BUDGET = 1500.0  # 25 minutes total ceiling

    for n in [13, 14, 15, 16]:
        remaining = HARD_TOTAL_BUDGET - cumulative_used
        if remaining < 60.0:
            print(f"n={n}: insufficient remaining budget ({remaining:.0f}s), skipping.")
            break
        budget = min(budgets[n], remaining)
        print(f"n={n}: starting (budget {budget:.0f}s, used {cumulative_used:.0f}s)")
        t0 = time.time()
        r = M_n_bb(n, time_budget_s=budget)
        dt = time.time() - t0
        cumulative_used += dt
        N = 2 * n
        A_vec = np.zeros(N, dtype=np.int8)
        if r["A_star"]:
            for p in r["A_star"]:
                A_vec[p - 1] = 1
        nb = num_blocks(A_vec)
        r["num_blocks"] = nb
        print(f"  -> M={r['M']}  blocks={nb}  t={dt:.1f}s  timed_out={r['timed_out']}")
        print(f"     A*={r['A_star']}")
        results.append(r)
        if r["timed_out"]:
            print("  (timed out; result may be only an upper bound)")
            # Continue to next n only if we have a lot of budget; usually stop
            if cumulative_used > HARD_TOTAL_BUDGET * 0.7:
                print("  stopping further searches due to budget.")
                break

    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
