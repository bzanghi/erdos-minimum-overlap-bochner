"""
Phase 1.1 — Brute-force M(n) for the Erdős minimum overlap problem.
Optimized version using numpy correlation.
"""
from itertools import combinations
import numpy as np
import time, json

def M_partition(a_vec, b_vec):
    """a_vec, b_vec are 0/1 numpy arrays of length 2n.
       M_k = sum_x a[x] b[x-k] = (a ⋆ b)[k]; we just need max of correlation."""
    # np.correlate(a, b, 'full') has length 4n-1, M_k for k=-(2n-1)..(2n-1)
    return int(np.correlate(a_vec, b_vec, mode='full').max())

def M_n_brute(n):
    N = 2 * n
    best = N  # trivial upper bound
    best_A = None
    base_A = np.zeros(N, dtype=np.int8)
    base_A[0] = 1  # 1 ∈ A (index 0 corresponds to integer 1)
    for A_rest in combinations(range(1, N), n - 1):
        A_vec = base_A.copy()
        for i in A_rest:
            A_vec[i] = 1
        B_vec = 1 - A_vec
        m = M_partition(A_vec, B_vec)
        if m < best:
            best = m
            best_A = (1,) + tuple(i + 1 for i in A_rest)
    return best, best_A

if __name__ == "__main__":
    table = []
    for n in range(2, 13):  # cap at n=12; C(23,11)~1.35M is manageable
        t0 = time.time()
        m, A = M_n_brute(n)
        dt = time.time() - t0
        ratio = m / n
        print(f"n={n:2d}  M(n)={m}  M(n)/n={ratio:.5f}   t={dt:6.2f}s   A*={A}")
        table.append({"n": n, "M": m, "ratio": ratio, "A_star": list(A)})
    with open("Mn_brute.json", "w") as f:
        json.dump(table, f, indent=2)
    print("saved Mn_brute.json")
