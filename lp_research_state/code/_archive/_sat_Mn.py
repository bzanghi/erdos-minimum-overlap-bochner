"""SAT-based exact M(n) solver for the Erdős minimum overlap problem.

M(n) = min over equipartitions A ⊔ B = [2n], |A|=|B|=n, of max_k |A ∩ (B+k)|.

Encoding
--------
Variables: for each i ∈ [2n], a boolean x_i = (i ∈ A). The complement is B.
Exactly-n constraint on x_1..x_{2n} (CardEnc.equals).
Symmetry: pin x_1 = 1 (A↔B complement symmetry).

For each shift k ∈ [1, 2n-1] (positive) and k ∈ [-(2n-1), -1] (negative):
  overlap_k = |{i : i, i-k ∈ [2n], i ∈ A, i-k ∈ B}|
            = Σ_{i: i, i-k ∈ [2n]} x_i ∧ ¬x_{i-k}
We introduce auxiliary y_i^k ↔ x_i ∧ ¬x_{i-k} via 3 clauses each, then assert
Σ y_i^k ≤ M (CardEnc.atmost, seqcounter).

Note: overlap_k is NOT symmetric in k for finite [2n]. E.g. A={1..n}, B={n+1..2n}
gives overlap_{-n} = n, overlap_{+n} = 0. So we must encode both directions.
For k > 0: i ∈ [k+1, 2n] (so that i-k ∈ [1, 2n-k]).
For k < 0: i ∈ [1, 2n+k] (so that i-k ∈ [1-k, 2n]). Let j = -k > 0; then
  overlap_{-j} = |{i : i ∈ [1, 2n-j], i ∈ A, i+j ∈ B}|.

Solver: Glucose4 / Cadical via pysat. Binary search on M.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from pysat.solvers import Glucose4, Cadical153
from pysat.card import CardEnc, EncType


def encode_Mn_at_most(n: int, M: int):
    """Build SAT formula: ∃ A ⊂ [2n] with |A|=n and max overlap ≤ M.

    Returns: (clauses, n_vars).
    """
    N = 2 * n
    next_var = N + 1  # auxiliary vars start here

    clauses = []

    # Cardinality: exactly n of x_1..x_N are true
    card_cnf = CardEnc.equals(
        lits=list(range(1, N + 1)), bound=n,
        top_id=next_var - 1, encoding=EncType.seqcounter,
    )
    clauses.extend(card_cnf.clauses)
    if card_cnf.nv >= next_var:
        next_var = card_cnf.nv + 1

    # Pin x_1 = 1 (break A↔B symmetry)
    clauses.append([1])

    # For each nonzero shift k in [-(N-1), N-1]
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        # i ranges so that both i and i-k are in [1, N]
        # i ∈ [max(1, 1+k), min(N, N+k)]
        i_lo = max(1, 1 + k)
        i_hi = min(N, N + k)
        if i_lo > i_hi:
            continue
        y_lits = []
        for i in range(i_lo, i_hi + 1):
            j = i - k  # j ∈ [1, N]
            y = next_var
            next_var += 1
            # y ↔ x_i ∧ ¬x_j
            clauses.append([-y, i])     # y → x_i
            clauses.append([-y, -j])    # y → ¬x_j
            clauses.append([y, -i, j])  # x_i ∧ ¬x_j → y
            y_lits.append(y)
        if M < len(y_lits):
            atmost = CardEnc.atmost(
                lits=y_lits, bound=M,
                top_id=next_var - 1, encoding=EncType.seqcounter,
            )
            clauses.extend(atmost.clauses)
            if atmost.nv >= next_var:
                next_var = atmost.nv + 1
    return clauses, next_var - 1


def verify_overlap(A: tuple, n: int) -> int:
    """Compute max overlap directly from A subset of [2n]."""
    N = 2 * n
    A_set = set(A)
    B_set = set(range(1, N + 1)) - A_set
    best = 0
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        count = sum(1 for i in A_set if (i - k) in B_set)
        if count > best:
            best = count
    return best


def solve_Mn_at_most(n: int, M: int, solver_cls=Cadical153, time_budget_sec: float = 600):
    """Single SAT call: is there A with max overlap ≤ M? Returns (sat, A or None, elapsed)."""
    t0 = time.time()
    clauses, _ = encode_Mn_at_most(n, M)
    enc_time = time.time() - t0
    N = 2 * n
    with solver_cls(bootstrap_with=clauses) as s:
        sat = s.solve()
        if sat:
            model = s.get_model()
            A = tuple(i for i in range(1, N + 1) if model[i - 1] > 0)
        else:
            A = None
    elapsed = time.time() - t0
    return sat, A, elapsed, enc_time


def solve_Mn_sat(n: int, time_budget_sec: float = 600, M_lo_init: int = None, M_hi_init: int = None):
    """Binary-search M(n). Returns dict with M, A_star, total_time, per-call times."""
    if M_lo_init is None:
        M_lo = max(1, int(0.36 * n))
    else:
        M_lo = M_lo_init
    if M_hi_init is None:
        M_hi = int(0.5 * n) + 2
    else:
        M_hi = M_hi_init
    M_best = None
    A_star = None
    t_start = time.time()
    calls = []

    while M_lo <= M_hi:
        if time.time() - t_start > time_budget_sec:
            break
        M = (M_lo + M_hi) // 2
        remaining = time_budget_sec - (time.time() - t_start)
        sat, A, elapsed, enc = solve_Mn_at_most(n, M, time_budget_sec=remaining)
        calls.append({"M": M, "sat": sat, "elapsed": elapsed, "enc_time": enc})
        print(f"    n={n} M={M}: {'SAT' if sat else 'UNSAT'} in {elapsed:.2f}s (enc {enc:.2f}s)", flush=True)
        if sat:
            # Verify
            obs = verify_overlap(A, n)
            assert obs <= M, f"Encoding bug: requested ≤{M}, got actual {obs}"
            M_best = M
            A_star = A
            M_hi = M - 1
        else:
            M_lo = M + 1

    total = time.time() - t_start
    return {
        "n": n,
        "M": M_best,
        "A_star": list(A_star) if A_star else None,
        "ratio": M_best / n if M_best is not None else None,
        "total_time": total,
        "calls": calls,
    }


def main():
    out_path = Path(__file__).parent.parent / "data" / "Mn_sat_results.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Known values (sanity checks)
    known_M = {
        13: 6, 14: 6, 15: 6, 16: 7, 17: 7, 18: 8,
    }

    results = []

    # Quick sanity: n=10, 12 (small)
    print("=== Sanity checks ===", flush=True)
    for n_sanity, expected in [(10, 5), (12, 5), (15, 6)]:
        r = solve_Mn_sat(n_sanity, time_budget_sec=120)
        ok = (r["M"] == expected)
        print(f"  n={n_sanity}: M={r['M']} (expected {expected}) ok={ok} t={r['total_time']:.2f}s", flush=True)
        if not ok:
            print("ABORT: sanity check failed.")
            return
        # Cross-check overlap value
        obs = verify_overlap(tuple(r["A_star"]), n_sanity)
        print(f"           verified max overlap = {obs}", flush=True)
        results.append({**r, "verified_overlap": obs, "expected_M": expected})

    # Extension run — start at known n and push
    total_budget = 25 * 60.0  # 25 minutes for the extension push
    start_extension = time.time()
    extension_ns = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

    for n in extension_ns:
        if time.time() - start_extension > total_budget:
            print(f"Total extension budget exhausted, stopping.", flush=True)
            break

        # Tight bracket based on prior n (M nondecreasing, roughly M(n) ≈ 0.4..0.5 n)
        if results:
            prior_Ms = [r["M"] for r in results if r["M"] is not None]
            M_lo_init = max(1, max(prior_Ms) - 1) if prior_Ms else max(1, int(0.36 * n))
        else:
            M_lo_init = max(1, int(0.36 * n))
        M_hi_init = int(0.55 * n) + 2

        budget_this = min(15 * 60.0, total_budget - (time.time() - start_extension))
        print(f"\n--- n={n} (bracket [{M_lo_init},{M_hi_init}], budget {budget_this:.0f}s) ---", flush=True)
        r = solve_Mn_sat(n, time_budget_sec=budget_this,
                         M_lo_init=M_lo_init, M_hi_init=M_hi_init)
        expected = known_M.get(n)
        match = (expected is None) or (r["M"] == expected)
        print(f"  -> M({n}) = {r['M']}  ratio = {r['ratio']:.6f}  t = {r['total_time']:.1f}s  expected={expected} match={match}", flush=True)
        if r["A_star"]:
            obs = verify_overlap(tuple(r["A_star"]), n)
            r["verified_overlap"] = obs
            print(f"     verified overlap = {obs}, A* = {r['A_star']}", flush=True)
        r["expected_M"] = expected
        r["match_known"] = match
        results.append(r)

        with open(out_path, "w") as f:
            json.dump(results, f, indent=2)

        # Stop if a single solve took too long
        if r["total_time"] > 10 * 60.0:
            print(f"  Single n took >10 min, stopping.", flush=True)
            break
        if r["M"] is None:
            print(f"  Did not finish binary search within budget, stopping.", flush=True)
            break

    print(f"\nSaved {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
