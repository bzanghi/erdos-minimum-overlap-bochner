"""SAT decision procedure for exact discrete M(n) — Erdős minimum overlap.

Independent re-implementation of the SAME problem solved by `discrete_M_ilp.py`
(McCormick-linearized HiGHS ILP). The two encodings must agree (project's
"independent re-implementation to 10+ digits" convention; here the value is an
integer so agreement is exact).

Convention (identical to discrete_M_ilp.py and known_Mn_values.json)
-------------------------------------------------------------------
    Full set [1..N], N = 2n. A ⊆ [N], |A| = n, B = [N] \\ A.
    overlap_k(A) = |{ i : i ∈ A, i-k ∈ B, i ∈ [1,N], i-k ∈ [1,N] }|.
    M(n)         = min_A  max_{k ∈ Z\\{0}}  overlap_k(A).
    µ            = lim M(n)/n.

SAT decision procedure
----------------------
M(n) is an integer.  Define the feasibility predicate

    FEAS(v)  ≡  "∃ A ⊆ [N], |A| = n, with max_k overlap_k(A) ≤ v".

FEAS is monotone in v (feasible at v ⇒ feasible at v+1).  So

    M(n) = min { v : FEAS(v) holds }.

We find it by binary search on v.  A solved n yields a CERTIFIED M(n) iff we
have BOTH a SAT certificate at v = M(n) (a witness A) AND an UNSAT certificate
at v = M(n) − 1.  We record the existence of both.

SAT encoding of FEAS(v)
-----------------------
Booleans x_1..x_N  (x_i ⇔ i ∈ A).
  * Cardinality  Σ x_i = n         (CardEnc.equals, totalizer).
  * Symmetry break 1 (complement A↔B): pin x_1 = 1.
        overlap_k(B) = overlap_{-k}(A) exactly, so M is invariant under A↔B and
        we may fix 1 ∈ A WLOG.
  * Symmetry break 2 (reflection r: i ↦ N+1−i): lex-leader  x ≤_lex x∘r.
        overlap_k(r(A)) = overlap_k(A) (reversal preserves the multiset of
        gaps), so M is invariant under reflection; we keep the lexicographically
        smaller of {A, r(A)}.  Encoded with a standard lex chain.
  * For each nonzero shift k with nonempty overlap domain I_k:
        aux y_{i,k} ⇔ x_i ∧ ¬x_{i-k}   (3 clauses each),
        then  Σ_{i∈I_k} y_{i,k} ≤ v    (CardEnc.atmost, totalizer).

Solver: Cadical195 (strongest bundled CDCL) by default; kissat404 selectable.
"""
from __future__ import annotations

import time
from typing import Optional

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool


_SOLVERS = {}


def _get_solver_cls(name: str):
    """Lazy import of a PySAT solver class by short name."""
    if name in _SOLVERS:
        return _SOLVERS[name]
    import pysat.solvers as S
    table = {
        "cadical195": S.Cadical195,
        "cadical153": S.Cadical153,
        "cadical103": S.Cadical103,
        "glucose4": S.Glucose4,
        "glucose42": S.Glucose42,
        "kissat404": getattr(S, "Kissat404", None),
        "minisat22": S.Minisat22,
    }
    cls = table.get(name)
    if cls is None:
        raise ValueError(f"unknown/unavailable solver {name!r}; have {list(table)}")
    _SOLVERS[name] = cls
    return cls


# ---------------------------------------------------------------------------
# Direct verification (no SAT) — the ground truth the encoding must match.
# ---------------------------------------------------------------------------
def max_overlap(A, N: int) -> int:
    """max over nonzero k of overlap_k(A), A ⊆ [1,N]."""
    Aset = set(A)
    Bset = set(range(1, N + 1)) - Aset
    best = 0
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        c = sum(1 for i in Aset if (i - k) in Bset)
        if c > best:
            best = c
    return best


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------
def _lex_le_clauses(a, b, pool: IDPool):
    """CNF for the lexicographic constraint  a ≤_lex b, a/b equal-length lists
    of literals (positive var ids).  Standard chain with equality-prefix aux:

        e_0 = True
        for t: ( e_t ∧ a_t ) → b_t           [if prefix equal, a_t ⇒ b_t]
               e_{t+1} ↔ e_t ∧ (a_t ↔ b_t)

    We only need the forward implication e_{t+1} → (e_t ∧ a_t=b_t) for
    correctness of the "≤" direction together with the per-position rule.
    Implemented via Tseitin on g_t = (a_t ↔ b_t).
    """
    clauses = []
    assert len(a) == len(b)
    m = len(a)
    # e_t : prefix [0..t-1] equal.  e_0 ≡ True (no aux var; treat as constant).
    e_prev = None  # None encodes the constant-True e_0
    for t in range(m):
        at, bt = a[t], b[t]
        # rule:  e_t ∧ a_t → b_t
        if e_prev is None:
            clauses.append([-at, bt])           # a_t → b_t  (e_0=True)
        else:
            clauses.append([-e_prev, -at, bt])  # e_t ∧ a_t → b_t
        if t == m - 1:
            break
        # g_t ↔ (a_t ↔ b_t)
        g = pool.id(("g", id(a), t))
        # g → (a_t→b_t) and (b_t→a_t):
        clauses.append([-g, -at, bt])
        clauses.append([-g, -bt, at])
        # (a_t↔b_t) → g  :  two clauses
        clauses.append([g, at, bt])
        clauses.append([g, -at, -bt])
        # e_{t+1} ↔ e_t ∧ g
        e_next = pool.id(("e", id(a), t + 1))
        if e_prev is None:
            # e_{t+1} ↔ g
            clauses.append([-e_next, g])
            clauses.append([e_next, -g])
        else:
            clauses.append([-e_next, e_prev])
            clauses.append([-e_next, g])
            clauses.append([e_next, -e_prev, -g])
        e_prev = e_next
    return clauses


def build_feasibility_cnf(
    n: int,
    v: int,
    use_reflection_sym: bool = True,
    card_enc: int = EncType.totalizer,
):
    """CNF for FEAS(v):  ∃ A⊆[2n], |A|=n, max overlap ≤ v.

    Returns (clauses, pool) where x_i has var id i (1..N).
    """
    N = 2 * n
    pool = IDPool(start_from=N + 1)  # x_1..x_N occupy 1..N
    clauses = []

    # |A| = n
    card = CardEnc.equals(
        lits=list(range(1, N + 1)), bound=n,
        vpool=pool, encoding=card_enc,
    )
    clauses.extend(card.clauses)

    # symmetry break 1: pin x_1 = 1
    clauses.append([1])

    # symmetry break 2: x ≤_lex reflect(x), reflect: i ↦ N+1-i
    if use_reflection_sym:
        a = list(range(1, N + 1))                 # x_1..x_N
        b = [N + 1 - i for i in range(1, N + 1)]  # x_N..x_1
        clauses.extend(_lex_le_clauses(a, b, pool))

    # overlap constraints
    for k in range(-(N - 1), N):
        if k == 0:
            continue
        i_lo = max(1, 1 + k)
        i_hi = min(N, N + k)
        if i_lo > i_hi:
            continue
        y_lits = []
        for i in range(i_lo, i_hi + 1):
            j = i - k
            y = pool.id(("y", k, i))
            # y ⇔ x_i ∧ ¬x_j
            clauses.append([-y, i])      # y → x_i
            clauses.append([-y, -j])     # y → ¬x_j
            clauses.append([y, -i, j])   # x_i ∧ ¬x_j → y
            y_lits.append(y)
        ndom = len(y_lits)
        if v < ndom:  # else constraint is vacuous
            atmost = CardEnc.atmost(
                lits=y_lits, bound=v, vpool=pool, encoding=card_enc,
            )
            clauses.extend(atmost.clauses)

    return clauses, pool


def feasible_at(
    n: int, v: int, solver: str = "cadical195",
    use_reflection_sym: bool = True, card_enc: int = EncType.totalizer,
):
    """Single decision: is FEAS(v) satisfiable?

    Returns (sat: bool, A: tuple|None, elapsed: float, enc_time: float).
    """
    t0 = time.time()
    clauses, _ = build_feasibility_cnf(
        n, v, use_reflection_sym=use_reflection_sym, card_enc=card_enc)
    enc_time = time.time() - t0
    N = 2 * n
    cls = _get_solver_cls(solver)
    s = cls(bootstrap_with=clauses)
    try:
        sat = s.solve()
        A = None
        if sat:
            model = s.get_model()
            mset = set(model)
            A = tuple(i for i in range(1, N + 1) if i in mset)
    finally:
        s.delete()
    return sat, A, time.time() - t0, enc_time


def solve_Mn(
    n: int,
    M_lo: Optional[int] = None,
    M_hi: Optional[int] = None,
    solver: str = "cadical195",
    time_budget_sec: float = 1800.0,
    use_reflection_sym: bool = True,
    card_enc: int = EncType.totalizer,
    verbose: bool = True,
) -> dict:
    """Binary-search M(n) with a per-n wall-clock budget.

    A returned result has  optimal=True  iff we obtained a SAT witness at
    M = M_value AND an UNSAT proof at M_value − 1 within budget.  Otherwise:
      * if we have a witness but no matching UNSAT  →  M_value is an UPPER bound
        on M(n) (feasible-only), optimal=False.
      * if budget hit mid-search with no witness    →  M_value None.
    """
    N = 2 * n
    if M_lo is None:
        M_lo = max(0, int(0.30 * n))      # safely below any M(n)
    if M_hi is None:
        M_hi = n                           # |A|=n ⇒ overlap ≤ n trivially
    lo, hi = M_lo, M_hi
    best_M = None
    best_A = None
    unsat_at = {}   # v -> True if proven UNSAT
    sat_at = {}     # v -> True if proven SAT (feasible)
    calls = []
    t_start = time.time()

    def log(msg):
        if verbose:
            print(msg, flush=True)

    log(f"n={n}: binary search M in [{lo},{hi}], solver={solver}")
    while lo <= hi:
        if time.time() - t_start > time_budget_sec:
            log(f"  budget exhausted ({time_budget_sec}s)")
            break
        v = (lo + hi) // 2
        sat, A, elapsed, enc = feasible_at(
            n, v, solver=solver,
            use_reflection_sym=use_reflection_sym, card_enc=card_enc)
        calls.append({"v": v, "sat": bool(sat), "elapsed": elapsed, "enc": enc})
        if sat:
            obs = max_overlap(A, N)
            assert obs <= v, f"ENCODING BUG n={n}: asked ≤{v}, witness has {obs}"
            sat_at[v] = True
            best_M, best_A = v, A
            hi = v - 1
            log(f"  v={v}: SAT  (witness overlap={obs}) {elapsed:.2f}s")
        else:
            unsat_at[v] = True
            lo = v + 1
            log(f"  v={v}: UNSAT  {elapsed:.2f}s")

    total = time.time() - t_start
    optimal = (
        best_M is not None
        and best_M in sat_at
        and (best_M - 1 in unsat_at or best_M == 0)
    )
    return {
        "n": n,
        "N": N,
        "M": best_M,
        "A_star": list(best_A) if best_A else None,
        "verified_overlap": max_overlap(best_A, N) if best_A else None,
        "ratio": best_M / n if best_M is not None else None,
        "optimal": optimal,
        "sat_certificate_at_M": (best_M in sat_at) if best_M is not None else False,
        "unsat_certificate_at_M_minus_1": (
            (best_M - 1) in unsat_at if best_M is not None else False),
        "wall_time_sec": total,
        "solver": solver,
        "calls": calls,
    }


if __name__ == "__main__":
    # Self-test against published Haugland values (including n=18 the ILP missed).
    KNOWN = {10: 5, 12: 5, 15: 6, 16: 7, 17: 7, 18: 8}
    for n in sorted(KNOWN):
        r = solve_Mn(n, time_budget_sec=300, verbose=False)
        exp = KNOWN[n]
        ok = r["M"] == exp and r["optimal"]
        print(f"n={n}: M={r['M']} (expect {exp}) opt={r['optimal']} "
              f"t={r['wall_time_sec']:.2f}s {'OK' if ok else 'FAIL'}")
