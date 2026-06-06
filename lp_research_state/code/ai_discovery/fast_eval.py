"""Fast SDP evaluator at small N for constraint-discovery screening.

Solves White's SDP at small (N, T, R, bochner_n) in <2s, with optional
extra constraint family added. Returns Omega = optimal SDP value.

The ΔΩ between baseline and a candidate constraint at fast-eval scale
is a screening signal: ΔΩ ≥ 5×10⁻⁵ at N=200 typically corresponds to
ΔΩ ≥ 1×10⁻⁴ at Phase 5 (N=10000) scale, based on empirical scaling.

Usage:
    from fast_eval import baseline_solve, solve_with_extra
    base = baseline_solve()
    delta = solve_with_extra(my_constraint_fn) - base
"""
from __future__ import annotations
import sys
from pathlib import Path
import cvxpy as cp
import numpy as np
import warnings

CODE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE))
from white_full_convex import build_problem  # noqa: E402


# Fast-eval defaults: row 4, small N, fast solve. Empirically:
#   - Solve time: ~0.1s
#   - Status: optimal (CLARABEL converges cleanly at small N)
#   - bochner_n=6 gives ΔΩ ~ +1e-5 over bn=2, so this resolves >= 1e-5
FAST_CONFIG = {
    "N": 200,
    "T": 100,
    "R": 5,
    "h1": 0.004,
    "h2": 0.004,
    "p1": 0.3875,
    "p2": 0.3875,
    "q1": -0.02,
    "q2": 0.02,
    "bochner_n": 6,
}


def build_baseline(**overrides):
    """Build the baseline SDP at fast-eval scale. Returns (Omega, w, v, c, d, eps, dlt, cons)."""
    cfg = {**FAST_CONFIG, **overrides}
    return build_problem(
        cfg["N"], cfg["T"], cfg["R"],
        cfg["h1"], cfg["h2"], cfg["p1"], cfg["p2"], cfg["q1"], cfg["q2"],
        bochner_n=cfg["bochner_n"],
    )


def baseline_solve(**overrides):
    """Solve baseline. Returns Omega value (float)."""
    Omega, *_, cons = build_baseline(**overrides)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        prob.solve(solver=cp.CLARABEL, verbose=False)
    return float(prob.value), prob.status


def solve_with_extra(extra_constraint_fn, **overrides):
    """Solve with one or more extra constraints.

    extra_constraint_fn(Omega, w, v, c, d, eps, dlt, cfg) -> list of cvxpy constraints

    Returns: (Omega_value, status)
    """
    cfg = {**FAST_CONFIG, **overrides}
    Omega, w, v, c, d, eps, dlt, cons = build_baseline(**overrides)
    extra = extra_constraint_fn(Omega, w, v, c, d, eps, dlt, cfg)
    if not isinstance(extra, list):
        extra = [extra]
    cons_all = cons + extra
    prob = cp.Problem(cp.Minimize(Omega), cons_all)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            prob.solve(solver=cp.CLARABEL, verbose=False)
        except cp.SolverError as e:
            return None, f"solver_error: {e}"
    return float(prob.value) if prob.value is not None else None, prob.status


def measure_delta(extra_constraint_fn, **overrides):
    """Return (Omega_baseline, Omega_with_extra, delta, status).

    delta > 0 means the new constraint TIGHTENS the SDP (pushes LB UP), which
    is what we want for lever discovery.
    """
    base, base_status = baseline_solve(**overrides)
    val, val_status = solve_with_extra(extra_constraint_fn, **overrides)
    if val is None:
        return base, None, None, f"{base_status} / {val_status}"
    return base, val, val - base, f"{base_status} / {val_status}"


# Cache baseline value (it's deterministic given config)
_baseline_cache: dict = {}


def baseline_cached(**overrides):
    key = tuple(sorted({**FAST_CONFIG, **overrides}.items()))
    if key not in _baseline_cache:
        _baseline_cache[key] = baseline_solve(**overrides)
    return _baseline_cache[key]


if __name__ == "__main__":
    import time
    print("Smoke test: baseline solve at FAST_CONFIG")
    t0 = time.time()
    val, status = baseline_solve()
    dt = time.time() - t0
    print(f"  baseline: Omega = {val:.7f}, status = {status}, time = {dt:.3f}s")

    # Sanity: adding `Omega >= val + 1e-3` should make problem infeasible
    def infeasible(Omega, w, v, c, d, eps, dlt, cfg):
        return [Omega >= val + 1e-3]

    val2, status2 = solve_with_extra(infeasible)
    print(f"  with bad constraint: Omega = {val2}, status = {status2}")

    # Adding a vacuous tightening: |c[0]| <= 100 (already implied) should give same value
    def vacuous(Omega, w, v, c, d, eps, dlt, cfg):
        return [cp.abs(c[0]) <= 100]

    val3, status3 = solve_with_extra(vacuous)
    print(f"  with vacuous extra: Omega = {val3:.7f}, status = {status3}, ΔΩ = {val3 - val:+.2e}")
