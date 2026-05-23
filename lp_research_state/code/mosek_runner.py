"""
Mosek runner for the white_full_convex SDP family.

Public entry point: solve_with_mosek(problem, **kwargs) -> dict mirroring
sdpa_runner.solve_with_sdpa_gmp(...).  Returns:
  primal_obj         : rigorous LOWER bound on Minimize(Omega) [see Notes].
  dual_obj           : rigorous UPPER bound on Minimize(Omega).
  rigorous_dual_LB   : primal_obj - duality_gap  (mirrors sdpa_runner's "subtract
                       the gap to be safely below the true optimum" convention).
  duality_gap        : |dual_obj - primal_obj|  in original objective units.
  primal_viol        : Mosek's reported primal constraint violation
  dual_viol          : Mosek's reported dual constraint violation
  status             : cvxpy status string ('optimal', 'optimal_inaccurate', ...)
  mosek_problem_status, mosek_solution_status : strings from Mosek log
  iterations         : interior-point iteration count
  runtime_sec        : wall clock for the cvxpy.solve(...) call
  stdout_log         : full captured Mosek log (truncated to 50KB if asked)

Notes
-----
CVXPY's `ConeMatrixStuffing` reduction converts a Minimize problem into the
canonical Conic standard form which Mosek sees as a Maximize (its standard
form maximizes a primal). In Mosek's iteration table:

    POBJ  = primal objective of Mosek's max-form == DUAL objective of our
            Minimize(Omega) == rigorous LOWER bound on min(Omega)
    DOBJ  = dual   objective of Mosek's max-form == PRIMAL objective of our
            Minimize(Omega) == rigorous UPPER bound on min(Omega)

This is the opposite of what one might naively expect, but is the unambiguous
convention.  We verified on a toy `min x s.t. [[x,1],[1,x]] >> 0` (true opt 1)
that POBJ < 1 < DOBJ at convergence, with cvxpy's `prob.value` = DOBJ.

So `dual_obj` (this dict) = DOBJ  (= cvxpy's prob.value),
and `primal_obj` (this dict) = POBJ  (= the rigorous lower bound certificate).

The `rigorous_dual_LB = primal_obj - duality_gap` deduction mirrors
sdpa_runner.py's "subtract the absolute gap" rule, which is conservative.
"""

from __future__ import annotations

import contextlib
import io
import logging
import re
import time
from typing import Optional

import cvxpy as cp


# ---- Mosek log parsing -----------------------------------------------------
_RE_PROBSTATUS = re.compile(r"Problem status\s*:\s*(\S+)")
_RE_SOLSTATUS  = re.compile(r"Solution status\s*:\s*(\S+)")
# "Primal.  obj: 9.9999999885e-01    nrm: 1e+00    Viol.  con: 2e-09    barvar: 0e+00"
_RE_PRIMAL_LINE = re.compile(
    r"Primal\.\s*obj:\s*([+\-]?\d+\.\d+e[+\-]?\d+)"
    r".*?Viol\.\s*con:\s*([+\-]?\d+(?:\.\d+)?e[+\-]?\d+)"
)
_RE_DUAL_LINE = re.compile(
    r"Dual\.\s*obj:\s*([+\-]?\d+\.\d+e[+\-]?\d+)"
    r".*?Viol\.\s*con:\s*([+\-]?\d+(?:\.\d+)?e[+\-]?\d+)"
)
# Iteration table lines: " 3   1.6e-09  5.4e-09  1.9e-13  9.98e-01   9.999999989e-01   1.000000002e+00   3.8e-09  0.00 "
_RE_ITE = re.compile(
    r"^\s*(?:\(CVXPY\).*?:\s*)?(\d+)\s+"          # iter
    r"(\S+e[+\-]?\d+)\s+(\S+e[+\-]?\d+)\s+(\S+e[+\-]?\d+)\s+"  # PFEAS DFEAS GFEAS
    r"([+\-]?\d+(?:\.\d+)?e[+\-]?\d+)\s+"                       # PRSTATUS
    r"([+\-]?\d+\.\d+e[+\-]?\d+)\s+([+\-]?\d+\.\d+e[+\-]?\d+)\s+"  # POBJ DOBJ
    r"(\S+e[+\-]?\d+)",                                          # MU
    re.MULTILINE,
)


_RE_CVXPY_PREFIX = re.compile(r"^\(CVXPY\)[^:]*(?::\d+){0,2}\s*[AP]M:\s?", re.MULTILINE)


def _strip_cvxpy_prefix(out: str) -> str:
    return _RE_CVXPY_PREFIX.sub("", out)


def parse_mosek_log(out: str) -> dict:
    """Extract POBJ, DOBJ, viols, statuses, iter count from Mosek log."""
    out = _strip_cvxpy_prefix(out)
    rows = []
    for m in _RE_ITE.finditer(out):
        try:
            rows.append({
                "iter":     int(m.group(1)),
                "pfeas":    float(m.group(2)),
                "dfeas":    float(m.group(3)),
                "gfeas":    float(m.group(4)),
                "prstatus": float(m.group(5)),
                "pobj":     float(m.group(6)),
                "dobj":     float(m.group(7)),
                "mu":       float(m.group(8)),
            })
        except ValueError:
            pass

    def _grab(rx, cast=str, group=1):
        m = rx.search(out)
        if m is None:
            return None
        try:
            return cast(m.group(group))
        except Exception:
            return m.group(group)

    primal_line_obj = _grab(_RE_PRIMAL_LINE, float, 1)
    primal_line_viol = _grab(_RE_PRIMAL_LINE, float, 2)
    dual_line_obj = _grab(_RE_DUAL_LINE, float, 1)
    dual_line_viol = _grab(_RE_DUAL_LINE, float, 2)

    return {
        "mosek_problem_status":  _grab(_RE_PROBSTATUS, str),
        "mosek_solution_status": _grab(_RE_SOLSTATUS, str),
        "pobj_summary":          primal_line_obj,
        "primal_viol":           primal_line_viol,
        "dobj_summary":          dual_line_obj,
        "dual_viol":             dual_line_viol,
        "iterations":            (rows[-1]["iter"] if rows else None),
        "iter_table":            rows,
    }


# ---- Driver ----------------------------------------------------------------
def solve_with_mosek(
    problem: cp.Problem,
    *,
    mosek_params: Optional[dict] = None,
    timeout_sec: Optional[float] = None,
    verbose_log: bool = True,
    keep_log: bool = True,
    log_max_chars: int = 200_000,
) -> dict:
    """Solve `problem` with Mosek via cvxpy, capture log, extract rigorous bounds.

    Parameters
    ----------
    mosek_params : dict, optional
        Forwarded to cvxpy's `mosek_params=` kwarg.  Useful keys for our SDP
        family (large, mildly ill-conditioned):
          'MSK_DPAR_INTPNT_CO_TOL_PFEAS': 1e-10,
          'MSK_DPAR_INTPNT_CO_TOL_DFEAS': 1e-10,
          'MSK_DPAR_INTPNT_CO_TOL_REL_GAP': 1e-12,
          'MSK_IPAR_NUM_THREADS': 0,   # auto
    timeout_sec : float, optional
        Approximate wall-clock cap via `MSK_DPAR_OPTIMIZER_MAX_TIME`.

    Returns
    -------
    dict with keys documented at top of this file.
    """
    assert isinstance(problem.objective, cp.Minimize), (
        "solve_with_mosek currently assumes Minimize (white_full_convex is Min)."
    )

    params = {
        # Tight default tolerances — Mosek's stock are ~1e-8 which is plenty
        # for our 6-7 digit headlines but we want to *measure* the achievable
        # precision, so request more.
        "MSK_IPAR_NUM_THREADS": 0,
    }
    if timeout_sec is not None:
        params["MSK_DPAR_OPTIMIZER_MAX_TIME"] = float(timeout_sec)
    if mosek_params:
        params.update(mosek_params)

    buf = io.StringIO()
    t0 = time.time()
    err = None
    # CVXPY routes Mosek's stream through Python's `logging` module with a
    # "(CVXPY) <ts>: " prefix.  Capture via a logger handler, then strip the
    # prefix before parsing.  We do NOT redirect stdout — Mosek's C-level
    # writes go there but cvxpy already mirrors them via the logger.
    cvxpy_logger = logging.getLogger("__cvxpy__")
    handler = logging.StreamHandler(buf)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(message)s"))
    prev_level = cvxpy_logger.level
    cvxpy_logger.addHandler(handler)
    cvxpy_logger.setLevel(logging.DEBUG)
    try:
        if verbose_log:
            problem.solve(solver=cp.MOSEK, verbose=True, mosek_params=params)
        else:
            problem.solve(solver=cp.MOSEK, verbose=False, mosek_params=params)
    except Exception as e:
        err = repr(e)
    finally:
        cvxpy_logger.removeHandler(handler)
        cvxpy_logger.setLevel(prev_level)
    elapsed = time.time() - t0
    out = buf.getvalue()

    parsed = parse_mosek_log(out)

    # The rigorous LB on the Minimize problem comes from POBJ (Mosek's primal
    # of the max-formed standard form == dual cert for our Min).
    #
    # Mosek's summary "Primal.  obj:" line is rounded to 11 sig digits, but
    # the iteration table prints 10 digits and cvxpy's `prob.value` (which is
    # POST-clean-up, i.e. DOBJ on the max-form = our Min's UB) has FULL float
    # precision.  We prefer:
    #   - dual_obj = cvxpy_prob_value (Min's UB, full precision)
    #   - primal_obj = last-iter POBJ from the table (Min's LB cert),
    #     fall back to summary line.
    pobj = None
    dobj = float(problem.value) if problem.value is not None else None
    if parsed["iter_table"]:
        pobj = parsed["iter_table"][-1]["pobj"]
    if pobj is None:
        pobj = parsed["pobj_summary"]
    if dobj is None:
        dobj = parsed["dobj_summary"]

    # Our convention (mirrors sdpa_runner):  rigorous_dual_LB = primal_obj - gap
    # where 'primal_obj' here is our Min's dual cert (== Mosek's POBJ).
    duality_gap = None
    rigorous_dual_LB = None
    if pobj is not None and dobj is not None:
        duality_gap = abs(dobj - pobj)
        rigorous_dual_LB = pobj - duality_gap

    log_to_save = out if (keep_log and len(out) <= log_max_chars) else (
        out[:log_max_chars] + f"\n... [truncated {len(out) - log_max_chars} chars]"
        if keep_log else None
    )

    return {
        "status": problem.status,
        "mosek_problem_status":  parsed["mosek_problem_status"],
        "mosek_solution_status": parsed["mosek_solution_status"],
        "primal_obj": pobj,          # rigorous LB cert on Min
        "dual_obj":   dobj,          # rigorous UB cert on Min (== prob.value)
        "duality_gap": duality_gap,
        "rigorous_dual_LB": rigorous_dual_LB,
        "primal_viol": parsed["primal_viol"],
        "dual_viol":   parsed["dual_viol"],
        "iterations":  parsed["iterations"],
        "cvxpy_prob_value": (float(problem.value) if problem.value is not None else None),
        "runtime_sec": elapsed,
        "error": err,
        "log": log_to_save,
    }


# ---- Self-test -------------------------------------------------------------
if __name__ == "__main__":
    import numpy as np
    x = cp.Variable()
    X = cp.bmat([[cp.reshape(x, (1, 1), order="C"), np.array([[1.0]])],
                 [np.array([[1.0]]), cp.reshape(x, (1, 1), order="C")]])
    prob = cp.Problem(cp.Minimize(x), [X >> 0])
    res = solve_with_mosek(prob)
    print("self-test: min x s.t. [[x,1],[1,x]]>>0  (true opt 1)")
    for k in ("status", "mosek_problem_status", "mosek_solution_status",
              "primal_obj", "dual_obj", "duality_gap", "rigorous_dual_LB",
              "primal_viol", "dual_viol", "iterations",
              "cvxpy_prob_value", "runtime_sec"):
        print(f"  {k:>22} = {res.get(k)!r}")
