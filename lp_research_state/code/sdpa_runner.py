"""
Build -> serialize -> exec sdpa_gmp -> parse pipeline.

Public entry point: solve_with_sdpa_gmp(problem, **kwargs) -> dict with
  primal_obj         (SDPA's objValPrimal, in original problem units)
  dual_obj           (SDPA's objValDual,   in original problem units)
  rigorous_dual_LB   (= dual_obj - duality_gap, preserving the project's
                      "subtract the gap to make it strictly rigorous" rule)
  duality_gap        (SDPA's reported absolute gap)
  precision_digits   (SDPA's "digits" line)
  phase              ('pdOPT', 'pdFEAS', 'dUNBD', ...)
  runtime_sec        (wall clock incl. file write + sdpa exec + parse)
  stdout_path, sdpa_dat_s_path  (kept for inspection)

Implementation notes
--------------------
- sdpa_gmp is invoked with -ds <data> -o <output> -p <param>. Output goes to a
  file we then parse for the summary block ("objValPrimal", "objValDual",
  "relative gap", "gap", "digits", "phase.value").
- For a MIN problem CVXPY/SDPA both report objValPrimal close to objValDual at
  convergence. The dual is the LOWER bound; we report it as rigorous_dual_LB
  minus the absolute gap to be safely below the true optimum, matching the
  convention in dual_extractor.py.
- For a MAX problem we DO NOT support here yet; the white_full_convex.py
  program is always a Minimize.

"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional

import cvxpy as cp

from sdpa_serializer import cvxpy_to_sdpa_s


_THIS_DIR = Path(__file__).resolve().parent
_REPO_BIN = (_THIS_DIR.parent / "bin").resolve()

DEFAULT_SDPA_BIN = str(_REPO_BIN / "sdpa_gmp")
# Custom param tuned for the white_full_convex problem family:
#   lambdaStar = 1.0    (default 1e4 violates Omega<=1 box bounds wildly)
#   epsilonStar = 1e-25 (default 1e-30 unreachable past ~25 digits at our scales)
#   betaBar = 0.2       (default 0.3 risks too-aggressive corrector steps)
# Falls back to the stock SDPA param if the white-tuned file is missing.
_PARAM_WHITE = _REPO_BIN / "param.sdpa.white"
_PARAM_STOCK = _REPO_BIN / "param.sdpa"
DEFAULT_PARAM_FILE = str(_PARAM_WHITE if _PARAM_WHITE.exists() else _PARAM_STOCK)


# --------------------------------------------------------------------------- #
# Output parsing
# --------------------------------------------------------------------------- #
_RE_PHASE  = re.compile(r"phase\.value\s*=\s*(\S+)")
_RE_PRIMAL = re.compile(r"objValPrimal\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_DUAL   = re.compile(r"objValDual\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_GAP    = re.compile(r"^\s*gap\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)", re.MULTILINE)
_RE_RELGAP = re.compile(r"relative gap\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_DIGITS = re.compile(r"digits\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_PFEAS  = re.compile(r"p\.feas\.error\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_DFEAS  = re.compile(r"d\.feas\.error\s*=\s*([+\-]?\d+\.\d+e[+\-]?\d+)")
_RE_ITER   = re.compile(r"Iteration\s*=\s*(\d+)")


def parse_sdpa_output(text: str) -> dict:
    """Pull summary fields out of sdpa_gmp's stdout/output file."""
    def _grab(rx, cast=float):
        m = rx.search(text)
        if m is None:
            return None
        try:
            return cast(m.group(1))
        except Exception:
            return m.group(1)
    out = {
        "phase": _grab(_RE_PHASE, cast=str),
        "primal_obj_raw": _grab(_RE_PRIMAL),
        "dual_obj_raw":   _grab(_RE_DUAL),
        "duality_gap":    _grab(_RE_GAP),
        "relative_gap":   _grab(_RE_RELGAP),
        "precision_digits": _grab(_RE_DIGITS),
        "p_feas_error":   _grab(_RE_PFEAS),
        "d_feas_error":   _grab(_RE_DFEAS),
        "iterations":     _grab(_RE_ITER, cast=int),
    }
    return out


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #
def solve_with_sdpa_gmp(
    problem: cp.Problem,
    *,
    sdpa_bin: str = DEFAULT_SDPA_BIN,
    param_file: Optional[str] = DEFAULT_PARAM_FILE,
    work_dir: Optional[str] = None,
    keep_files: bool = True,
    timeout_sec: Optional[float] = None,
    verbose: bool = False,
) -> dict:
    """Serialize `problem` to SDPA-S, run sdpa_gmp, parse the result.

    Returns a dict including primal_obj, dual_obj, rigorous_dual_LB (all in the
    ORIGINAL CVXPY problem's objective units), runtime_sec, precision_digits,
    phase, and paths to the kept artifacts.
    """
    assert isinstance(problem.objective, cp.Minimize), (
        "Only Minimize problems are supported (white_full_convex is Minimize)."
    )

    if not os.path.exists(sdpa_bin):
        raise FileNotFoundError(f"sdpa_gmp binary not found at {sdpa_bin}")
    if param_file is not None and not os.path.exists(param_file):
        raise FileNotFoundError(f"param file not found at {param_file}")

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="sdpa_gmp_")
    else:
        os.makedirs(work_dir, exist_ok=True)
    work_dir = str(Path(work_dir).resolve())

    dat_s = os.path.join(work_dir, "problem.dat-s")
    out_f = os.path.join(work_dir, "sdpa.out")
    log_f = os.path.join(work_dir, "sdpa.log")

    t0 = time.time()
    meta = cvxpy_to_sdpa_s(problem, dat_s)
    t_serialize = time.time() - t0

    cmd = [sdpa_bin, "-ds", dat_s, "-o", out_f]
    if param_file is not None:
        cmd += ["-p", param_file]
    if verbose:
        print("  [sdpa_runner] cmd:", " ".join(cmd))

    t1 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_sec,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "status": "timeout",
            "phase": None,
            "primal_obj": None, "dual_obj": None, "duality_gap": None,
            "relative_gap": None, "rigorous_dual_LB": None,
            "precision_digits": None, "iterations": None,
            "p_feas_error": None, "d_feas_error": None,
            "runtime_sec": time.time() - t0,
            "serialize_sec": t_serialize, "solve_sec": None,
            "sdpa_dat_s_path": dat_s, "sdpa_out_path": None,
            "stdout_path": log_f, "m": None, "block_structure": None,
            "block_kinds": None, "returncode": None,
            "error": f"sdpa_gmp timed out after {timeout_sec}s",
        }
    t_solve = time.time() - t1

    with open(log_f, "w") as fh:
        fh.write("=== STDOUT ===\n")
        fh.write(proc.stdout or "")
        fh.write("\n=== STDERR ===\n")
        fh.write(proc.stderr or "")

    sdpa_text = ""
    if os.path.exists(out_f):
        with open(out_f) as fh:
            sdpa_text = fh.read()
    # also include stdout — sdpa_gmp prints the same summary to both
    if not sdpa_text:
        sdpa_text = proc.stdout or ""

    parsed = parse_sdpa_output(sdpa_text)

    # The objective vector we wrote IS CVXPY's canonical c-vector. CVXPY's
    # Minimize objective in canonical form may have an additive offset; for
    # white_full_convex.build_problem the objective is literally the Variable
    # `Omega`, so the offset is 0. We accept SDPA's reported objValPrimal /
    # objValDual at face value.
    pv = parsed["primal_obj_raw"]
    dv = parsed["dual_obj_raw"]
    gap = parsed["duality_gap"]
    rig = None
    if dv is not None and gap is not None:
        rig = dv - abs(gap)
    elif dv is not None:
        rig = dv

    runtime = time.time() - t0
    result = {
        "status": "ok" if (parsed["phase"] in ("pdOPT", "pdFEAS")) else (parsed["phase"] or "unknown"),
        "phase": parsed["phase"],
        "primal_obj": pv,
        "dual_obj": dv,
        "duality_gap": gap,
        "relative_gap": parsed["relative_gap"],
        "rigorous_dual_LB": rig,
        "precision_digits": parsed["precision_digits"],
        "iterations": parsed["iterations"],
        "p_feas_error": parsed["p_feas_error"],
        "d_feas_error": parsed["d_feas_error"],
        "runtime_sec": runtime,
        "serialize_sec": t_serialize,
        "solve_sec": t_solve,
        "sdpa_dat_s_path": dat_s,
        "sdpa_out_path": out_f,
        "stdout_path": log_f,
        "m": meta["m"],
        "block_structure": meta["block_structure"],
        "block_kinds": meta["block_kinds"],
        "returncode": proc.returncode,
    }
    if not keep_files:
        shutil.rmtree(work_dir, ignore_errors=True)
    return result


if __name__ == "__main__":
    # Tiny self-test: a 2x2 LMI we know the answer to.
    import numpy as np
    x = cp.Variable()
    # min x  s.t.  [[x, 1], [1, x]] >= 0   ->  x >= 1
    X = cp.bmat([[cp.reshape(x, (1, 1), order="C"), np.array([[1.0]])],
                 [np.array([[1.0]]), cp.reshape(x, (1, 1), order="C")]])
    prob = cp.Problem(cp.Minimize(x), [X >> 0])
    res = solve_with_sdpa_gmp(prob, keep_files=True, verbose=True)
    print("self-test:", {k: v for k, v in res.items() if k not in ("block_kinds",)})
