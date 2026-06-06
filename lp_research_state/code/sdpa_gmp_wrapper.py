"""Subprocess wrapper for the SDPA-GMP arbitrary-precision SDP solver.

SDPA-GMP is built locally from
https://github.com/sdpa-python/sdpa-multiprecision (a fork of Nakata's
original SDPA-GMP). The binary lives at
`lp_research_state/bin/sdpa_gmp`.

This module:
1. Locates the binary.
2. Accepts an SDP problem in standard SDPA-S sparse format (.dat-s) and
   runs sdpa_gmp on it, returning parsed primal/dual objective values and
   feasibility errors at GMP precision (typically ~10⁻⁷⁵ or better).
3. (Future) provides a builder that emits .dat-s from a cvxpy problem so
   the SDPs in white_full_convex.py can be cross-verified.

Usage
-----
>>> from sdpa_gmp_wrapper import solve_sdpa_file
>>> result = solve_sdpa_file('/tmp/myproblem.dat-s')
>>> print(result['objValPrimal'], result['objValDual'])

The wrapper is INTENTIONALLY thin — for any non-trivial use we'd build a
proper cvxpy→SDPA-S serializer (deferred to a follow-up; the gating value
is having the binary, not the interface).
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_BIN = _PROJECT_ROOT / "lp_research_state" / "bin" / "sdpa_gmp"
_PARAM = _PROJECT_ROOT / "lp_research_state" / "bin" / "param.sdpa"


def binary_path() -> Path:
    """Return path to the sdpa_gmp binary; raise if not present."""
    if not _BIN.exists():
        raise FileNotFoundError(
            f"sdpa_gmp not found at {_BIN}. Build it via "
            "/tmp/sdpa_build/build_all.sh and copy to lp_research_state/bin/."
        )
    return _BIN


def solve_sdpa_file(dat_s: str | Path,
                    output_path: str | Path | None = None,
                    param_path: str | Path | None = None,
                    pt: int = 0,
                    timeout: float = 3600.0) -> dict[str, Any]:
    """Run sdpa_gmp on a .dat-s file; return parsed result dict.

    Parameters
    ----------
    dat_s : path to SDPA-S sparse input file.
    output_path : where to write sdpa_gmp's output (default: input + '.out').
    param_path : path to parameter file (default: bundled param.sdpa).
    pt : print level (0=quiet, 1=normal, 2=verbose). Default 0.
    timeout : seconds before subprocess kill.

    Returns dict with keys:
      'status' : 'SDP_SUCCESS' / 'pdOPT' / etc., as reported by SDPA
      'objValPrimal', 'objValDual'  (strings retaining full GMP precision)
      'pFeasError', 'dFeasError'    (GMP-precision residuals)
      'totalTime', 'iterations'
      'raw_output' : the full sdpa_gmp output for further inspection
    """
    dat_s = Path(dat_s).resolve()
    if not dat_s.exists():
        raise FileNotFoundError(dat_s)
    output_path = Path(output_path or (str(dat_s) + ".out")).resolve()
    param_path = Path(param_path or _PARAM).resolve()

    cmd = [
        str(binary_path()),
        str(dat_s),
        str(output_path),
        "-p", str(param_path),
        "-pt", str(pt),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    raw = proc.stdout + ("\n--- stderr ---\n" + proc.stderr if proc.stderr else "")

    result = {"raw_output": raw, "returncode": proc.returncode}
    # Parse standard SDPA outputs from stdout
    patterns = {
        "status": r"phase\.value\s*=\s*(\S+)",
        "objValPrimal": r"objValPrimal\s*=\s*(\S+)",
        "objValDual": r"objValDual\s*=\s*(\S+)",
        "pFeasError": r"p\.feas\.error\s*=\s*(\S+)",
        "dFeasError": r"d\.feas\.error\s*=\s*(\S+)",
        "totalTime": r"total time\s*=\s*(\S+)",
        "iterations": r"iteration\s*=\s*(\d+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, raw)
        result[key] = m.group(1) if m else None

    return result


def selftest() -> None:
    """Solve the bundled example1.dat-s and report the result."""
    example = Path(__file__).resolve().parents[2] / "lp_research_state" / "bin" / "example1.dat-s"
    if not example.exists():
        # Try the build tree
        example = Path("/tmp/sdpa_build/sdpa-multiprecision/example1.dat-s")
    if not example.exists():
        raise FileNotFoundError(
            "example1.dat-s not found; copy it from the SDPA-multiprecision "
            "source tree to lp_research_state/bin/ first."
        )
    print(f"Solving {example} ...")
    r = solve_sdpa_file(example, pt=0)
    print(f"  status        = {r['status']}")
    print(f"  objValPrimal  = {r['objValPrimal']}")
    print(f"  objValDual    = {r['objValDual']}")
    print(f"  pFeasError    = {r['pFeasError']}")
    print(f"  dFeasError    = {r['dFeasError']}")
    print(f"  iterations    = {r['iterations']}")
    print(f"  totalTime     = {r['totalTime']}")


if __name__ == "__main__":
    selftest()
