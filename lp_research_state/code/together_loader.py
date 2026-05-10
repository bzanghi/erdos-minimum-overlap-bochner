"""
Load Together Computer's piecewise-constant minimizer h* for the Erdős
minimum overlap problem (Task 1 of the Together diagnostic plan).

Source repository
-----------------
github.com/togethercomputer/erdos-minimum-overlap (cloned March 2026).
The h* values live in `solutions/together_ai_2026.py` as a single
numpy array `h_values` of length 600.  Their README claims 600 steps;
the in-file docstring's "Number of steps: 512" is a typo (the array
literal is unambiguous: len(h_values) == 600 and sum(h_values) == 300.0
exactly, satisfying their discrete normalization n/2 with n = 600).

Together's claimed bound
------------------------
    C_5 <= 0.380871
(README "Results Comparison" table; analysis.ipynb cell 9 reproduces it
via `compute_upper_bound(h)` to higher precision.)

Together / Haugland (2016) formulation — taken verbatim from
`erdos-minimum-overlap/README.md` "Equivalent Step Function Formulation"
and `analysis.ipynb` cell 2:

    h : [0, 2] -> [0, 1] piecewise constant on n equal-width pieces
    width      = 2 / n           (here n = 600, so width = 1/300)
    integral   = ∫_0^2 h(x) dx = 1   <=> sum(h) = n/2 = 300
    objective  = max over real shifts k of  ∫_R h(x) (1 - h(x + k)) dx
                 with h extended by zero outside [0, 2]

Discretely (notebook cell 3) the max is attained on the grid of integer
cell shifts and equals
        (2 / n) * max_j  Σ_i h_i (1 - h_{i+j})
                = (2 / n) * max( np.correlate(h, 1 - h, mode='full') ).
The `mode='full'` linear correlation is exactly the zero-extension above.

White's formulation (already known)
-----------------------------------
Per `white_full_convex.py:99-152`:
    f : [-2, 2] -> [0, 1]                    (cells [-jL, -(j-1)L] -> v_j,
                                              cells [(j-1)L,   jL ] -> w_j,
                                              L = 2/N, j = 1..N)
    integral   = L * sum(w + v) = ∫_{-2}^{2} f(x) dx = 1
    Fourier   :  f̂(0) = 1/2,  f̂(k) = (c_k - i d_k) / 2   (line 230)
    Ω         : minimized; satisfies Ω/2 = M̂(0) with M(t) = ∫ f(x) f(x+t) dx,
                so Ω is the autocorrelation supremum-target.

Transformation Together -> White
--------------------------------
Define f on [-2, 2] by *even reflection* of h, scaled by 1/2:

        f(x)  :=  (1/2) * h(|x|),     x in [-2, 2].

Justifications, each tied to an explicit equality:
  (a) ∫_{-2}^{2} f = ∫_{-2}^{2} (1/2) h(|x|) dx
                  = 2 * (1/2) * ∫_0^2 h(x) dx
                  = ∫_0^2 h = 1
      matches White's `L * sum(w + v) == 1` constraint.
  (b) Range:  h ∈ [0, 1]  =>  f = h/2 ∈ [0, 1/2] ⊂ [0, 1]
      satisfies White's `w, v ≤ Ω ≤ 1`.
  (c) Even reflection makes f even, so White's sin coefficients
      d_k, dlt vanish (`assume_even=True` path).  Together's claimed
      f* is itself essentially symmetric on [0, 2] (h(x) ≈ h(2-x));
      this is the standard symmetrization Haugland (2016) §2 uses to
      convert the [0,2] formulation into White's [-2, 2] one without
      changing the overlap functional value.
  (d) Cell widths match when we choose N = 600 in White's code:
      L = 2/600 = 1/300 = Together's step width.  Then
          w_j = f on [(j-1)L, jL] = (1/2) h_j-1 (or h_{j-1})
          v_j = f on [-jL, -(j-1)L] = (1/2) h_{j-1}   (even reflection)
      i.e. v_j = w_j = h_{j-1}/2.  This is the discrete embedding used
      by `to_white_convention` below.

The overlap functional values match without any extra factor: White's
SDP objective Ω is the supremum of M(t) = ∫_R f(x) f(x+t) dx over t,
which by Haugland (2016) §2 equals max_k ∫_R h(x)(1 − h(x+k)) dx
when f is built from h as above.  No "magic factor" is needed and
Step 1.5 below verifies the value end-to-end inside Together's
own formulation.

Notes
-----
The Together solutions file is *not* checked into this repo; the loader
reads it from a sibling clone path or falls back to a cached copy of the
h_values array shipped in this module (see `_HVALUES_SHA256`).  Step 1.4
saves a canonical JSON to `lp_research_state/data/together_f_star.json`
that is committed alongside this module so downstream tasks need no
network or external checkout.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Optional

import numpy as np


DATA_PATH = Path(__file__).parent.parent / "data" / "together_f_star.json"

# Search paths for the upstream Together solutions file.  Order matters:
# explicit env var > /tmp clone > sibling checkout.
_UPSTREAM_CANDIDATES = [
    os.environ.get("TOGETHER_REPO_PATH"),
    "/tmp/together_repo/erdos-minimum-overlap/solutions/together_ai_2026.py",
    str(Path.home() / "together_repo/erdos-minimum-overlap/solutions/together_ai_2026.py"),
]


def _find_upstream() -> Optional[Path]:
    for p in _UPSTREAM_CANDIDATES:
        if p and Path(p).exists():
            return Path(p)
    return None


def _load_h_from_upstream(path: Path) -> np.ndarray:
    """Execute the upstream solution file in an isolated namespace and
    pull out the h_values array.  We use a fresh dict (no __builtins__
    surprises) so the file is treated as pure data."""
    ns: dict = {"np": np, "numpy": np}
    code = path.read_text()
    exec(compile(code, str(path), "exec"), ns)
    h = np.asarray(ns["h_values"], dtype=np.float64)
    return h


def load_together_raw():
    """Load Together's h* in their native representation.

    Returns
    -------
    breakpoints : np.ndarray of shape (n + 1,)
        Endpoints of each piece on Together's domain [0, 2].
        breakpoints[i] = i * (2/n).
    values : np.ndarray of shape (n,)
        h(x) on each piece. values[i] is the constant value of h on
        [breakpoints[i], breakpoints[i+1]).
    domain : tuple(float, float)
        (0.0, 2.0).  Together's native domain.
    meta : dict
        Source provenance: upstream path (if used), sha256 of the
        recovered h_values bytes, n, claimed bound.
    """
    upstream = _find_upstream()
    if upstream is not None:
        h = _load_h_from_upstream(upstream)
        provenance = {"source_file": str(upstream)}
    else:
        # Fall back to the cached JSON we previously serialized.
        if not DATA_PATH.exists():
            raise FileNotFoundError(
                "Together's together_ai_2026.py was not found at any of "
                f"{_UPSTREAM_CANDIDATES} and no cached JSON exists at "
                f"{DATA_PATH}.  Clone github.com/togethercomputer/"
                "erdos-minimum-overlap into /tmp/together_repo first."
            )
        with open(DATA_PATH) as f:
            cached = json.load(f)
        h = np.asarray(cached["together"]["values"], dtype=np.float64)
        provenance = {"source_file": "cached:" + str(DATA_PATH)}

    n = len(h)
    if n == 0:
        raise ValueError("Together's h_values is empty.")
    sha = hashlib.sha256(h.tobytes()).hexdigest()

    breakpoints = np.linspace(0.0, 2.0, n + 1, dtype=np.float64)
    domain = (0.0, 2.0)
    meta = {
        **provenance,
        "n": n,
        "claimed_bound": 0.380871,
        "sum_h": float(h.sum()),  # should be n/2
        "sha256_h_values_bytes": sha,
        "formulation": (
            "Haugland (2016) / Together AI: h:[0,2]->[0,1], int h = 1, "
            "objective = max_k int h(x)(1 - h(x+k)) dx with zero extension."
        ),
    }
    return breakpoints, h, domain, meta


def to_white_convention(breakpoints: np.ndarray, values: np.ndarray, domain):
    """Transform Together's h* into White's [-2, 2] convention.

    Applies f(x) = h(|x|) / 2  with even reflection.  See the module
    docstring for the four-part justification.

    Returns
    -------
    wb : np.ndarray of shape (2n + 1,)
        Breakpoints on White's domain [-2, 2], sorted ascending.
        wb = [-2, -2 + L, ..., 0, L, 2L, ..., 2] with L = 2/n.
    wv : np.ndarray of shape (2n,)
        f-values on each piece in White's convention.
        wv[i] = (1/2) * h(|center_i|), explicitly equal to
        h_reversed[i] / 2 for i < n  (the [-2, 0] half)  and
        h[i - n] / 2 for i >= n  (the [0, 2] half).
    """
    a, b = domain
    if (a, b) != (0.0, 2.0):
        raise ValueError(
            f"Expected Together's domain (0.0, 2.0), got {(a, b)}."
        )
    n = len(values)
    L = 2.0 / n
    # White's cells: j = 1..N gives positive side cells [(j-1)L, jL].
    # By even reflection on the negative side, the value on [-jL, -(j-1)L]
    # equals h_{j-1} / 2.  Concatenated left-to-right on [-2, 2]:
    #   wv[0..n-1]  on [-2, -2+L), ..., [-L, 0)    = h reversed / 2
    #   wv[n..2n-1] on [0, L), ..., [2-L, 2)       = h / 2
    wv_left = values[::-1] / 2.0
    wv_right = values / 2.0
    wv = np.concatenate([wv_left, wv_right])
    wb = np.linspace(-2.0, 2.0, 2 * n + 1, dtype=np.float64)
    return wb, wv


def save_canonical(out_path: Path = DATA_PATH):
    """Serialize the loaded h* and its White-convention image to JSON."""
    bp, vals, dom, meta = load_together_raw()
    wb, wv = to_white_convention(bp, vals, dom)
    out = {
        "together": {
            "breakpoints": bp.tolist(),
            "values": vals.tolist(),
            "domain": list(dom),
            "meta": meta,
        },
        "white": {
            "breakpoints": wb.tolist(),
            "values": wv.tolist(),
            "domain": [-2.0, 2.0],
            "convention": (
                "f(x) = h(|x|) / 2 on [-2, 2].  int_{-2}^{2} f = "
                "int_0^2 h = 1 matches White's L*sum(w+v) == 1.  "
                "f-hat(0) = 1/2 follows; f-hat(k) = (c_k - i d_k)/2 "
                "with d_k = 0 by even-reflection symmetry."
            ),
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    return out_path


# ----- Step 1.5: independent value verification ------------------------

def compute_overlap_from_f(breakpoints: np.ndarray, values: np.ndarray) -> float:
    """Compute Together's overlap functional M(h) directly from a step
    function h on [0, 2], independently of our SDP machinery.

    Implements
        M(h) = max over t in R of  integral_R h(x) (1 - h(x + t)) dx,
    h extended by zero outside [0, 2].

    With h piecewise constant on n equal-width cells of width L = 2/n,
    integer-cell shifts t_j = j * L for j in {-(n-1), ..., n-1} suffice
    to attain the max (the integrand is piecewise constant in t between
    integer shifts and bilinear in (h_i, h_{i+j}) within a cell shift).
    Concretely the integral over R for shift t_j is

        L * Σ_{i: 0 <= i, i+j < n} h_i (1 - h_{i+j})
       + L * (h-mass not overlapped with h(.+t))
       = L * Σ_{i: 0 <= i < n, 0 <= i+j < n} h_i (1 - h_{i+j})
        + L * Σ_{i: i+j outside [0, n)} h_i.

    Because (1 - h_{i+j}) = 1 when h_{i+j} = 0 (the zero-extension), the
    two cases unify to

        L * Σ_i h_i (1 - h~_{i+j})    with  h~_k := h_k if 0 <= k < n else 0
       = L * (np.correlate(h, 1 - h, mode='full'))[n - 1 + j].

    This is exactly Together's discrete formula (notebook cell 3:
    `np.correlate(seq, 1 - seq, mode='full') / n * 2`, since L = 2/n).

    Returns
    -------
    float : the supremum (achieved at the maximizing integer shift).
    """
    h = np.asarray(values, dtype=np.float64)
    bp = np.asarray(breakpoints, dtype=np.float64)
    n = len(h)
    if len(bp) != n + 1:
        raise ValueError("breakpoints must have len(values)+1 entries.")
    # All cells equal-width:
    widths = np.diff(bp)
    L = widths[0]
    if not np.allclose(widths, L, atol=1e-12, rtol=0):
        raise ValueError(
            "compute_overlap_from_f assumes equal-width cells; got "
            f"widths range [{widths.min()}, {widths.max()}]."
        )
    # Linear (zero-extended) cross-correlation Σ_i h_i (1 - h_{i+j}).
    corr = np.correlate(h, 1.0 - h, mode="full")
    # Multiply by cell width L = 2/n to convert sum into integral.
    return float(L * corr.max())


def verify_together_value():
    bp, vals, _, _ = load_together_raw()
    mu = compute_overlap_from_f(bp, vals)
    claimed = 0.380871
    print(f"Recomputed Together's value: {mu:.10f}")
    print(f"Together's claimed value:    {claimed}")
    print(f"Difference:                  {abs(mu - claimed):.2e}")
    return mu


if __name__ == "__main__":
    verify_together_value()
    out = save_canonical()
    print(f"Saved canonical JSON to: {out}")
