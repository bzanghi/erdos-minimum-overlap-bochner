"""
Load Together Computer's piecewise-constant minimizer h* for the Erdős
minimum overlap problem (Task 1 of the Together diagnostic plan).

Source repository
-----------------
github.com/togethercomputer/erdos-minimum-overlap (cloned March 2026).
The h* values live in `solutions/together_ai_2026.py` as a single
numpy array `h_values`.  Their in-file docstring claims "Number of
steps: 512" but the array literal is unambiguously of length 600 with
sum 300.0 (= 600/2, matching their discrete normalization n/2).  This
discrepancy is asserted in `load_together_raw` so a future upstream
change cannot silently produce a different `wv`.

Together's claimed bound
------------------------
    µ <= 0.380871
(README "Results Comparison" table; analysis.ipynb cell 9 reproduces it
via `compute_upper_bound(h)` to higher precision.)

Together / Haugland (2016) formulation — taken verbatim from
`erdos-minimum-overlap/README.md` "Equivalent Step Function Formulation"
and `analysis.ipynb` cell 2:

    h : [0, 2] -> [0, 1] piecewise constant on n equal-width pieces
    width      = 2 / n           (here n = 600, so width = 1/300)
    integral   = ∫_0^2 h(x) dx = 1   <=> sum(h) = n/2 = 300
    objective  = M(h) := max over real shifts k of
                 ∫_R h(x) (1 - h(x + k)) dx,  h zero-extended.

Discretely (notebook cell 3) the max is attained on the grid of integer
cell shifts and equals
        (2 / n) * max_j  Σ_i h_i (1 - h_{i+j})
                = (2 / n) * max( np.correlate(h, 1 - h, mode='full') ).

White's formulation (already known)
-----------------------------------
Per `white_full_convex.py:99-152`:
    f : [-2, 2] -> [0, 1]                    (cells [(j-1)L,   jL ] -> w_j,
                                              cells [-jL, -(j-1)L] -> v_j,
                                              L = 2/N, j = 1..N)
    ∫f = 1    via  L * sum(w + v) == 1           (line 142)
    0 ≤ f     via  w >= 0, v >= 0                (line 141)
    f ≤ Ω ≤ 1 via  w <= Ω, v <= Ω, Ω <= 1        (line 141)
    f̂(0) = 1/2  follows from ∫f = 1 over [-2,2] when f is identified
                with its periodic extension of period 4 (line 230 of
                `white_full_convex.py` writes "f̂(0)=1/2").
    Ω         : minimized; satisfies Ω/2 = M̂(0) with
                M(t) = ∫ f(x) f(x+t) dx (autocorrelation of f over [-2,2]).

The non-`assume_even=True` path (line 149) does NOT impose v == w, so
White's variable space naturally accommodates non-even f.  We use this
to offer two embeddings (see below).

CRITICAL: White's Ω and Together's M are DIFFERENT functionals
---------------------------------------------------------------
White's program minimizes Ω, which majorizes the SUPREMUM over t of
M_W(t) := ∫_{-2}^{2} f(x) f(x+t) dx     (autocorrelation of f).

Together's M majorizes the SUPREMUM over k of
M_T(k) := ∫_R h(x) (1 - h(x+k)) dx       (Haugland 2016 formulation,
                                          zero extension).

These two suprema both equal the Erdős minimum-overlap constant µ
in the limit of optimal f / optimal h, but they are NOT pointwise
equal on a given input.  Empirically for Together's h*:

    M_T(h*) = 0.38087  (Together's claimed bound)
    Ω(f)    = 0.38734  using f(x) = h*(|x|)/2 on [-2, 2]
                       (computed via verify_white_embedding below)

The 0.38734 vs 0.38087 gap is not a bug in either project; it reflects
the gap between sup_t M_W(t) and sup_k M_T(k) at this specific h*.

Two embeddings supplied
-----------------------
We expose BOTH embeddings.  Downstream tasks decide which to use.

  1) `to_white_convention_even(bp, vals, dom)`
        f_even(x) := (1/2) * h(|x|),   x in [-2, 2]
     - even by construction, so White's `assume_even=True` path is valid;
     - ∫ f = ∫_0^2 h = 1  (matches line 142);
     - f in [0, 1/2] ⊂ [0, 1]  (matches line 141);
     - LOSES information: symmetrizes Together's asymmetric h*.
       Empirically max |h_i - h_{n-1-i}| ≈ 0.53, so Together's h* is
       materially asymmetric and the symmetrized image is a strictly
       different function.

  2) `to_white_convention_direct(bp, vals, dom)`
        f_dir on [-2, 2] defined by
            f_dir(x) := (1/2) * h(x)       for x in [0, 2]
            f_dir(x) := (1/2) * h(-x)      for x in [-2, 0]
        i.e. the same even reflection BUT we also offer the variant
            f_asym(x) := h(x) / 2          for x in [0, 2]
            f_asym(x) := 0                  for x in [-2, 0)
       (zero-extension to the negative half-domain).
       - NOT even; must be used with `assume_even=False`;
       - ∫ f_asym = ∫_0^2 (h/2) = 1/2 ≠ 1  → does NOT satisfy line 142.
         So this is NOT admissible as-is.  To make it admissible, one
         can use  g(x) := h(x) on [0, 2], g(x) := 0 elsewhere; then
         ∫ g = 1 but g ∈ [0, 1] still, however g̅(0) = 1/4 not 1/2.
         The natural way to retain admissibility AND asymmetry is the
         "shifted" embedding  f(x) := h(x+2)/2 + h(-x+2)/2 type
         reflections — but those re-symmetrize again.  In short, every
         even, [0,2]-zero-extended embedding into White's domain that
         preserves ∫f = 1 and 0 ≤ f ≤ 1 forces symmetrization.

     We therefore expose, as the second embedding, the truly direct
     no-rescale variant:
        f_direct(x) := h(x)  for x in [0, 2], 0 elsewhere on [-2, 2].
     This is admissible (∫f = 1, 0 ≤ f ≤ 1, since ∫_0^2 h = 1 and
     h ∈ [0,1]).  It is NOT even, so it must be supplied with
     `assume_even=False`.  In White's discretization (line 134),
     v_j = 0 for all j and w_j = h_{j-1} for j = 1..N.

The two embeddings (`_even` and `_direct`) give different White-side
diagnostics; downstream tasks may want both.  Verification is provided
by `verify_white_embedding(kind=...)`.

Notes
-----
The Together solutions file is *not* checked into this repo.  We:
  - parse it via `ast.literal_eval` (NOT `exec`) for safety;
  - assert n == 600 and sum(h) == 300.0 exactly;
  - cache a SHA-256 of the raw float64 bytes for provenance.
Step 1.4 saves a canonical JSON to `lp_research_state/data/together_f_star.json`.
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
from pathlib import Path
from typing import Optional, Tuple

import numpy as np


DATA_PATH = Path(__file__).parent.parent / "data" / "together_f_star.json"

# Search paths for the upstream Together solutions file.  Order matters:
# explicit env var > /tmp clone > sibling checkout.
_UPSTREAM_CANDIDATES = [
    os.environ.get("TOGETHER_REPO_PATH"),
    "/tmp/together_repo/erdos-minimum-overlap/solutions/together_ai_2026.py",
    str(Path.home() / "together_repo/erdos-minimum-overlap/solutions/together_ai_2026.py"),
]

# Hard invariants of Together's file (asserted on every load).
_EXPECTED_N = 600
_EXPECTED_SUM = 300.0
_CLAIMED_BOUND = 0.380871


def _find_upstream() -> Optional[Path]:
    for p in _UPSTREAM_CANDIDATES:
        if p and Path(p).exists():
            return Path(p)
    return None


def _load_h_from_upstream(path: Path) -> np.ndarray:
    """Parse the upstream solution file SAFELY (no `exec`).

    The file contains exactly one top-level assignment of the form
    `h_values = np.array([...literal floats...])`.  We walk the AST,
    locate that assignment, confirm the RHS is a `np.array(...)` call
    with a single list-of-numbers argument, and use `ast.literal_eval`
    on that list.  No arbitrary code is executed.
    """
    code = path.read_text()
    tree = ast.parse(code, filename=str(path))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if len(node.targets) != 1:
            continue
        tgt = node.targets[0]
        if not (isinstance(tgt, ast.Name) and tgt.id == "h_values"):
            continue
        # Expect: np.array([...]) call.
        rhs = node.value
        if not isinstance(rhs, ast.Call):
            raise ValueError(
                f"Expected h_values RHS to be a np.array(...) call in {path}; "
                f"got {ast.dump(rhs)[:80]}..."
            )
        # Confirm the callee is np.array (or numpy.array).
        callee = rhs.func
        callee_ok = False
        if isinstance(callee, ast.Attribute) and callee.attr == "array":
            callee_ok = isinstance(callee.value, ast.Name) and callee.value.id in {"np", "numpy"}
        if not callee_ok:
            raise ValueError(
                f"Expected h_values = np.array(...) in {path}; "
                f"got callee {ast.dump(callee)}."
            )
        if not rhs.args:
            raise ValueError(f"np.array(...) call in {path} has no arguments.")
        # The first positional arg should be a list literal of numbers.
        list_node = rhs.args[0]
        # literal_eval handles list of int/float and unary minus uniformly.
        try:
            data = ast.literal_eval(list_node)
        except (ValueError, SyntaxError) as e:
            raise ValueError(
                f"Could not literal-eval the np.array list in {path}: {e}"
            ) from e
        return np.asarray(data, dtype=np.float64)
    raise ValueError(f"No `h_values = ...` assignment found in {path}.")


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

    Asserts
    -------
    - len(h) == 600 exactly (Together's in-file docstring says 512 but
      the array literal is length 600; assert protects against silent
      upstream change).
    - h.sum() == 300.0 exactly (Together's discrete normalization n/2).
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
    # Hard invariants — fail loudly on upstream changes.
    assert n == _EXPECTED_N, (
        f"Expected {_EXPECTED_N} cells (per Together's array literal), got {n}."
    )
    sum_h = float(h.sum())
    # Exact equality holds for the upstream file; allow ~1e-12 for the
    # cached-JSON round-trip via base-10 text serialization.
    assert abs(sum_h - _EXPECTED_SUM) < 1e-9, (
        f"Expected sum(h) == {_EXPECTED_SUM} (Together's n/2 normalization), "
        f"got {sum_h}."
    )
    sha = hashlib.sha256(h.tobytes()).hexdigest()

    breakpoints = np.linspace(0.0, 2.0, n + 1, dtype=np.float64)
    domain = (0.0, 2.0)
    meta = {
        **provenance,
        "n": n,
        "claimed_bound": _CLAIMED_BOUND,
        "sum_h": sum_h,  # should be n/2 = 300.0
        "sha256_h_values_bytes": sha,
        "formulation": (
            "Haugland (2016) / Together AI: h:[0,2]->[0,1], int h = 1, "
            "objective = max_k int h(x)(1 - h(x+k)) dx with zero extension."
        ),
    }
    return breakpoints, h, domain, meta


def _check_domain(domain) -> None:
    a, b = domain
    if (a, b) != (0.0, 2.0):
        raise ValueError(
            f"Expected Together's domain (0.0, 2.0), got {(a, b)}."
        )


def to_white_convention_even(
    breakpoints: np.ndarray, values: np.ndarray, domain
) -> Tuple[np.ndarray, np.ndarray]:
    """Even-reflection embedding into White's [-2, 2] convention.

    f_even(x) := (1/2) * h(|x|)     for x in [-2, 2].

    Derivation (each step tied to a `white_full_convex.py` line):
      (a) ∫_{-2}^{2} f_even = 2 * (1/2) * ∫_0^2 h = 1
          matches line 142  (L * cp.sum(w + v) == 1).
      (b) h ∈ [0, 1]  ⇒  f_even ∈ [0, 1/2] ⊂ [0, 1]
          satisfies line 141  (w, v >= 0;  w, v <= Ω <= 1).
      (c) f̂(0) = (1/4) ∫_{-2}^{2} f = 1/4 in the convention of
          period-4 Fourier; line 230 documents f̂(0) = 1/2 in White's
          ½-renormalized convention, so consistency follows because
          line 230's "f̂(0) = 1/2" is a *statement of the constraint*
          line 142, not an independent claim.

    Discretization (N = len(values) = n):
      L = 2/n  matches Together's step width 1/300.
      For j = 1..N:
        w_j = average of f_even on [(j-1)L, jL] = h_{j-1} / 2
        v_j = average of f_even on [-jL, -(j-1)L] = h_{j-1} / 2  (even reflection)
      So v_j = w_j (the `assume_even=True` invariant on line 150 is
      automatically satisfied here, though we don't ENFORCE it — caller
      decides whether to set `assume_even=True`).

    Caveat
    ------
    Together's h* is empirically asymmetric: max |h_i - h_{n-1-i}| ≈ 0.53.
    This embedding symmetrizes it, so f_even ≠ "Together's function";
    it is the symmetrized version.  Use `to_white_convention_direct`
    if you want to preserve Together's actual shape.

    Returns
    -------
    wb : np.ndarray of shape (2n + 1,)  — breakpoints on [-2, 2].
    wv : np.ndarray of shape (2n,)      — f_even values on each cell,
                                          ordered left-to-right.
    """
    _check_domain(domain)
    wv_left = values[::-1] / 2.0
    wv_right = values / 2.0
    wv = np.concatenate([wv_left, wv_right])
    wb = np.linspace(-2.0, 2.0, 2 * len(values) + 1, dtype=np.float64)
    return wb, wv


def to_white_convention_direct(
    breakpoints: np.ndarray, values: np.ndarray, domain
) -> Tuple[np.ndarray, np.ndarray]:
    """Direct (asymmetric, zero-extended) embedding.

    f_direct(x) := h(x)     for x in [0, 2]
                := 0         for x in [-2, 0)

    Derivation:
      (a) ∫_{-2}^{2} f_direct = ∫_0^2 h = 1
          matches line 142.
      (b) h ∈ [0, 1]  ⇒  f_direct ∈ [0, 1]  (and = 0 on [-2, 0))
          satisfies line 141.
      (c) f_direct is NOT even — d_k, dlt in White's program will be
          nonzero in general.  Must NOT be combined with
          `assume_even=True` (line 149).

    Discretization (N = len(values) = n):
      L = 2/n. For j = 1..N:
        w_j = h_{j-1}     (positive side carries all of Together's h)
        v_j = 0           (negative side is zero-extended)

    Preserves Together's h* shape verbatim on the positive half.
    """
    _check_domain(domain)
    n = len(values)
    wv_left = np.zeros(n, dtype=np.float64)
    wv_right = np.asarray(values, dtype=np.float64).copy()
    wv = np.concatenate([wv_left, wv_right])
    wb = np.linspace(-2.0, 2.0, 2 * n + 1, dtype=np.float64)
    return wb, wv


# Back-compat alias: the original Task-1 function name.
def to_white_convention(breakpoints, values, domain):
    """Deprecated name for `to_white_convention_even`.

    Retained so existing callers (if any) keep working.  New code should
    pick `to_white_convention_even` or `to_white_convention_direct`
    explicitly to make the embedding choice visible.
    """
    return to_white_convention_even(breakpoints, values, domain)


def save_canonical(out_path: Path = DATA_PATH):
    """Serialize the loaded h* and its White-convention images to JSON.

    Saves both the even and direct embeddings under `white_even` and
    `white_direct` keys, plus the raw `together` data.
    """
    bp, vals, dom, meta = load_together_raw()
    wb_e, wv_e = to_white_convention_even(bp, vals, dom)
    wb_d, wv_d = to_white_convention_direct(bp, vals, dom)
    out = {
        "together": {
            "breakpoints": bp.tolist(),
            "values": vals.tolist(),
            "domain": list(dom),
            "meta": meta,
        },
        "white_even": {
            "breakpoints": wb_e.tolist(),
            "values": wv_e.tolist(),
            "domain": [-2.0, 2.0],
            "convention": (
                "f_even(x) = h(|x|) / 2 on [-2, 2].  int f = 1 matches "
                "white_full_convex.py:142.  0 <= f <= 1/2 matches line 141. "
                "Even by construction (d_k = dlt = 0; v_j = w_j); "
                "compatible with assume_even=True."
            ),
        },
        "white_direct": {
            "breakpoints": wb_d.tolist(),
            "values": wv_d.tolist(),
            "domain": [-2.0, 2.0],
            "convention": (
                "f_direct(x) = h(x) on [0, 2], 0 on [-2, 0).  int f = 1 "
                "matches line 142.  0 <= f <= 1 matches line 141. "
                "NOT even — must use assume_even=False."
            ),
        },
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    return out_path


# ----- Step 1.5: independent value verification ------------------------

def compute_overlap_from_f(breakpoints: np.ndarray, values: np.ndarray) -> float:
    """Compute Together's overlap functional M_T(h) directly from a step
    function h on [0, 2], independently of our SDP machinery.

    Implements
        M_T(h) = max over t in R of  integral_R h(x) (1 - h(x + t)) dx,
    h extended by zero outside [0, 2].

    With h piecewise constant on n equal-width cells of width L = 2/n,
    integer-cell shifts t_j = j * L for j in {-(n-1), ..., n-1} suffice
    to attain the max (the integrand is piecewise constant in t between
    integer shifts; the candidate value at each j is bilinear in (h_i,
    h_{i+j})).  Equals  L * max(np.correlate(h, 1 - h, mode='full')).
    """
    h = np.asarray(values, dtype=np.float64)
    bp = np.asarray(breakpoints, dtype=np.float64)
    n = len(h)
    if len(bp) != n + 1:
        raise ValueError("breakpoints must have len(values)+1 entries.")
    widths = np.diff(bp)
    L = widths[0]
    if not np.allclose(widths, L, atol=1e-12, rtol=0):
        raise ValueError(
            "compute_overlap_from_f assumes equal-width cells; got "
            f"widths range [{widths.min()}, {widths.max()}]."
        )
    corr = np.correlate(h, 1.0 - h, mode="full")
    return float(L * corr.max())


def compute_white_omega_from_f(
    breakpoints: np.ndarray, values: np.ndarray
) -> float:
    """Compute White's Ω(f) = sup_t ∫_{-2}^{2} f(x) f(x+t) dx for a step
    function f on [-2, 2] (zero-extended).

    Discretely, with f piecewise constant on 2n equal-width cells of width
    L = 2/n:
        Ω(f) = L * max_j Σ_i f_i f_{i+j}
             = L * max(np.correlate(f, f, mode='full')).
    """
    f = np.asarray(values, dtype=np.float64)
    bp = np.asarray(breakpoints, dtype=np.float64)
    if len(bp) != len(f) + 1:
        raise ValueError("breakpoints must have len(values)+1 entries.")
    widths = np.diff(bp)
    L = widths[0]
    if not np.allclose(widths, L, atol=1e-12, rtol=0):
        raise ValueError("compute_white_omega_from_f assumes equal-width cells.")
    corr = np.correlate(f, f, mode="full")
    return float(L * corr.max())


def verify_together_value():
    bp, vals, _, _ = load_together_raw()
    mu = compute_overlap_from_f(bp, vals)
    claimed = _CLAIMED_BOUND
    print(f"Recomputed Together's value: {mu:.10f}")
    print(f"Together's claimed value:    {claimed}")
    print(f"Difference:                  {abs(mu - claimed):.2e}")
    return mu


def verify_white_embedding(kind: str = "even") -> dict:
    """Empirically verify that the White-convention embedding lands in
    White's feasible set, and report Ω(f) alongside M_T(h).

    kind : "even" or "direct"
        Which embedding to verify.

    Checks (each tied to a `white_full_convex.py` line):
      - L * sum(wv) == 1                         (line 142,  ∫f = 1)
      - 0 <= wv <= 1                              (line 141,  0 ≤ f ≤ Ω ≤ 1)
      - 2 * mean(wv) == 1  (i.e. f̂(0) = 1/2 in White's convention,
        which is the same statement as ∫f = 1 over [-2, 2]) — line 230.

    Prints, for the same h*:
      M_T(h)  — Together's overlap functional   (sup_k ∫ h(x)(1-h(x+k)) dx)
      Ω(f)    — White's autocorrelation         (sup_t ∫ f(x) f(x+t) dx)
    These are DIFFERENT functionals (see module docstring).  The expected
    qualitative relation is  Ω(f) >= 2 * (∫f)^2 / |support| - O(.)  — we
    do NOT assert a sharp inequality, only that both values are positive
    and at most 1 (since 0 ≤ f ≤ 1 forces ∫ f f(.+t) ≤ ∫ f = 1).
    """
    bp, vals, dom, _ = load_together_raw()
    if kind == "even":
        wb, wv = to_white_convention_even(bp, vals, dom)
    elif kind == "direct":
        wb, wv = to_white_convention_direct(bp, vals, dom)
    else:
        raise ValueError(f"kind must be 'even' or 'direct'; got {kind!r}.")

    L = (wb[-1] - wb[0]) / len(wv)
    int_f = L * float(np.sum(wv))
    f_min = float(np.min(wv))
    f_max = float(np.max(wv))
    fhat0 = int_f / 2.0  # f̂(0) in period-4 Fourier convention.

    m_together = compute_overlap_from_f(bp, vals)
    omega_white = compute_white_omega_from_f(wb, wv)

    print(f"--- verify_white_embedding(kind={kind!r}) ---")
    print(f"  L = 2/N    = {L:.12f}")
    print(f"  int f      = {int_f:.12f}   (expect 1; line 142)")
    print(f"  min f      = {f_min:.6e}    (expect >= 0; line 141)")
    print(f"  max f      = {f_max:.6f}    (expect <= 1; line 141)")
    print(f"  f-hat(0)   = {fhat0:.12f}   (expect 1/2; line 230)")
    print(f"  M_T(h*)    = {m_together:.10f}  (Together's functional)")
    print(f"  Omega(f)   = {omega_white:.10f}  (White's autocorrelation)")
    print(
        "  Note: M_T and Omega are DIFFERENT functionals; they coincide "
        "only in the µ-limit, not pointwise on a given input."
    )

    # Hard checks for admissibility into White's feasible set.
    tol = 1e-9
    ok_int = abs(int_f - 1.0) < tol
    ok_min = f_min >= -tol
    ok_max = f_max <= 1.0 + tol
    ok_fhat = abs(fhat0 - 0.5) < tol
    all_ok = ok_int and ok_min and ok_max and ok_fhat
    print(
        f"  Admissible into White (int f=1, 0<=f<=1, f-hat(0)=1/2): "
        f"{all_ok} "
        f"(int_f={ok_int}, min>=0={ok_min}, max<=1={ok_max}, fhat0={ok_fhat})"
    )

    return {
        "kind": kind,
        "L": L,
        "int_f": int_f,
        "f_min": f_min,
        "f_max": f_max,
        "fhat0": fhat0,
        "M_together": m_together,
        "Omega_white": omega_white,
        "admissible": all_ok,
    }


if __name__ == "__main__":
    verify_together_value()
    print()
    verify_white_embedding(kind="even")
    print()
    verify_white_embedding(kind="direct")
    print()
    out = save_canonical()
    print(f"Saved canonical JSON to: {out}")
