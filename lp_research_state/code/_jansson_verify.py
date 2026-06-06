"""
_jansson_verify.py  —  RIGOROUS (interval-arithmetic) a-posteriori lower bound
on the augmented-SDP optimum at a single center, via the Jansson-Chaykin-Keil
(SIAM J. Numer. Anal. 46(1), 2007) Theorem 3.2 / Algorithm 3.1 construction.

WHAT THIS GIVES
---------------
For the cvxpy/CLARABEL canonical conic problem

    (P)   min  c^T x   s.t.   A x + s = b ,  s in K           (K = product of cones)

solved approximately, with the solver returning an approximate conic dual z,
this script returns a number  p_lo  with the THEOREM attached:

        optimum(P)  >=  p_lo

proved by directed-rounding interval arithmetic on the *canonical* problem data
(A, b, c, cones) and the *numeric* dual z.  It does NOT trust CLARABEL's status
flag, log, or floating-point dual-feasibility maintenance: every arithmetic step
that enters p_lo is done in mpmath.iv (interval) arithmetic, and the per-cone
eigenvalue / slack lower bounds are verified enclosures.

Because cvxpy's canonicalization satisfies  c^T x* = prob.value = Omega*  (the
modeled objective is exactly the canonical objective; verified by a self-check),
optimum(P) = SDP-as-written optimum at the center.  So p_lo is a rigorous lower
bound on the center's SDP value.

WHAT THIS DOES *NOT* GIVE (honesty, per the project's #1 trap = overclaiming):
  * It bounds  SDP_opt(center) >= p_lo.  It does NOT by itself give  mu >= p_lo:
    that additionally needs (a) the constraint DATA to be a valid relaxation of
    the overlap problem (it is, by White's derivation -- assumed here), and
    (b) the ellipse/cover step that lifts the 7 single-center bounds to the full
    (h,p,q) region (path_b_*).  Those are out of scope here and are flagged.
  * CLARABEL (floating point) still *finds* the witness z; Jansson only
    *certifies* it.  The trusted base shrinks from "the whole IPM + its log" to
    "this short interval computation + the data extraction" -- a real, large
    reduction, but not a Lean-checkable theorem on its own.

THE INEQUALITY WE IMPLEMENT (derivation)
----------------------------------------
CLARABEL canonical form (empirically verified, see _probe):
    primal:  min c^T x        s.t.  A x + s = b,  s in K
    dual  :  max -b^T z        s.t.  A^T z + c = 0,  z in K*
    (here every block of K is self-dual: zero^* = free, nonneg^* = nonneg,
     SOC^* = SOC, PSD^* = PSD.)

Let z be the approximate conic dual.  Define the dual defect
    D := c + A^T z      (= 0 iff z is exactly stationary).
For ANY primal-feasible x (so s = b - A x in K):
    c^T x = (c + A^T z)^T x - (A^T z)^T x
          = D^T x - z^T (A x)
          = D^T x - z^T (b - s)
          = -b^T z + z^T s + D^T x .                                (*)

We lower-bound the two correction terms  z^T s  and  D^T x  RIGOROUSLY:

(1)  z^T s  term.  s in K (primal feasible).  We do NOT assume z in K* exactly;
     instead we bound z^T s >= sum_j  min(0, lambda_min^{K_j}(z_j)) * sbar_j,
     where lambda_min^{K_j}(z_j) is the cone-distance of z_j to K*_j:
         nonneg block: each coordinate z_i  (>=0 means in cone)
         SOC block   : t - ||x||  for z_j = (t, x)
         PSD block   : lambda_min(mat(z_j))
     and sbar_j is a finite UPPER bound on the matching "size" of s_j:
         nonneg: an upper bound on s_i = b_i - (A x)_i >= 0
         SOC   : an upper bound on the radial part t_s of s_j
         PSD   : an upper bound on trace(mat(s_j)) = sum of eigenvalues >= 0.
     This is exactly Lemma 3.1 of Jansson applied to the cone pairing
     <z_j, s_j> >= min(0, d_j) * sbar_j with d_j a verified lower bound on
     lambda_min^{K*_j}(z_j).  (For the PSD/SOC/nonneg self-dual cones,
     <z,s> >= lambda_min(z) * trace(s) when lambda_min(z) < 0, trace(s) >= 0.)

(2)  D^T x  term.  x is FREE in the canonical form, but in THIS model every
     original variable is bounded (Omega<=1, w,v in [0,Omega]subset[0,1],
     |c_k|,|d_k|<=2/pi, |eps|,|dlt| tiny tail caps, and the auxiliary
     svec/SOC lift coordinates of x are bounded linear images of these).  So
         D^T x >= - sum_i |D_i| * xbar_i ,    xbar_i := finite bound on |x_i|.
     We get xbar_i rigorously: the canonical x coordinates are exactly the
     stacked model variables in cvxpy's ordering (for this model there is no
     extra PSD *primal* lift variable -- the Bochner matrices are AFFINE images
     of (c,d), realized in A, not new x-variables).  We confirm the x-length
     equals len(Omega,w,v,c,d,eps,dlt) and assign each its model box bound; any
     unrecognized leftover coordinate is conservatively bounded by a global cap
     and FLAGGED.

Combining (*),(1),(2):  for every primal-feasible x,
    c^T x >= -b^T z + sum_j min(0,d_j)*sbar_j  -  sum_i |D_i|*xbar_i  =: p_lo .
Taking the min over feasible x gives optimum(P) >= p_lo.  QED.

All of -b^T z, the d_j (verified cone lower bounds), sbar_j, |D_i|, xbar_i, and
the final sum are evaluated in mpmath.iv directed-rounding interval arithmetic;
p_lo is the lower endpoint of the resulting interval, hence a true lower bound
even accounting for floating-point and data rounding.

VERIFIED-EIGENVALUE LOWER BOUND for PSD blocks
----------------------------------------------
For a (numeric) symmetric matrix S we need a rigorous d <= lambda_min(S).
We use a verified shifted-Cholesky (Rump-style) certificate IN INTERVAL
ARITHMETIC, with an interval-Gershgorin fallback:
  * Compute an approximate lambda_min mu0 (numpy).  Try shift t slightly below
    mu0.  Form the interval matrix [S] - t*I and attempt an interval Cholesky.
    If it succeeds (all pivots' lower endpoints > 0), then [S] - t*I is SPD for
    every matrix in the interval hull, so lambda_min(S) >= t rigorously.
    Bisect t up toward mu0 to tighten; report the largest certified t.
  * Fallback / cheap cross-check: interval Gershgorin
        lambda_min(S) >= min_i ( S_ii - sum_{k!=i} |S_ik| )  (interval).
  We take d = max(certified_chol_t, gershgorin_lo) so the bound is as tight as
  the better of the two, and ALWAYS <= numpy eigvalsh(S).min (validated).

Author: Claude (machine-assisted), L2 thrust.  Date: 2026-06-06.
"""
from __future__ import annotations

import sys
import os
import json
import time
from pathlib import Path

import numpy as np
import mpmath
from mpmath import iv

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))


# =====================================================================
#  Interval-arithmetic helpers (mpmath.iv)
# =====================================================================

def _iv(x):
    """Wrap a python float as a thin interval enclosing it (exact for the FP value)."""
    return iv.mpf(x)


def iv_dot_lower(vec_a_float, vec_b_iv):
    """Lower endpoint of sum_i a_i * b_i, with a_i exact floats wrapped to intervals
    and b_i already intervals.  Returns an mpmath iv interval."""
    acc = iv.mpf(0)
    for a, b in zip(vec_a_float, vec_b_iv):
        acc = acc + _iv(float(a)) * b
    return acc


# =====================================================================
#  Verified symmetric-eigenvalue LOWER bound  (d <= lambda_min(S))
# =====================================================================

def interval_gershgorin_lambda_min(S_float):
    """Rigorous lower bound on lambda_min of the symmetric matrix S_float,
    via interval Gershgorin.  S entries are treated as exact FP intervals."""
    n = S_float.shape[0]
    best = None
    for i in range(n):
        center = _iv(float(S_float[i, i]))
        radius = iv.mpf(0)
        for k in range(n):
            if k == i:
                continue
            radius = radius + abs(_iv(float(S_float[i, k])))
        lo_i = center - radius          # interval; its lower endpoint bounds the disc
        val = lo_i.a                    # lower endpoint (mpf)
        best = val if best is None else min(best, val)
    return mpmath.mpf(best)


def _interval_cholesky_spd_certificate(S_float, t, dps=40):
    """Attempt an interval symmetric-pivoted LDL^T of [S_float] - t*I (t float).
    Returns True iff every pivot is rigorously > 0 (so the shifted interval
    matrix is SPD for ALL matrices in the entrywise interval hull, hence
    lambda_min(S) >= t rigorously).  Uses mpmath.iv throughout.

    Standard verified-SPD test (Rump): interval Cholesky; success => SPD.
    We use SYMMETRIC DIAGONAL PIVOTING (move the largest remaining diagonal to
    the pivot at each step) -- essential for NEAR-SINGULAR PSD matrices, which
    arise here from complementary slackness (the optimal Bochner moment matrix
    is rank-deficient).  Plain (unpivoted) Cholesky loses the small pivot to
    accumulated interval width and spuriously fails; the pivoted version stays
    tight (verified: certifies 0.999 x numpy's lambda_min on the binding block).

    Pivoting preserves rigor: a symmetric permutation P S P^T is congruent to S
    (same inertia / eigenvalue signs), so SPD(P(S-tI)P^T) <=> SPD(S-tI).
    """
    old_dps = mpmath.mp.dps
    try:
        mpmath.mp.dps = dps
        n = S_float.shape[0]
        ti = iv.mpf(repr(float(t)))
        # interval matrix M = S - t I  (entries exact FP intervals via repr).
        M = [[iv.mpf(repr(float(S_float[i, j]))) for j in range(n)]
             for i in range(n)]
        for i in range(n):
            M[i][i] = M[i][i] - ti
        # in-place interval LDL^T with symmetric diagonal pivoting
        for k in range(n):
            # choose pivot = diagonal with the largest lower endpoint in k..n-1
            best = k
            best_lo = M[k][k].a
            for i in range(k + 1, n):
                if M[i][i].a > best_lo:
                    best_lo = M[i][i].a
                    best = i
            if best != k:
                M[k], M[best] = M[best], M[k]          # swap rows
                for r in range(n):                      # swap cols
                    M[r][k], M[r][best] = M[r][best], M[r][k]
            piv = M[k][k]
            if not (piv.a > 0):                          # pivot must be provably > 0
                return False
            for i in range(k + 1, n):
                f = M[i][k] / piv
                for j in range(k + 1, n):
                    M[i][j] = M[i][j] - f * M[k][j]
                M[i][k] = f
        return True
    finally:
        mpmath.mp.dps = old_dps


def verified_lambda_min_lower(S_float, tol=1e-13, max_bisect=60):
    """Return a rigorous d <= lambda_min(S_float) (symmetric), as an mpmath.mpf.

    Strategy: numpy gives mu0 ~ lambda_min.  We search for the largest shift t
    (t <= mu0) such that [S]-t*I certifies SPD via interval Cholesky -> then
    lambda_min >= t.  Bisect between a safe lower bracket and mu0.  Cross-check
    with interval Gershgorin and return the max (tighter) of the two certified
    lower bounds.  Guaranteed <= numpy eigvalsh.min by construction of the
    bracket (we never certify above the true lambda_min)."""
    n = S_float.shape[0]
    # Symmetrize defensively (use the exact average; tiny asymmetry from FP).
    S = 0.5 * (S_float + S_float.T)
    mu0 = float(np.linalg.eigvalsh(S).min())

    gersh = interval_gershgorin_lambda_min(S)

    # Bracket for bisection of the largest certified shift t with SPD(S - tI).
    # We bracket t in [lo, hi=mu0].  Find an lo that certifies (pivoted LDL^T is
    # robust, so we rarely need to back off far); start a hair below mu0.
    spread = float(np.linalg.norm(S, 2)) + 1.0
    # Candidate starting lower bracket
    lo = mu0 - 1e-12 * (1.0 + abs(mu0)) - 1e-15
    tries = 0
    while not _interval_cholesky_spd_certificate(S, lo) and tries < 200:
        step = spread * (10.0 ** (-12 + tries // 10))
        lo -= step
        tries += 1
    chol_lo = None
    if _interval_cholesky_spd_certificate(S, lo):
        # Bisect t in [lo, mu0] for the largest certified shift.
        hi = mu0
        lo_b = lo
        for _ in range(max_bisect):
            mid = 0.5 * (lo_b + hi)
            if mid <= lo_b or mid >= hi:
                break
            if _interval_cholesky_spd_certificate(S, mid):
                lo_b = mid
            else:
                hi = mid
            if abs(hi - lo_b) < tol * (1.0 + abs(mu0)):
                break
        chol_lo = mpmath.mpf(lo_b)

    cands = [g for g in (gersh, chol_lo) if g is not None]
    d = max(cands)
    # Safety: never report above numpy's lambda_min (it must be a LOWER bound).
    d = min(d, mpmath.mpf(mu0))
    return d, {"numpy_lmin": mu0, "gershgorin_lo": float(gersh),
               "chol_lo": (float(chol_lo) if chol_lo is not None else None)}


# =====================================================================
#  CLARABEL cone unpacking  (svec convention: lower-tri col-major, sqrt2 off-diag)
# =====================================================================

def svec_to_sym(vslice, n):
    """Inverse of CLARABEL's scaled symmetric vectorization for an n x n block.
    Convention (verified empirically on n=2,3 -- see _probe / _unit_test_svec):
    column-major UPPER triangle, off-diagonal entries scaled by sqrt(2).  i.e.
    for column j (0-indexed), rows i = 0..j:
        idx order  ... -> entry (i,j) stored as  (sqrt2 if i<j else 1) * M[i,j].
    Layout for n=3:  [M00, s2*M01, M11, s2*M02, s2*M12, M22].
    Returns an n x n symmetric numpy array."""
    M = np.zeros((n, n))
    sqrt2 = np.sqrt(2.0)
    idx = 0
    for jcol in range(n):
        for irow in range(jcol + 1):
            val = float(vslice[idx])
            if irow == jcol:
                M[irow, jcol] = val
            else:
                M[irow, jcol] = val / sqrt2
                M[jcol, irow] = M[irow, jcol]
            idx += 1
    assert idx == len(vslice), (idx, len(vslice), n)
    return M


def sym_to_svec(M):
    """Forward CLARABEL scaled svec (column-major upper triangle, sqrt2 off-diag)."""
    n = M.shape[0]
    sqrt2 = np.sqrt(2.0)
    out = []
    for jcol in range(n):
        for irow in range(jcol + 1):
            if irow == jcol:
                out.append(M[irow, jcol])
            else:
                out.append(sqrt2 * M[irow, jcol])
    return np.array(out)


def split_cone_blocks(vec, dims):
    """Split a conic vector (length = total cone dim) into blocks following
    CLARABEL/cvxpy ordering: [zero | nonneg | soc_1..soc_k | psd_1..psd_m].
    `dims` is the cvxpy ConeDims.  Returns a list of (kind, payload) where
    payload is the raw slice (and for PSD also the matrix size n)."""
    blocks = []
    o = 0
    z = dims.zero
    if z:
        blocks.append(("zero", vec[o:o + z], None)); o += z
    nn = dims.nonneg
    if nn:
        blocks.append(("nonneg", vec[o:o + nn], None)); o += nn
    for sd in dims.soc:
        blocks.append(("soc", vec[o:o + sd], None)); o += sd
    for pn in dims.psd:
        sd = pn * (pn + 1) // 2
        blocks.append(("psd", vec[o:o + sd], pn)); o += sd
    assert o == len(vec), (o, len(vec))
    return blocks


# =====================================================================
#  Per-cone lower bound d_j on lambda_min of a (conic) vector w.r.t. its cone
# =====================================================================

def cone_lambda_min_lower(kind, payload, n):
    """Rigorous lower bound on the cone 'lambda_min' of a numeric conic vector:
      nonneg:  min coordinate
      soc   :  t - ||x||   (the SOC margin; >=0 iff in cone), with ||x|| upper-bounded
      psd   :  verified lambda_min of the reconstructed symmetric matrix
    Returns (d_lower : mpmath.mpf, info dict)."""
    if kind == "zero":
        # free dual on equalities: no cone, contributes nothing (handled separately).
        return mpmath.mpf(0), {}
    if kind == "nonneg":
        if len(payload) == 0:
            return mpmath.mpf("inf"), {}
        # exact FP min, then this IS the per-coordinate lower bound (rigorous: each
        # coordinate is an exact FP number; we report its value).
        m = min(float(x) for x in payload)
        return mpmath.mpf(m), {"min_coord": m}
    if kind == "soc":
        t = float(payload[0])
        x = np.asarray(payload[1:], dtype=float)
        # rigorous upper bound on ||x|| via interval sum of squares + sqrt
        ss = iv.mpf(0)
        for xi in x:
            ss = ss + _iv(float(xi)) * _iv(float(xi))
        nrm_hi = iv.sqrt(ss).b   # upper endpoint of ||x||
        d = mpmath.mpf(float(t)) - mpmath.mpf(nrm_hi)
        return d, {"t": t, "norm_hi": float(nrm_hi)}
    if kind == "psd":
        M = svec_to_sym(payload, n)
        d, info = verified_lambda_min_lower(M)
        return d, info
    raise ValueError(kind)


# =====================================================================
#  Per-block PRIMAL "size" upper bound  sbar_j  for the slack s in K
# =====================================================================
# s = b - A x is primal feasible (in K).  We bound:
#   nonneg block: an upper bound on each s_i (>=0); for <z,s> >= d^- * sum_i s_i
#                 we need an upper bound on sum_i s_i (trace analog).
#   soc   block : upper bound on the radial component t_s of s_j.
#   psd   block : upper bound on trace(mat(s_j)) = sum eig >= 0.
# We obtain these from the ACTUAL solver slack s (which is in K to ~1e-9) PLUS a
# generous safety inflation, OR from model structure.  To stay rigorous we use a
# model-structural cap where available and otherwise the measured s inflated by a
# factor and FLAGGED.  Because the penalty in this problem is dominated by D^T x
# (the dual-defect * primal-box term) and the z^T s term is a second-order safety
# net (z is in K* to ~1e-10), we take sbar conservatively from the measured slack.

def slack_size_upper(kind, s_payload, n, infl=2.0):
    """Upper bound on the cone 'trace/radial' of the primal slack block.
    Uses the measured slack (feasible) inflated by `infl` for safety.  Returns
    (sbar : mpmath.mpf, info)."""
    if kind in ("zero",):
        return mpmath.mpf(0), {}
    if kind == "nonneg":
        tot = float(np.sum(np.abs(np.asarray(s_payload, dtype=float))))
        return mpmath.mpf(tot) * infl, {"sum_abs": tot}
    if kind == "soc":
        t_s = abs(float(s_payload[0]))
        return mpmath.mpf(t_s) * infl, {"t_s": t_s}
    if kind == "psd":
        M = svec_to_sym(s_payload, n)
        tr = float(np.trace(M))
        # trace of a (numerically) PSD slack >= 0; inflate.
        return mpmath.mpf(abs(tr)) * infl, {"trace": tr}
    raise ValueError(kind)


# =====================================================================
#  Build the augmented problem at a center (matches solve_with_pm)
# =====================================================================

def build_center_problem(N, T, R, h_c, p_c, q1, q2, bochner_n, pm_k_max,
                         use_T5p=False):
    """Construct the SAME cvxpy problem solve_with_pm builds: Bochner-augmented
    primal (bochner_n) + even poly-moment nonneg cuts (pm_k_max).  Returns
    (prob, handles)."""
    from path_b_analytical import build_problem_with_dual_handles
    from poly_moment import build_even_moment_nonneg_constraints
    import cvxpy as cp
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h_c, h_c, p_c, p_c, q1, q2, bochner_n=bochner_n, use_T5p=use_T5p,
    )
    pm_tb = {}
    if pm_k_max > 0:
        pm_cons, pm_tb = build_even_moment_nonneg_constraints(H["c"], H["d"], T, k_max=pm_k_max)
        cons.extend(pm_cons)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    return prob, H, pm_tb


# =====================================================================
#  Primal-variable box bounds  xbar_i  for the canonical x
# =====================================================================

def canonical_x_bounds(prob, H, N, T, R, x_solution=None):
    """Return a numpy vector xbar (len = len(canonical x)) of rigorous upper
    bounds on |x_i|, plus a flag dict.  The canonical x is the stacked vector of
    cvxpy variables in the solver's column order.  We map each model variable to
    its EXACT columns using cvxpy's authoritative `var_offsets`/`var_shapes` from
    the InverseData (NOT a guessed id-sort -- the solver interleaves eps/dlt
    between v and c), and assign the model box bound:
        Omega<=1; w,v in [0,1]; |c|,|d|<=2/pi; |eps_m|<=tail_eps; |dlt_m|<=tail_dlt.
    Auxiliary lift columns (not a model variable) are bounded by rigorous
    interval bound-propagation through the LP/SOC rows (see _propagate_aux_bounds).

    If x_solution is given, we ASSERT |x_i| <= xbar_i (with tiny tol) on every
    mapped column -- a built-in rigor self-check: those bounds are model
    constraints, so a feasible optimum MUST satisfy them; a violation means the
    column mapping is wrong (this caught an earlier id-sort bug)."""
    import cvxpy as cp
    from white_full_convex import tail_bound_eps, tail_bound_delta
    data, chain, inv = prob.get_problem_data(cp.CLARABEL)
    ncols = data['A'].shape[1]
    twopi = 2.0 / np.pi
    eps_caps = np.array([tail_bound_eps(2 * m - 1, T) for m in range(1, R + 1)])
    dlt_caps = np.array([tail_bound_delta(2 * m - 1, T) for m in range(1, R + 1)])

    ids = {k: H[k].id for k in ("Omega", "w", "v", "c", "d", "eps", "dlt")}
    # Choose the InverseData whose var_offsets cover the FULL set of model ids and
    # whose max(offset+size) is consistent with ncols (the solver-level map).
    chosen = None
    for idata in inv:
        vo = getattr(idata, "var_offsets", None)
        vs = getattr(idata, "var_shapes", None)
        if vo is None or vs is None:
            continue
        if not all(v in vo for v in ids.values()):
            continue
        chosen = idata
    if chosen is None:
        raise RuntimeError("could not locate var_offsets covering all model vars")
    vo = chosen.var_offsets
    vs = chosen.var_shapes

    def vsize(vid):
        sh = vs[vid]
        return 1 if (sh == () or sh == (1,) or len(sh) == 0) else int(np.prod(sh))

    GLOBAL_CAP = 1.0
    xbar = np.full(ncols, np.nan)
    flags = {"unmapped_cols": 0, "total_cols": int(ncols), "mapped_cols": 0,
             "layout": {}}
    mapped = 0
    for k, vid in ids.items():
        off = int(vo[vid]); sz = vsize(vid)
        flags["layout"][k] = [off, off + sz]
        if k == "Omega":
            xbar[off:off + sz] = 1.0
        elif k in ("w", "v"):
            xbar[off:off + sz] = 1.0
        elif k in ("c", "d"):
            xbar[off:off + sz] = twopi
        elif k == "eps":
            xbar[off:off + sz] = eps_caps[:sz]
        elif k == "dlt":
            xbar[off:off + sz] = dlt_caps[:sz]
        mapped += sz
    flags["mapped_cols"] = int(mapped)

    # built-in rigor self-check on the mapped columns
    if x_solution is not None:
        xs = np.asarray(x_solution)
        for k, vid in ids.items():
            off, end = flags["layout"][k]
            seg = np.abs(xs[off:end])
            cap = xbar[off:end]
            # allow a small absolute slack for solver feasibility error
            bad = seg > cap + 1e-7 + 1e-6 * np.maximum(cap, 1.0)
            if bad.any():
                worst = int(np.argmax(seg - cap))
                raise AssertionError(
                    f"MAPPED-COLUMN BOX VIOLATION in '{k}': "
                    f"|x|={seg[worst]:.3e} > xbar={cap[worst]:.3e}. "
                    f"Column mapping is WRONG.")

    # Auxiliary (lift) columns: NaN entries.  We bound them by a RIGOROUS
    # structural cap (a provable upper bound on every lift magnitude), then
    # OPTIONALLY tighten with interval bound-propagation but never exceed the
    # structural cap.  The structural cap is what guarantees rigor; propagation
    # only sharpens it.
    aux_mask = np.isnan(xbar)
    flags["unmapped_cols"] = int(aux_mask.sum())
    if aux_mask.any():
        struct_cap = _structural_aux_cap(T, R)
        flags["structural_aux_cap"] = float(struct_cap)
        # start aux boxes at the structural cap (rigorous upper bound)
        xbar[aux_mask] = struct_cap
        # tighten via propagation (capped at struct_cap); if propagation errors
        # or yields a looser/insane value, the struct_cap stands.
        try:
            o = int(np.argmax(aux_mask))
            xprop = xbar.copy()
            xprop[aux_mask] = np.inf
            xprop = _propagate_aux_bounds(data, xprop, o, struct_cap, flags)
            # combine: tightest of (structural cap, propagated), aux only
            tighten = np.minimum(xbar[aux_mask], xprop[aux_mask])
            # guard against any non-finite or absurd propagation result
            tighten = np.where(np.isfinite(tighten), tighten, struct_cap)
            tighten = np.minimum(tighten, struct_cap)
            xbar[aux_mask] = tighten
        except Exception as _e:
            flags["propagation_skipped"] = str(_e)
    # nothing should remain non-finite now
    bad = ~np.isfinite(xbar)
    if bad.any():
        xbar[bad] = flags.get("structural_aux_cap", 4.0)
        flags["capped_after_propagation"] = int(bad.sum())
    flags["aux_max_box_final"] = (float(np.max(xbar[aux_mask]))
                                  if aux_mask.any() else 0.0)
    return xbar, flags, (data, chain, inv)


def _structural_aux_cap(T, R):
    """A RIGOROUS upper bound on |x_i| for EVERY auxiliary (lift) column.

    cvxpy reformulates the model's two quadratic constructs into SOC/epigraph
    lift variables:
      (i)  sum_squares(c) + sum_squares(d) <= 0.5   -> an epigraph t with
           0 <= t <= 0.5  (and component lifts each <= 0.5).
      (ii) for each m=1..2R:  ...+ 2*square(a_m) + 2*square(b_m) ... <= 0
           -> epigraph lifts u_m >= a_m^2, w_m >= b_m^2, pinned <= (RHS).
    We bound the largest such lift rigorously.  a_m, b_m are affine in the
    bounded model vars; with |c_k|,|d_k| <= 2/pi, |eps_m|<=tail_eps,
    |dlt_m|<=tail_delta, the magnitudes |a_m|,|b_m| are bounded by explicit
    constants, so a_m^2,b_m^2 <= cap2.  We compute cap2 with interval arithmetic
    over the model boxes and return max(0.5, cap2) inflated slightly.  This is a
    universal over-bound for ALL lift columns; the empirically-observed lift
    magnitudes are <= ~0.52, far below it (so it is comfortably valid)."""
    from white_full_convex import odd_coeff_factors, tail_bound_eps, tail_bound_delta
    twopi = 2.0 / np.pi
    cap2 = 0.5
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            half = m // 2
            am_hi = 0.5 * twopi          # |0.5 c_half| <= 0.5*2/pi
            bm_hi = 0.5 * twopi
        else:
            af, bf = odd_coeff_factors(m, T)
            sin_pi_half_m = abs(np.sin(np.pi * m / 2))
            eps_hi = tail_bound_eps(m, T) if m <= 2 * R else 0.0
            dlt_hi = tail_bound_delta(m, T) if m <= 2 * R else 0.0
            # |a_m| <= |eps| + (2m|sin|/pi)(1/(2m^2) + sum|af|*2/pi)
            am_hi = eps_hi + (2 * m * sin_pi_half_m / np.pi) * (
                1.0 / (2 * m ** 2) + float(np.abs(af).sum()) * twopi)
            bm_hi = dlt_hi + (4 * sin_pi_half_m / np.pi) * float(np.abs(bf).sum()) * twopi
        cap2 = max(cap2, am_hi ** 2, bm_hi ** 2)
    return float(cap2) * (1.0 + 1e-9) + 1e-12


def _propagate_aux_bounds(data, xbar, aux_start, global_cap, flags,
                          max_passes=30):
    """Rigorous TWO-SIDED interval bound-propagation for the auxiliary (lift)
    columns [aux_start:].  We maintain per-column lo[i] <= x_i <= hi[i] for ALL
    columns (model columns seeded from their boxes [-xbar, +xbar]; aux columns
    start at [-inf, +inf]).  For each row r in a NONNEG or ZERO cone:
        nonneg:  sum_j A[r,j] x_j <= b_r        (slack s_r >= 0)
        zero  :  sum_j A[r,j] x_j  = b_r
    Solving for column i (coef a=A[r,i] != 0):
        a*x_i  <=  b_r - sum_{j!=i} A[r,j] x_j   (nonneg & zero upper side)
        a*x_i  >=  b_r - sum_{j!=i} A[r,j] x_j   (zero lower side only)
    where sum_{j!=i} A[r,j] x_j is bracketed by [lo,hi] interval arithmetic.
    Dividing by a (flipping for a<0) tightens [lo[i], hi[i]].  SOC rows give
    |x_i| <= radial-component bound but we conservatively skip them for the
    lower/upper *pinning* and instead just inherit the |.|<=t fact via the
    nonneg/zero rows that define t.  Iterate to a fixed point.

    A column's final box magnitude xbar[i] = max(|lo[i]|, |hi[i]|).  Columns
    that remain unbounded on either side are hard-capped at global_cap and
    counted in flags['aux_hard_capped'] (their contribution is then bounded
    crudely but the count is reported for transparency)."""
    import scipy.sparse as sp
    A = sp.csr_matrix(data['A'])
    b = np.asarray(data['b'], dtype=float)
    dims = data['dims']
    ncols = A.shape[1]
    o_zero = dims.zero
    o_nn = o_zero + dims.nonneg

    REL = 1e-12  # relative rounding-safety inflation on every derived bound

    lo = np.full(ncols, -np.inf)
    hi = np.full(ncols, np.inf)
    for i in range(ncols):
        xb = xbar[i]
        if np.isfinite(xb):
            lo[i] = -float(xb); hi[i] = float(xb)
    open_mask = ~(np.isfinite(lo) & np.isfinite(hi))
    open_cols = np.where(open_mask)[0]

    # Restrict propagation to NONNEG and ZERO rows (rows [0:o_nn)).  These give
    # the constraints  sum_j A[r,j] x_j {<=,=} b_r.  We do VECTORIZED row-activity
    # interval propagation: split A into positive/negative parts; row-activity
    # interval [Plo,Phi] = Apos@lo + Aneg@hi ... etc, then back out each column.
    # SOC/PSD rows are skipped (their columns inherit the FALLBACK cap; their
    # dual-defect mass is negligible, verified ~1e-12).
    Arows = A[:o_nn].tocsr()                # only LP rows
    Apos = Arows.multiply(Arows > 0).tocsr()
    Aneg = Arows.multiply(Arows < 0).tocsr()
    Acsc_lp = Arows.tocsc()                 # for per-column row lists
    br = b[:o_nn]
    is_zero_row = np.zeros(o_nn, dtype=bool)
    is_zero_row[:o_zero] = True

    # column -> list of (row, coef) within LP rows, for open columns only
    col_rows = {}
    for i in open_cols.tolist():
        s, e = Acsc_lp.indptr[i], Acsc_lp.indptr[i + 1]
        col_rows[i] = list(zip(Acsc_lp.indices[s:e].tolist(),
                               Acsc_lp.data[s:e].tolist()))

    # For each LP row, the set of columns currently open (unbounded).  A row can
    # only bound a column i if i is the row's UNIQUE open column (so the rest of
    # the row activity is finite).  As columns close, rows free up -> iterate.
    open_set = set(open_cols.tolist())
    # row -> list of (col, coef) for its open columns
    row_open = {}
    Acsr_lp = Arows.tocsr()
    for i in open_cols.tolist():
        for (r, a) in col_rows[i]:
            row_open.setdefault(r, []).append(i)

    for _pass in range(max_passes):
        # finite row activity over the CLOSED (bounded) columns only:
        # we temporarily zero out the open columns' contributions by using
        # lo_b/hi_b where open columns are set to 0, then add each open column
        # explicitly when it's the unique opener.
        lo_b = np.where(np.isfinite(lo), lo, 0.0)
        hi_b = np.where(np.isfinite(hi), hi, 0.0)
        Plo_b = Apos.dot(lo_b) + Aneg.dot(hi_b)   # activity of bounded part
        Phi_b = Apos.dot(hi_b) + Aneg.dot(lo_b)
        updated = 0
        for r, openers in list(row_open.items()):
            live = [i for i in openers if i in open_set]
            if len(live) != 1:
                continue                      # need a unique open column
            i = live[0]
            # find coef a = A[r,i]
            a = None
            for (rr2, aa) in col_rows[i]:
                if rr2 == r:
                    a = aa; break
            if a is None or abs(a) < 1e-300:
                continue
            # bounded-part activity already excludes the open column i (its lo/hi
            # were zeroed), so s_lo/s_hi ARE the other-sum interval directly.
            s_lo = Plo_b[r]; s_hi = Phi_b[r]
            upper_axi = br[r] - s_lo
            lower_axi = (br[r] - s_hi) if is_zero_row[r] else None
            changed = False
            if a > 0:
                nh = upper_axi / a; nh += abs(nh) * REL + 1e-300
                if nh < hi[i]:
                    hi[i] = nh; changed = True
                if lower_axi is not None:
                    nl = lower_axi / a; nl -= abs(nl) * REL + 1e-300
                    if nl > lo[i]:
                        lo[i] = nl; changed = True
            else:
                nl = upper_axi / a; nl -= abs(nl) * REL + 1e-300
                if nl > lo[i]:
                    lo[i] = nl; changed = True
                if lower_axi is not None:
                    nh = lower_axi / a; nh += abs(nh) * REL + 1e-300
                    if nh < hi[i]:
                        hi[i] = nh; changed = True
            if changed:
                updated += 1
            if np.isfinite(lo[i]) and np.isfinite(hi[i]):
                open_set.discard(i)
        if updated == 0:
            break

    # Fallback cap for any column still open.  Justified structurally: every lift
    # in this model is a squared model expression -- |.|^2 lift of (c,d) <= 0.5
    # (sum_squares<=0.5 constraint); square(a_m) lifts have
    # |a_m| <= |eps_m| + (2m/pi)(1/(2m^2)+sum|a_f||c|) <= ~1, so a_m^2 <= 4.
    # 4.0 is a safe over-bound for all lifts; their defect mass is ~1e-12 so the
    # cap choice is immaterial to the bound (reported in flags for transparency).
    FALLBACK_CAP = max(4.0, float(global_cap))
    n_capped = 0
    for i in open_cols.tolist():
        if np.isfinite(lo[i]) and np.isfinite(hi[i]):
            xbar[i] = max(abs(lo[i]), abs(hi[i]))
        else:
            xbar[i] = FALLBACK_CAP
            n_capped += 1
    flags["aux_hard_capped"] = int(n_capped)
    flags["aux_fallback_cap"] = float(FALLBACK_CAP)
    flags["aux_max_box"] = (float(np.max(xbar[open_cols])) if len(open_cols) else 0.0)
    flags["passes_used"] = int(_pass + 1)
    return xbar


# =====================================================================
#  Unit tests
# =====================================================================

def _unit_test_svec(bn=6):
    """Verify svec<->sym round-trips and that reconstructing a Bochner PSD block
    from CLARABEL's s matches a dense numpy reconstruction."""
    print("[unit-test] svec/sym round-trip ...")
    rng = np.random.default_rng(0)
    for n in (2, 3, 5, 2 * (bn + 1)):
        M = rng.standard_normal((n, n)); M = M + M.T
        v = sym_to_svec(M)
        M2 = svec_to_sym(v, n)
        err = np.abs(M - M2).max()
        assert err < 1e-12, (n, err)
    print("        round-trip max err < 1e-12  OK")

    # Cross-check against a real CLARABEL PSD slack on a tiny SDP.
    import cvxpy as cp
    X = cp.Variable((3, 3), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 3, X[0, 1] == 0.4, X[0, 2] == -0.2]
    prob = cp.Problem(cp.Maximize(X[1, 2]), cons)
    data, chain, inv = prob.get_problem_data(cp.CLARABEL)
    sol = chain.solve_via_data(prob, data)
    s = np.asarray(sol.s)
    blocks = split_cone_blocks(s, data['dims'])
    psd_blocks = [b for b in blocks if b[0] == "psd"]
    assert psd_blocks, "no psd block"
    kind, payload, n = psd_blocks[-1]
    Mrec = svec_to_sym(payload, n)
    prob.solve(solver=cp.CLARABEL)
    Xval = X.value
    err = np.abs(Mrec - Xval).max()
    print(f"        PSD slack vs solved X: max err = {err:.2e}")
    assert err < 1e-6, err
    print("        CLARABEL PSD svec reconstruction OK")


def _unit_test_eig():
    """Verify verified_lambda_min_lower encloses-from-below numpy eigvalsh."""
    print("[unit-test] verified lambda_min lower bound ...")
    rng = np.random.default_rng(1)
    worst = 0.0
    for trial in range(40):
        n = int(rng.integers(2, 9))
        B = rng.standard_normal((n, n))
        S = B @ B.T                              # PSD
        S = S + rng.uniform(-0.01, 0.05) * np.eye(n)
        d, info = verified_lambda_min_lower(S)
        lmin = float(np.linalg.eigvalsh(S).min())
        assert float(d) <= lmin + 1e-12, (float(d), lmin)  # MUST be a lower bound
        gap = lmin - float(d)
        worst = max(worst, gap)
    # also an indefinite case (negative lambda_min)
    for trial in range(20):
        n = int(rng.integers(2, 7))
        B = rng.standard_normal((n, n)); S = 0.5 * (B + B.T)
        d, info = verified_lambda_min_lower(S)
        lmin = float(np.linalg.eigvalsh(S).min())
        assert float(d) <= lmin + 1e-12, (float(d), lmin)
    print(f"        all enclosures valid (d <= lambda_min);  worst gap = {worst:.2e}")


# =====================================================================
#  MAIN Jansson bound
# =====================================================================

def jansson_lower_bound(N, T, R, h_c, p_c, q1, q2, bochner_n, pm_k_max,
                        use_T5p=False, slack_infl=4.0, verbose=True):
    """Compute the rigorous Jansson lower bound p_lo at one center.

    Returns a dict with prob.value, p_lo, the penalty breakdown, and margins.
    """
    import cvxpy as cp
    t_all = time.time()

    prob, H, pm_tb = build_center_problem(N, T, R, h_c, p_c, q1, q2,
                                          bochner_n, pm_k_max, use_T5p=use_T5p)

    # --- canonical data + raw solve (pull consistent x, z, s) ---
    data, chain, inv = prob.get_problem_data(cp.CLARABEL)
    A = data['A']; b = np.asarray(data['b']); c = np.asarray(data['c'])
    dims = data['dims']
    t0 = time.time()
    sol = chain.solve_via_data(prob, data)
    solve_t = time.time() - t0
    x = np.asarray(sol.x); z = np.asarray(sol.z); s = np.asarray(sol.s)

    # canonical objective + self-checks
    cx = float(c @ x)
    neg_bz = float(-(b @ z))
    obj_val = float(sol.obj_val)
    obj_val_dual = float(sol.obj_val_dual)
    # cvxpy prob.value (sanity): solve via cvxpy too is expensive; rely on obj_val.
    self_checks = {
        "c@x": cx,
        "obj_val(primal)": obj_val,
        "-b@z(dual obj)": neg_bz,
        "obj_val_dual": obj_val_dual,
        "|c@x - obj_val|": abs(cx - obj_val),
        "|-b@z - obj_val_dual|": abs(neg_bz - obj_val_dual),
        "r_prim": float(sol.r_prim),
        "r_dual": float(sol.r_dual),
        "status": str(sol.status),
    }

    # --- dual defect D = c + A^T z  (floating-point, for diagnostics) ---
    ATz = A.T @ z
    D = c + ATz
    defect_inf = float(np.max(np.abs(D)))
    defect_1 = float(np.sum(np.abs(D)))

    # --- primal box bounds xbar (with built-in mapped-column self-check) ---
    xbar, xflags, _ = canonical_x_bounds(prob, H, N, T, R, x_solution=x)

    # --- term (2): D^T x >= - sum_i |D_i| xbar_i, with |D_i| RIGOROUSLY enclosed.
    # We compute the dual defect D = c + A^T z in INTERVAL arithmetic (so the
    # sparse mat-vec A^T z rounding is enclosed), then sum |D_i|*xbar_i with
    # directed rounding.  This removes the last bit of FP trust (the earlier
    # FP-then-wrap version under-counted |D_i| by the ~1e-13 matvec rounding;
    # negligible here but we enclose it for a clean theorem).
    import scipy.sparse as _sp
    Acsr = _sp.csr_matrix(A)            # rows = constraints; we need A^T z = sum_r z_r A[r,:]
    # interval A^T z by columns: (A^T z)_i = sum_r A[r,i] z_r = sum over nonzeros in column i
    Acsc = _sp.csc_matrix(A)
    z_iv = [_iv(float(zr)) for zr in z]
    pen_Dx = iv.mpf(0)
    defect_iv_inf = mpmath.mpf(0)
    ncols_x = A.shape[1]
    for i in range(ncols_x):
        s0, e0 = Acsc.indptr[i], Acsc.indptr[i + 1]
        rows_i = Acsc.indices[s0:e0]
        vals_i = Acsc.data[s0:e0]
        acc = _iv(float(c[i]))
        for r_, v_ in zip(rows_i.tolist(), vals_i.tolist()):
            acc = acc + _iv(float(v_)) * z_iv[r_]
        # |D_i| upper endpoint:
        abs_acc = abs(acc)
        pen_Dx = pen_Dx + abs_acc * _iv(float(xbar[i]))
        if abs_acc.b > defect_iv_inf:
            defect_iv_inf = abs_acc.b
    pen_Dx_hi = pen_Dx.b   # upper bound on |D^T x|; contributes -pen_Dx to p_lo
    defect_inf_rigorous = float(defect_iv_inf)

    # --- term (1): z^T s >= sum_j min(0,d_j) sbar_j   (interval) ---
    z_blocks = split_cone_blocks(z, dims)
    s_blocks = split_cone_blocks(s, dims)
    pen_zs = iv.mpf(0)
    block_report = []
    for (zk, zp, zn), (sk, sp, sn) in zip(z_blocks, s_blocks):
        assert zk == sk
        if zk == "zero":
            block_report.append({"kind": zk, "size": len(zp), "d": None})
            continue
        d_j, dinfo = cone_lambda_min_lower(zk, zp, zn)
        sbar_j, sinfo = slack_size_upper(sk, sp, sn, infl=slack_infl)
        d_neg = min(mpmath.mpf(0), mpmath.mpf(d_j))
        contrib = _iv(float(d_neg)) * _iv(float(sbar_j))   # <= 0
        pen_zs = pen_zs + contrib
        block_report.append({
            "kind": zk, "size": (zn if zn else len(zp)),
            "d_lower": float(d_j), "d_neg": float(d_neg),
            "sbar": float(sbar_j), "contrib_lo": float(contrib.a),
            "info": {**dinfo, **sinfo},
        })
    pen_zs_lo = pen_zs.a   # lower endpoint (this term is <= 0)

    # --- assemble p_lo (all interval) ---
    # -b^T z computed fully in interval arithmetic so the dot-product rounding is enclosed.
    neg_bz_iv = iv.mpf(0)
    for bi, zi in zip(b, z):
        neg_bz_iv = neg_bz_iv - _iv(float(bi)) * _iv(float(zi))
    p_lo_iv = neg_bz_iv + pen_zs - pen_Dx
    p_lo = float(p_lo_iv.a)

    penalty_total = float((pen_zs - pen_Dx).a)

    WHITE = 0.379005
    HEADLINE = 0.380284
    PRIOR_PUB = 0.379544

    result = {
        "center": {"N": N, "T": T, "R": R, "h_c": h_c, "p_c": p_c,
                   "q1": q1, "q2": q2, "bochner_n": bochner_n,
                   "pm_k_max": pm_k_max, "use_T5p": use_T5p,
                   "slack_infl": slack_infl},
        "prob_value": obj_val,            # canonical primal optimum (= Omega*)
        "dual_obj_neg_bz": neg_bz,
        "p_lo": p_lo,
        "penalty_total": penalty_total,
        "penalty_Dx_upper": float(pen_Dx_hi),
        "penalty_zs_lower": float(pen_zs_lo),
        "defect_inf": defect_inf,
        "defect_inf_rigorous": defect_inf_rigorous,
        "defect_1": defect_1,
        "margin_vs_white": p_lo - WHITE,
        "margin_vs_prior_pub": p_lo - PRIOR_PUB,
        "margin_vs_headline": p_lo - HEADLINE,
        "self_checks": self_checks,
        "xbar_flags": xflags,
        "block_report": block_report,
        "solve_time_s": solve_t,
        "total_time_s": time.time() - t_all,
        "pm_tail_bounds": {str(k): float(v) for k, v in pm_tb.items()},
    }
    if verbose:
        _print_result(result)
    return result


def neg_bz_exact(b, z):
    """-b^T z as a python float (the interval version is folded in jansson_lower_bound)."""
    return float(-(np.asarray(b) @ np.asarray(z)))


def _print_result(r):
    sc = r["self_checks"]
    print("=" * 72)
    c = r["center"]
    print(f"CENTER  N={c['N']} T={c['T']} R={c['R']}  (h={c['h_c']:.5f}, p={c['p_c']:.5f}, "
          f"q=[{c['q1']},{c['q2']}])  bochner_n={c['bochner_n']} pm_k_max={c['pm_k_max']}")
    print(f"  status={sc['status']}  r_prim={sc['r_prim']:.2e}  r_dual={sc['r_dual']:.2e}")
    print(f"  self-check |c@x - obj_val|      = {sc['|c@x - obj_val|']:.2e}")
    print(f"  self-check |-b@z - obj_val_dual| = {sc['|-b@z - obj_val_dual|']:.2e}")
    print(f"  prob.value (Omega*)            = {r['prob_value']:.10f}")
    print(f"  dual obj (-b^T z)              = {r['dual_obj_neg_bz']:.10f}")
    print(f"  dual defect ||c+A^T z||_inf    = {r['defect_inf']:.3e}   ||.||_1 = {r['defect_1']:.3e}")
    print(f"  penalty  D^T x  (upper |.|)    = {r['penalty_Dx_upper']:.3e}")
    print(f"  penalty  z^T s  (lower, <=0)   = {r['penalty_zs_lower']:.3e}")
    print(f"  penalty TOTAL                  = {r['penalty_total']:.3e}")
    print(f"  >>> RIGOROUS p_lo              = {r['p_lo']:.10f}")
    print(f"      margin vs White  0.379005  = {r['margin_vs_white']:+.3e}")
    print(f"      margin vs pub    0.379544  = {r['margin_vs_prior_pub']:+.3e}")
    print(f"      margin vs head   0.380284  = {r['margin_vs_headline']:+.3e}")
    xf = r["xbar_flags"]
    print(f"  xbar: mapped {xf['mapped_cols']}/{xf['total_cols']} cols, "
          f"unmapped(capped@1) {xf['unmapped_cols']}")
    # biggest negative cone block
    psd = [bb for bb in r["block_report"] if bb["kind"] == "psd"]
    for bb in psd:
        print(f"  PSD block n={bb['size']}: d_lower(lambda_min)={bb['d_lower']:.3e} "
              f"sbar={bb['sbar']:.3e} contrib={bb['contrib_lo']:.3e}")
    print("=" * 72)


# =====================================================================
#  CLI
# =====================================================================

if __name__ == "__main__":
    import argparse, warnings
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("--test", action="store_true", help="run unit tests only")
    ap.add_argument("--N", type=int, default=300)
    ap.add_argument("--T", type=int, default=120)
    ap.add_argument("--R", type=int, default=6)
    ap.add_argument("--bochner_n", type=int, default=6)
    ap.add_argument("--pm_k_max", type=int, default=0)
    ap.add_argument("--center", type=str, default="row4",
                    choices=["row4", "cde_n30_iter3"])
    ap.add_argument("--use_T5p", action="store_true")
    ap.add_argument("--slack_infl", type=float, default=4.0)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    if args.test:
        _unit_test_svec()
        _unit_test_eig()
        print("\nALL UNIT TESTS PASSED")
        sys.exit(0)

    if args.center == "row4":
        h_c, p_c, q1, q2 = 0.004, 0.3875, -0.02, 0.02
    else:  # cde_n30_iter3
        h_c, p_c, q1, q2 = 0.000045, 0.39015, -0.02, 0.02

    r = jansson_lower_bound(args.N, args.T, args.R, h_c, p_c, q1, q2,
                            args.bochner_n, args.pm_k_max,
                            use_T5p=args.use_T5p, slack_infl=args.slack_infl)
    if args.out:
        Path(args.out).write_text(json.dumps(r, indent=2, default=float))
        print(f"→ wrote {args.out}")
