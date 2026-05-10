"""
Lasserre level-3 (Toeplitz-reduced) augmentation for the Erdős
minimum-overlap problem.

THEORETICAL OBSERVATION (justifies what is implemented)
=======================================================
In White (2023)'s Fourier formulation, the constraint f² ≤ f is a
function-positivity constraint on the circle. By Fejér-Riesz, every
nonneg trig polynomial of degree 2n equals |p|² with deg p ≤ n. Hence
the Bochner-style localizing matrix
    Loc[j, k] := (f - f²)̂(j - k),   j, k = 0..T_loc
of order T_loc captures ALL nonneg degree-2T_loc trig polynomial test
functions. Climbing Lasserre levels in the standard sense gives no
new constraints "on the circle" — the only meaningful strengthening
comes from MOMENT lifts of the auxiliary variables (c_k, d_k).

CONCRETE LEVEL-3 LIFT IMPLEMENTED HERE
======================================
We add a TRILINEAR moment lift on a small inner basis
    ξ^lift = (1, c_1, ..., c_{T3_lift}, d_1, ..., d_{T3_lift})
of length L_lift = 2 T3_lift + 1, plus a 4th-order COUPLING moment
matrix Q_m for each outer coordinate c_m, m = 1..M3_max. Specifically
we introduce two auxiliary cvxpy variables:

  W_m  of shape (L_lift, L_lift), interpreted as
          W_m[a, b] = E[c_m · ξ^lift_a · ξ^lift_b]                (3rd-order)

  Q_m  of shape (L_lift, L_lift), interpreted as
          Q_m[a, b] = E[c_m² · ξ^lift_a · ξ^lift_b]               (4th-order)

These are tied to the level-2 lift M_top via:
  (1) W_m[0, b] = M_top[m, big_idx[b]]                           (link)
  (2) Q_m[0, 0] = M_top[m, m]                                    (link)
  (3) The bordered matrix
          [[ M_top sub-block on (1, ξ^lift)    | W_m ]
           [ W_m^T                              | Q_m ]] ⪰ 0     (PSD lift)

Constraint (3) is the proper Lasserre LMI for a 4-th order moment matrix
of (1, ξ^lift, c_m·ξ^lift), which generalizes the bilinear PSD lift to
include the trilinear/quadrilinear moments.

This is a STRICT TIGHTENING IF (and only if) Q_m enters another
constraint. We add the 4th-order LOCALIZING matrix for c_m² (f - f²) ≥ 0:
since c_m² ≥ 0 always, this is a valid constraint:
    Loc^{(m)}_3[j, k] := (c_m² · (f - f²))̂(j - k),  j, k = 0..T_loc3
which expands using Q_m + the lift.

(c_m² · (f - f²))̂(ℓ)
   = c_m² · f̂(ℓ) - Σ_n c_m² · f̂(n) · f̂(ℓ - n)
The first piece is c_m² · (linear in c, d) — a 3rd-order moment
captured by the Q_m link Q_m[0, b] = "c_m² · ξ^lift_b" (extended).
The second piece is 4th-order, captured by Q_m via Q_m[a, b] = c_m² · ξ^lift_a · ξ^lift_b
where (a, b) encode the f̂(n), f̂(ℓ-n) factors.

Toeplitz block decomposition
----------------------------
For each outer m (m = 1..M3_max), Loc^{(m)}_3 is a Hermitian Toeplitz
matrix of size (T_loc3 + 1) — independent of the others. So the full
level-3 SDP decomposes into M3_max independent (T_loc3+1)-dim Hermitian
PSD constraints, plus M3_max bordered LMIs of size (1+L_lift+L_lift) =
2 L_lift + 1 each. Total cost: O(M3_max · (T_loc3³ + L_lift³)) — tiny
compared to a flat level-3 lift O(T_max⁶).
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


# ---------------------------------------------------------------------------
# Helpers (mirrored from lasserre.py)
# ---------------------------------------------------------------------------
def _f_hat_re_im(c, d, k):
    if k == 0:
        return cp.Constant(0.5), cp.Constant(0.0)
    elif k > 0:
        return 0.5 * c[k - 1], -0.5 * d[k - 1]
    else:
        return 0.5 * c[-k - 1], 0.5 * d[-k - 1]


def _bilinear_in_M(M, T_max, kind1, k1, kind2, k2):
    def idx(kind, k):
        if kind == "1": return 0
        elif kind == "c": return k
        elif kind == "d": return T_max + k
        raise ValueError
    i = idx(kind1, k1); j = idx(kind2, k2)
    if i == 0 and j == 0: return cp.Constant(1.0)
    if i == 0: return M[0, j]
    if j == 0: return M[i, 0]
    return M[i, j]


def _f2_hat_re_im(M, T_max, m):
    n_lo = max(-T_max, m - T_max); n_hi = min(T_max, m + T_max)
    re = cp.Constant(0.0); im = cp.Constant(0.0)
    for n in range(n_lo, n_hi + 1):
        n2 = m - n
        if n == 0: re_a_kind="1"; re_a_k=0; w_re_a=0.5; im_a_kind="1"; im_a_k=0; w_im_a=0.0
        elif n > 0: re_a_kind="c"; re_a_k=n; w_re_a=0.5; im_a_kind="d"; im_a_k=n; w_im_a=-0.5
        else: re_a_kind="c"; re_a_k=-n; w_re_a=0.5; im_a_kind="d"; im_a_k=-n; w_im_a=+0.5
        if n2 == 0: re_b_kind="1"; re_b_k=0; w_re_b=0.5; im_b_kind="1"; im_b_k=0; w_im_b=0.0
        elif n2 > 0: re_b_kind="c"; re_b_k=n2; w_re_b=0.5; im_b_kind="d"; im_b_k=n2; w_im_b=-0.5
        else: re_b_kind="c"; re_b_k=-n2; w_re_b=0.5; im_b_kind="d"; im_b_k=-n2; w_im_b=+0.5
        if w_re_a and w_re_b:
            re = re + (w_re_a*w_re_b) * _bilinear_in_M(M, T_max, re_a_kind, re_a_k, re_b_kind, re_b_k)
        if w_im_a and w_im_b:
            re = re - (w_im_a*w_im_b) * _bilinear_in_M(M, T_max, im_a_kind, im_a_k, im_b_kind, im_b_k)
        if w_re_a and w_im_b:
            im = im + (w_re_a*w_im_b) * _bilinear_in_M(M, T_max, re_a_kind, re_a_k, im_b_kind, im_b_k)
        if w_im_a and w_re_b:
            im = im + (w_im_a*w_re_b) * _bilinear_in_M(M, T_max, im_a_kind, im_a_k, re_b_kind, re_b_k)
    return re, im


def _Q_hat_re_im_for_cm2_f2(Q, T_max, T3_lift, m):
    """(c_m² · f²)̂(m_freq) using the 4th-order lift Q.
    Q[a, b] = E[c_m² · ξ^lift_a · ξ^lift_b].
    (c_m² · f²)̂(ℓ) = Σ_n c_m² · f̂(n) · f̂(ℓ - n).
    For each (n, ℓ-n) with |n|, |ℓ-n| ≤ T3_lift, the summand has the form
    c_m² · (kind_a ± i d_kind_a)/2 · (kind_b ± i d_kind_b)/2, expanding
    to bilinear-in-(c, d, Q) terms via Q[a, b].
    """
    # ξ^lift index lookup
    def idx(kind, k):
        if kind == "1": return 0
        elif kind == "c":
            if not (1 <= k <= T3_lift): return None
            return k
        elif kind == "d":
            if not (1 <= k <= T3_lift): return None
            return T3_lift + k
        return None

    n_lo = max(-T3_lift, m - T3_lift); n_hi = min(T3_lift, m + T3_lift)
    re = cp.Constant(0.0); im = cp.Constant(0.0)
    for n in range(n_lo, n_hi + 1):
        n2 = m - n
        if abs(n) > T3_lift or abs(n2) > T3_lift:
            continue
        if n == 0: re_a_kind="1"; re_a_k=0; w_re_a=0.5; im_a_kind="1"; im_a_k=0; w_im_a=0.0
        elif n > 0: re_a_kind="c"; re_a_k=n; w_re_a=0.5; im_a_kind="d"; im_a_k=n; w_im_a=-0.5
        else: re_a_kind="c"; re_a_k=-n; w_re_a=0.5; im_a_kind="d"; im_a_k=-n; w_im_a=+0.5
        if n2 == 0: re_b_kind="1"; re_b_k=0; w_re_b=0.5; im_b_kind="1"; im_b_k=0; w_im_b=0.0
        elif n2 > 0: re_b_kind="c"; re_b_k=n2; w_re_b=0.5; im_b_kind="d"; im_b_k=n2; w_im_b=-0.5
        else: re_b_kind="c"; re_b_k=-n2; w_re_b=0.5; im_b_kind="d"; im_b_k=-n2; w_im_b=+0.5

        def Q_at(ka, kk_a, kb, kk_b):
            ia = idx(ka, kk_a); ib = idx(kb, kk_b)
            if ia is None or ib is None:
                return None
            return Q[ia, ib]

        if w_re_a and w_re_b:
            t = Q_at(re_a_kind, re_a_k, re_b_kind, re_b_k)
            if t is not None: re = re + (w_re_a*w_re_b) * t
        if w_im_a and w_im_b:
            t = Q_at(im_a_kind, im_a_k, im_b_kind, im_b_k)
            if t is not None: re = re - (w_im_a*w_im_b) * t
        if w_re_a and w_im_b:
            t = Q_at(re_a_kind, re_a_k, im_b_kind, im_b_k)
            if t is not None: im = im + (w_re_a*w_im_b) * t
        if w_im_a and w_re_b:
            t = Q_at(im_a_kind, im_a_k, re_b_kind, re_b_k)
            if t is not None: im = im + (w_im_a*w_re_b) * t
    return re, im


def _cm2_f_hat_re_im(Q, T_max, T3_lift, m):
    """(c_m² · f)̂(ℓ).
    f̂(0) = 1/2 → contributes c_m²/2 if ℓ = 0.
    f̂(±k) = (c_k ∓ i d_k)/2 → contributes c_m²·c_k/2 ∓ i·c_m²·d_k/2.
    All three terms are 3rd-order moments, captured by Q[0, idx_k] (since
    Q[0, b] = E[c_m² · ξ^lift_b]).
    """
    def idx(kind, k):
        if kind == "1": return 0
        elif kind == "c":
            if not (1 <= k <= T3_lift): return None
            return k
        elif kind == "d":
            if not (1 <= k <= T3_lift): return None
            return T3_lift + k
        return None

    if m == 0:
        return 0.5 * Q[0, 0], cp.Constant(0.0)
    elif m > 0:
        ic = idx("c", m); id_ = idx("d", m)
        if ic is None or id_ is None:
            return cp.Constant(0.0), cp.Constant(0.0)   # truncated
        return 0.5 * Q[0, ic], -0.5 * Q[0, id_]
    else:
        ic = idx("c", -m); id_ = idx("d", -m)
        if ic is None or id_ is None:
            return cp.Constant(0.0), cp.Constant(0.0)
        return 0.5 * Q[0, ic], 0.5 * Q[0, id_]


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def add_lasserre3_toeplitz_constraint(
    cons, c, d, M_top, T_max, T_loc3, T3_lift=1, M3_max=None
):
    """Add the level-3 trilinear/quadrilinear lift and the
    Loc^{(m)}_3[j, k] = (c_m² · (f - f²))̂(j-k) Hermitian Toeplitz PSD
    constraint for each m = 1..M3_max."""
    if M3_max is None:
        M3_max = min(T_max, 3)
    L_lift = 2 * T3_lift + 1

    # ξ^lift index → big-ξ index in M_top
    big_idx = [0]
    for a in range(1, T3_lift + 1):
        big_idx.append(a)        # c_a
    for a in range(1, T3_lift + 1):
        big_idx.append(T_max + a)  # d_a

    W = {}; Q = {}
    for m_full in range(1, M3_max + 1):
        # Trilinear lift W_m[a, b] = E[c_m · ξ^lift_a · ξ^lift_b]
        Wm = cp.Variable((L_lift, L_lift), symmetric=True)
        # 4th-order lift Q_m[a, b] = E[c_m² · ξ^lift_a · ξ^lift_b]
        Qm = cp.Variable((L_lift, L_lift), symmetric=True)
        W[m_full] = Wm; Q[m_full] = Qm

        # Linking: W_m[0, b] = E[c_m · ξ^lift_b] = M_top[m_full, big_idx[b]].
        for b in range(L_lift):
            cons.append(Wm[0, b] == M_top[m_full, big_idx[b]])
        # Linking: Q_m[0, 0] = E[c_m²] = M_top[m_full, m_full].
        cons.append(Qm[0, 0] == M_top[m_full, m_full])
        # Q_m[0, b] = E[c_m² · ξ^lift_b]: free trilinear moment (NOT pinned).

        # Compose the 4th-order moment matrix on (1, ξ^lift, c_m·ξ^lift):
        #   [[ M_top_sub        |  W_m   ],
        #    [ W_m^T            |  Q_m   ]]   ⪰ 0
        # M_top_sub is the (L_lift x L_lift) sub-block of M_top on big_idx.
        Mtop_sub = cp.bmat([
            [M_top[big_idx[i], big_idx[j]] if (big_idx[i] != 0 or big_idx[j] != 0)
             else cp.Constant(1.0)
             for j in range(L_lift)] for i in range(L_lift)
        ])
        big = cp.bmat([
            [Mtop_sub, Wm],
            [Wm.T,    Qm],
        ])
        cons.append(big >> 0)

        # Localizing matrix Loc^{(m)}_3[j, k] = (c_m² · (f - f²))̂(j - k),
        # j, k = 0..T_loc3, Hermitian Toeplitz.  Encoded as 2(T_loc3+1) PSD.
        Re_rows = []; Im_rows = []
        for j in range(T_loc3 + 1):
            re_row = []; im_row = []
            for k in range(T_loc3 + 1):
                ell = j - k
                # (c_m² · f)̂(ell)
                f_re, f_im = _cm2_f_hat_re_im(Qm, T_max, T3_lift, ell)
                # (c_m² · f²)̂(ell)
                f2_re, f2_im = _Q_hat_re_im_for_cm2_f2(Qm, T_max, T3_lift, ell)
                re_row.append(f_re - f2_re)
                im_row.append(f_im - f2_im)
            Re_rows.append(re_row); Im_rows.append(im_row)
        Re_M = cp.bmat(Re_rows); Im_M = cp.bmat(Im_rows)
        cons.append(cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]]) >> 0)

    return {"W": W, "Q": Q, "big_idx": big_idx, "M3_max": M3_max,
            "T3_lift": T3_lift}


# ---------------------------------------------------------------------------
# Public hook
# ---------------------------------------------------------------------------
def add_lasserre3(cons, c, d, T_max, T_loc, T3_lift=1, M3_max=None):
    """Add level-2 + level-3 lift to `cons`. Returns (M_top, level3_dict)."""
    import importlib.util as _ilu
    from pathlib import Path as _Path
    here = _Path(__file__).resolve().parent
    spec = _ilu.spec_from_file_location("lasserre", here / "lasserre.py")
    mod = _ilu.module_from_spec(spec); spec.loader.exec_module(mod)
    M_top = mod.add_lasserre2_constraint(cons, c, d, T_max=T_max, T_loc=T_loc)
    info = add_lasserre3_toeplitz_constraint(
        cons, c, d, M_top, T_max, T_loc3=T_loc, T3_lift=T3_lift, M3_max=M3_max
    )
    return M_top, info


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    print("=== Lasserre-3 self-test 1: feasible f = 1/2 + 0.4 cos(πx) ===")
    T = 5; T_max = 5; T_loc = 3
    c = cp.Variable(T); d = cp.Variable(T)
    cons = []
    M_top, info = add_lasserre3(cons, c, d, T_max=T_max, T_loc=T_loc, T3_lift=1, M3_max=2)
    cons += [c[0] == 0.4, c[1] == 0.0, c[2] == 0.0, c[3] == 0.0, c[4] == 0.0]
    cons += [d == 0]
    xi = np.zeros(2*T_max+1); xi[0] = 1.0; xi[1] = 0.4
    cons.append(M_top == np.outer(xi, xi))
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  Feasible f: status = {prob.status}")

    print("=== Lasserre-3 self-test 2: infeasible f = 1/2 + 0.7 cos(πx) ===")
    c = cp.Variable(T); d = cp.Variable(T)
    cons = []
    M_top, info = add_lasserre3(cons, c, d, T_max=T_max, T_loc=T_loc, T3_lift=1, M3_max=2)
    cons += [c[0] == 0.7, c[1] == 0, c[2] == 0, c[3] == 0, c[4] == 0]
    cons += [d == 0]
    xi = np.zeros(2*T_max+1); xi[0] = 1.0; xi[1] = 0.7
    cons.append(M_top == np.outer(xi, xi))
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  Infeasible f: status = {prob.status} (expect: infeasible*)")
