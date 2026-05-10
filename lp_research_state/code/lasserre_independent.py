"""
Independent re-implementation of Lasserre level-2 SDP encoding for the
polynomial-optimization relaxation of f^2 <= f, on top of the
Section-5 program of White (arXiv:2201.05704).

Conventions
-----------
Fourier expansion of the (real, period-2) density f(x):
    f(x) = 1/2 + sum_{k>=1} ( c_k cos(k*pi*x) + d_k sin(k*pi*x) ).
Hence
    f̂(0)  = 1/2
    f̂(k)  = (c_k - i d_k)/2          for k >= 1
    f̂(-k) = conj(f̂(k)) = (c_k + i d_k)/2.

The Lasserre level-2 lift (re-derived below) lifts the bilinear products
f̂(j)·f̂(k) into NEW DECISION VARIABLES Y[j,k], j,k in {-T_max, ..., T_max},
forming a (2*T_max+1)x(2*T_max+1) Hermitian matrix Y with the constraint
    [[1, ξ*],
     [ξ, Y ]] ⪰ 0      (Schur complement form)
where ξ = (f̂(-T_max), ..., f̂(T_max))^T is the linear-feature vector.

Conjugate-pair reductions enforce
    Y[j,k] = conj(Y[k,j]),      Y[-j,-k] = conj(Y[j,k]),
which collapses the 'free' real parameters to roughly half.

Localizing constraint for f^2 <= f
----------------------------------
The constraint f - f^2 >= 0 (a.e.) lifted via Lasserre level-2 is encoded
via PSD on the Hermitian Toeplitz-like matrix
    Loc[j,k] = (f - f^2)̂(j-k)
            = f̂(j-k)  -  sum_{m=-T_max}^{T_max}  f̂(m) * f̂(j-k-m)
            = f̂(j-k)  -  sum_{m}  Y[m, j-k-m]
indexed by j,k in {0, 1, ..., T_loc}.  Loc must be PSD (it is automatically
Hermitian).  This tests f - f^2 >= 0 against ALL nonnegative degree-T_loc
trigonometric polynomials of the form |p|^2, which strictly contains the
Bochner-on-f constraint (which is only the Toeplitz of f̂, not f̂ - (f^2)̂).

Author: independent verification (no peeking at lasserre.py).
"""

from __future__ import annotations
import numpy as np
import cvxpy as cp


def _fhat_re_im_from_cd(c: cp.Variable, d: cp.Variable, k: int):
    """
    Return Re(f̂(k)), Im(f̂(k)) as cvxpy expressions.  k may be negative.

    Convention: f̂(0) = 1/2; f̂(k) = (c_k - i d_k)/2 for k>=1.
    """
    if k == 0:
        return cp.Constant(0.5), cp.Constant(0.0)
    if k > 0:
        return 0.5 * c[k - 1], -0.5 * d[k - 1]
    # k < 0: conjugate of f̂(|k|)
    return 0.5 * c[-k - 1], 0.5 * d[-k - 1]


def add_lasserre2_constraint_indep(
    cons: list,
    c: cp.Variable,
    d: cp.Variable,
    T_max: int,
    T_loc: int,
):
    """
    Append the Lasserre level-2 SDP constraints to `cons`.

    Parameters
    ----------
    cons   : list of cvxpy constraints to extend.
    c, d   : cvxpy.Variable of length T (the cosine and sine Fourier
             coefficients of f, indexed 1..T).
    T_max  : moment-matrix order. The lifted variables Y[j,k] are introduced
             for j,k in {-T_max, ..., T_max}.
    T_loc  : localizing-matrix order. The constraint Loc[j,k] = (f-f^2)̂(j-k)
             ⪰ 0 is imposed for j,k in {0, ..., T_loc}.

    Notes
    -----
    For Y to support all entries of the localizing matrix, we need entries
    (f^2)̂(ℓ) for |ℓ| <= T_loc, i.e. Y[m, ℓ-m] for |m|, |ℓ-m| <= T_max.  In
    particular, taking m = -T_max and ℓ = T_loc requires j-k-m = T_loc + T_max
    <= T_max in absolute value of the f̂ index, but here we are summing over
    m so the index j-k-m ranges over [ℓ-T_max, ℓ+T_max] which can exceed
    T_max.  We therefore restrict the inner sum to those m with |m| <= T_max
    AND |ℓ-m| <= T_max -- i.e. we only use legal Y entries.  This is a
    valid relaxation (we ARE truncating Lasserre, after all).
    """
    T = c.shape[0]
    assert T_max >= 1, "T_max must be >= 1"
    assert T_max <= T, f"T_max={T_max} must be <= T={T}"
    assert T_loc >= 1, "T_loc must be >= 1"
    assert T_loc <= T, f"T_loc={T_loc} must be <= T={T}"

    n = 2 * T_max + 1  # size of moment matrix
    # Index map: matrix-row index r in {0, ..., n-1} <-> Fourier index j = r - T_max.
    def jof(r):
        return r - T_max

    # --- Step 1: declare the lifted Y variables ----------------------------
    # Y is Hermitian, size n x n.  We store Y_re (n x n symmetric) and
    # Y_im (n x n antisymmetric).
    Y_re = cp.Variable((n, n), symmetric=True)
    Y_im = cp.Variable((n, n))
    # antisymmetric: Y_im^T = -Y_im
    cons.append(Y_im + Y_im.T == 0)

    # --- Step 2: conjugate-pair symmetry --------------------------------------
    # Y[j,k] is supposed to represent f̂(j)*f̂(k). With f real, f̂(-j) = conj(f̂(j)),
    # so f̂(-j)*f̂(-k) = conj(f̂(j)*f̂(k)) = conj(Y[j,k]).
    # In matrix entries: Y[r_neg_j, r_neg_k] = conj(Y[r_j, r_k]).
    # Equivalently Y_re[r_neg_j, r_neg_k] = Y_re[r_j, r_k]
    #              Y_im[r_neg_j, r_neg_k] = -Y_im[r_j, r_k].
    for r1 in range(n):
        for r2 in range(r1, n):  # r2 >= r1 (symmetry already enforced for re/im as above)
            j = jof(r1); k = jof(r2)
            r1n = -j + T_max  # index for -j
            r2n = -k + T_max  # index for -k
            if (r1n, r2n) == (r1, r2):
                continue
            cons.append(Y_re[r1, r2] == Y_re[r1n, r2n])
            cons.append(Y_im[r1, r2] == -Y_im[r1n, r2n])

    # --- Step 3: build the linear-feature vector ξ and the Schur block -------
    # ξ has length n, components ξ[r] = f̂(jof(r)).  Build Re/Im as cvxpy expressions.
    xi_re_list = []
    xi_im_list = []
    for r in range(n):
        re_, im_ = _fhat_re_im_from_cd(c, d, jof(r))
        xi_re_list.append(re_)
        xi_im_list.append(im_)
    # Stack as column vectors of cvxpy expressions
    xi_re = cp.reshape(cp.hstack(xi_re_list), (n, 1), order="C")
    xi_im = cp.reshape(cp.hstack(xi_im_list), (n, 1), order="C")

    # The bordered matrix
    #   [[ 1     ξ* ],
    #    [ ξ     Y  ]]
    # is Hermitian PSD iff its real-form representation
    #   [[ Re_block, -Im_block ],
    #    [ Im_block,  Re_block ]]
    # is real symmetric PSD.
    # With first row/column "real scalar 1", the block sizes are:
    #   Re_block = [[1, ξ_reᵀ], [ξ_re, Y_re]]   (size (n+1)x(n+1))
    #   Im_block = [[0, -ξ_imᵀ], [ξ_im, Y_im]]  (size (n+1)x(n+1), antisymm)
    # NOTE: the (1,1) entry of the bordered Hermitian matrix is the constant 1.
    # The first row of the bordered matrix is (1, ξ*), i.e. (1, ξ_re^T - i ξ_im^T).
    # Its real part is (1, ξ_re^T), its imag part is (0, -ξ_im^T).
    # The first column is (1, ξ), i.e. (1, ξ_re + i ξ_im); real part (1, ξ_re),
    # imag part (0, ξ_im).  The off-diagonal Im block must be antisymmetric:
    # (Im_block)[0, k] = -ξ_im[k-1]; (Im_block)[k, 0] = +ξ_im[k-1]; consistent with antisymm.
    one = cp.Constant(np.array([[1.0]]))
    zero_row = cp.Constant(np.zeros((1, n)))
    zero_col = cp.Constant(np.zeros((n, 1)))
    zero_11 = cp.Constant(np.array([[0.0]]))

    Re_block = cp.bmat([
        [one,         cp.reshape(xi_re, (1, n), order="C")],
        [xi_re,       Y_re],
    ])
    Im_block = cp.bmat([
        [zero_11,                                cp.reshape(-xi_im, (1, n), order="C")],
        [xi_im,                                  Y_im],
    ])
    real_form = cp.bmat([[Re_block, -Im_block], [Im_block, Re_block]])
    cons.append(real_form >> 0)

    # --- Step 4: localizing matrix for f^2 <= f -------------------------------
    # Loc[r1, r2] = f̂(jof(r1) - jof(r2)) - sum_m f̂(m) f̂( (jof(r1) - jof(r2)) - m )
    #            = f̂(ℓ) - (f^2)̂(ℓ) where ℓ = jof(r1) - jof(r2).
    # We build Loc as an (T_loc+1) x (T_loc+1) Hermitian matrix indexed by
    # r1, r2 in {0, ..., T_loc} corresponding to "Fourier index = r1, r2"
    # in the localizing test polynomial.  Then Loc[r1, r2] = (f-f^2)̂(r1 - r2).
    nL = T_loc + 1
    Loc_re_rows = []
    Loc_im_rows = []
    for r1 in range(nL):
        re_row = []
        im_row = []
        for r2 in range(nL):
            ell = r1 - r2  # Fourier lag
            # First term: f̂(ell) - linear in c, d.
            fhat_re, fhat_im = _fhat_re_im_from_cd(c, d, ell)
            # Second term: (f^2)̂(ell) = sum_m f̂(m) f̂(ell - m).
            # In Y entries: Y[m, ell-m] (interpreting Y[j,k] = f̂(j) f̂(k)).
            # We sum only over m for which BOTH |m| <= T_max AND |ell-m| <= T_max.
            f2_re_terms = []
            f2_im_terms = []
            for m in range(-T_max, T_max + 1):
                k = ell - m
                if abs(k) > T_max:
                    continue
                # Y[m, k] in matrix indices:
                rm = m + T_max
                rk = k + T_max
                f2_re_terms.append(Y_re[rm, rk])
                f2_im_terms.append(Y_im[rm, rk])
            if f2_re_terms:
                f2_re = cp.sum(cp.hstack(f2_re_terms))
                f2_im = cp.sum(cp.hstack(f2_im_terms))
            else:
                f2_re = cp.Constant(0.0)
                f2_im = cp.Constant(0.0)
            re_row.append(fhat_re - f2_re)
            im_row.append(fhat_im - f2_im)
        Loc_re_rows.append(re_row)
        Loc_im_rows.append(im_row)
    Loc_Re = cp.bmat(Loc_re_rows)
    Loc_Im = cp.bmat(Loc_im_rows)
    Loc_real_form = cp.bmat([[Loc_Re, -Loc_Im], [Loc_Im, Loc_Re]])
    cons.append(Loc_real_form >> 0)


# =============================================================================
# Self-test
# =============================================================================
if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    print("=" * 70)
    print("SELF-TEST 1: f(x) = 1/2 + a cos(πx) for a = 0.4 (nonneg)")
    print("=" * 70)
    a_ = 0.4
    T = 5
    T_max = 3
    T_loc = 3
    # Pin c_1 = a, others = 0; d = 0; check whether feasibility region admits
    # this (it must, since f is nonneg and bounded above by 1).
    c = cp.Variable(T)
    d = cp.Variable(T)
    cons = [c[0] == a_, d == 0]
    for k in range(2, T + 1):
        cons.append(c[k - 1] == 0)
    add_lasserre2_constraint_indep(cons, c, d, T_max=T_max, T_loc=T_loc)
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL")
    print(f"  status: {prob.status}    (expected: optimal — f is nonneg AND f^2 <= f a.e.)")

    print()
    print("=" * 70)
    print("SELF-TEST 2: f(x) = 1/2 + a cos(πx) for a = 0.6 (NEGATIVE at x=1)")
    print("=" * 70)
    a_ = 0.6
    c = cp.Variable(T)
    d = cp.Variable(T)
    cons = [c[0] == a_, d == 0]
    for k in range(2, T + 1):
        cons.append(c[k - 1] == 0)
    add_lasserre2_constraint_indep(cons, c, d, T_max=T_max, T_loc=T_loc)
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL")
    print(f"  status: {prob.status}    (expected: infeasible — f<0 at x=1)")

    print()
    print("=" * 70)
    print("SELF-TEST 3: Y[1,1] = a^2/4 for the trivial case f = 1/2 + a cos(πx)")
    print("=" * 70)
    # We re-solve and inspect Y_re[T_max+1, T_max+1].
    a_ = 0.4
    c = cp.Variable(T)
    d = cp.Variable(T)
    cons = [c[0] == a_, d == 0]
    for k in range(2, T + 1):
        cons.append(c[k - 1] == 0)
    # Custom: extract Y_re directly.  Re-call but capture the variable.
    # Quick way: instantiate variables manually and re-add.
    # Simpler: inspect by minimizing slack.
    # We minimize Y_re[T_max+1, T_max+1] (which = f̂(1)^2 - slack = a^2/4 - slack >= a^2/4
    # at the optimum, since the Schur block requires Y >= ξ ξ*).
    # Hack: temporarily build the Lasserre matrix WITH access to Y_re.
    n = 2 * T_max + 1
    Y_re = cp.Variable((n, n), symmetric=True)
    Y_im = cp.Variable((n, n))
    cons2 = [c[0] == a_, d == 0]
    for k in range(2, T + 1):
        cons2.append(c[k - 1] == 0)
    cons2.append(Y_im + Y_im.T == 0)
    # Pure replication of the encoding above, but with our local Y_re/Y_im so
    # we can read out the value.  We just call the function on a NEW cons list,
    # then stash a reference.  Simpler: parse out from the dual variables. But
    # easier still — solve with objective Y_re[T_max+1, T_max+1] minimized:
    # Re-build encoding from scratch (manually, dup of the function)
    def jof(r): return r - T_max
    def fhat(k_):
        if k_ == 0: return cp.Constant(0.5), cp.Constant(0.0)
        if k_ > 0:  return 0.5 * c[k_-1], -0.5 * d[k_-1]
        return 0.5 * c[-k_-1], 0.5 * d[-k_-1]
    for r1 in range(n):
        for r2 in range(r1, n):
            j = jof(r1); k = jof(r2)
            r1n = -j + T_max
            r2n = -k + T_max
            if (r1n, r2n) == (r1, r2):
                continue
            cons2.append(Y_re[r1, r2] == Y_re[r1n, r2n])
            cons2.append(Y_im[r1, r2] == -Y_im[r1n, r2n])
    xi_re_list = []; xi_im_list = []
    for r in range(n):
        re_, im_ = fhat(jof(r))
        xi_re_list.append(re_); xi_im_list.append(im_)
    xi_re = cp.reshape(cp.hstack(xi_re_list), (n, 1), order="C")
    xi_im = cp.reshape(cp.hstack(xi_im_list), (n, 1), order="C")
    one = cp.Constant(np.array([[1.0]]))
    zero_11 = cp.Constant(np.array([[0.0]]))
    Re_block = cp.bmat([[one, cp.reshape(xi_re, (1, n), order="C")],
                        [xi_re, Y_re]])
    Im_block = cp.bmat([[zero_11, cp.reshape(-xi_im, (1, n), order="C")],
                        [xi_im, Y_im]])
    real_form = cp.bmat([[Re_block, -Im_block], [Im_block, Re_block]])
    cons2.append(real_form >> 0)
    # Localizing
    Loc_re_rows = []; Loc_im_rows = []
    for r1 in range(T_loc + 1):
        rr_re = []; rr_im = []
        for r2 in range(T_loc + 1):
            ell = r1 - r2
            fr, fi = fhat(ell)
            f2_re = []; f2_im = []
            for m in range(-T_max, T_max+1):
                kk = ell - m
                if abs(kk) > T_max: continue
                rm = m + T_max; rk = kk + T_max
                f2_re.append(Y_re[rm, rk])
                f2_im.append(Y_im[rm, rk])
            f2_re_s = cp.sum(cp.hstack(f2_re)) if f2_re else cp.Constant(0.0)
            f2_im_s = cp.sum(cp.hstack(f2_im)) if f2_im else cp.Constant(0.0)
            rr_re.append(fr - f2_re_s)
            rr_im.append(fi - f2_im_s)
        Loc_re_rows.append(rr_re); Loc_im_rows.append(rr_im)
    Loc_Re = cp.bmat(Loc_re_rows); Loc_Im = cp.bmat(Loc_im_rows)
    Loc_real_form = cp.bmat([[Loc_Re, -Loc_Im], [Loc_Im, Loc_Re]])
    cons2.append(Loc_real_form >> 0)
    # Minimize Y_re[T_max+1, T_max+1] (the (j=1, k=1) entry).
    prob2 = cp.Problem(cp.Minimize(Y_re[T_max + 1, T_max + 1]), cons2)
    prob2.solve(solver="CLARABEL")
    print(f"  status: {prob2.status}")
    print(f"  Y[1,1] minimized: {Y_re.value[T_max+1, T_max+1]:.8f}")
    print(f"  Expected lower bound: a^2/4 = {a_**2/4:.8f}")
    print(f"  Match: {abs(Y_re.value[T_max+1, T_max+1] - a_**2/4) < 1e-6}")
