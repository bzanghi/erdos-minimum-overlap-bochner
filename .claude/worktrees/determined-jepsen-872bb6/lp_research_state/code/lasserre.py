"""
Lasserre level-2 augmentation of White's Fourier convex program for the
Erdős minimum overlap problem.

Goal:  enforce  f^2 <= f  pointwise (since f in [0,1])  against ALL
nonneg degree-2 trigonometric polynomials simultaneously, by lifting the
bilinear products of Fourier coefficients into a moment matrix and adding
PSD localizing-matrix constraints.

White's existing (5.5) tests f^2 <= f against single trig polynomials per m;
Bochner-on-f tests against the family of |p|^2 with deg(p) <= n_b (i.e.,
all nonneg deg-2n_b trig polys factoring as a single square).  Lasserre
level-2 tests against ALL nonneg degree-2 trig polynomials, which is the
full SOS cone at degree-2 — a strictly stronger family.

Lift used here (real basis):
    xi = (1, c_1, ..., c_{T_max}, d_1, ..., d_{T_max}),       length L = 2T_max+1.
    M  = (L+1)x(L+1) symmetric SDP variable encoding [[ xi xi^T, xi ],
                                                      [ xi^T,    1  ]] ⪰ 0
    via Schur complement equivalent to ξ ξ^T ⪯ M_top, with M_top := M[0:L,0:L].

Then every degree-2 polynomial in (c, d) (in particular every |\hat f(k)|^2,
\hat f(j)\hat f(k), etc.) becomes LINEAR in M_top.

For real f, with the conventions
    f̂(0)=1/2,   f̂(k)=(c_k - i d_k)/2  (k>=1),   f̂(-k) = conj(f̂(k)),
the convolution coefficients of f^2 are
    (f^2)̂(m) = sum_n f̂(n) f̂(m-n)        (n ∈ Z),
each summand a constant + linear-in-(c,d) + bilinear-in-(c,d) expression
that becomes linear in (c, d, M_top) after the lift.

The Lasserre level-2 localizing constraint for (f - f^2) >= 0 is the
Hermitian matrix
    Loc[j,k] := (f - f^2)̂(j - k),  j,k = 0..T_loc,
must be PSD (Bochner). This matrix is linear in (c, d, M_top), so the
constraint is a (T_loc+1)-dim Hermitian PSD constraint encoded as a real
2(T_loc+1)x2(T_loc+1) PSD constraint.

We also include the symmetric Lasserre constraint for (f^2 - f) >= 0  is NOT
a true constraint (we only have f^2 <= f, not equality), so only one side.
However, Bochner-on-f and Bochner-on-(1-f) (already in white_full_convex.py)
are included via existing flags; they remain compatible.

Key implementation caveats verified below:
  * Hermitian symmetry of Loc:     Loc[k,j] = conj(Loc[j,k])  ✓ since
    (f-f^2) is real, so its Fourier coeffs satisfy ĝ(-m) = conj(ĝ(m)).
  * Truncation:  in (f^2)̂(m) = sum_n f̂(n)f̂(m-n) we truncate to
    |n|, |m-n| <= T_max.  The dropped terms have product of two coefficients
    each at high frequency; their omission is conservative when handled via
    a tail bound, but for an SDP UPPER bound on Ω we can include only the
    truncated sum (giving a TIGHTER constraint, hence a LOOSER bound on Ω).
    To remain a valid LOWER bound on µ we'd need to upper-bound the missing
    tail; here we instead include the constraint with the truncated formula,
    treating the result as a heuristic test (matching the spirit of how
    bochner_n is used in the existing code).
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def _fourier_real_indexed(c: cp.Variable, d: cp.Variable, k: int):
    """Return real and imaginary parts of f̂(k) as cvxpy expressions.

    f̂(0) = 1/2 (real).
    f̂(+k) = (c_k - i d_k)/2.
    f̂(-k) = (c_k + i d_k)/2.
    """
    if k == 0:
        return cp.Constant(0.5), cp.Constant(0.0)
    elif k > 0:
        return 0.5 * c[k - 1], -0.5 * d[k - 1]
    else:
        return 0.5 * c[-k - 1], 0.5 * d[-k - 1]


def _index_in_xi(coord: str, k: int) -> int:
    """Index of (1 / c_k / d_k) in xi = (1, c_1..c_T, d_1..d_T)."""
    if coord == "1":
        return 0
    raise NotImplementedError


def _bilinear_in_M(M: cp.Variable, T_max: int, kind1: str, k1: int,
                   kind2: str, k2: int):
    """Return the cvxpy expression for the product (kind1_k1 * kind2_k2)
    using the lifted (L+1)x(L+1) matrix M, where
        xi = (1, c_1..c_{T_max}, d_1..d_{T_max}),  L = 2T_max+1,
        M[0:L, 0:L] represents xi xi^T,
        M[0:L, L]  represents xi  (the link row),
        M[L, L]    = 1.

    kind in {"1","c","d"}; k=0 only valid for kind="1".
    """
    def idx(kind, k):
        if kind == "1":
            return 0
        elif kind == "c":
            assert 1 <= k <= T_max
            return k       # 1..T_max
        elif kind == "d":
            assert 1 <= k <= T_max
            return T_max + k   # T_max+1 .. 2T_max
        raise ValueError

    i = idx(kind1, k1)
    j = idx(kind2, k2)
    if i == 0 and j == 0:
        return cp.Constant(1.0)
    # Use the linking row to keep linear pieces honest:
    # If one of the indices is 0 (the constant), we want the LINEAR
    # part of xi (i.e., c_k or d_k itself), which by the lift is M[0, j]
    # via M[0, j] = xi_0 * xi_j = 1 * xi_j = xi_j.
    if i == 0:
        return M[0, j]
    if j == 0:
        return M[i, 0]
    return M[i, j]


def _f_hat_re_im(c: cp.Variable, d: cp.Variable, k: int):
    """Real / imaginary parts of f̂(k) as cvxpy LINEAR expressions in (c,d)."""
    return _fourier_real_indexed(c, d, k)


def _f2_hat_re_im(M: cp.Variable, T_max: int, m: int):
    """Real / imaginary parts of (f^2)̂(m) as cvxpy LINEAR expressions in M.

    (f^2)̂(m) = sum_{n: |n|<=T_max and |m-n|<=T_max} f̂(n) f̂(m-n).

    For each n, f̂(n) and f̂(m-n) each have the form (1/2 if index 0,
    else (c_|k| ∓ i d_|k|)/2), so the product is a sum of (1/4)-scaled
    bilinear monomials in (1, c_*, d_*).
    """
    n_lo = max(-T_max, m - T_max)
    n_hi = min(T_max, m + T_max)
    re = cp.Constant(0.0)
    im = cp.Constant(0.0)
    for n in range(n_lo, n_hi + 1):
        # f̂(n) = (Re_n + i Im_n) where:
        #   n=0:  Re=1/2,  Im=0
        #   n>0:  Re=c_n/2, Im=-d_n/2
        #   n<0:  Re=c_|n|/2, Im=+d_|n|/2
        # Same for n2 = m-n.
        n2 = m - n
        # We need (Re_n + i Im_n)(Re_{n2} + i Im_{n2})
        # = (Re_n Re_{n2} - Im_n Im_{n2}) + i(Re_n Im_{n2} + Im_n Re_{n2}).
        # Each of these four products is a scaled (1/4) bilinear:
        # determined by the (kind, k) for each factor.
        def factor(idx):
            if idx == 0:
                return ("1", 0, +1)   # (kind, k, im_sign)
            elif idx > 0:
                return ("pos", idx, -1)
            else:
                return ("neg", -idx, +1)

        kind_a, k_a, sgn_a = factor(n)
        kind_b, k_b, sgn_b = factor(n2)

        # Real factor of f̂(n):  Re_n = 1/2 if kind=="1" else c_{k_a}/2.
        # Imag factor of f̂(n):  Im_n = 0 if kind=="1" else sgn_a * d_{k_a}/2.

        def re_im(kind, k, sgn):
            if kind == "1":
                return ("1", 0, "1", 0, +1.0, +0.0)  # (re_kind, re_k, im_kind, im_k, re_coef, im_coef)
            else:
                return ("c", k, "d", k, +0.5, sgn * 0.5)  # f̂ = c/2 + i sgn d/2

        re_kind_a, re_k_a, im_kind_a, im_k_a, re_coef_a, im_coef_a = re_im(kind_a, k_a, sgn_a)
        re_kind_b, re_k_b, im_kind_b, im_k_b, re_coef_b, im_coef_b = re_im(kind_b, k_b, sgn_b)

        # Each factor's real part = re_coef * (1 or c_{k}); we need it
        # multiplied by the OTHER factor's real or imag part.
        # We'll construct as bilinear products in M.
        # Product Re = Re_a Re_b - Im_a Im_b
        # Product Im = Re_a Im_b + Im_a Re_b.

        def get_lifted(re_kind, re_k):
            return ("1", 0) if re_kind == "1" else ("c", re_k)
        # Real piece of f̂(n): kind=re_kind_a (either "1" or "c"); coef = re_coef_a if not "1", else 0.5.
        # We'll use the convention: the "Re factor" is always either constant 0.5 (if the index is 0)
        # or 0.5*c_{k}.

        # Re part of f̂(n) = kind RE_a, weight w_re_a:
        if kind_a == "1":
            re_a_kind = "1"; re_a_k = 0; w_re_a = 0.5
            im_a_kind = "1"; im_a_k = 0; w_im_a = 0.0
        else:
            re_a_kind = "c"; re_a_k = k_a; w_re_a = 0.5
            im_a_kind = "d"; im_a_k = k_a; w_im_a = sgn_a * 0.5

        if kind_b == "1":
            re_b_kind = "1"; re_b_k = 0; w_re_b = 0.5
            im_b_kind = "1"; im_b_k = 0; w_im_b = 0.0
        else:
            re_b_kind = "c"; re_b_k = k_b; w_re_b = 0.5
            im_b_kind = "d"; im_b_k = k_b; w_im_b = sgn_b * 0.5

        # Real of product: w_re_a * w_re_b * <re_a, re_b>  -  w_im_a * w_im_b * <im_a, im_b>
        # Imag of product: w_re_a * w_im_b * <re_a, im_b>  +  w_im_a * w_re_b * <im_a, re_b>

        # Each <X, Y> is _bilinear_in_M for the appropriate (kind, k).
        if w_re_a != 0.0 and w_re_b != 0.0:
            re = re + (w_re_a * w_re_b) * _bilinear_in_M(M, T_max, re_a_kind, re_a_k, re_b_kind, re_b_k)
        if w_im_a != 0.0 and w_im_b != 0.0:
            re = re - (w_im_a * w_im_b) * _bilinear_in_M(M, T_max, im_a_kind, im_a_k, im_b_kind, im_b_k)
        if w_re_a != 0.0 and w_im_b != 0.0:
            im = im + (w_re_a * w_im_b) * _bilinear_in_M(M, T_max, re_a_kind, re_a_k, im_b_kind, im_b_k)
        if w_im_a != 0.0 and w_re_b != 0.0:
            im = im + (w_im_a * w_re_b) * _bilinear_in_M(M, T_max, im_a_kind, im_a_k, re_b_kind, re_b_k)

    return re, im


def add_lasserre2_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                              T_max: int, T_loc: int = None):
    """Add the Lasserre level-2 augmentation for f^2 <= f.

    Parameters
    ----------
    cons : list
        cvxpy constraint list to extend in place.
    c, d : cvxpy.Variable
        Real Fourier coefficient variables of length T (T >= T_max).
    T_max : int
        Lasserre lifting cutoff: bilinear products f̂(j)f̂(k) considered
        for |j|,|k| <= T_max.  The lift moment matrix has size 2T_max+2
        (in the Schur-complement form).
    T_loc : int, optional
        Order of the Hermitian localizing matrix Loc[j,k] = (f-f^2)̂(j-k)
        for j,k = 0..T_loc.  Defaults to T_max.

    Returns
    -------
    M_top : cvxpy.Variable
        The (2T_max+1)x(2T_max+1) lifted moment matrix (for diagnostics).
    """
    if T_loc is None:
        T_loc = T_max
    assert T_loc <= T_max, "T_loc must be <= T_max for the truncation to be consistent."
    T = c.shape[0]
    assert T_max <= T, f"T_max={T_max} must be <= T={T}."

    L = 2 * T_max + 1   # length of xi = (1, c_1..c_T_max, d_1..d_T_max)

    # Build xi as a cvxpy expression of length L.
    xi = cp.hstack([cp.Constant(1.0), c[:T_max], d[:T_max]])

    # Lifted symmetric moment matrix M_top of size LxL, plus the bordered
    # full matrix of size (L+1)x(L+1).
    # Schur complement form: [[ M_top, xi ], [ xi^T, 1 ]] ⪰ 0   ⇔   M_top ⪰ xi xi^T.
    M_top = cp.Variable((L, L), symmetric=True)

    # Bordered matrix block.
    border = cp.bmat([
        [M_top, cp.reshape(xi, (L, 1), order="C")],
        [cp.reshape(xi, (1, L), order="C"), cp.Constant(np.array([[1.0]]))],
    ])
    cons.append(border >> 0)

    # Pin M_top[0, 0] = 1 (since xi_0 = 1).  This is implied by the above PSD,
    # but adding it as an equality tightens the SDP.
    cons.append(M_top[0, 0] == 1.0)
    # Pin M_top[0, i] = xi_i for all i (link row), again implied but useful.
    for i in range(1, T_max + 1):
        cons.append(M_top[0, i] == c[i - 1])
    for i in range(1, T_max + 1):
        cons.append(M_top[0, T_max + i] == d[i - 1])

    # ----- Localizing matrix for (f - f^2) >= 0 ---------------------------
    # Loc[j, k] = (f - f^2)̂(j - k) for j, k = 0..T_loc.
    # f̂(j-k) is linear in (c, d).  (f^2)̂(j-k) is bilinear in (c, d), thus
    # linear in M_top after the lift.  Build the Hermitian PSD constraint
    # via real-form 2(T_loc+1) x 2(T_loc+1) block matrix.
    Re_rows = []
    Im_rows = []
    for j in range(T_loc + 1):
        re_row = []
        im_row = []
        for k in range(T_loc + 1):
            ell = j - k
            f_re, f_im = _f_hat_re_im(c, d, ell)
            f2_re, f2_im = _f2_hat_re_im(M_top, T_max, ell)
            re_row.append(f_re - f2_re)
            im_row.append(f_im - f2_im)
        Re_rows.append(re_row)
        Im_rows.append(im_row)
    Re_M = cp.bmat(Re_rows)
    Im_M = cp.bmat(Im_rows)
    real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
    cons.append(real_form >> 0)

    return M_top


# ----- self-test --------------------------------------------------------------
if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")

    # Sanity check 1: For f(x) = 1/2 + a cos(πx), test that the lift gives the
    # correct value of (f^2)̂(0) = 1/4 + a^2/2.
    print("=== Self-test 1: f = 1/2 + a cos(πx),  a = 0.4 ===")
    T = 3
    T_max = 3
    c = cp.Variable(T)
    d = cp.Variable(T)

    # Hand-compute (f^2)̂(0) for a=0.4: (1/4) + (0.4)^2/2 = 0.25 + 0.08 = 0.33.
    # Build the lift programmatically and check.
    cons = []
    M_top = add_lasserre2_constraint(cons, c, d, T_max=T_max, T_loc=2)
    # Pin c, d to a=0.4, b=0:
    cons.append(c[0] == 0.4)
    cons.append(c[1] == 0.0)
    cons.append(c[2] == 0.0)
    cons.append(d[0] == 0.0)
    cons.append(d[1] == 0.0)
    cons.append(d[2] == 0.0)
    # And M_top to the rank-1 lift of xi = (1, 0.4, 0, 0, 0, 0, 0):
    xi_val = np.array([1.0, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0])
    cons.append(M_top == np.outer(xi_val, xi_val))

    # Objective: minimize the (0,0) entry of (f - f^2) localizing matrix
    # = f̂(0) - (f^2)̂(0) = 1/2 - (1/4 + 0.08) = 0.17.
    f_re, _ = _f_hat_re_im(c, d, 0)
    f2_re, _ = _f2_hat_re_im(M_top, T_max, 0)
    obj = cp.Minimize(f_re - f2_re)
    prob = cp.Problem(obj, cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  f̂(0) - (f^2)̂(0) computed = {prob.value:.6f}   "
          f"expected = {0.5 - 0.25 - 0.08:.6f}")
    print(f"  status = {prob.status}")

    # Sanity check 2: Same but for a=0.6 (f<0 at x=1, so f^2 > f at x=1).
    # The localizing matrix Loc should NOT be PSD, and so the Lasserre-2
    # constraint should make the program infeasible if we ALSO require
    # the lift to match (c_1=0.6).  Test that.
    print()
    print("=== Self-test 2: f = 1/2 + 0.6 cos(πx)  → NEG at x=1 ===")
    T = 3
    T_max = 3
    c = cp.Variable(T)
    d = cp.Variable(T)
    cons = []
    M_top = add_lasserre2_constraint(cons, c, d, T_max=T_max, T_loc=3)
    cons.append(c[0] == 0.6)
    cons.append(c[1] == 0.0)
    cons.append(c[2] == 0.0)
    cons.append(d[0] == 0.0)
    cons.append(d[1] == 0.0)
    cons.append(d[2] == 0.0)
    xi_val = np.array([1.0, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0])
    cons.append(M_top == np.outer(xi_val, xi_val))

    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  Solver status with f = 1/2 + 0.6 cos(πx) (NOT a true density): {prob.status}")
    print(f"  Expected: 'infeasible' or 'infeasible_inaccurate'.")
