r"""Half-integer (period-4) Bochner constraint  —  transplant from Kim & Pilanci
(arXiv:2606.31182, Appendix A / Algorithm 1 lines 20-21), plus a strict improvement.

WHAT THIS IS
------------
White's program already carries, for m = 1..2R, the quantities `a_m`, `b_m`
(white_full_convex.py lines ~157-176):

    a_{2k} = c_k/2,   b_{2k} = d_k/2                                  (even m)
    a_m    = eps_r + (2m sin(pi m/2)/pi)(1/(2m^2) + sum_k (-1)^k c_k/(m^2-4k^2))
    b_m    = dlt_r + (4 sin(pi m/2)/pi) sum_k k(-1)^k d_k/(m^2-4k^2)  (odd m)
    |eps_r| <= tail_bound_eps,  |dlt_r| <= tail_bound_delta

VERIFIED numerically (10 digits, m = 1..7, both parities) that these are

    a_m = (1/2) \int_{-1}^{1} f(x) cos(m pi x/2) dx
    b_m = (1/2) \int_{-1}^{1} f(x) sin(m pi x/2) dx

i.e. the HALF-INTEGER frequency samples of f-hat: the period-4 Fourier
coefficients of F := f * 1_{[-1,1]}, extended by zero to [-2,2] and then
4-periodically.  (The repo's c_k, d_k are the period-2 coefficients,
c_k = \int_{-1}^1 f cos(k pi x) dx; this is confirmed independently by the
program's own bounds |c_k| <= 2/pi = max_{|S|=1} \int_S cos(pi x) dx and
sum(c^2+d^2) <= 1/2 = Parseval with \int f^2 <= \int f = 1.  NOTE:
`_fourier_convention_notes.md` section 3 states a different, inconsistent
convention and is wrong on this point.)

Write A_0 := 1/4 and A_m := (a_m - i b_m)/2 for m >= 1.

THE CONSTRAINT.  For any n <= 2R, with

    H     := [A_{j-k}]_{j,k=0..n}            (Hermitian Toeplitz, diagonal 1/4)
    Theta := [theta_{j-k}],  theta_0 = 1/2,  theta_l = sin(l pi/2)/(l pi)

(Theta = the Toeplitz matrix of the indicator 1_{[-1,1]} on the period-4 torus;
it vanishes at every even lag >= 2), every admissible f satisfies

    (i)   H         >= 0
    (ii)  Theta - H >= 0

PROOF.  (i) F >= 0 pointwise on R/4Z, so Bochner gives H >= 0.
(ii) 1_{[-1,1]} - F >= 0 pointwise (f <= 1 on [-1,1], both vanish off it);
the moment map is linear, so its Toeplitz matrix is Theta - H.  QED.

The tail slacks eps, dlt do not threaten validity: the true f's exact a_m
satisfy the transfer relation with |eps| inside the program's bound, so the
true point remains feasible.  Slack only weakens the cut.

RELATION TO WHAT THE REPO ALREADY HAS
-------------------------------------
`bochner.py` builds the Toeplitz matrix indexed by c_l, d_l at INTEGER lag l.
Verified exactly (max abs diff 0.0 and 4e-17):

    repo's M(f)   == 2 * (even-index principal submatrix of H)
    repo's M(1-f) == 2 * (even-index principal submatrix of Theta - H)

the second identity working precisely because theta vanishes at even lags >= 2.
So the repo has been imposing Bochner only on the EVEN frequency sublattice.
Every odd row/column of H is new.

RELATION TO KIM & PILANCI
-------------------------
Their Eq. (15)+(16) is exactly (i) [their T_f = H].  Their second constraint is
I - T_f >= 0, i.e. Theta replaced by I.  Since I - Theta = Toeplitz(1_{1<|x|<2})
>= 0, constraint (ii) here is STRICTLY STRONGER than theirs: they relax
"f <= 1 on [-1,1]" to "F <= 1 on [-2,2]", which is pure slack on half the torus.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp

from white_full_convex import odd_coeff_factors


def a_b_exprs(c, d, eps, dlt, T: int, R: int):
    """Rebuild White's half-integer coefficients a_m, b_m for m = 1..2R.

    Mirrors white_full_convex.build_problem lines ~157-176 exactly.
    """
    A, B = [], []
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            A.append(0.5 * c[m // 2 - 1])
            B.append(0.5 * d[m // 2 - 1])
        else:
            af, bf = odd_coeff_factors(m, T)
            s = np.sin(np.pi * m / 2)
            A.append(eps[(m - 1) // 2]
                     + (2 * m * s / np.pi) * (1.0 / (2 * m ** 2)
                                              + cp.sum(cp.multiply(af, c))))
            B.append(dlt[(m - 1) // 2] + (4 * s / np.pi) * cp.sum(cp.multiply(bf, d)))
    return A, B


def theta(ell: int) -> float:
    """Period-4 Fourier coefficient of 1_{[-1,1]} on [-2,2]."""
    return 0.5 if ell == 0 else float(np.sin(ell * np.pi / 2) / (ell * np.pi))


def add_halfint_bochner(cons: list, A, B, n: int, complement: bool = True):
    """Append H >= 0 and (optionally) Theta - H >= 0, as real-form PSD blocks.

    `A`, `B` come from `a_b_exprs`; requires n <= len(A) = 2R.
    Set complement=False to impose only H >= 0 (the literal Kim-Pilanci cut,
    minus their weaker I - H).
    """
    assert n <= len(A), f"n={n} exceeds 2R={len(A)}"
    modes = ["f", "comp"] if complement else ["f"]
    for which in modes:
        Re_rows, Im_rows = [], []
        for j in range(n + 1):
            re_row, im_row = [], []
            for k in range(n + 1):
                ell = j - k
                al = abs(ell)
                re = cp.Constant(0.25) if al == 0 else 0.5 * A[al - 1]
                im = (cp.Constant(0.0) if al == 0
                      else (-0.5 * B[al - 1] if ell > 0 else 0.5 * B[al - 1]))
                if which == "comp":
                    re = cp.Constant(theta(al)) - re
                    im = -im
                re_row.append(re)
                im_row.append(im)
            Re_rows.append(re_row)
            Im_rows.append(im_row)
        Re_M, Im_M = cp.bmat(Re_rows), cp.bmat(Im_rows)
        cons.append(cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]]) >> 0)


# ----- validity self-test: exact 50-digit check on admissible extreme points ----
if __name__ == "__main__":
    import mpmath as mp
    mp.mp.dps = 50
    pi = mp.pi

    def coeffs(ivals, n):
        a, b = [], []
        for m in range(1, n + 1):
            s = t = mp.mpf(0)
            for (al, be) in ivals:
                s += (mp.sin(m * pi * be / 2) - mp.sin(m * pi * al / 2)) / (m * pi)
                t += (mp.cos(m * pi * al / 2) - mp.cos(m * pi * be / 2)) / (m * pi)
            a.append(s); b.append(t)
        return a, b

    def mineig(n, a, b, comp):
        M = mp.zeros(2 * (n + 1))
        for j in range(n + 1):
            for k in range(n + 1):
                l = j - k; al = abs(l)
                re = mp.mpf(1) / 4 if al == 0 else a[al - 1] / 2
                im = mp.mpf(0) if al == 0 else (-b[al - 1] / 2 if l > 0 else b[al - 1] / 2)
                if comp:
                    re = (mp.mpf(1) / 2 if al == 0 else mp.sin(al * pi / 2) / (al * pi)) - re
                    im = -im
                M[j, k] = re; M[n + 1 + j, n + 1 + k] = re
                M[j, n + 1 + k] = -im; M[n + 1 + j, k] = im
        return min(mp.eigsy(M, eigvals_only=True))

    M = mp.mpf
    CASES = {
        "[-1/2,1/2]":        [(M(-1) / 2, M(1) / 2)],
        "[0,1]":             [(M(0), M(1))],
        "[-.9,-.4]u[.1,.6]": [(M('-0.9'), M('-0.4')), (M('0.1'), M('0.6'))],
        "3 pieces asym":     [(M('-1'), M('-0.7')), (M('-0.2'), M('0.25')),
                              (M('0.55'), M('0.8'))],
        "5 pieces":          [(M('-0.95'), M('-0.8')), (M('-0.6'), M('-0.35')),
                              (M('-0.1'), M('0.15')), (M('0.4'), M('0.6')),
                              (M('0.75'), M('0.9'))],
    }
    n = 24
    print(f"50-digit validity test, matrix size {n+1}. Both columns must be >= 0.")
    print(f"{'S with |S|=1':24s} {'min eig H':>16s} {'min eig Theta-H':>18s}")
    for name, iv in CASES.items():
        a, b = coeffs(iv, n)
        print(f"{name:24s} {mp.nstr(mineig(n,a,b,False),8):>16s} "
              f"{mp.nstr(mineig(n,a,b,True),8):>18s}")
