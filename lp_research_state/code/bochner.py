"""
Bochner-PSD constraints for f ≥ 0 and 1-f ≥ 0 (i.e., f ∈ [0,1]).

Theorem (Bochner). For a real measurable f: [-1,1] → R extended periodically,
f ≥ 0 a.e. iff for every n ≥ 0, the Hermitian moment matrix
    Toeplitz_n(f) := [ f̂(j-k) ]_{j,k=0..n}
is positive semidefinite. With our Fourier convention
    f̂(0) = 1/2,  f̂(k) = (c_k - i d_k)/2,  f̂(-k) = (c_k + i d_k)/2,
the matrix is Hermitian Toeplitz with first row (1/2, (c_1+id_1)/2, ...,
(c_n+id_n)/2).

Symmetric application: 1-f ≥ 0 gives the analogous PSD constraint with the
off-diagonal entries negated (since (1-f)̂(0) = 1/2, (1-f)̂(k) = -(c_k-id_k)/2).

This is a MUCH stronger constraint than just f² ≤ f tested against single
trig polynomials (the latter is implied by Bochner).

We encode the (n+1)x(n+1) complex Hermitian PSD constraint as a (2n+2)x(2n+2)
real symmetric PSD constraint:
    M_real = [[Re(M), -Im(M)], [Im(M), Re(M)]]   ⪰ 0

Inputs to add_bochner_constraints:
    cons     — the cvxpy constraint list
    c, d     — Fourier coefficient variables (length T)
    n        — Bochner order (matrix size n+1)
    sign     — +1 for f≥0, -1 for 1-f≥0
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def add_bochner_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                           n: int, sign: int = +1):
    """Append the (n+1)x(n+1) Hermitian PSD constraint for f or 1-f."""
    T = c.shape[0]
    assert n <= T, f"Bochner order n={n} cannot exceed T={T}"

    # Build Re(M) (symmetric) and Im(M) (antisymmetric) of size (n+1)×(n+1).
    # Convention: M[j, k] = f̂(j - k).
    # f̂(0) = 1/2.
    # For ℓ > 0: f̂(ℓ) = (c_ℓ - i d_ℓ)/2,   f̂(-ℓ) = (c_ℓ + i d_ℓ)/2.
    # So Re M[j,k] = (1/2 if j==k else c_{|j-k|}/2);
    #    Im M[j,k] = (0 if j==k else (k-j > 0 ? +d_{|j-k|}/2 : -d_{|j-k|}/2)).
    #    Sign of imaginary part: M[j,k] = f̂(j-k); j-k > 0 → -i d_{|j-k|}/2; j-k < 0 → +i d_{|j-k|}/2.

    Re_blocks = []
    Im_blocks = []

    # We construct row by row; entries of Re M are linear in (c) (and a constant 1/2 on diag);
    # entries of Im M are linear in (d).
    # Use cvxpy bmat-like construction: each entry is a (real or zero) cvxpy expression.

    Re_rows = []
    Im_rows = []
    half = 0.5
    for j in range(n + 1):
        re_row = []
        im_row = []
        for k in range(n + 1):
            ell = j - k
            if ell == 0:
                re_row.append(cp.Constant(half))
                im_row.append(cp.Constant(0.0))
            else:
                aell = abs(ell)
                # Re part = sign * c_{aell}/2  (off-diagonals get sign factor for 1-f case).
                re_row.append(cp.Constant(sign * 0.5) * c[aell - 1])
                # Im part: ell > 0 gives -d/2; ell < 0 gives +d/2; with sign factor.
                if ell > 0:
                    im_row.append(cp.Constant(-sign * 0.5) * d[aell - 1])
                else:
                    im_row.append(cp.Constant(+sign * 0.5) * d[aell - 1])
        Re_rows.append(re_row)
        Im_rows.append(im_row)

    Re_M = cp.bmat(Re_rows)
    Im_M = cp.bmat(Im_rows)

    # Real-form 2(n+1) x 2(n+1) PSD constraint:
    #     [[ Re_M, -Im_M ],
    #      [ Im_M,  Re_M ]]
    # is real symmetric PSD iff complex Hermitian Re_M + i Im_M is PSD.
    real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
    cons.append(real_form >> 0)


# ----- self-test ---------------------------------------------------------
if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    # Sanity: a real f(x) = 1/2 + a cos(πx) + b sin(πx). Take (a, b) = (0.4, 0).
    # The Bochner constraint at level n=2 should be satisfied if this is a valid
    # nonneg density; not necessarily satisfied otherwise.

    # Build a small test: c, d fixed, check PSD by hand.
    a_, b_ = 0.4, 0.0
    n = 2
    half = 0.5
    Re = np.array([
        [half, a_ / 2, 0.0],     # j=0; k=0,1,2; |ell|=0,1,2
        [a_ / 2, half, a_ / 2],
        [0.0, a_ / 2, half],
    ])
    Im = np.array([
        [0.0, +b_ / 2, +0.0 / 2],
        [-b_ / 2, 0.0, +b_ / 2],
        [-0.0 / 2, -b_ / 2, 0.0],
    ])
    real_form = np.block([[Re, -Im], [Im, Re]])
    print("Test: f(x) = 1/2 + 0.4 cos(πx)")
    print(f"  6x6 real form eigenvalues: {np.sort(np.linalg.eigvalsh(real_form))}")
    print(f"  → PSD: {(np.linalg.eigvalsh(real_form) >= -1e-12).all()}")
    print(f"  (note: f(x) at x = 1: 1/2 - 0.4 = 0.1 ≥ 0; at x = 0: 0.9 ≥ 0; nonneg ✓)")

    # Now test a NEGATIVE-violating example: a = 0.6 (f(1) = 1/2 - 0.6 = -0.1 < 0).
    a_ = 0.6
    Re[0, 1] = Re[1, 0] = Re[1, 2] = Re[2, 1] = a_ / 2
    real_form = np.block([[Re, -Im], [Im, Re]])
    eigs = np.linalg.eigvalsh(real_form)
    print()
    print("Test: f(x) = 1/2 + 0.6 cos(πx)  → NEGATIVE at x=1")
    print(f"  6x6 real form eigenvalues: {np.sort(eigs)}")
    print(f"  → PSD: {(eigs >= -1e-12).all()}  (expected False, since f<0 at x=1)")
