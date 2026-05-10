"""
M-side Bochner-PSD constraint via EXACT Schur-complement lifting.

Lemma 2 (White 2023, arXiv:2201.05704): for the autocorrelation-like nonneg
measure M arising from f, M̂(0) = Ω/2 and
    M̂(m) = a_m f̂(m) − 4 |f̂(m)|²    for m ≥ 1,
where a_m = (4/(mπ)) sin(mπ/2)  (= 0 for even m, ±4/(mπ) for odd m).

Bochner ⇒ M ≥ 0 ⟺ ∀n the (n+1)×(n+1) Hermitian Toeplitz matrix
    T_M(c, d) := [M̂(j-k)]_{j,k=0..n}
is positive semidefinite. Because M̂(m) is QUADRATIC in (c, d) through the
−4|f̂(m)|² term, T_M ⪰ 0 is non-convex in the LP variables.

EXACT SCHUR LIFTING.
For each m=1..n_M introduce a real scalar slack s_m and impose the 3×3 real-
symmetric PSD Schur complement
    [[ s_m,           Re f̂(m),  −Im f̂(m) ],
     [ Re f̂(m),       1,        0        ],
     [ −Im f̂(m),      0,        1        ]] ⪰ 0
which (taking 2×2 Schur on the bottom-right identity block) is equivalent to
    s_m ≥ (Re f̂(m))² + (Im f̂(m))² = |f̂(m)|².
With f̂(m) = (c_m − i d_m)/2, Re f̂(m) = c_m/2, Im f̂(m) = −d_m/2, this is
    s_m ≥ (c_m² + d_m²)/4 = |f̂(m)|².
Mathematically this defines the SAME convex set as the SOC inequality used in
mside_bochner.py; the Schur form is simply written as one PSD constraint.

SIGN/RELAXATION DIRECTION.
Replace |f̂(m)|² with the over-estimate s_m in the M-Toeplitz off-diagonals:
    M̂_relax(m) = a_m f̂(m) − 4 s_m,    s_m ≥ |f̂(m)|².
Define the convex feasible set
    F_2 = F_0 ∩ {(c,d,s) : s_m ≥ |f̂(m)|² ∀m, T_relax(c,d,s) ⪰ 0}.
Take any (c,d) in the TRUE non-convex feasible set
    F_1 = F_0 ∩ {(c,d) : T_M(c,d) ⪰ 0}.
Set s_m = |f̂(m)|² (smallest feasible). Then T_relax = T_M ⪰ 0, so
(c,d,s) ∈ F_2; hence F_1 ⊆ F_2. Together with F_2 ⊆ F_0 this gives a VALID
lower bound:  min_F0 Ω ≤ min_F2 Ω ≤ min_F1 Ω ≤ µ.

Imposing T_relax ⪰ 0 with s_m ≥ |f̂(m)|² is one-sided in the right direction
(the LP optimum will WANT to push s_m down to |f̂(m)|² to relax the PSD
constraint, so the binding case is exactly s_m = |f̂(m)|² which is what we
want).

Encoding analogous to bochner.py: the complex Hermitian (n+1)×(n+1) PSD
constraint becomes a real-symmetric (2n+2)×(2n+2) PSD constraint via
    [[ Re_M, -Im_M ],
     [ Im_M,  Re_M ]] ⪰ 0.

Inputs to add_mside_bochner_schur_constraint:
    cons   — cvxpy constraint list to extend in place
    c, d   — Fourier-coefficient cvxpy.Variable arrays of length T
    Omega  — LP objective scalar variable (M̂(0) = Ω/2)
    n_M    — Toeplitz order (matrix size n_M+1; 0 = no constraint)
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def a_m(m: int) -> float:
    """Lemma-2 coefficient a_m = (4/(mπ)) sin(mπ/2).

    Vanishes for even m; equals (4/(mπ))·(−1)^((m−1)/2) for odd m.
    """
    return 4.0 / (m * np.pi) * np.sin(m * np.pi / 2.0)


def add_mside_bochner_schur_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                                       Omega: cp.Variable, n_M: int):
    """Append the Schur-lifted M-side Bochner PSD constraint of order n_M.

    Adds:
      • n_M scalar slack variables s_m,
      • n_M small 3x3 PSD Schur constraints enforcing s_m ≥ |f̂(m)|²,
      • one (2n_M+2)×(2n_M+2) real-form PSD constraint encoding the relaxed
        Hermitian-Toeplitz PSD condition T_relax(c,d,s) ⪰ 0.
    """
    T = c.shape[0]
    assert n_M <= T, f"M-side Bochner order n_M={n_M} cannot exceed T={T}"
    if n_M <= 0:
        return

    # Schur complement slack variables.
    # s_m ≥ |f̂(m)|² = (c_m² + d_m²)/4.
    # We write the Schur block as a 3x3 real-symmetric PSD matrix
    #   [[ s_m,    c_m/2,   -(-d_m/2) ],
    #    [ c_m/2,  1,        0        ],
    #    [ -(-d_m/2), 0,     1        ]]
    # = [[s_m, c_m/2, d_m/2], [c_m/2, 1, 0], [d_m/2, 0, 1]] ⪰ 0
    # whose Schur complement on the lower-right 2x2 identity gives
    #   s_m - [c_m/2, d_m/2] · [c_m/2, d_m/2]^T = s_m - (c_m²+d_m²)/4 ≥ 0,
    # i.e., s_m ≥ |f̂(m)|². ✓
    s = cp.Variable(n_M, nonneg=True)
    for m in range(1, n_M + 1):
        block = cp.bmat([
            [cp.reshape(s[m - 1], (1, 1), order='C'),
             cp.reshape(c[m - 1] / 2.0, (1, 1), order='C'),
             cp.reshape(d[m - 1] / 2.0, (1, 1), order='C')],
            [cp.reshape(c[m - 1] / 2.0, (1, 1), order='C'),
             cp.Constant(np.array([[1.0]])),
             cp.Constant(np.array([[0.0]]))],
            [cp.reshape(d[m - 1] / 2.0, (1, 1), order='C'),
             cp.Constant(np.array([[0.0]])),
             cp.Constant(np.array([[1.0]]))],
        ])
        cons.append(block >> 0)

    # Build the relaxed Hermitian Toeplitz M̂_relax(j-k):
    #   M̂_relax(0)  = Ω/2
    #   M̂_relax(+m) = (a_m/2) c_{m-1} - 4 s_{m-1}  - i (a_m/2) d_{m-1}    (m > 0)
    #   M̂_relax(-m) = conj(M̂_relax(+m))                                  (m > 0)
    Re_rows = []
    Im_rows = []
    for j in range(n_M + 1):
        re_row = []
        im_row = []
        for k in range(n_M + 1):
            ell = j - k
            if ell == 0:
                re_row.append(0.5 * Omega)
                im_row.append(cp.Constant(0.0))
            else:
                m = abs(ell)
                am = a_m(m)
                re_entry = (am / 2.0) * c[m - 1] - 4.0 * s[m - 1]
                if ell > 0:
                    im_entry = -(am / 2.0) * d[m - 1]
                else:
                    im_entry = (am / 2.0) * d[m - 1]
                re_row.append(re_entry)
                im_row.append(im_entry)
        Re_rows.append(re_row)
        Im_rows.append(im_row)

    Re_M = cp.bmat(Re_rows)
    Im_M = cp.bmat(Im_rows)
    # Hermitian PSD ⟺ real-symmetric block matrix PSD.
    real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
    cons.append(real_form >> 0)


# ----- self-test (math sanity, no solve) ---------------------------------
if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    print("a_m for m=1..6:")
    for m in range(1, 7):
        print(f"  a_{m} = {a_m(m): .6f}")
    # Sanity: verify the Schur block has the intended Schur complement.
    # Pick (c_1, d_1) = (0.6, 0.2). |f̂(1)|² = (0.36 + 0.04)/4 = 0.1.
    # Set s_1 = 0.1 (boundary). Block should be PSD with one zero eigenvalue.
    s_val = 0.1
    c1, d1 = 0.6, 0.2
    block = np.array([
        [s_val, c1/2, d1/2],
        [c1/2,  1.0,  0.0],
        [d1/2,  0.0,  1.0],
    ])
    eigs = np.linalg.eigvalsh(block)
    print(f"\nSchur block at s=|f̂|²: eigenvalues = {eigs}")
    print(f"  smallest = {eigs[0]:.3e}  (expected ≈ 0, s_1=|f̂(1)|²=0.1)")

    # Now s_val = 0.05 (s < |f̂|²): block should NOT be PSD.
    block[0, 0] = 0.05
    eigs = np.linalg.eigvalsh(block)
    print(f"Schur block at s=0.05 < |f̂|²=0.1: eigenvalues = {eigs}")
    print(f"  smallest = {eigs[0]:.3e}  (expected NEGATIVE)")

    # cvxpy variable signature smoke test.
    print("\ncvxpy build smoke test:")
    T = 5
    c = cp.Variable(T)
    d = cp.Variable(T)
    Omega = cp.Variable()
    cons = []
    add_mside_bochner_schur_constraint(cons, c, d, Omega, n_M=3)
    print(f"  added {len(cons)} constraints "
          f"(3 Schur PSD + 1 Toeplitz PSD = 4 expected): "
          f"{'OK' if len(cons) == 4 else 'MISMATCH'}")
