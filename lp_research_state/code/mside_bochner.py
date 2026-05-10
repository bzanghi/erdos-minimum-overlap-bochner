"""
M-side Bochner-PSD (SOC-relaxed) constraint for M ≥ 0.

Lemma 2 (White 2023, arXiv:2201.05704): for the autocorrelation-like nonneg
measure M arising from f, M̂(0) = Ω/2 and
    M̂(m) = a_m f̂(m) − 4 |f̂(m)|²    for m ≥ 1,
where a_m = (4/(mπ)) sin(mπ/2)  (= 0 for even m, ±4/(mπ) for odd m).

Bochner ⇒ M ≥ 0 ⟺ ∀n the (n+1)×(n+1) Hermitian Toeplitz matrix
    T_M(c, d) := [M̂(j-k)]_{j,k=0..n}
is positive semidefinite. Because M̂(m) is QUADRATIC in (c, d) through the
−4|f̂(m)|² term, T_M ⪰ 0 is non-convex in the LP variables.

SOC RELAXATION (validity proof in findings.md, dated 2026-05-10).
Introduce SOC slack U_m ≥ |f̂(m)|² = (c_m² + d_m²)/4 (a single second-order
cone constraint per m). Replace |f̂(m)|² with U_m to get a linear-in-(c,d,U)
Toeplitz entry, and impose the resulting Hermitian-Toeplitz PSD constraint.

Let
    F_0 = White's LP feasible set,
    F_1 = F_0 ∩ {(c,d) : T_M(c, d) ⪰ 0}                    (true tightening),
    F_2 = F_0 ∩ {(c,d) : ∃ U ≥ |f̂|² with T_relax(c,d,U) ⪰ 0}  (SOC relax).
For any (c,d) ∈ F_1, set U_m = |f̂(m)|² to get T_relax = T_M ⪰ 0, so
F_1 ⊆ F_2; F_2 ⊆ F_0 by construction. Hence
    min_F0 Ω ≤ min_F2 Ω ≤ min_F1 Ω.
The middle quantity is a RIGOROUS lower bound on µ that is no smaller than the
plain White LP optimum. F_2 may be strictly larger than F_1 (and the bound
correspondingly looser) when, at the F_2 optimum, U_m > |f̂(m)|² is needed.

Encoded analogously to bochner.py: a complex Hermitian (n+1)×(n+1) PSD
constraint becomes a real-symmetric (2n+2)×(2n+2) PSD constraint via
    [[Re_M, -Im_M],
     [Im_M,  Re_M]] ⪰ 0.

Inputs to add_mside_bochner_constraint:
    cons   — the cvxpy constraint list to extend (in place)
    c, d   — Fourier-coefficient cvxpy.Variable arrays of length T
    Omega  — the LP objective scalar variable (M̂(0) = Ω/2)
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


def add_mside_bochner_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                                 Omega: cp.Variable, n_M: int):
    """Append the SOC-relaxed M-side Bochner PSD constraint of order n_M.

    Adds n_M SOC slack variables U_m ≥ |f̂(m)|², n_M scalar SOC constraints,
    and one (2n_M+2)×(2n_M+2) real-form PSD constraint encoding the relaxed
    Hermitian-Toeplitz PSD condition on T_relax(c, d, U).
    """
    T = c.shape[0]
    assert n_M <= T, f"M-side Bochner order n_M={n_M} cannot exceed T={T}"
    if n_M <= 0:
        return

    # SOC slack variables: U_m ≥ |f̂(m)|² = (c_m² + d_m²)/4
    # i.e., cp.sum_squares([c[m-1], d[m-1]]) ≤ 4 * U[m-1].
    U = cp.Variable(n_M, nonneg=True)
    for m in range(1, n_M + 1):
        cons.append(
            cp.sum_squares(cp.hstack([c[m - 1], d[m - 1]])) <= 4.0 * U[m - 1]
        )

    # Build the relaxed Hermitian Toeplitz M̂_relax(j-k):
    #   M̂_relax(0)  = Ω/2
    #   M̂_relax(+m) = (a_m/2) c_{m-1} - 4 U_{m-1}  - i (a_m/2) d_{m-1}    (m > 0)
    #   M̂_relax(-m) = conj(M̂_relax(+m))                                  (m > 0)
    # In code the m-th Fourier coefficient is c[m-1] (zero-based), and U[m-1]
    # is the SOC slack for that lag.
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
                # Real part is identical for ell > 0 and ell < 0 (Hermitian).
                re_entry = (am / 2.0) * c[m - 1] - 4.0 * U[m - 1]
                # Imaginary part flips sign with ell.
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


# ----- self-test (math sanity, no solve) ----------------------------------
if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    # Print a_m for m=1..6 to sanity-check the ±4/(mπ)/0 pattern.
    print("a_m for m=1..6:")
    for m in range(1, 7):
        print(f"  a_{m} = {a_m(m): .6f}   (expected: "
              f"{4.0/(m*np.pi)*np.sin(m*np.pi/2):+.6f})")
    # m=1: +4/π = +1.2732
    # m=2: 0
    # m=3: -4/(3π) = -0.4244
    # m=4: 0
    # m=5: +4/(5π) = +0.2546
    # m=6: 0
    print("\nVariable signature check:")
    T = 5
    c = cp.Variable(T)
    d = cp.Variable(T)
    Omega = cp.Variable()
    cons = []
    add_mside_bochner_constraint(cons, c, d, Omega, n_M=3)
    print(f"  added {len(cons)} constraints "
          f"({3} SOC + 1 PSD = 4 expected): "
          f"{'OK' if len(cons) == 4 else 'MISMATCH'}")
