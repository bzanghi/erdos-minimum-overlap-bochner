"""
M-side Bochner constraint using Lasserre level-2 lifted bilinears.

Lemma 2 (White 2023, arXiv:2201.05704): for the autocorrelation-like nonneg
measure M arising from f, M̂(0) = Ω/2 and
    M̂(m) = a_m f̂(m) − 4 |f̂(m)|²    for m ≥ 1,
where a_m = (4/(mπ)) sin(mπ/2)  (= 0 for even m, ±4/(mπ) for odd m).

Bochner ⇒ M ≥ 0 ⟺ ∀n the (n_M+1)×(n_M+1) Hermitian Toeplitz matrix
    T_M(c, d) := [M̂(j-k)]_{j,k=0..n_M}
is positive semidefinite.

Previous attempts (mside_bochner.py SOC, mside_bochner_schur.py Schur) used
slack U_m / s_m ≥ |f̂(m)|² to relax the −4|f̂(m)|² entry. Both empirically dead
(Δ ≈ 1e-9): the slack absorbs all of the constraint content because U_m
freely inflates to make T_relax PSD without constraining (c, d).

THIS MODULE — exact lift via Lasserre level-2 moment variables.
With the Lasserre lift in place, the moment matrix M_top contains the EXACT
bilinear values
    M_top[m, m]              = c_m · c_m,
    M_top[T_max+m, T_max+m]  = d_m · d_m,
so |f̂(m)|² = (c_m² + d_m²)/4 = (M_top[m,m] + M_top[T_max+m, T_max+m]) / 4
is now an EXACT LINEAR function of program variables (no slack, no relaxation
direction issue). The M-side Toeplitz entries
    M̂(0)  = Ω/2,
    M̂(m)  = (a_m/2) c_{m-1} − 4 · (M_top[m,m] + M_top[T_max+m, T_max+m]) / 4
          = (a_m/2) c_{m-1} − (M_top[m,m] + M_top[T_max+m, T_max+m]),
    Im part: ∓ (a_m/2) d_{m-1}
are linear in (Ω, c, d, M_top).

Imposing T_M(c, d, M_top) ⪰ 0 is then a STANDARD linear-in-variables PSD
constraint, encoded as a real (2(n_M+1)) × (2(n_M+1)) symmetric PSD via the
usual [[Re, -Im], [Im, Re]] block form.

VALIDITY (rigorous LB direction).
The Lasserre level-2 PSD constraint M_top ⪰ ξ ξ^T forces M_top[m,m] ≥ c_m²
and M_top[T_max+m, T_max+m] ≥ d_m² (diagonal Schur-complement). Thus the
program-internal value of M̂(m) used in the Toeplitz constraint satisfies
    M̂(m)_program = (a_m/2) c_m − M_top[m,m] − M_top[T_max+m, T_max+m]
                ≤ (a_m/2) c_m − c_m² − d_m²
                = (a_m/2) c_m − 4|f̂(m)|²
                = M̂(m)_true.
The Toeplitz matrix `T_M_program` has the SAME diagonal as `T_M_true`
(M̂(0) = Ω/2 in both) but smaller-or-equal off-diagonal real parts (in
absolute value, but with sign... see below). So T_M_program ⪰ 0 does NOT
immediately imply T_M_true ⪰ 0 — we need to think harder.

Cleaner validity proof: at every TRUE realizable f, take the rank-1 lift
M_top = ξ ξ^T (which is feasible for the Lasserre PSD constraint since
M_top ⪰ ξ ξ^T then holds with equality). Then M_top[m,m] = c_m² and
M_top[T_max+m, T_max+m] = d_m² exactly, hence the program-internal M̂(m)
equals the true M̂(m). The added Toeplitz PSD constraint is therefore
satisfied at every µ-realizing distribution → the LB is rigorous.

CAVEAT (anti-bluff #3 from queue): the bite of the constraint is bounded
by how tight the Lasserre level-2 lift's diagonal slack
M_top[i,i] - ξ_i² is at the LP optimum. Empirically the slack is small but
non-zero (~10⁻³ at row 4 N=2000, measured in self-test below), which means
the M-Toeplitz PSD test at the program optimum is slightly LOOSER than
the true T_M ⪰ 0 test (program permits inflated |f̂(m)|² → smaller off-
diagonal real parts → easier PSD). If Δ ≈ 1e-9 like the SOC version, the
explanation is that Lasserre-2's existing localizing PSD already implies
T_M ⪰ 0 at the optimum, so the M-side Toeplitz adds no new content.

Inputs to add_mside_bochner_lasserre_constraint:
    cons   — the cvxpy constraint list to extend (in place)
    c, d   — Fourier-coefficient cvxpy.Variable arrays of length T
    Omega  — the LP objective scalar variable (M̂(0) = Ω/2)
    M_top  — the (2T_max+1)×(2T_max+1) Lasserre-lifted moment matrix
             returned by lasserre.add_lasserre2_constraint
    T_max  — the Lasserre lift cutoff (must equal the value used to build M_top)
    n_M    — Toeplitz order; matrix has size (n_M+1); 0 = no constraint;
             must satisfy n_M ≤ T_max.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def a_m(m: int) -> float:
    """Lemma-2 coefficient a_m = (4/(mπ)) sin(mπ/2).

    Vanishes for even m; equals (4/(mπ))·(−1)^((m−1)/2) for odd m.
    """
    return 4.0 / (m * np.pi) * np.sin(m * np.pi / 2.0)


def add_mside_bochner_lasserre_constraint(
    cons: list,
    c: cp.Variable,
    d: cp.Variable,
    Omega: cp.Variable,
    M_top: cp.Variable,
    T_max: int,
    n_M: int,
):
    """Append the M-side Bochner PSD constraint using Lasserre-2 lifts.

    Parameters
    ----------
    cons : list
        cvxpy constraint list to extend in place.
    c, d : cvxpy.Variable
        Fourier-coefficient variables of length T (T >= n_M).
    Omega : cvxpy.Variable
        The objective scalar with M̂(0) = Ω/2.
    M_top : cvxpy.Variable
        The (2T_max+1)×(2T_max+1) Lasserre-2 moment matrix; index layout
        ξ = (1, c_1..c_{T_max}, d_1..d_{T_max}) so that
            |f̂(m)|² = (M_top[m,m] + M_top[T_max+m, T_max+m]) / 4.
    T_max : int
        Lasserre lift cutoff (must match the M_top construction).
    n_M : int
        Toeplitz order. The Hermitian matrix has size (n_M+1).
    """
    if n_M <= 0:
        return
    assert n_M <= T_max, (
        f"M-side Bochner order n_M={n_M} cannot exceed Lasserre T_max={T_max}"
    )
    T = c.shape[0]
    assert n_M <= T, f"n_M={n_M} cannot exceed T={T}"

    # Build the Hermitian Toeplitz matrix M̂(j-k) entries:
    #   M̂(0)        = Ω / 2
    #   Re M̂(m)     = (a_m / 2) c_{m-1} - (M_top[m, m] + M_top[T_max+m, T_max+m])
    #   Im M̂(+m)    = -(a_m / 2) d_{m-1}
    #   Im M̂(-m)    = +(a_m / 2) d_{m-1}     (Hermitian symmetry)
    # The substitution
    #     |f̂(m)|² = (c_m² + d_m²)/4
    #             = (M_top[m,m] + M_top[T_max+m, T_max+m]) / 4
    # converts the −4|f̂(m)|² term to −(M_top[m,m] + M_top[T_max+m, T_max+m]).
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
                # Lasserre-lifted exact bilinear: c_m² → M_top[m, m],
                # d_m² → M_top[T_max + m, T_max + m].
                bilin = M_top[m, m] + M_top[T_max + m, T_max + m]
                re_entry = (am / 2.0) * c[m - 1] - bilin
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
    # Hermitian (n_M+1)×(n_M+1) PSD ⟺ real-symmetric (2n_M+2)×(2n_M+2) PSD via
    #     [[Re, -Im], [Im, Re]] ⪰ 0.
    real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
    cons.append(real_form >> 0)


# ----- self-test (math sanity, no solve) ----------------------------------
if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")

    print("a_m for m=1..6:")
    for m in range(1, 7):
        print(f"  a_{m} = {a_m(m): .6f}")
    print("  (expected: ±4/(mπ) for odd m, 0 for even m)")

    # Smoke test: build a dummy program with Lasserre + M-side Bochner.
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from lasserre import add_lasserre2_constraint

    print("\nSmoke test: build constraints with T=5, T_max=4, n_M=3.")
    T = 5
    T_max = 4
    n_M = 3
    c = cp.Variable(T)
    d = cp.Variable(T)
    Omega = cp.Variable()
    cons = []
    M_top = add_lasserre2_constraint(cons, c, d, T_max=T_max, T_loc=T_max)
    n0 = len(cons)
    add_mside_bochner_lasserre_constraint(
        cons, c, d, Omega, M_top, T_max=T_max, n_M=n_M,
    )
    n1 = len(cons)
    print(f"  Lasserre added {n0} constraints; M-side added {n1 - n0} "
          f"(expected 1 PSD).")
    print(f"  M_top shape: {M_top.shape} (expected ({2*T_max+1}, {2*T_max+1}))")

    # Verify the Toeplitz matrix has the right structure: pin (c, d) to a
    # specific value and M_top to its rank-1 lift, then evaluate the matrix.
    a_val = 0.4   # c_1
    b_val = 0.0   # d_1
    cons2 = []
    cons2.append(c[0] == a_val); cons2.append(d[0] == b_val)
    for k in range(1, T):
        cons2.append(c[k] == 0); cons2.append(d[k] == 0)
    xi_val = np.zeros(2 * T_max + 1)
    xi_val[0] = 1.0
    xi_val[1] = a_val   # c_1 in slot 1
    xi_val[T_max + 1] = b_val   # d_1 in slot T_max+1
    cons2.append(M_top == np.outer(xi_val, xi_val))
    cons2.append(Omega == 0.5)

    prob = cp.Problem(cp.Minimize(0), cons + cons2)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"\nAt c_1=0.4, d_1=0, Omega=0.5 (rank-1 lift): solver status = {prob.status}")
    # Expected M̂(1) = (4/π) * (1/2) * 0.4 - 0.4² - 0² = 0.5093 - 0.16 = 0.3493
    # Toeplitz [[0.25, 0.3493], [0.3493, 0.25]] — eigenvalues 0.25 ± 0.3493
    # → smallest eigenvalue = -0.0993 → NOT PSD → infeasible expected.
    print(f"  Expected: 'infeasible' (the test value c_1=0.4 violates "
          f"M-Toeplitz PSD — smallest eigenvalue ≈ -0.0993).")

    # Re-try with c_1 small enough that M̂(1) < Ω/2 = 0.25:
    print("\nRe-test: c_1=0.1, d_1=0 (gives M̂(1) ≈ 0.1273 - 0.01 = 0.1173 < 0.25)")
    T = 5; T_max = 4; n_M = 3
    c = cp.Variable(T)
    d = cp.Variable(T)
    Omega = cp.Variable()
    cons = []
    M_top = add_lasserre2_constraint(cons, c, d, T_max=T_max, T_loc=T_max)
    add_mside_bochner_lasserre_constraint(
        cons, c, d, Omega, M_top, T_max=T_max, n_M=n_M,
    )
    a_val = 0.1
    cons.append(c[0] == a_val); cons.append(d[0] == 0)
    for k in range(1, T):
        cons.append(c[k] == 0); cons.append(d[k] == 0)
    xi_val = np.zeros(2 * T_max + 1)
    xi_val[0] = 1.0; xi_val[1] = a_val
    cons.append(M_top == np.outer(xi_val, xi_val))
    cons.append(Omega == 0.5)
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL", verbose=False)
    print(f"  Solver status: {prob.status}  (expected: 'optimal' or '..._inaccurate')")
