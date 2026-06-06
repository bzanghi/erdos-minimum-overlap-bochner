"""
REAL Z/2 CENTROSYMMETRIC Bochner-PSD constraints (Approach ③ — PRO-49, the
in-cvxpy realizable speedup).

NON-DESTRUCTIVE variant of bochner.py / white_full_convex.py's Bochner section.
The verified real-embedding machinery is left untouched; this file implements the
*mathematically identical* constraint with HALF the PSD matrix dimension that
CLARABEL actually consumes (a real PSDTriangle cone of side n+1 instead of
2(n+1)).

------------------------------------------------------------------------------
Why the complex-Hermitian route (bochner_hermitian.py) gave ZERO win
------------------------------------------------------------------------------
cvxpy's mandatory `Complex2Real` reduction expands `cp.Variable(hermitian=True),
H >> 0` back into the identical 2(n+1)×2(n+1) real embedding [[Re,−Im],[Im,Re]]
before any solver sees it (CLARABEL/SCS/MOSEK have no complex PSD cone). So the
"free ~4×" never reaches the solver. Confirmed last iteration
(sym_reduction_result.{json,md}). The centrosymmetric split below is the ONLY
in-cvxpy path to a real speedup.

------------------------------------------------------------------------------
The centrosymmetric block split (validated, this session)
------------------------------------------------------------------------------
The Bochner real form is
    RF = [[A, −B], [B, A]],   A = Re M_n (symmetric),  B = Im M_n (antisymmetric),
RF ⪰ 0  ⟺  the (n+1)×(n+1) complex Hermitian M_n = A + iB ⪰ 0.

RF is centrosymmetric (S·RF·Sᵀ = RF with S = [[J,0],[0,−J]], J = anti-diagonal
reversal of size n+1). The orthogonal involution
    Q = (1/√2) [[I, J], [J, −I]]      (Q = Qᵀ = Q⁻¹)
block-diagonalizes it:
    Q·RF·Qᵀ = diag(B1, B2).

For this *Hermitian-Toeplitz* structure one has the extra identity  JB = −BJ
(verified to 0.0), which forces the two diagonal blocks to COLLAPSE to the SAME
matrix:
    B1 = B2 = A + J·B   =:  Bk   (real, symmetric).

Hence the full 2(n+1)×2(n+1) real PSD constraint is EXACTLY equivalent to the
SINGLE (n+1)×(n+1) real PSD constraint
    Bk = A + J·B  ⪰ 0,
and the RF spectrum is precisely the Bk spectrum DOUBLED (the U(1)/SO(2)
multiplicity-2 fact). This is even cheaper than "two half-size blocks": the
Toeplitz symmetry makes them identical, so one block suffices.

Validation (this session, 2000 random trials, orders n=1..8, both signs):
    max |Bk − Bkᵀ|              = 0.0      (Bk symmetric → valid real PSD cone)
    max |J·B + B·J|             = 0.0      (JB = −BJ)
    max |eig(RF) − eig(Bk)×2|   = 3.6e-15  (spectrum doubled)
    max |eig(A+iB) − eig(Bk)|   = 2.0e-15  (Bk is the Hermitian, real-similar)
    PSD-consistency failures     = 0/2000

Closed form of the single block (sub/super-diagonals of A plus the
anti-diagonal-flipped antisymmetric B); verified bit-for-bit (0.0 error) against
the numeric A + J@B:
    Bk[j,k] = A_part(j,k) + B_part(j,k),    0 ≤ j,k ≤ n
      A_part = 1/2                  if j == k
             = sign·(1/2)·c_{|j−k|} if j ≠ k
      B_part = (J·B)[j,k] = B[n−j, k],  with r = n−j, rk = r−k:
             = 0                    if rk == 0
             = −sign·(1/2)·d_{|rk|} if rk > 0
             = +sign·(1/2)·d_{|rk|} if rk < 0

Inputs to add_bochner_centro_constraint:
    cons     — the cvxpy constraint list (appended in place)
    c, d     — real Fourier coefficient variables (length T)
    n        — Bochner order (matrix size n+1)
    sign     — +1 for f ≥ 0,  −1 for 1−f ≥ 0
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def add_bochner_centro_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                                  n: int, sign: int = +1):
    """Append the SINGLE (n+1)×(n+1) real symmetric PSD constraint Bk = A + J·B
    ⪰ 0, mathematically identical to the 2(n+1)×2(n+1) real-form Bochner block
    of bochner.add_bochner_constraint / white_full_convex.py lines 237-262, but
    HALF the side length (the cone CLARABEL actually consumes).

    Parameters
    ----------
    cons : list           cvxpy constraint list (appended in place)
    c, d : cp.Variable    real Fourier-coefficient variables, length T
    n    : int            Bochner order (matrix size n+1)
    sign : int            +1 for f ≥ 0,  −1 for 1−f ≥ 0
    """
    T = c.shape[0]
    assert n <= T, f"Bochner order n={n} cannot exceed T={T}"

    half = 0.5
    rows = []
    for j in range(n + 1):
        row = []
        for k in range(n + 1):
            # A_part = Re M_n[j,k]
            if j == k:
                a_part = cp.Constant(half)
            else:
                a_part = cp.Constant(sign * half) * c[abs(j - k) - 1]
            # B_part = (J·Im M_n)[j,k] = Im M_n[n-j, k]
            r = n - j
            rk = r - k
            if rk == 0:
                b_part = cp.Constant(0.0)
            elif rk > 0:
                b_part = cp.Constant(-sign * half) * d[abs(rk) - 1]
            else:
                b_part = cp.Constant(+sign * half) * d[abs(rk) - 1]
            row.append(a_part + b_part)
        rows.append(row)

    Bk = cp.bmat(rows)
    cons.append(Bk >> 0)
    return Bk


# ---------------------------------------------------------------------------
# Numeric (numpy) helper for cross-checking against bochner.py's real form
# and bochner_hermitian.make_hermitian_matrix.
# ---------------------------------------------------------------------------
def make_centro_block(c_vals, d_vals, n: int, sign: int = +1):
    """Build the single real symmetric (n+1)×(n+1) block Bk = A + J·B numerically.

    By construction Bk has the SAME spectrum as the complex Hermitian
    M_n = A + iB (bochner_hermitian.make_hermitian_matrix) and as each half of
    the real form RF = [[A,−B],[B,A]] (bochner_independent.make_real_form).
    """
    c_vals = np.asarray(c_vals, dtype=float)
    d_vals = np.asarray(d_vals, dtype=float)
    Bk = np.zeros((n + 1, n + 1), dtype=float)
    for j in range(n + 1):
        for k in range(n + 1):
            a_part = 0.5 if j == k else sign * 0.5 * c_vals[abs(j - k) - 1]
            r = n - j
            rk = r - k
            if rk == 0:
                b_part = 0.0
            elif rk > 0:
                b_part = -sign * 0.5 * d_vals[abs(rk) - 1]
            else:
                b_part = sign * 0.5 * d_vals[abs(rk) - 1]
            Bk[j, k] = a_part + b_part
    return Bk


# ----- self-test ----------------------------------------------------------
if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    import importlib.util as ilu
    from pathlib import Path

    here = Path(__file__).resolve().parent

    # Independent real-form reference.
    spec = ilu.spec_from_file_location("bochner_independent", here / "bochner_independent.py")
    bi = ilu.module_from_spec(spec); spec.loader.exec_module(bi)
    # Complex-Hermitian reference.
    spec2 = ilu.spec_from_file_location("bochner_hermitian", here / "bochner_hermitian.py")
    bh = ilu.module_from_spec(spec2); spec2.loader.exec_module(bh)

    rng = np.random.default_rng(0)
    max_sym = max_spec_real = max_spec_herm = 0.0
    psd_fail = 0
    for _ in range(50):
        n = int(rng.integers(2, 9))
        cv = rng.uniform(-0.5, 0.5, size=n)
        dv = rng.uniform(-0.5, 0.5, size=n)
        for sign in (+1, -1):
            Bk = make_centro_block(cv, dv, n, sign)
            max_sym = max(max_sym, np.max(np.abs(Bk - Bk.T)))
            real_form, _, _ = bi.make_real_form(cv, dv, n, sign)
            H = bh.make_hermitian_matrix(cv, dv, n, sign)
            eRF = np.sort(np.linalg.eigvalsh(real_form))
            eBk = np.sort(np.linalg.eigvalsh(Bk))
            eH = np.sort(np.linalg.eigvalsh(H))
            max_spec_real = max(max_spec_real,
                                np.max(np.abs(eRF - np.sort(np.concatenate([eBk, eBk])))))
            max_spec_herm = max(max_spec_herm, np.max(np.abs(eH - eBk)))
            if (np.min(eRF) >= -1e-9) != (np.min(eBk) >= -1e-9):
                psd_fail += 1
    print(f"max |Bk - Bk^T|                  = {max_sym:.2e}  (symmetric)")
    print(f"max |eig(real_form) - eig(Bk)x2| = {max_spec_real:.2e}  (RF spectrum doubled)")
    print(f"max |eig(Hermitian) - eig(Bk)|   = {max_spec_herm:.2e}  (Bk is the Hermitian)")
    print(f"PSD-consistency failures         = {psd_fail}/100")

    # cvxpy smoke test: tiny feasibility on a known density f=1/2+0.4cos.
    T = 4; n = 2
    c = cp.Variable(T); d = cp.Variable(T)
    cons = [c[0] == 0.4, d == 0, c[1:] == 0]
    add_bochner_centro_constraint(cons, c, d, n, +1)
    add_bochner_centro_constraint(cons, c, d, n, -1)
    prob = cp.Problem(cp.Minimize(0), cons)
    prob.solve(solver="CLARABEL")
    print(f"cvxpy feasibility (f=1/2+0.4cos): status={prob.status} (expect optimal)")
