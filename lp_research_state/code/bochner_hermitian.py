"""
COMPLEX-HERMITIAN Bochner-PSD constraints (Approach ③ — representation-theoretic
symmetry reduction, the "free ~4×" U(1) factor).

This is a NON-DESTRUCTIVE variant of bochner.py / white_full_convex.py's Bochner
section. The verified real-embedding machinery is left untouched; this file
implements the *mathematically identical* constraint with HALF the PSD matrix
dimension.

------------------------------------------------------------------------------
Background (real embedding, the current verified encoding)
------------------------------------------------------------------------------
Bochner: real f ≥ 0 a.e. (extended 2-periodically) ⟺ for every n ≥ 0 the
(n+1)×(n+1) Hermitian Toeplitz moment matrix
    M_n[j,k] = f̂(j−k),   0 ≤ j,k ≤ n
is PSD. With White's convention
    f̂(0)   = 1/2,
    f̂(k)   = (c_k − i d_k)/2     (k ≥ 1),
    f̂(−k)  = (c_k + i d_k)/2 = conj(f̂(k)).
For 1−f ≥ 0:  (1−f)̂(0)=1/2, (1−f)̂(k) = −f̂(k) ⇒ off-diagonals negated.

The verified code encodes the complex Hermitian PSD constraint M_n ⪰ 0 as a
2(n+1)×2(n+1) REAL symmetric PSD constraint via the standard embedding
    real_form = [[Re M, −Im M], [Im M, Re M]]  ⪰ 0   ⟺   M ⪰ 0.

------------------------------------------------------------------------------
This file: direct complex-Hermitian PSD
------------------------------------------------------------------------------
cvxpy ≥ 1.8 supports complex PSD natively. We declare a complex Hermitian
variable H of size (n+1)×(n+1), impose H ⪰ 0, and tie its entries to the
Fourier coefficient variables so that H = M_n(f) (resp. M_n(1−f)).

WHY THIS IS EXACTLY EQUIVALENT (and not a relaxation):
  The map  M ↦ [[Re M, −Im M], [Im M, Re M]]  is an injective *-algebra
  homomorphism C^{(n+1)×(n+1)} → R^{2(n+1)×2(n+1)} whose image is exactly the
  matrices commuting with J=[[0,−I],[I,0]]. It preserves the PSD cone
  bijectively: M ⪰ 0 (complex) ⟺ real_form ⪰ 0 (real). Every real eigenvalue
  of M appears with multiplicity 2 in real_form (the U(1)/SO(2) action), so the
  real form carries no extra/fewer constraints — only twice the dimension. The
  complex form is the U(1)-isotypic reduction: the documented free ~4× on PSD
  flops/memory (2·(½)³ = ¼), obtained with zero representation theory.

Ties imposed on the complex variable H (size (n+1)×(n+1)):
  * H Hermitian:                 enforced by hermitian=True (so only the upper
                                 triangle + real diagonal are free).
  * Toeplitz + coefficient tie:  H[j,k] = f̂(j−k) for the chosen sign.
      diagonal           H[j,j] = 1/2
      super-diagonal ℓ=k−j>0:    H[j,k] = f̂(−ℓ) = (c_ℓ + i d_ℓ)/2 · sign
      sub-diagonal   ℓ=j−k>0:    H[j,k] = f̂(ℓ)  = (c_ℓ − i d_ℓ)/2 · sign
    (the sub-diagonal tie is automatically the conjugate of the super-diagonal
     one because H is Hermitian; we tie the super-diagonal explicitly and let
     Hermiticity propagate, which is the cheaper/cleaner encoding.)

CONVENTION-TRAP CHECKLIST (the classic svec/√2/sign/conjugation traps):
  * cvxpy's complex `>> 0` canonicalizes Hermitian H via the SAME real embedding
    [[Re,−Im],[Im,Re]] internally — so the dual it returns lives on the identical
    cone; rigorous_dual_LB is directly comparable. (Verified empirically by the
    10-digit cross-check driver _herm_equiv_check.py.)
  * Sign of the imaginary part: f̂(j−k) with j−k>0 (sub-diagonal) is c−id, i.e.
    NEGATIVE imaginary; with j−k<0 (super-diagonal) it is c+id. We tie the
    SUPER-diagonal (k>j) to (c+id)/2 — matching Im M[j,k] = +sign·d/2 there,
    exactly the real-form Im_M of bochner.py.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def add_bochner_hermitian_constraint(cons: list, c: cp.Variable, d: cp.Variable,
                                     n: int, sign: int = +1):
    """Append the (n+1)×(n+1) COMPLEX Hermitian PSD constraint for f or 1−f.

    Mathematically identical to bochner.add_bochner_constraint / the
    white_full_convex.py real-form block, but the PSD cone is half the size.

    Parameters
    ----------
    cons : list           cvxpy constraint list (appended in place)
    c, d : cp.Variable    real Fourier-coefficient variables, length T
    n    : int            Bochner order (matrix size n+1)
    sign : int            +1 for f ≥ 0,  −1 for 1−f ≥ 0
    """
    T = c.shape[0]
    assert n <= T, f"Bochner order n={n} cannot exceed T={T}"

    H = cp.Variable((n + 1, n + 1), hermitian=True)

    ties = []
    # Diagonal: H[j,j] = 1/2 (real).
    for j in range(n + 1):
        ties.append(H[j, j] == 0.5)
    # Upper triangle (k > j): ℓ = k − j > 0, entry = f̂(−ℓ) = (c_ℓ + i d_ℓ)/2,
    # scaled by `sign` for the 1−f case. Hermiticity fixes the lower triangle as
    # the conjugate automatically, so we do NOT tie it again (would be redundant
    # but harmless; we keep the encoding minimal).
    for j in range(n + 1):
        for k in range(j + 1, n + 1):
            ell = k - j  # > 0
            re = cp.Constant(sign * 0.5) * c[ell - 1]
            im = cp.Constant(sign * 0.5) * d[ell - 1]
            # cvxpy complex literal: re + 1j*im
            ties.append(H[j, k] == re + 1j * im)

    cons.extend(ties)
    cons.append(H >> 0)
    return H


# ---------------------------------------------------------------------------
# Numeric (numpy) helper for cross-checking against bochner_independent.py.
# ---------------------------------------------------------------------------
def make_hermitian_matrix(c_vals, d_vals, n: int, sign: int = +1):
    """Build the complex Hermitian (n+1)×(n+1) M_n given numeric c, d.

    Returns the complex matrix H. By construction
        H = (Re_form_top_left) + i (Im_form_bottom_left)
    of bochner_independent.make_real_form, i.e. H is exactly A + iB there.
    """
    c_vals = np.asarray(c_vals, dtype=float)
    d_vals = np.asarray(d_vals, dtype=float)
    H = np.zeros((n + 1, n + 1), dtype=complex)
    for j in range(n + 1):
        for k in range(n + 1):
            ell = j - k
            if ell == 0:
                H[j, k] = 0.5
            else:
                aell = abs(ell)
                # f̂(ℓ) = (c−id)/2 for ℓ>0 (sub-diag); (c+id)/2 for ℓ<0 (super-diag)
                im_sign = -1.0 if ell > 0 else +1.0
                H[j, k] = sign * (c_vals[aell - 1] + 1j * im_sign * d_vals[aell - 1]) / 2.0
    return H


# ----- self-test ----------------------------------------------------------
if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    import importlib.util as ilu
    from pathlib import Path

    # 1) Numeric equivalence to the real embedding (bochner_independent.make_real_form).
    here = Path(__file__).resolve().parent
    spec = ilu.spec_from_file_location("bochner_independent", here / "bochner_independent.py")
    bi = ilu.module_from_spec(spec); spec.loader.exec_module(bi)

    rng = np.random.default_rng(0)
    max_err = 0.0
    for trial in range(20):
        n = int(rng.integers(2, 8))
        cv = rng.uniform(-0.3, 0.3, size=n)
        dv = rng.uniform(-0.3, 0.3, size=n)
        for sign in (+1, -1):
            H = make_hermitian_matrix(cv, dv, n, sign)
            real_form, A, B = bi.make_real_form(cv, dv, n, sign)
            # Real embedding of H should equal real_form bit-for-bit.
            emb = np.block([[H.real, -H.imag], [H.imag, H.real]])
            max_err = max(max_err, np.max(np.abs(emb - real_form)))
            # Eigenvalues: real_form eigs = H eigs each doubled.
            eH = np.sort(np.linalg.eigvalsh(H))
            eR = np.sort(np.linalg.eigvalsh(real_form))
            eH2 = np.sort(np.concatenate([eH, eH]))
            max_err = max(max_err, np.max(np.abs(eH2 - eR)))
    print(f"Numeric Hermitian↔real-embedding max abs error over 20 trials: {max_err:.2e}")
    print("  (expect < 1e-13: the complex form IS the real form, halved.)")

    # 2) cvxpy smoke test: tiny feasibility solve agrees on a known density.
    import cvxpy as cp
    T = 4; n = 2
    c = cp.Variable(T); d = cp.Variable(T)
    cons = [c[0] == 0.4, d == 0]  # f = 1/2 + 0.4 cos(πx), nonneg
    add_bochner_hermitian_constraint(cons, c, d, n, +1)
    add_bochner_hermitian_constraint(cons, c, d, n, -1)
    prob = cp.Problem(cp.Minimize(0), cons + [c[1:] == 0])
    prob.solve(solver="CLARABEL")
    print(f"cvxpy feasibility (f=1/2+0.4cos): status={prob.status} (expect optimal)")
