"""
NON-DESTRUCTIVE complex-Hermitian variant of white_full_convex.build_problem.

Approach ③ (representation-theoretic symmetry reduction). The ONLY difference
from the verified program is the Bochner f≥0 / 1−f≥0 block: instead of the
2(n+1)×2(n+1) real embedding [[Re,−Im],[Im,Re]], we use a direct (n+1)×(n+1)
complex Hermitian PSD constraint (the U(1)-isotypic reduction → free ~4× on the
PSD blocks). Everything else (cell envelope, M-side SOC, tail bounds,
normalization, objective) is reused verbatim by delegating to the verified
build_problem with bochner_n=0 and bolting the Hermitian block onto the returned
constraint list.

This guarantees the rest of the SDP is byte-for-byte the verified program; the
equivalence cross-check (_herm_equiv_check.py) then only has to validate the
Bochner swap itself.
"""
from __future__ import annotations
import cvxpy as cp

from white_full_convex import build_problem as _build_problem_real
from bochner_hermitian import add_bochner_hermitian_constraint


def build_problem_hermitian(
    N, T, R, h1, h2, p1, p2, q1, q2,
    cell_mode="exact",
    use_T3=False, use_T5=False, use_T5p=False,
    mside_sin_coeff=4.0,
    bochner_n=0,
    mside_bochner_n=0,
    mside_bochner_schur_n=0,
    assume_even=False,
):
    """Same signature/semantics as build_problem, but the f-side Bochner block
    is encoded as a COMPLEX Hermitian PSD instead of the real embedding.

    Returns (Omega, w, v, c, d, eps, dlt, cons) exactly like build_problem, plus
    the Hermitian-block list is appended to cons.
    """
    # Delegate to the verified program with the real-form Bochner DISABLED.
    Omega, w, v, c, d, eps, dlt, cons = _build_problem_real(
        N, T, R, h1, h2, p1, p2, q1, q2, cell_mode,
        use_T3=use_T3, use_T5=use_T5, use_T5p=use_T5p,
        mside_sin_coeff=mside_sin_coeff,
        bochner_n=0,                      # <-- key: no real-form Bochner
        mside_bochner_n=mside_bochner_n,
        mside_bochner_schur_n=mside_bochner_schur_n,
        assume_even=assume_even,
    )

    # Bolt on the complex-Hermitian f-side Bochner block (both f≥0 and 1−f≥0),
    # mirroring white_full_convex.py lines 237-262 exactly but in complex form.
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in (+1, -1):
            add_bochner_hermitian_constraint(cons, c, d, n_b, sign)

    return Omega, w, v, c, d, eps, dlt, cons
