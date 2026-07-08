"""
Reproduce + extend White (Acta Arith. 2023, arXiv:2201.05704) Section 5.

Two sets of cell bounds α^±, β^± available:
  mode='lipschitz'  — White's Lipschitz envelope (α± = cos(midpoint) ± πmL/4)
  mode='exact'      — exact min/max of cos/sin on each cell (tighter)

Other knobs:
  N, T, R           — paper's parameters
  parameter ranges  — (h1, h2), (p1, p2), (q1, q2)

Returns the convex problem and primal solution.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp
import time
from typing import Optional, Literal


# ----- cell bounds: tight versions ------------------------------------------------
def cos_cell_bounds_exact(j_arr, m, L):
    """min, max of cos(πmx/2) on each cell [(j-1)L, jL]."""
    x_lo = (j_arr - 1) * L
    x_hi = j_arr * L
    cl = np.cos(np.pi * m * x_lo / 2)
    cr = np.cos(np.pi * m * x_hi / 2)
    lo = np.minimum(cl, cr)
    hi = np.maximum(cl, cr)
    # Critical points of cos(πmx/2) are at x = 2n/m  (cos(nπ) = (-1)^n)
    # For each cell, find integer n with x_lo < 2n/m < x_hi.
    # ⇔ m*x_lo/2 < n < m*x_hi/2
    n_lo = np.ceil(m * x_lo / 2 + 1e-15).astype(int)
    n_hi = np.floor(m * x_hi / 2 - 1e-15).astype(int)
    # Each cell gets at most a few critical points; iterate.
    for k in range(int(np.min(n_lo)), int(np.max(n_hi)) + 1):
        mask = (n_lo <= k) & (k <= n_hi)
        if not mask.any(): continue
        v = (-1) ** k
        if v == -1:
            lo = np.where(mask, np.minimum(lo, -1.0), lo)
        else:
            hi = np.where(mask, np.maximum(hi, 1.0), hi)
    return lo, hi


def sin_cell_bounds_exact(j_arr, m, L):
    """min, max of sin(πmx/2) on each cell."""
    x_lo = (j_arr - 1) * L
    x_hi = j_arr * L
    sl = np.sin(np.pi * m * x_lo / 2)
    sr = np.sin(np.pi * m * x_hi / 2)
    lo = np.minimum(sl, sr)
    hi = np.maximum(sl, sr)
    # Critical points of sin(πmx/2) are at x = (2n+1)/m  (sin = (-1)^n)
    # ⇔ m*x_lo/2 < (2n+1)π/2 / (π/1)... let me redo: sin(πmx/2) has critical points where
    # cos(πmx/2) = 0, i.e., πmx/2 = π/2 + nπ, i.e., x = (2n+1)/m.
    # So we need n with x_lo < (2n+1)/m < x_hi  ⇔  (m*x_lo - 1)/2 < n < (m*x_hi - 1)/2.
    n_lo = np.ceil((m * x_lo - 1) / 2 + 1e-15).astype(int)
    n_hi = np.floor((m * x_hi - 1) / 2 - 1e-15).astype(int)
    for k in range(int(np.min(n_lo)), int(np.max(n_hi)) + 1):
        mask = (n_lo <= k) & (k <= n_hi)
        if not mask.any(): continue
        v = (-1) ** k
        if v == -1:
            lo = np.where(mask, np.minimum(lo, -1.0), lo)
        else:
            hi = np.where(mask, np.maximum(hi, 1.0), hi)
    return lo, hi


def cos_cell_bounds_lipschitz(j_arr, m, L):
    x_mid = (j_arr - 0.5) * L
    cm = np.cos(np.pi * m * x_mid / 2)
    return cm - np.pi * m * L / 4, cm + np.pi * m * L / 4


def sin_cell_bounds_lipschitz(j_arr, m, L):
    x_mid = (j_arr - 0.5) * L
    sm = np.sin(np.pi * m * x_mid / 2)
    return sm - np.pi * m * L / 4, sm + np.pi * m * L / 4


def odd_coeff_factors(m: int, T: int):
    k = np.arange(1, T + 1)
    denom = m * m - 4 * k * k
    sgn = (-1) ** k
    return sgn / denom, k * sgn / denom


def tail_bound_eps(m_odd: int, T: int) -> float:
    return (1.0 / (4 - m_odd ** 2 / T ** 2)) * (2 * m_odd / (np.pi * np.sqrt(6 * T ** 3)))


def tail_bound_delta(m_odd: int, T: int) -> float:
    return (1.0 / (4 - m_odd ** 2 / T ** 2)) * (4 / (np.pi * np.sqrt(2 * T)))


def build_problem(
    N: int, T: int, R: int,
    h1: float, h2: float,
    p1: float, p2: float,
    q1: float, q2: float,
    cell_mode: Literal["lipschitz", "exact"] = "exact",
    use_T3: bool = False,        # Tightening 3: L Σ (w² + v²) ≤ Ω    (NEW)
    use_T5: bool = False,        # Tightening 5: φ = 1+cos(πx)
    use_T5p: bool = False,       # Tightening 5': φ = 1-cos(πx) — VIOLATED at White's optimum (NEW)
    mside_sin_coeff: float = 4.0,  # RHS coeff of imag/sine cell-consistency constraint 5.6/5.7.
                                 # White's 2026-05-31 email correction: this was 8, should be 4
                                 # (default now 4.0; pass mside_sin_coeff=8.0 to reproduce old behavior).
    bochner_n: int = 0,          # Bochner level (0 = off). Adds PSD constraints for f≥0 and 1-f≥0
    mside_bochner_n: int = 0,    # M-side Bochner (SOC-relaxed) level (0 = off). Adds PSD on T_relax(c,d,U)
    mside_bochner_schur_n: int = 0,  # M-side Bochner via EXACT Schur lifting (0 = off). Adds PSD on T_relax(c,d,s)
    assume_even: bool = False,   # NEW: if True, enforce d_k = 0, dlt = 0, and v_j = w_j (even f*)
    lasserre_T_max: int = 0,     # Lasserre level-2 cutoff (0 = off). Adds the moment-matrix lift
                                 # plus localizing PSD constraint for f^2 <= f against ALL nonneg
                                 # degree-2 trig polynomials (a strict tightening of Bochner-on-f).
    lasserre_T_loc: int = 0,     # Localizing-matrix order (0 = same as lasserre_T_max).
    mside_bochner_lasserre_n: int = 0,  # M-side Bochner using Lasserre-lifted bilinears (NEW).
                                        # Requires lasserre_T_max > 0. Replaces the SOC slack
                                        # |f̂(m)|² with the EXACT lifted bilinear M_top diagonal
                                        # entries — turns M-side Toeplitz PSD into a standard
                                        # linear-in-variables PSD constraint (no slack absorption).
):
    L = 2.0 / N
    j = np.arange(1, N + 1)

    if cell_mode == "lipschitz":
        cos_bnds = cos_cell_bounds_lipschitz
        sin_bnds = sin_cell_bounds_lipschitz
    else:
        cos_bnds = cos_cell_bounds_exact
        sin_bnds = sin_cell_bounds_exact

    Omega = cp.Variable()
    w = cp.Variable(N)
    v = cp.Variable(N)
    c = cp.Variable(T)
    d = cp.Variable(T)
    eps = cp.Variable(R)
    dlt = cp.Variable(R)

    cons = []
    cons += [w >= 0, v >= 0, w <= Omega, v <= Omega, Omega <= 1]
    cons.append(L * cp.sum(w + v) == 1)

    # ----- Even-f assumption (NEW): d_k = 0 ∀k, dlt = 0, and v_j = w_j ∀j.
    # If f* is even (f(x) = f(-x)) then sin Fourier coeffs vanish (d=0, dlt=0)
    # and the discretized cell-averages agree on positive and negative sides
    # (v = w). This halves the LP/SDP variable count. Resulting bound is
    # CONDITIONAL on the (open) even-f conjecture (cf. White 2023, §4 & §6).
    if assume_even:
        cons += [d == 0, dlt == 0, v == w]
    cons.append(L ** 2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h1)
    cons.append(L ** 3 * cp.sum(cp.multiply((j - 1) ** 2, (w + v))) <= 2.0 / 3 + h2 ** 2 / 2)

    a_expr = []
    b_expr = []
    for m in range(1, 2 * R + 1):
        if m % 2 == 0:
            half = m // 2
            a_expr.append(0.5 * c[half - 1])
            b_expr.append(0.5 * d[half - 1])
        else:
            af, bf = odd_coeff_factors(m, T)
            sin_pi_half_m = np.sin(np.pi * m / 2)
            am = (
                eps[(m - 1) // 2]
                + (2 * m * sin_pi_half_m / np.pi)
                * (1.0 / (2 * m ** 2) + cp.sum(cp.multiply(af, c)))
            )
            bm = (
                dlt[(m - 1) // 2]
                + (4 * sin_pi_half_m / np.pi) * cp.sum(cp.multiply(bf, d))
            )
            a_expr.append(am)
            b_expr.append(bm)

    for m in range(1, 2 * R + 1):
        am = a_expr[m - 1]; bm = b_expr[m - 1]
        a_minus, _ = cos_bnds(j, m, L)
        lhs = (L / 2) * (a_minus @ (w + v))
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
        cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)

    for m in range(1, 2 * R + 1):
        bm = b_expr[m - 1]
        b_minus, b_plus = sin_bnds(j, m, L)
        sin_pi_half_m = np.sin(np.pi * m / 2)
        # White's 2026-05-31 email correction: constraints 5.6/5.7 had an 8 in the RHS numerator, should be 4 (default now 4.0; pass mside_sin_coeff=8.0 to reproduce old behavior).
        rhs = -(mside_sin_coeff / (m * np.pi)) * sin_pi_half_m * bm
        cons.append((L / 2) * (b_minus @ w - b_plus @ v) <= rhs)
        cons.append((L / 2) * (b_plus @ w - b_minus @ v) >= rhs)

    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [
            cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
            cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T),
        ]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)
    cons += [c[0] >= p1, c[0] <= p2, d[0] >= q1, d[0] <= q2]

    _, a_plus_2 = cos_bnds(j, 2, L)
    rhs_513 = -0.5 * (max(p1 ** 2, p2 ** 2) + max(q1 ** 2, q2 ** 2))
    cons.append((L / 2) * (a_plus_2 @ (w + v)) >= rhs_513)

    # ----- Tightening 3 (NEW): ∫ M² ≤ Ω · ∫ M = Ω,
    # combined with cell-wise Cauchy-Schwarz L · wⱼ² ≤ ∫_cell M²
    # gives L · Σ(wⱼ² + vⱼ²) ≤ Ω.  Convex.
    if use_T3:
        cons.append(L * (cp.sum_squares(w) + cp.sum_squares(v)) <= Omega)

    # ----- Tightening 5 (NEW): from f² ≤ f tested against φ(x) = 1+cos(πx) ≥ 0
    # gives  Σ(c_k² + d_k²) + Σ(c_k c_{k+1} + d_k d_{k+1}) ≤ 1/2.
    # Convex form:  c.T Q c + d.T Q d ≤ 1/2  where Q is tridiagonal-toeplitz
    # with diagonal 1 and off-diagonal 1/2  (PSD, eigenvalues 1+cos(kπ/(T+1))).
    if use_T5:
        Q = np.eye(T) + 0.5 * np.eye(T, k=1) + 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Q)) + cp.quad_form(d, cp.psd_wrap(Q)) <= 0.5)

    # ----- Tightening 5' (NEW, biggest impact): from f² ≤ f tested against
    # φ(x) = 1 - cos(πx) ≥ 0.  Gives Σ(c² + d²) - Σ(c_k c_{k+1} + d_k d_{k+1}) ≤ 1/2.
    # Q' is PSD with eigenvalues 1 - cos(kπ/(T+1)) ∈ (0, 2].
    if use_T5p:
        Qp = np.eye(T) - 0.5 * np.eye(T, k=1) - 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Qp)) + cp.quad_form(d, cp.psd_wrap(Qp)) <= 0.5)

    # ----- Bochner moment-matrix PSD constraints (NEW)
    # Theorem (Bochner). f ≥ 0 ⟺ ∀n, the Hermitian Toeplitz moment matrix
    # M_n(f) := [f̂(j-k)]_{j,k=0..n} is PSD.  With f̂(0)=1/2, f̂(k)=(c_k - id_k)/2.
    # Encoded as 2(n+1)x2(n+1) real symmetric PSD via [[Re,-Im],[Im,Re]].
    # Both f ≥ 0 and 1-f ≥ 0 are encoded (sign flip on off-diagonal entries).
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            # Build Re_M and Im_M as (n+1)x(n+1) cvxpy expressions
            half = 0.5
            Re_rows, Im_rows = [], []
            for j in range(n_b + 1):
                re_row, im_row = [], []
                for k in range(n_b + 1):
                    ell = j - k
                    if ell == 0:
                        re_row.append(cp.Constant(half))
                        im_row.append(cp.Constant(0.0))
                    else:
                        aell = abs(ell)
                        re_row.append(cp.Constant(sign * 0.5) * c[aell - 1])
                        if ell > 0:
                            im_row.append(cp.Constant(-sign * 0.5) * d[aell - 1])
                        else:
                            im_row.append(cp.Constant(+sign * 0.5) * d[aell - 1])
                Re_rows.append(re_row)
                Im_rows.append(im_row)
            Re_M = cp.bmat(Re_rows)
            Im_M = cp.bmat(Im_rows)
            real_form = cp.bmat([[Re_M, -Im_M], [Im_M, Re_M]])
            cons.append(real_form >> 0)

    # ----- M-side Bochner moment-matrix PSD (SOC-relaxed) constraint (NEW)
    # Lemma 2: M̂(0) = Ω/2; M̂(m) = a_m f̂(m) − 4|f̂(m)|² for m ≥ 1, with
    # a_m = (4/(mπ)) sin(mπ/2). M ≥ 0 ⟺ ∀n: Hermitian Toeplitz T_M ⪰ 0.
    # Non-convex through the −4|f̂(m)|² term; we relax via SOC slack
    # U_m ≥ |f̂(m)|² and impose PSD on the resulting linear-in-(c,d,U)
    # T_relax. Validity (F_1 ⊆ F_2 ⊆ F_0 ⇒ valid LB) proved in findings.md.
    if mside_bochner_n > 0:
        # Lazy import: keep base build_problem usable even if mside_bochner.py
        # is missing (e.g. at very old code-mount paths).
        import importlib.util as _ilu
        from pathlib import Path as _Path
        _here = _Path(__file__).resolve().parent
        _spec = _ilu.spec_from_file_location("mside_bochner", _here / "mside_bochner.py")
        _mb = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_mb)
        _mb.add_mside_bochner_constraint(cons, c, d, Omega, mside_bochner_n)

    # ----- M-side Bochner via EXACT Schur lifting (NEW)
    # Same convex set as the SOC version, but written as a 3x3 Schur block
    # [[s, c/2, d/2], [c/2, 1, 0], [d/2, 0, 1]] ⪰ 0  ⟺  s ≥ |f̂|²
    # for each lag, plus the Hermitian-Toeplitz PSD on T_relax(c,d,s).
    if mside_bochner_schur_n > 0:
        import importlib.util as _ilu
        from pathlib import Path as _Path
        _here = _Path(__file__).resolve().parent
        _spec = _ilu.spec_from_file_location(
            "mside_bochner_schur", _here / "mside_bochner_schur.py"
        )
        _mbs = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_mbs)
        _mbs.add_mside_bochner_schur_constraint(
            cons, c, d, Omega, mside_bochner_schur_n
        )

    # ----- Lasserre level-2 augmentation (NEW)
    # Lift the bilinear products f̂(j) f̂(k) into a moment matrix M_top of
    # size (2T_max+1)x(2T_max+1) and Schur-bordered to PSD-constrain
    # M_top ⪰ xi xi^T. Then enforce f^2 <= f via the Hermitian PSD
    # localizing-matrix constraint Loc[j,k] = (f - f^2)̂(j-k) ⪰ 0.
    # Tests f^2 <= f against ALL nonneg deg-2 trig polynomials
    # (strict tightening of Bochner-on-f, which only tests against |p|^2).
    _M_top_lifted = None  # captures the Lasserre M_top for downstream use
    if lasserre_T_max > 0:
        import importlib.util as _ilu
        from pathlib import Path as _Path
        _here = _Path(__file__).resolve().parent
        _spec = _ilu.spec_from_file_location("lasserre", _here / "lasserre.py")
        _las = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_las)
        _T_loc = lasserre_T_loc if lasserre_T_loc > 0 else lasserre_T_max
        _M_top_lifted = _las.add_lasserre2_constraint(
            cons, c, d, T_max=lasserre_T_max, T_loc=_T_loc
        )

    # ----- M-side Bochner via Lasserre-lifted bilinears (NEW)
    # Reuses the Lasserre level-2 moment-matrix M_top to substitute
    # |f̂(m)|² = (c_m² + d_m²)/4 = (M_top[m,m] + M_top[T_max+m, T_max+m]) / 4
    # EXACTLY (no slack), turning the M-side Toeplitz PSD into a standard
    # linear-in-(Ω, c, d, M_top) PSD constraint.
    if mside_bochner_lasserre_n > 0:
        if _M_top_lifted is None:
            raise ValueError(
                "mside_bochner_lasserre_n > 0 requires lasserre_T_max > 0 "
                "(the Lasserre level-2 moment matrix M_top is needed for the "
                "exact-bilinear M-side encoding)."
            )
        import importlib.util as _ilu
        from pathlib import Path as _Path
        _here = _Path(__file__).resolve().parent
        _spec = _ilu.spec_from_file_location(
            "mside_via_lasserre", _here / "mside_via_lasserre.py"
        )
        _mvl = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_mvl)
        _mvl.add_mside_bochner_lasserre_constraint(
            cons, c, d, Omega, _M_top_lifted,
            T_max=lasserre_T_max, n_M=mside_bochner_lasserre_n,
        )

    return Omega, w, v, c, d, eps, dlt, cons


def solve_full_program(
    N, T, R, h1, h2, p1, p2, q1, q2,
    cell_mode="exact", solver="CLARABEL", verbose=False,
    use_T3=False, use_T5=False, use_T5p=False, bochner_n=0, mside_bochner_n=0,
    mside_bochner_schur_n=0,
    assume_even=False,
    lasserre_T_max=0, lasserre_T_loc=0,
    mside_bochner_lasserre_n=0,
):
    Omega, w, v, c, d, eps, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, q1, q2, cell_mode,
        use_T3=use_T3, use_T5=use_T5, use_T5p=use_T5p, bochner_n=bochner_n,
        mside_bochner_n=mside_bochner_n,
        mside_bochner_schur_n=mside_bochner_schur_n,
        assume_even=assume_even,
        lasserre_T_max=lasserre_T_max, lasserre_T_loc=lasserre_T_loc,
        mside_bochner_lasserre_n=mside_bochner_lasserre_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=solver, verbose=verbose)
    return {
        "value": prob.value, "status": prob.status, "time": time.time() - t0,
        "Omega": Omega.value, "w": w.value, "v": v.value,
        "c": c.value, "d": d.value, "eps": eps.value, "dlt": dlt.value,
    }


WHITE_TABLE3 = [
    # (h, p, q1, q2, label)
    (0.015, 0.381,   -0.02, 0.02, "row1"),
    (0.015, 0.385,   -0.02, 0.02, "row2"),
    (0.020, 0.375,   -0.02, 0.02, "row3"),
    (0.004, 0.3875,  -0.02, 0.02, "row4"),
    (0.000, 0.4,     -0.02, 0.02, "row5"),
    (0.000, 0.381,   -0.02, 0.02, "row6"),
    (0.030, 0.375,   -0.02, 0.02, "row7"),
]


def sweep(N, T, R, mode="exact", solver="CLARABEL"):
    print(f"--- N={N}, T={T}, R={R}, mode={mode} ---")
    bounds = []
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        res = solve_full_program(N, T, R, h, h, p, p, qm, qp, cell_mode=mode, solver=solver)
        print(f"  {label}: h={h:.3f} p={p:.4f} q∈[{qm:.2f},{qp:.2f}]: "
              f"Ω*={res['value']:.7f}  ({res['status']}, {res['time']:.1f}s)")
        bounds.append(res["value"])
    print(f"  → MIN over Table 3 rows: {min(bounds):.7f}  (White's reported: 0.379005)")
    return bounds


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    sweep(N=5000, T=2000, R=10, mode="exact")
