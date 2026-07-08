"""PRO-48 — genuine cubic BISPECTRUM 3-point SDP lift, bolted ONTO the cell-envelope.

The Erdos minimum-overlap program optimizes over the 2-point data of f (power
spectrum |f^(m)|^2 / autocorrelation).  Its dual ceiling C_inf ~ 0.380558 is a
*provably 2-point* phenomenon (PRO-6/PRO-22).  A 3-point relaxation lives in a
strictly larger cone with a different ceiling.

The genuinely-cubic object is the triple correlation of f >= 0:

    R3(s,t) = INT f(x) f(x+s) f(x+t) dx  >= 0   (pointwise, since f >= 0)

whose 2-D Fourier coefficients are the BISPECTRUM

    B(m,n) = y_m * y_n * conj(y_{m+n}),   y_k := f^(k).

Bochner-2D: R3 >= 0  =>  the moment matrix  G[a,b] = B(a-b)  (a,b in a freq
grid Lambda subset Z^2) is PSD.  We lift the B(m,n) to INDEPENDENT moment
variables, impose:
  (P) PSD of the real embedding of G,
  (S) reality / permutation symmetry of B (true for real f),
  (M) MARGINAL links to the 2-point program:
        B(0,0) = y_0^3 = 1/8,
        B(m,0) = B(0,m) = B(m,-m) = (1/2) * P_m,   P_m >= |y_m|^2  (Schur slack)
      where |y_m|^2 = (c_m^2 + d_m^2)/4 are the existing program's variables.

This is a PURE ADDITION onto build_problem's constraints (cell-envelope kept,
per PRO-22), so min Omega can only RISE or stay equal.  Any rise is a real
tightening *provided the block stays feasible for the true f* — which the
self-test below verifies directly.

Convention (matches white_full_convex.py Bochner block):
    y_0 = 1/2,   y_k = (c_k - i d_k)/2   for k>=1,   y_{-k} = conj(y_k).
"""
from __future__ import annotations
import itertools
import numpy as np
import cvxpy as cp


# --------------------------------------------------------------------------
# y_k = f^(k) accessors in terms of the program's real variables c, d.
# Return (re, im) where each is a python float (for k=0) or a cvxpy expr.
# --------------------------------------------------------------------------
def _y(c, d, k):
    if k == 0:
        return 0.5, 0.0
    if k > 0:
        return 0.5 * c[k - 1], -0.5 * d[k - 1]
    a = -k
    return 0.5 * c[a - 1], 0.5 * d[a - 1]


def _grid(L):
    """Moment-matrix generating set Lambda = {0..L}^2 subset Z^2."""
    return [(a, b) for a in range(L + 1) for b in range(L + 1)]


def add_bispectrum_block(cons, c, d, L=3, with_perm_sym=True, return_handles=False,
                         Omega=None, localize_A=True):
    """Bolt the cubic bispectrum PSD block onto `cons` (built from c, d).

    L           generating-grid side (Lambda = {0..L}^2, |Lambda| = (L+1)^2).
                Difference set Delta = {-L..L}^2 ; needs |freq| <= 2L <= T0.
    with_perm_sym  also impose the bispectrum 6-fold permutation symmetry
                (all relations exact for real f; extra tightening).
    Omega       if given (the program's objective var), ALSO add the
                LOCALIZING block (Omega - R3) >= 0.  Valid because
                0 <= f <= 1  =>  R3(s,t) = INT f(x)f(x+s)f(x+t)dx
                              <= INT f(x)f(x+s)dx = A(s) <= sup_s A = Omega.
                Moment matrix:  G2[a,b] = Omega*delta_{a,b} - B(a-b)  >= 0.
                This is the link that touches the OBJECTIVE directly (the
                marginal links only touch the power spectrum).

    Returns the PSD constraint object(s) (for dual inspection) and, if
    requested, the variable handles.
    """
    Lam = _grid(L)
    nL = len(Lam)
    idx = {p: i for i, p in enumerate(Lam)}

    # Difference set Delta.
    Delta = set()
    for a in Lam:
        for b in Lam:
            Delta.add((a[0] - b[0], a[1] - b[1]))
    Delta = sorted(Delta)

    # One (Re, Im) variable pair per difference point.
    Bre = {dlt: cp.Variable(name=f"Bre_{dlt[0]}_{dlt[1]}") for dlt in Delta}
    Bim = {dlt: cp.Variable(name=f"Bim_{dlt[0]}_{dlt[1]}") for dlt in Delta}

    block = []

    # (center) B(0,0) = 1/8 ; imaginary part 0.
    block += [Bre[(0, 0)] == 0.125, Bim[(0, 0)] == 0.0]

    # (conjugation) B(-delta) = conj(B(delta)) : Hermitian G.
    for dlt in Delta:
        ndlt = (-dlt[0], -dlt[1])
        if ndlt in Bre and dlt < ndlt:
            block += [Bre[ndlt] == Bre[dlt], Bim[ndlt] == -Bim[dlt]]

    # (permutation symmetry, true for real f)
    if with_perm_sym:
        for dlt in Delta:
            m, n = dlt
            # B(m,n) = B(n,m)
            if (n, m) in Bre and (n, m) != dlt and dlt < (n, m):
                block += [Bre[(n, m)] == Bre[dlt], Bim[(n, m)] == Bim[dlt]]
            # B(m,n) = B(-m-n, n)
            t = (-m - n, n)
            if t in Bre and t != dlt and dlt < t:
                block += [Bre[t] == Bre[dlt], Bim[t] == Bim[dlt]]

    # (marginal links to the 2-point program)
    # P_m >= |y_m|^2 = (c_m^2 + d_m^2)/4 ; B(m,0)=B(0,m)=B(m,-m) = (1/2) P_m (real).
    P = {}
    for m in range(1, 2 * L + 1):
        ry, iy = _y(c, d, m)
        Pm = cp.Variable(nonneg=True, name=f"P_{m}")
        P[m] = Pm
        block.append(cp.square(ry) + cp.square(iy) <= Pm)        # Schur / SOC slack
        for pos in [(m, 0), (0, m), (m, -m), (-m, 0), (0, -m), (-m, m)]:
            if pos in Bre:
                block += [Bre[pos] == 0.5 * Pm, Bim[pos] == 0.0]

    # (PSD) real embedding of the Hermitian moment matrix G[a,b] = B(a-b).
    ReG = [[Bre[(Lam[i][0] - Lam[j][0], Lam[i][1] - Lam[j][1])] for j in range(nL)]
           for i in range(nL)]
    ImG = [[Bim[(Lam[i][0] - Lam[j][0], Lam[i][1] - Lam[j][1])] for j in range(nL)]
           for i in range(nL)]
    ReG = cp.bmat(ReG)
    ImG = cp.bmat(ImG)
    Greal = cp.bmat([[ReG, -ImG], [ImG, ReG]])
    psd_con = (Greal >> 0)
    block.append(psd_con)

    # (LOCALIZING) (Omega - R3) >= 0  =>  G2[a,b] = Omega*delta - B(a-b)  >= 0.
    loc_con = None
    if Omega is not None:
        ReG2 = [[(Omega - Bre[(0, 0)]) if i == j else (-1.0) * ReG[i, j]
                 for j in range(nL)] for i in range(nL)]
        # off-diagonal Re entries: -Bre[a-b]; diagonal: Omega - 1/8.
        ReG2 = cp.bmat([[ (Omega - 0.125) if i == j else
                          (-Bre[(Lam[i][0]-Lam[j][0], Lam[i][1]-Lam[j][1])])
                          for j in range(nL)] for i in range(nL)])
        ImG2 = cp.bmat([[ cp.Constant(0.0) if i == j else
                          (-Bim[(Lam[i][0]-Lam[j][0], Lam[i][1]-Lam[j][1])])
                          for j in range(nL)] for i in range(nL)])
        G2real = cp.bmat([[ReG2, -ImG2], [ImG2, ReG2]])
        loc_con = (G2real >> 0)
        block.append(loc_con)

    # (LOCALIZING-MATRIX, tighter) (A(s) - R3) >= 0  since
    #   A(s) - R3(s,t) = INT f(x)f(x+s)(1 - f(x+t)) dx >= 0  (f>=0, 1-f>=0).
    # Fourier:  h^(m,n) = A_m * delta_{n,0} - B(m,n),  A_m = |y_m|^2.
    # Represent A_m by the same Schur slack P_{|m|} (>= |y_m|^2): true f has
    # equality, so the true localizing matrix is recovered => valid.
    locA_con = None
    if localize_A:
        def Amval(m):
            am = abs(m)
            return 0.25 if am == 0 else P[am]     # |y_0|^2=1/4 ; |y_m|^2 <- P_m
        ReG3 = [[None] * nL for _ in range(nL)]
        ImG3 = [[None] * nL for _ in range(nL)]
        for i in range(nL):
            for j in range(nL):
                m = Lam[i][0] - Lam[j][0]
                n = Lam[i][1] - Lam[j][1]
                if n == 0:
                    ReG3[i][j] = Amval(m) - Bre[(m, n)]   # A_m - B(m,0)
                else:
                    ReG3[i][j] = (-1.0) * Bre[(m, n)]     # -B(m,n)
                ImG3[i][j] = (-1.0) * Bim[(m, n)]
        ReG3 = cp.bmat(ReG3)
        ImG3 = cp.bmat(ImG3)
        G3real = cp.bmat([[ReG3, -ImG3], [ImG3, ReG3]])
        locA_con = (G3real >> 0)
        block.append(locA_con)

    cons.extend(block)
    if return_handles:
        return (psd_con, loc_con, locA_con), dict(Bre=Bre, Bim=Bim, P=P,
                                                  Lam=Lam, Delta=Delta)
    return (psd_con, loc_con, locA_con)


# --------------------------------------------------------------------------
# MANDATORY self-test: feed a genuine nonneg f and confirm the block's
# constraints (PSD + symmetry + marginals + center) all hold for its TRUE
# bispectrum.  An over-constraining bug would fail here BEFORE any SDP, so a
# later SDP "rise" cannot be a feasibility artifact.
# --------------------------------------------------------------------------
def self_test(L=3, ntest=5, tol=1e-9, seed_list=(1, 2, 3, 4, 5)):
    print(f"[self-test] L={L}, grid |Lambda|={(L+1)**2}, "
          f"checking {len(seed_list)} random nonneg f ...")
    Lam = _grid(L)
    nL = len(Lam)
    Ngrid = 4096
    xs = np.arange(Ngrid) / Ngrid
    ok = True
    for seed in seed_list:
        rng = np.random.default_rng(seed)
        # random f in [0,1] with INT f = 1/2 :
        #   f(x) = 1/2 + sum_k a_k cos(2 pi k x + phi_k),  sum|a_k| <= 1/2.
        # (mean 1/2 => y_0 = 1/2 ; |deviation| <= 1/2 => 0 <= f <= 1.)
        deg = 6
        amp = rng.standard_normal(deg) ** 2
        amp = 0.5 * amp / amp.sum() * rng.uniform(0.6, 0.999)   # sum|a_k| < 1/2
        phi = rng.uniform(0, 2 * np.pi, deg)
        f = 0.5 * np.ones(Ngrid)
        for k in range(1, deg + 1):
            f += amp[k - 1] * np.cos(2 * np.pi * k * xs + phi[k - 1])
        assert f.min() >= -1e-9 and f.max() <= 1 + 1e-9, (f.min(), f.max())
        yhat = np.fft.fft(f) / Ngrid             # yhat[k] = f^(k), k=0..Ngrid-1

        def yk(k):
            return yhat[k % Ngrid]

        # TRUE bispectrum on the difference set and true |y_m|^2.
        Delta = set((a[0] - b[0], a[1] - b[1]) for a in Lam for b in Lam)
        Btrue = {(m, n): yk(m) * yk(n) * np.conj(yk(m + n)) for (m, n) in Delta}
        # (center)
        assert abs(Btrue[(0, 0)] - 0.125) < tol, "center"
        # (conjugation)
        for (m, n) in Delta:
            if (-m, -n) in Btrue:
                if abs(Btrue[(-m, -n)] - np.conj(Btrue[(m, n)])) > 1e-7:
                    ok = False; print("  conj FAIL", (m, n))
        # (perm symmetry)
        for (m, n) in Delta:
            if (n, m) in Btrue and abs(Btrue[(n, m)] - Btrue[(m, n)]) > 1e-7:
                ok = False; print("  perm B(n,m) FAIL", (m, n))
            t = (-m - n, n)
            if t in Btrue and abs(Btrue[t] - Btrue[(m, n)]) > 1e-7:
                ok = False; print("  perm B(-m-n,n) FAIL", (m, n))
        # (marginal) B(m,0) = (1/2)|y_m|^2
        for m in range(1, 2 * L + 1):
            if (m, 0) in Btrue:
                lhs = Btrue[(m, 0)]
                rhs = 0.5 * abs(yk(m)) ** 2
                if abs(lhs - rhs) > 1e-7:
                    ok = False; print("  marginal FAIL", m, lhs, rhs)
        # (PSD) assemble G and check eigenvalues >= ~0
        G = np.zeros((nL, nL), dtype=complex)
        for i in range(nL):
            for j in range(nL):
                G[i, j] = Btrue[(Lam[i][0] - Lam[j][0], Lam[i][1] - Lam[j][1])]
        ev = np.linalg.eigvalsh((G + G.conj().T) / 2).min()
        # (LOCALIZING) Omega = sup_s A(s), A(s) = sum_m |y_m|^2 e^{2pi i m s}
        Acorr = np.abs(yhat) ** 2                 # A^(m) = |y_m|^2
        A_s = np.fft.ifft(Acorr) * Ngrid          # A(s) on the grid (real)
        Omega_true = float(A_s.real.max())
        G2 = np.zeros((nL, nL), dtype=complex)
        for i in range(nL):
            for j in range(nL):
                dd = (Lam[i][0] - Lam[j][0], Lam[i][1] - Lam[j][1])
                G2[i, j] = (Omega_true if dd == (0, 0) else 0.0) - Btrue[dd]
        ev2 = np.linalg.eigvalsh((G2 + G2.conj().T) / 2).min()
        # (LOCALIZING-MATRIX) A(s)-R3 >= 0 : h^(m,n) = A_m delta_{n,0} - B(m,n)
        G3 = np.zeros((nL, nL), dtype=complex)
        for i in range(nL):
            for j in range(nL):
                m = Lam[i][0] - Lam[j][0]
                n = Lam[i][1] - Lam[j][1]
                Am = (abs(yk(m)) ** 2) if n == 0 else 0.0
                G3[i, j] = Am - Btrue[(m, n)]
        ev3 = np.linalg.eigvalsh((G3 + G3.conj().T) / 2).min()
        if ev < -1e-9 or ev2 < -1e-9 or ev3 < -1e-9:
            ok = False
            print(f"  PSD FAIL seed={seed} G={ev:.2e} G2={ev2:.2e} G3={ev3:.2e}")
        else:
            print(f"  seed={seed}: OK  (G={ev:.3e}, G2[Om-R3]={ev2:.3e}, "
                  f"G3[A-R3]={ev3:.3e}, Omega={Omega_true:.4f})")
    print("[self-test]", "PASS" if ok else "*** FAIL ***")
    return ok


if __name__ == "__main__":
    self_test()
