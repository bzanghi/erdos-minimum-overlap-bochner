"""
INDEPENDENT re-derivation of White (2023) Section 5.1 / Appendix II
ellipse-extension argument applied to the Bochner-augmented program.

Author: independent verifier (no peeking at path_b_analytical.py).

White's argument (Lemma 10 / §8.3):
  The dual SOCP's CONSTRAINTS do not depend on (h1, h2, p1, p2). Hence any
  feasible dual point gives a valid LB for any (h1, h2, p1, p2) just by
  re-evaluating the dual OBJECTIVE.

Equivalently, by the envelope theorem:  d/dθ [primal_opt] = ± λ * d/dθ [rhs(θ)]
for each parameter-bearing constraint with multiplier λ.  A FIRST-ORDER
re-evaluation around (h_c, p_c) gives a valid LB function Φ_row(h, p, q):

  For each constraint of form "LHS ≥ rhs(θ)" with cvxpy dual λ ≥ 0:
      contribution = + λ * (rhs(θ) - rhs(θ_c))
  For each constraint of form "LHS ≤ rhs(θ)" with cvxpy dual λ ≥ 0:
      contribution = - λ * (rhs(θ) - rhs(θ_c))

This is EXACT, not first-order, when only the rhs depends on θ (the dual
function inherits the same affine/quadratic dependence on θ as the rhs).

Parameter-bearing constraints (cvxpy primal in white_full_convex.py, with h=h1=h2):
  (5.3)    L^2 sum(jw - (j-1)v) >= h          rhs = h, linear in h
  (5.4)    L^3 sum((j-1)^2(w+v)) <= 2/3 + h^2/2   rhs depends on h^2
  (5.12pL) c[0] >= p                          rhs = p
  (5.12pU) c[0] <= p                          rhs = p
  (5.12qL) d[0] >= q                          rhs = q   (q1 = q used at q_c1 = -0.02 originally)
  (5.12qU) d[0] <= q                          rhs = q   (q2 = q used at q_c2 = +0.02 originally)
  (5.13)   (L/2)(a_plus_2 (w+v)) >= -0.5*(p^2 + max(q1^2,q2^2))
                                              rhs depends on p^2 (and q^2 via max)

So Φ_row(h, p, q1, q2) = V_c
   + λ_53 * (h - h_c)
   - λ_54 * (h^2/2 - h_c^2/2)
   + λ_pL * (p - p_c)
   - λ_pU * (p - p_c)
   + λ_qL * (q1 - q1_c)
   - λ_qU * (q2 - q2_c)
   - λ_513 * (p^2/2 - p_c^2/2)
   - λ_513 * (max(q1^2,q2^2)/2 - max(q1_c^2, q2_c^2)/2)

For the residual region (5.16):  h ∈ [0, 0.06], p = c1 ∈ [0.35, 0.45],
q = d1 ∈ [-0.02, 0.02]. To get a valid LB on f* with d1* = q, we need a dual
that is feasible for an SDP with q1 = q2 = q (a single-point d1 constraint).
For that dual evaluation, set q1 = q2 = q in Lemma 10:
   contribution_q = λ_qL*(q - q1_c) - λ_qU*(q - q2_c) - λ_513*(q^2/2 - max(q1_c^2,q2_c^2)/2)

The cvxpy SDP is solved with q1_c = -0.02, q2_c = +0.02 (so the dual is feasible
for THAT q-interval). Lemma 10 says re-evaluating the SAME dual variables with a
NARROWER interval q1' = q2' = q is always a valid LB (constraints don't change;
only the objective contribution does, monotonically by the formula since it
shrinks the feasible set, but the dual objective is recomputed). The Φ formula
above gives that re-evaluation.
"""
from __future__ import annotations
import os, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp
from white_full_convex import build_problem


# White's Table 3 — 7 ellipse centers
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


def solve_one_center(N, T, R, h_c, p_c, q1_c, q2_c, bochner_n=20, solver="CLARABEL"):
    """Solve the Bochner-augmented SDP at a single (h_c, p_c, q1_c, q2_c) center,
    extract dual values for each parameter-bearing constraint."""
    L = 2.0 / N
    j = np.arange(1, N + 1)

    # Build everything from scratch but using the same constraint set as build_problem,
    # explicitly tracking which constraint object corresponds to which dual.
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

    # (5.3) L^2 sum(jw - (j-1)v) >= h1
    con_53 = L ** 2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h_c
    cons.append(con_53)

    # (5.4) L^3 sum((j-1)^2(w+v)) <= 2/3 + h2^2/2
    con_54 = L ** 3 * cp.sum(cp.multiply((j - 1) ** 2, (w + v))) <= 2.0 / 3 + h_c ** 2 / 2
    cons.append(con_54)

    # cell bounds (exact, from white_full_convex)
    from white_full_convex import cos_cell_bounds_exact, sin_cell_bounds_exact, odd_coeff_factors, tail_bound_eps, tail_bound_delta

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
        a_minus, _ = cos_cell_bounds_exact(j, m, L)
        lhs = (L / 2) * (a_minus @ (w + v))
        sin_pi_half_m = np.sin(np.pi * m / 2)
        rhs_lin = (4 * sin_pi_half_m / (m * np.pi)) * am
        cons.append(lhs + 2 * cp.square(am) + 2 * cp.square(bm) - rhs_lin <= 0)

    for m in range(1, 2 * R + 1):
        bm = b_expr[m - 1]
        b_minus, b_plus = sin_cell_bounds_exact(j, m, L)
        sin_pi_half_m = np.sin(np.pi * m / 2)
        # White's 2026-05-31 email correction: constraints 5.6/5.7 had an 8 in the RHS numerator, should be 4.
        rhs = -(4.0 / (m * np.pi)) * sin_pi_half_m * bm
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

    # (5.12) p, q box constraints — track each individually
    con_512_pL = c[0] >= p_c   # c1 >= p1
    con_512_pU = c[0] <= p_c   # c1 <= p2
    con_512_qL = d[0] >= q1_c  # d1 >= q1
    con_512_qU = d[0] <= q2_c  # d1 <= q2
    cons += [con_512_pL, con_512_pU, con_512_qL, con_512_qU]

    # (5.13) (L/2)(a_plus_2 (w+v)) >= -0.5*(max(p1^2,p2^2) + max(q1^2,q2^2))
    _, a_plus_2 = cos_cell_bounds_exact(j, 2, L)
    rhs_513 = -0.5 * (max(p_c ** 2, p_c ** 2) + max(q1_c ** 2, q2_c ** 2))
    con_513 = (L / 2) * (a_plus_2 @ (w + v)) >= rhs_513
    cons.append(con_513)

    # ----- Bochner moment matrix PSD (n_b)
    n_b = min(bochner_n, T)
    for sign in [+1, -1]:
        half_val = 0.5
        Re_rows, Im_rows = [], []
        for jj in range(n_b + 1):
            re_row, im_row = [], []
            for kk in range(n_b + 1):
                ell = jj - kk
                if ell == 0:
                    re_row.append(cp.Constant(half_val))
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

    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=solver, verbose=False)
    elapsed = time.time() - t0

    return {
        "status": prob.status,
        "value": float(prob.value),
        "time": elapsed,
        "duals": {
            "lam_53":  float(con_53.dual_value),
            "lam_54":  float(con_54.dual_value),
            "lam_pL":  float(con_512_pL.dual_value),
            "lam_pU":  float(con_512_pU.dual_value),
            "lam_qL":  float(con_512_qL.dual_value),
            "lam_qU":  float(con_512_qU.dual_value),
            "lam_513": float(con_513.dual_value),
        },
        "center": {"h_c": h_c, "p_c": p_c, "q1_c": q1_c, "q2_c": q2_c},
    }


def Phi_row(rec, h, p, q):
    """Reconstructed dual lower-bound as function of (h, p, q).
    Uses q1' = q2' = q (single-point d1 constraint, Lemma 10 valid)."""
    duals = rec["duals"]
    cen = rec["center"]
    Vc = rec["value"]

    L53 = duals["lam_53"]
    L54 = duals["lam_54"]
    LpL = duals["lam_pL"]
    LpU = duals["lam_pU"]
    LqL = duals["lam_qL"]
    LqU = duals["lam_qU"]
    L513 = duals["lam_513"]

    h_c = cen["h_c"]; p_c = cen["p_c"]; q1_c = cen["q1_c"]; q2_c = cen["q2_c"]
    max_q2_c = max(q1_c ** 2, q2_c ** 2)

    # Convention: cvxpy dual >= 0 for inequality. For (LHS >= rhs(θ)), envelope says
    # d val/d θ = + λ * (d rhs/d θ). For (LHS <= rhs(θ)), d val/dθ = - λ * (d rhs/d θ).
    # Wait — let me redo. cvxpy: for "expr <= 0" Lagrangian is +λ*expr. So:
    #   "LHS - rhs <= 0" (i.e. LHS <= rhs):  Lag += λ*(LHS - rhs).  d/dθ_rhs = -λ.
    #   "rhs - LHS <= 0" (i.e. LHS >= rhs):  Lag += λ*(rhs - LHS).  d/dθ_rhs = +λ.
    # By envelope: d primal_opt / d θ_rhs equals the d Lag / d θ_rhs at optimum:
    #   <= form: d/d rhs = -λ
    #   >= form: d/d rhs = +λ
    # So Φ(θ) = V_c + Σ (sign_i * λ_i) * (rhs_i(θ) - rhs_i(θ_c)).

    # (5.3) >= h:  +λ_53 * (h - h_c)
    delta = L53 * (h - h_c)
    # (5.4) <= 2/3 + h^2/2:  -λ_54 * ((2/3 + h^2/2) - (2/3 + h_c^2/2)) = -λ_54 * (h^2 - h_c^2)/2
    delta += -L54 * (h ** 2 / 2 - h_c ** 2 / 2)
    # (5.12pL) c1 >= p:  +λ_pL * (p - p_c)
    delta += LpL * (p - p_c)
    # (5.12pU) c1 <= p:  -λ_pU * (p - p_c)
    delta += -LpU * (p - p_c)
    # (5.12qL) d1 >= q:   +λ_qL * (q - q1_c)
    delta += LqL * (q - q1_c)
    # (5.12qU) d1 <= q:   -λ_qU * (q - q2_c)
    delta += -LqU * (q - q2_c)
    # (5.13) >= -0.5*(p^2 + max(q1^2,q2^2)):  +λ_513 * (rhs(p,q) - rhs(p_c,q_c))
    rhs_513_new = -0.5 * (p ** 2 + q ** 2)
    rhs_513_old = -0.5 * (p_c ** 2 + max_q2_c)
    delta += L513 * (rhs_513_new - rhs_513_old)

    return Vc + delta


def grid_min(records, h_grid, p_grid, q_grid):
    """Compute min over (h, p, q) of max over rows of Phi_row(...).
    Returns minimum value and its location."""
    best = float('inf')
    loc = None
    binding = None
    for h in h_grid:
        for p in p_grid:
            for q in q_grid:
                # max over rows = best LB available at this point
                best_at = -float('inf')
                best_row = None
                for rec in records:
                    val = Phi_row(rec, h, p, q)
                    if val > best_at:
                        best_at = val
                        best_row = rec.get("label", "?")
                if best_at < best:
                    best = best_at
                    loc = (h, p, q)
                    binding = best_row
    return best, loc, binding


def grid_min_vectorized(records, h_grid, p_grid, q_grid):
    """Vectorized over (p, q) per h slice — bounds memory.
    For each (h, p, q): compute max over rows of Phi_row, then min over the grid."""
    nrec = len(records)
    P, Q = np.meshgrid(p_grid, q_grid, indexing='ij')  # (np, nq)
    Pf = P.ravel(); Qf = Q.ravel()  # (np*nq,)

    # Precompute per-record params
    params = []
    for rec in records:
        d = rec["duals"]; cen = rec["center"]; Vc = rec["value"]
        h_c = cen["h_c"]; p_c = cen["p_c"]; q1_c = cen["q1_c"]; q2_c = cen["q2_c"]
        max_q2_c = max(q1_c ** 2, q2_c ** 2)
        L53, L54, LpL, LpU, LqL, LqU, L513 = (
            d["lam_53"], d["lam_54"], d["lam_pL"], d["lam_pU"],
            d["lam_qL"], d["lam_qU"], d["lam_513"],
        )
        params.append((Vc, h_c, p_c, q1_c, q2_c, max_q2_c, L53, L54, LpL, LpU, LqL, LqU, L513))

    best_min = float('inf')
    best_loc = None
    best_row = None
    for h in h_grid:
        # vals[r, idx] = Phi at (h, Pf[idx], Qf[idx]) for record r
        max_per_pt = np.full(Pf.shape, -np.inf)
        argmax_row = np.zeros(Pf.shape, dtype=int)
        for r, prm in enumerate(params):
            (Vc, h_c, p_c, q1_c, q2_c, max_q2_c,
             L53, L54, LpL, LpU, LqL, LqU, L513) = prm
            delta = (
                L53 * (h - h_c)
                - L54 * (h ** 2 / 2 - h_c ** 2 / 2)
                + LpL * (Pf - p_c) - LpU * (Pf - p_c)
                + LqL * (Qf - q1_c) - LqU * (Qf - q2_c)
                + L513 * (-0.5 * (Pf ** 2 + Qf ** 2) - (-0.5 * (p_c ** 2 + max_q2_c)))
            )
            vals_r = Vc + delta
            mask = vals_r > max_per_pt
            max_per_pt = np.where(mask, vals_r, max_per_pt)
            argmax_row = np.where(mask, r, argmax_row)
        imin = int(np.argmin(max_per_pt))
        if max_per_pt[imin] < best_min:
            best_min = float(max_per_pt[imin])
            best_loc = (float(h), float(Pf[imin]), float(Qf[imin]))
            best_row = records[int(argmax_row[imin])].get("label", "?")

    return best_min, best_loc, best_row


def main(N=10000, T=4000, R=10, bochner_n=20, n_grid=1001, out_path=None,
         checkpoint_dir=None, only_label=None):
    print(f"[indep] Solving 7 SDPs with N={N} T={T} R={R} bochner_n={bochner_n}")
    records = []
    t_start = time.time()
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        if only_label is not None and label != only_label:
            continue
        # Try to load from checkpoint
        ckpt = None
        if checkpoint_dir:
            ckpt = os.path.join(checkpoint_dir, f"{label}.json")
            if os.path.exists(ckpt):
                with open(ckpt) as f:
                    rec = json.load(f)
                rec["label"] = label
                records.append(rec)
                print(f"  [{label}] loaded from checkpoint: V_c={rec['value']:.10f}")
                continue
        print(f"  [{label}] center=(h={h}, p={p}, q1={qm}, q2={qp})  ...", flush=True)
        rec = solve_one_center(N, T, R, h, p, qm, qp, bochner_n=bochner_n)
        rec["label"] = label
        records.append(rec)
        print(f"    status={rec['status']}, V_c={rec['value']:.10f}, time={rec['time']:.1f}s")
        print(f"    duals: " + ", ".join(f"{k}={v:.4e}" for k, v in rec["duals"].items()))
        if ckpt:
            os.makedirs(checkpoint_dir, exist_ok=True)
            with open(ckpt, "w") as f:
                json.dump(rec, f, indent=2)
            print(f"    saved checkpoint: {ckpt}")
    if only_label is not None:
        return None

    # Now grid-search over the residual region (5.16)
    h_grid = np.linspace(0.0, 0.06, n_grid)
    p_grid = np.linspace(0.35, 0.45, n_grid)
    q_grid = np.linspace(-0.02, 0.02, 41)  # finer in (h,p), coarser in q
    print(f"[indep] Grid: {n_grid} x {n_grid} x {len(q_grid)} = {n_grid*n_grid*len(q_grid)} points")

    # vectorized search
    min_val, loc, binding = grid_min_vectorized(records, h_grid, p_grid, q_grid)
    print(f"[indep] MIN over (5.16) = {min_val:.10f}")
    print(f"        at (h, p, q) = {loc}")
    print(f"        binding row: {binding}")

    # also corner check
    corners = [(h, p, q) for h in [0.0, 0.06] for p in [0.35, 0.45] for q in [-0.02, 0.0, 0.02]]
    print("[indep] Corner check:")
    for (h, p, q) in corners:
        vmax = max(Phi_row(rec, h, p, q) for rec in records)
        print(f"   (h={h}, p={p}, q={q}): max_row = {vmax:.10f}")

    elapsed_total = time.time() - t_start
    print(f"[indep] Total time: {elapsed_total:.1f}s")

    out = {
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n,
                   "n_grid_h": n_grid, "n_grid_p": n_grid, "n_grid_q": len(q_grid)},
        "records": [
            {"label": r["label"], "value": r["value"], "status": r["status"],
             "center": r["center"], "duals": r["duals"], "time": r["time"]}
            for r in records
        ],
        "result": {
            "min_value": min_val,
            "loc_h": loc[0], "loc_p": loc[1], "loc_q": loc[2],
            "binding_row": binding,
            "vs_white_0p379005": min_val - 0.379005,
        },
        "total_time_s": elapsed_total,
    }
    if out_path:
        with open(out_path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"[indep] Wrote {out_path}")
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    # Defaults: a quick run with smaller N, then upgrade
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=10000)
    ap.add_argument("--T", type=int, default=4000)
    ap.add_argument("--R", type=int, default=10)
    ap.add_argument("--bochner_n", type=int, default=20)
    ap.add_argument("--n_grid", type=int, default=1001)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--ckpt", type=str, default=None)
    ap.add_argument("--only", type=str, default=None)
    args = ap.parse_args()
    main(N=args.N, T=args.T, R=args.R, bochner_n=args.bochner_n,
         n_grid=args.n_grid, out_path=args.out,
         checkpoint_dir=args.ckpt, only_label=args.only)
