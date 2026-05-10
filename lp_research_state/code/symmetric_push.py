"""
CONDITIONAL bound on the Erdős minimum overlap constant µ assuming f* is even.

Under the even-f assumption (d_k = 0, dlt_k = 0, v_j = w_j), the SDP variable
count is roughly halved AND the Bochner moment matrix collapses to a real
symmetric (n+1)x(n+1) matrix (the imaginary block vanishes since d=0). This
allows much higher Bochner level n and Lasserre level T_max than the
unconditional encoding.

Only White Table-3 rows 5 (h=0, p=0.4) and 6 (h=0, p=0.381) are feasible
under the even assumption (other rows have h>0, contradicting f even since
∫ x f(x) dx = 0). The MIN over rows {5, 6} is the conditional bound.

Path B ellipse extension under symmetry: parameter space collapses to the
p-axis (h=0, q=0). Verify that the row-5 and row-6 1-D ellipse intervals
jointly cover p in [0.35, 0.45].

Usage:
    python symmetric_push.py build <row> <bochner_n> <lasserre_T_max> <pickle_path>
    python symmetric_push.py solve <pickle_path> <out_json>
    python symmetric_push.py oneshot <row> <bochner_n> <lasserre_T_max> <out_json>
    python symmetric_push.py aggregate <out_dir>
"""
from __future__ import annotations
import sys
import os
import json
import time
import pickle
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp

from white_full_convex import build_problem, WHITE_TABLE3
from path_b_analytical import build_problem_with_dual_handles


UNCONDITIONAL_BOUND = 0.379828

# Only rows feasible under f even (h = 0).
SYMMETRIC_ROWS = [
    (0.000, 0.4, -0.02, 0.02, "row5"),
    (0.000, 0.381, -0.02, 0.02, "row6"),
]


# ---------------------------------------------------------------------------
# Build a problem with assume_even=True and DUAL HANDLES so we can do path B.
# We also expose Bochner under the even-f reduction (Im=0 block).
# ---------------------------------------------------------------------------
def build_problem_even_with_handles(N, T, R, h, p, q1, q2, bochner_n, lasserre_T_max):
    """Wrapper around build_problem that adds (h,p,q)-dependent constraint handles
    (for Path B) and forces assume_even=True."""
    from white_full_convex import (
        cos_cell_bounds_exact, sin_cell_bounds_exact, odd_coeff_factors,
        tail_bound_eps, tail_bound_delta,
    )
    L = 2.0 / N
    j = np.arange(1, N + 1)

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

    # Even-f: d=0, dlt=0, v=w.
    cons += [d == 0, dlt == 0, v == w]

    # 5.3 / 5.4
    con_53 = L**2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h
    con_54 = L**3 * cp.sum(cp.multiply((j - 1) ** 2, (w + v))) <= 2.0 / 3 + h ** 2 / 2
    cons += [con_53, con_54]

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
            am = (eps[(m - 1) // 2]
                  + (2 * m * sin_pi_half_m / np.pi)
                  * (1.0 / (2 * m ** 2) + cp.sum(cp.multiply(af, c))))
            bm = (dlt[(m - 1) // 2]
                  + (4 * sin_pi_half_m / np.pi) * cp.sum(cp.multiply(bf, d)))
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
        rhs = -(8.0 / (m * np.pi)) * sin_pi_half_m * bm
        cons.append((L / 2) * (b_minus @ w - b_plus @ v) <= rhs)
        cons.append((L / 2) * (b_plus @ w - b_minus @ v) >= rhs)

    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
                 cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T)]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)

    # 5.12
    con_512_pL = c[0] >= p
    con_512_pU = c[0] <= p
    con_512_qL = d[0] >= q1
    con_512_qU = d[0] <= q2
    cons += [con_512_pL, con_512_pU, con_512_qL, con_512_qU]

    # 5.13
    _, a_plus_2 = cos_bnds(j, 2, L)
    rhs_513 = -0.5 * (max(p ** 2, p ** 2) + max(q1 ** 2, q2 ** 2))
    con_513 = (L / 2) * (a_plus_2 @ (w + v)) >= rhs_513
    cons.append(con_513)

    # ----- Bochner (real (n+1)x(n+1) form because d=0).
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            half_v = 0.5
            Re_rows = []
            for jj in range(n_b + 1):
                row = []
                for kk in range(n_b + 1):
                    ell = jj - kk
                    if ell == 0:
                        row.append(cp.Constant(half_v))
                    else:
                        aell = abs(ell)
                        row.append(cp.Constant(sign * 0.5) * c[aell - 1])
                Re_rows.append(row)
            Re_M = cp.bmat(Re_rows)
            cons.append(Re_M >> 0)

    # ----- Lasserre level-2 (under d=0): the moment-matrix uses xi=(c_1,...,c_T);
    # localizing matrix uses Hermitian Toeplitz of (f - f^2)^.
    if lasserre_T_max > 0:
        import importlib.util as _ilu
        from pathlib import Path as _Path
        _here = _Path(__file__).resolve().parent
        _spec = _ilu.spec_from_file_location("lasserre", _here / "lasserre.py")
        _las = _ilu.module_from_spec(_spec)
        _spec.loader.exec_module(_las)
        # We pass the same c, d as in the unconditional case; d=0 is enforced
        # already by the cons above. The Lasserre moment matrix in (c,d) becomes
        # block-diagonal with the d-block trivially zero, but cvxpy will still
        # build the larger object. That's fine — variable count is O(T_max) anyway.
        _las.add_lasserre2_constraint(cons, c, d, T_max=lasserre_T_max,
                                      T_loc=lasserre_T_max)

    handles = {
        "con_53": con_53, "con_54": con_54,
        "con_512_pL": con_512_pL, "con_512_pU": con_512_pU,
        "con_512_qL": con_512_qL, "con_512_qU": con_512_qU,
        "con_513": con_513,
        "Omega": Omega, "w": w, "v": v, "c": c, "d": d, "eps": eps, "dlt": dlt,
    }
    return Omega, cons, handles


# ---------------------------------------------------------------------------
# Solve a row at given (bochner_n, lasserre_T_max). Records primal + duals.
# ---------------------------------------------------------------------------
def solve_row(row_label, h, p, q1, q2, N, T, R, bochner_n, lasserre_T_max,
              solver="CLARABEL", time_limit=None, verbose=False):
    Omega, cons, H = build_problem_even_with_handles(
        N, T, R, h, p, q1, q2, bochner_n, lasserre_T_max
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    solver_kwargs = {"verbose": verbose}
    if time_limit is not None:
        solver_kwargs["time_limit"] = time_limit
    try:
        prob.solve(solver=solver, **solver_kwargs)
    except Exception as e:
        return {
            "row": row_label, "h": h, "p": p, "q1": q1, "q2": q2,
            "N": N, "T": T, "R": R, "bochner_n": bochner_n,
            "lasserre_T_max": lasserre_T_max,
            "value": None, "status": f"exception: {type(e).__name__}: {e}",
            "time_seconds": time.time() - t0,
        }
    elapsed = time.time() - t0

    duals = {}
    for key in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                "con_512_qL", "con_512_qU", "con_513"):
        dv = H[key].dual_value
        duals[key] = float(dv) if dv is not None else 0.0

    # Verification: d should be ~0, v-w ~0
    d_max = float(np.max(np.abs(H["d"].value))) if H["d"].value is not None else None
    vmw_max = float(np.max(np.abs(H["v"].value - H["w"].value))) \
        if H["v"].value is not None else None

    return {
        "row": row_label, "h": h, "p": p, "q1": q1, "q2": q2,
        "N": N, "T": T, "R": R,
        "bochner_n": bochner_n, "lasserre_T_max": lasserre_T_max,
        "value": float(prob.value) if prob.value is not None else None,
        "status": prob.status,
        "time_seconds": elapsed,
        "duals": duals,
        "verified_d_zero_max_abs": d_max,
        "verified_v_minus_w_max_abs": vmw_max,
    }


# ---------------------------------------------------------------------------
# Path-B ellipse on the p-axis (h=0, q=0): compute 1-D interval where
# obj(p) = V_c + linear_p (p - p_c) + quad_p (p - p_c)^2 + (h shift terms at h=0)
# >= unconditional_bound.
# Under h=0 fixed and q=q_c fixed, the shift becomes ONLY a function of p:
#   shift(p) = (LpL - LpU) (p - p_c) - 0.5 L513 (p^2 - p_c^2)
# For an even f the Path-B reduces to a 1-D quadratic.
# ---------------------------------------------------------------------------
def find_p_interval(center_value, duals, h_c, p_c, q1, q2,
                    target=UNCONDITIONAL_BOUND, fixed_h=0.0, fixed_q1=-0.02,
                    fixed_q2=0.02):
    """Compute interval of p where obj(p) >= target, with h=fixed_h, q=fixed."""
    qm2 = max(fixed_q1 ** 2, fixed_q2 ** 2)
    qm2_c = max(q1 ** 2, q2 ** 2)
    # h shift (here fixed_h = h_c = 0 typically)
    Drhs_53 = fixed_h - h_c
    Drhs_54 = (fixed_h ** 2 - h_c ** 2) / 2
    h_shift = duals["con_53"] * Drhs_53 - duals["con_54"] * Drhs_54

    # q shift (we evaluate at the same q-range as center, so q-part = 0 unless
    # caller specifies otherwise)
    Drhs_qL = fixed_q1 - q1
    Drhs_qU = fixed_q2 - q2
    Drhs_513_q = -0.5 * (qm2 - qm2_c)
    q_shift = (duals["con_512_qL"] * Drhs_qL
               - duals["con_512_qU"] * Drhs_qU
               + duals["con_513"] * Drhs_513_q)

    # p shift: (LpL - LpU)(p - p_c) - 0.5 * L513 * (p^2 - p_c^2)
    LpL = duals["con_512_pL"]
    LpU = duals["con_512_pU"]
    L513 = duals["con_513"]

    # obj(p) = center_value + h_shift + q_shift
    #          + (LpL - LpU)(p - p_c) - 0.5 L513 (p^2 - p_c^2)
    # = A p^2 + B p + K
    A = -0.5 * L513
    B = (LpL - LpU)
    K0 = center_value + h_shift + q_shift \
         - (LpL - LpU) * p_c + 0.5 * L513 * p_c ** 2
    # obj(p) = A p^2 + B p + K0 >= target  =>  A p^2 + B p + (K0 - target) >= 0

    if abs(A) < 1e-15:
        # Linear case: B p >= target - K0
        if abs(B) < 1e-15:
            covers = K0 >= target
            return {
                "A": A, "B": B, "K0": K0, "p_lo": -np.inf if covers else np.nan,
                "p_hi": np.inf if covers else np.nan, "p_star": p_c,
                "V_max": K0, "concave": False, "covers_uncond": covers,
            }
        if B > 0:
            p_lo = (target - K0) / B; p_hi = np.inf
        else:
            p_lo = -np.inf; p_hi = (target - K0) / B
        return {
            "A": A, "B": B, "K0": K0, "p_lo": p_lo, "p_hi": p_hi,
            "p_star": p_c, "V_max": K0, "concave": False,
            "covers_uncond": True,
        }

    if A < -1e-15:
        # Concave: obj(p_star) = -B^2/(4A) + K0; interval is closed.
        p_star = -B / (2 * A)
        V_max = -B ** 2 / (4 * A) + K0
        if V_max < target:
            return {
                "A": A, "B": B, "K0": K0, "p_lo": np.nan, "p_hi": np.nan,
                "p_star": p_star, "V_max": V_max, "concave": True,
                "covers_uncond": False,
            }
        # Discriminant
        disc = B ** 2 - 4 * A * (K0 - target)
        sqd = np.sqrt(disc)
        # roots: (-B ± sqd) / (2A); since A < 0, larger A*p formula flips order.
        r1 = (-B - sqd) / (2 * A)
        r2 = (-B + sqd) / (2 * A)
        p_lo, p_hi = sorted((r1, r2))
        return {
            "A": A, "B": B, "K0": K0, "p_lo": float(p_lo), "p_hi": float(p_hi),
            "p_star": float(p_star), "V_max": float(V_max),
            "concave": True, "covers_uncond": True,
        }
    else:
        # Convex (A > 0): exterior is feasible. Treat as unbounded on tails.
        # We just check the box [0.35, 0.45].
        return {
            "A": A, "B": B, "K0": K0, "p_lo": -np.inf, "p_hi": np.inf,
            "p_star": p_c, "V_max": np.inf, "concave": False,
            "covers_uncond": True,
        }


# ---------------------------------------------------------------------------
# Aggregation: read all per-row JSONs in a dir, take MIN, run 1-D coverage
# ---------------------------------------------------------------------------
def aggregate(out_dir):
    rows = {}
    for fn in os.listdir(out_dir):
        if not fn.endswith(".json"): continue
        if not (fn.startswith("row5_") or fn.startswith("row6_")):
            continue
        with open(os.path.join(out_dir, fn)) as fh:
            data = json.load(fh)
        rows.setdefault(data["row"], []).append((data, fn))

    print("=== Per-row results ===")
    best = {}
    for label, lst in rows.items():
        # Choose the highest-quality (largest bochner_n + lasserre_T_max)
        # successful (non-None) value.
        ok = [(d, n) for d, n in lst if d.get("value") is not None]
        if not ok:
            print(f"  {label}: no successful runs!")
            continue
        # rank by bochner_n + lasserre_T_max
        ok.sort(key=lambda t: (t[0]["bochner_n"] + t[0]["lasserre_T_max"],
                                t[0]["value"]),
                reverse=True)
        best_d, best_fn = ok[0]
        best[label] = best_d
        print(f"  {label}: best from {best_fn}")
        print(f"    value = {best_d['value']:.10f}  status = {best_d['status']}")
        print(f"    bochner_n = {best_d['bochner_n']}, T_max = {best_d['lasserre_T_max']}")
        print(f"    time = {best_d['time_seconds']:.1f}s")

    if "row5" in best and "row6" in best:
        v5 = best["row5"]["value"]
        v6 = best["row6"]["value"]
        mu_even = min(v5, v6)
        print(f"\n=== µ_even >= MIN(row5, row6) = {mu_even:.10f}")
        print(f"  vs unconditional {UNCONDITIONAL_BOUND}: "
              f"diff = {mu_even - UNCONDITIONAL_BOUND:+.6e}")

        # 1-D Path B ellipse coverage on p in [0.35, 0.45], h=q=0
        print(f"\n=== Path-B 1-D ellipse coverage on p in [0.35, 0.45] ===")
        intervals = []
        for label in ("row5", "row6"):
            r = best[label]
            ell = find_p_interval(
                r["value"], r["duals"], r["h"], r["p"], r["q1"], r["q2"],
                target=UNCONDITIONAL_BOUND,
            )
            intervals.append((label, ell))
            print(f"  {label} (p_c={r['p']:.4f}): "
                  f"[p_lo, p_hi] = [{ell['p_lo']:.6f}, {ell['p_hi']:.6f}], "
                  f"V_max={ell['V_max']:.7f}")

        # Coverage of [0.35, 0.45]
        # Take union of the two intervals.
        ints = [(e["p_lo"], e["p_hi"]) for _, e in intervals]
        ints = [(max(a, -1e9), min(b, 1e9)) for a, b in ints]
        ints.sort()
        # Merge
        merged = []
        for a, b in ints:
            if not merged or a > merged[-1][1] + 1e-12:
                merged.append([a, b])
            else:
                merged[-1][1] = max(merged[-1][1], b)
        # Check coverage of [0.35, 0.45]
        target_lo, target_hi = 0.35, 0.45
        gap = []
        cur = target_lo
        for a, b in merged:
            if a > cur:
                gap.append((cur, min(a, target_hi)))
                cur = b
            else:
                cur = max(cur, b)
            if cur >= target_hi: break
        if cur < target_hi:
            gap.append((cur, target_hi))
        gap = [(a, b) for a, b in gap if b > a + 1e-12]
        print(f"  merged intervals: {merged}")
        print(f"  uncovered subintervals of [0.35, 0.45]: {gap}")

        # Worst-case lower bound on the box [0.35, 0.45] for h=0, q in [-0.02, 0.02]:
        # for each p, take MAX over rows of obj(p), then MIN over p in [0.35, 0.45].
        n_grid = 5001
        p_grid = np.linspace(0.35, 0.45, n_grid)
        envelope_min = np.inf
        argmin_p = None
        argmin_label = None
        for p in p_grid:
            best_val = -np.inf
            best_label = None
            for label, _ in intervals:
                r = best[label]
                ell = find_p_interval(
                    r["value"], r["duals"], r["h"], r["p"], r["q1"], r["q2"],
                    target=UNCONDITIONAL_BOUND,
                )
                A, B, K0 = ell["A"], ell["B"], ell["K0"]
                v = A * p ** 2 + B * p + K0
                if v > best_val:
                    best_val = v
                    best_label = label
            if best_val < envelope_min:
                envelope_min = best_val
                argmin_p = p
                argmin_label = best_label
        print(f"  envelope MIN over p in [0.35, 0.45]: {envelope_min:.10f} "
              f"at p={argmin_p:.4f} (witness {argmin_label})")
        print(f"  vs unconditional {UNCONDITIONAL_BOUND}: "
              f"diff = {envelope_min - UNCONDITIONAL_BOUND:+.6e}")

        return {
            "rows": best,
            "mu_even_min_over_rows": mu_even,
            "diff_vs_unconditional": mu_even - UNCONDITIONAL_BOUND,
            "intervals": [{"row": label,
                           "p_lo": e["p_lo"], "p_hi": e["p_hi"],
                           "p_star": e["p_star"], "V_max": e["V_max"],
                           "A": e["A"], "B": e["B"], "K0": e["K0"]}
                          for label, e in intervals],
            "merged_intervals": merged,
            "uncovered_in_035_045": gap,
            "envelope_min_grid": float(envelope_min),
            "envelope_argmin_p": float(argmin_p),
            "envelope_argmin_row": argmin_label,
        }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_and_dump(row, N, T, R, bochner_n, lasserre_T_max, pickle_path):
    """Step A of split-solve: build cvxpy problem, reduce to standard form,
    pickle (data, chain, inv_data, handles_meta) for the solver step."""
    import cvxpy as cp
    spec = next(s for s in SYMMETRIC_ROWS if s[4] == row)
    h, p, qm, qp, lbl = spec
    print(f"Building {lbl}: h={h} p={p} q=[{qm},{qp}] N={N} T={T} R={R} "
          f"bochner_n={bochner_n} T_max={lasserre_T_max}", flush=True)
    t0 = time.time()
    Omega, cons, H = build_problem_even_with_handles(
        N, T, R, h, p, qm, qp, bochner_n, lasserre_T_max
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    print(f"  built in {time.time()-t0:.1f}s, reducing to standard form...", flush=True)
    t1 = time.time()
    data, chain, inv_data = prob.get_problem_data(solver=cp.CLARABEL)
    print(f"  reduced in {time.time()-t1:.1f}s, A shape {data['A'].shape}", flush=True)
    # Pickle the entire problem too, plus handles for dual extraction.
    blob = {
        "row": row, "h": h, "p": p, "q1": qm, "q2": qp,
        "N": N, "T": T, "R": R,
        "bochner_n": bochner_n, "lasserre_T_max": lasserre_T_max,
        "data": data, "chain": chain, "inv_data": inv_data,
        "prob": prob,
        # Handle indices into prob.constraints for dual readout
        "handle_keys": list(H.keys()),
        # Save the position of each handle constraint within prob.constraints
        # by id matching
        "handle_pos": {
            k: prob.constraints.index(v) if isinstance(v, cp.constraints.constraint.Constraint) else None
            for k, v in H.items()
            if not isinstance(v, cp.expressions.expression.Expression)
        },
        "v_idx": [id(c) for c in prob.constraints],
    }
    os.makedirs(os.path.dirname(pickle_path), exist_ok=True)
    t2 = time.time()
    with open(pickle_path, "wb") as f:
        pickle.dump(blob, f)
    print(f"  pickled to {pickle_path} in {time.time()-t2:.1f}s "
          f"({os.path.getsize(pickle_path)/1e6:.1f}MB)", flush=True)


def solve_dumped(pickle_path, out_json, time_limit=None):
    """Step B of split-solve: load pickle, run CLARABEL via prob._solve from
    the cached compiled problem, extract dual values, write JSON."""
    import cvxpy as cp
    print(f"Loading {pickle_path}...", flush=True)
    t0 = time.time()
    with open(pickle_path, "rb") as f:
        blob = pickle.load(f)
    print(f"  loaded in {time.time()-t0:.1f}s", flush=True)
    prob = blob["prob"]
    t1 = time.time()
    solver_kwargs = {"verbose": False}
    if time_limit is not None:
        solver_kwargs["time_limit"] = time_limit
    try:
        prob.solve(solver=cp.CLARABEL, **solver_kwargs)
    except Exception as e:
        res = {
            "row": blob["row"], "h": blob["h"], "p": blob["p"],
            "q1": blob["q1"], "q2": blob["q2"],
            "N": blob["N"], "T": blob["T"], "R": blob["R"],
            "bochner_n": blob["bochner_n"],
            "lasserre_T_max": blob["lasserre_T_max"],
            "value": None, "status": f"exception: {type(e).__name__}: {e}",
            "time_seconds": time.time() - t1,
        }
        os.makedirs(os.path.dirname(out_json), exist_ok=True)
        with open(out_json, "w") as fh:
            json.dump(res, fh, indent=2)
        print(f"FAILED: {e}", flush=True)
        return
    elapsed = time.time() - t1
    print(f"  solved in {elapsed:.1f}s. status={prob.status}, "
          f"value={prob.value}", flush=True)
    # Read dual values
    duals = {}
    handle_keys = ("con_53", "con_54", "con_512_pL", "con_512_pU",
                   "con_512_qL", "con_512_qU", "con_513")
    handle_pos = blob["handle_pos"]
    for key in handle_keys:
        if key in handle_pos and handle_pos[key] is not None:
            con = prob.constraints[handle_pos[key]]
            dv = con.dual_value
            duals[key] = float(dv) if dv is not None else 0.0
        else:
            duals[key] = 0.0
    res = {
        "row": blob["row"], "h": blob["h"], "p": blob["p"],
        "q1": blob["q1"], "q2": blob["q2"],
        "N": blob["N"], "T": blob["T"], "R": blob["R"],
        "bochner_n": blob["bochner_n"],
        "lasserre_T_max": blob["lasserre_T_max"],
        "value": float(prob.value) if prob.value is not None else None,
        "status": prob.status,
        "time_seconds": elapsed,
        "duals": duals,
    }
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"  written to {out_json}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="cmd")

    p_build = sub.add_parser("build")
    p_build.add_argument("row")
    p_build.add_argument("bochner_n", type=int)
    p_build.add_argument("lasserre_T_max", type=int)
    p_build.add_argument("pickle_path")
    p_build.add_argument("--N", type=int, default=10000)
    p_build.add_argument("--T", type=int, default=4000)
    p_build.add_argument("--R", type=int, default=10)

    p_solve = sub.add_parser("solve")
    p_solve.add_argument("pickle_path")
    p_solve.add_argument("out_json")
    p_solve.add_argument("--time_limit", type=float, default=None)

    p_one = sub.add_parser("oneshot")
    p_one.add_argument("row")
    p_one.add_argument("bochner_n", type=int)
    p_one.add_argument("lasserre_T_max", type=int)
    p_one.add_argument("out_json")
    p_one.add_argument("--N", type=int, default=10000)
    p_one.add_argument("--T", type=int, default=4000)
    p_one.add_argument("--R", type=int, default=10)
    p_one.add_argument("--time_limit", type=float, default=None)

    p_agg = sub.add_parser("aggregate")
    p_agg.add_argument("out_dir")

    args = parser.parse_args()

    if args.cmd == "build":
        build_and_dump(args.row, args.N, args.T, args.R, args.bochner_n,
                        args.lasserre_T_max, args.pickle_path)
    elif args.cmd == "solve":
        solve_dumped(args.pickle_path, args.out_json,
                     time_limit=args.time_limit)
    elif args.cmd == "oneshot":
        # Find the row spec
        spec = None
        for (h, p, qm, qp, lbl) in SYMMETRIC_ROWS:
            if lbl == args.row:
                spec = (h, p, qm, qp, lbl)
                break
        if spec is None:
            print(f"ERROR: row {args.row} not in SYMMETRIC_ROWS")
            sys.exit(1)
        h, p, qm, qp, lbl = spec
        print(f"Solving {lbl}: h={h}, p={p}, q in [{qm},{qp}], "
              f"N={args.N}, T={args.T}, R={args.R}, "
              f"bochner_n={args.bochner_n}, T_max={args.lasserre_T_max}",
              flush=True)
        res = solve_row(lbl, h, p, qm, qp, args.N, args.T, args.R,
                        args.bochner_n, args.lasserre_T_max,
                        time_limit=args.time_limit, verbose=False)
        os.makedirs(os.path.dirname(args.out_json), exist_ok=True)
        with open(args.out_json, "w") as fh:
            json.dump(res, fh, indent=2)
        print(f"Result: value={res['value']} status={res['status']} "
              f"time={res['time_seconds']:.1f}s", flush=True)
        print(f"Written to {args.out_json}")
    elif args.cmd == "aggregate":
        result = aggregate(args.out_dir)
        out_file = os.path.join(args.out_dir, "..", "symmetric_high_n.json")
        out_file = os.path.normpath(out_file)
        # Add metadata
        full_result = {
            "description": "CONDITIONAL bound on µ assuming f* even, pushed to "
                           "high (bochner_n, lasserre_T_max) under the halved-variable "
                           "even-f encoding.",
            "interpretation": "This is a CONDITIONAL bound. µ >= reported only "
                              "if the f-even conjecture holds (White §6).",
            "unconditional_bound": UNCONDITIONAL_BOUND,
            "infeasible_rows_under_even": ["row1", "row2", "row3", "row4", "row7"],
            "result": result,
        }
        with open(out_file, "w") as fh:
            json.dump(full_result, fh, indent=2)
        print(f"\nFull JSON written to {out_file}")
    else:
        parser.print_help()


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()
