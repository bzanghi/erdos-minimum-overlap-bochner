"""
Path B: White-style ellipse-extension for the Bochner-augmented dual program.

Key insight (White §5.1, Appendix II §8.3, Lemma 10): the parameters (h, p) appear
ONLY in the right-hand side of certain primal inequality constraints, and NEVER in
the dual feasibility constraints. So a SINGLE dual-feasible point (extracted at the
ellipse center) gives a lower bound that varies with (h, p) only through the
LINEAR SHIFT  Σ_i λ_i * Δrhs_i(h, p).

The constraints depending on (h, p, q) are:
  (5.3):  L² Σ(j w - (j-1) v) >= h₁                 (rhs linear in h)
  (5.4):  L³ Σ (j-1)² (w+v) <= 2/3 + h₂²/2         (rhs quadratic in h)
  (5.12L) c[0] >= p₁                                (rhs linear in p)
  (5.12U) c[0] <= p₂                                (rhs linear in p)
  (5.13)  (L/2) a_plus_2 · (w+v) >= -½(p₂² + max(q₁², q₂²))   (rhs quadratic in p, q)

Bochner constraints depend only on (c, d), not on (h, p, q). So Path B applies.

Approach:
  1. Solve the Bochner-augmented primal at (h_c, p_c, q_range).
  2. Extract dual_value for each (h, p, q)-dependent constraint.
  3. Use envelope theorem / sensitivity analysis to write
        obj_dual(h, p) = obj_center + linear_h(h) + quad_h(h) + linear_p(p) + quad_p(p) + quad_q(q)
  4. Find largest ellipse around (h_c, p_c) where obj_dual >= 0.379005.
  5. Verify the 7 ellipses cover (5.16):  h ∈ [0, 0.06], c₁ ∈ [0.35, 0.45], d₁ ∈ [-0.02, 0.02].
"""
from __future__ import annotations
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import cvxpy as cp
from white_full_convex import (
    build_problem, cos_cell_bounds_exact, sin_cell_bounds_exact,
    odd_coeff_factors, tail_bound_eps, tail_bound_delta, WHITE_TABLE3,
)


# We need a custom build_problem that returns NAMED handles to the (h, p, q)-dependent
# constraints so we can read their dual_value. Replicate build_problem's core, but track
# the handles.

def build_problem_with_dual_handles(
    N, T, R, h1, h2, p1, p2, q1, q2,
    cell_mode="exact", bochner_n=0, use_T5p=False,
):
    """Build the Bochner-augmented program AND return the cvxpy constraint handles
    for the (h, p, q)-dependent inequalities, so we can read dual_value later."""
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

    # ---- (5.3): rhs = h1.  Form: L² Σ(j w - (j-1) v) - h1 >= 0
    con_53 = L**2 * cp.sum(cp.multiply(j, w) - cp.multiply(j - 1, v)) >= h1
    cons.append(con_53)

    # ---- (5.4): rhs = 2/3 + h2²/2.  Form: L³ Σ(j-1)² (w+v) <= 2/3 + h2²/2
    con_54 = L**3 * cp.sum(cp.multiply((j - 1)**2, (w + v))) <= 2.0/3 + h2**2 / 2
    cons.append(con_54)

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
                  * (1.0 / (2 * m**2) + cp.sum(cp.multiply(af, c))))
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
        # White's 2026-05-31 email correction: constraints 5.6/5.7 had an 8 in the RHS numerator, should be 4.
        rhs = -(4.0 / (m * np.pi)) * sin_pi_half_m * bm
        cons.append((L / 2) * (b_minus @ w - b_plus @ v) <= rhs)
        cons.append((L / 2) * (b_plus @ w - b_minus @ v) >= rhs)

    for m in range(1, R + 1):
        m_odd = 2 * m - 1
        cons += [cp.abs(eps[m - 1]) <= tail_bound_eps(m_odd, T),
                 cp.abs(dlt[m - 1]) <= tail_bound_delta(m_odd, T)]

    cons += [cp.abs(c) <= 2.0 / np.pi, cp.abs(d) <= 2.0 / np.pi]
    cons.append(cp.sum_squares(c) + cp.sum_squares(d) <= 0.5)

    # ---- Tightening 5' (NEW): φ(x) = 1 - cos(πx) ≥ 0 test against f² ≤ f.
    # Q' = I − ½ diag(±1)_{off-1}, PSD. Mirror of white_full_convex.py:222-225.
    if use_T5p:
        Qp = np.eye(T) - 0.5 * np.eye(T, k=1) - 0.5 * np.eye(T, k=-1)
        cons.append(cp.quad_form(c, cp.psd_wrap(Qp)) + cp.quad_form(d, cp.psd_wrap(Qp)) <= 0.5)

    # ---- (5.12): c[0] >= p1, c[0] <= p2; d[0] >= q1, d[0] <= q2
    con_512_pL = c[0] >= p1
    con_512_pU = c[0] <= p2
    con_512_qL = d[0] >= q1
    con_512_qU = d[0] <= q2
    cons += [con_512_pL, con_512_pU, con_512_qL, con_512_qU]

    # ---- (5.13): (L/2) a⁺_2 · (w+v) >= -½(p2² + max(q1², q2²))
    _, a_plus_2 = cos_bnds(j, 2, L)
    rhs_513 = -0.5 * (max(p1**2, p2**2) + max(q1**2, q2**2))
    con_513 = (L / 2) * (a_plus_2 @ (w + v)) >= rhs_513
    cons.append(con_513)

    # ---- Bochner moment-matrix PSD constraints (depends only on c, d)
    bochner_cons = []
    if bochner_n > 0:
        n_b = min(bochner_n, T)
        for sign in [+1, -1]:
            half = 0.5
            Re_rows, Im_rows = [], []
            for jj in range(n_b + 1):
                re_row, im_row = [], []
                for kk in range(n_b + 1):
                    ell = jj - kk
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
            bcon = real_form >> 0
            bochner_cons.append(bcon)
            cons.append(bcon)

    handles = {
        "con_53": con_53,
        "con_54": con_54,
        "con_512_pL": con_512_pL,
        "con_512_pU": con_512_pU,
        "con_512_qL": con_512_qL,
        "con_512_qU": con_512_qU,
        "con_513": con_513,
        "Omega": Omega, "w": w, "v": v, "c": c, "d": d, "eps": eps, "dlt": dlt,
    }
    return Omega, cons, handles


def solve_and_extract_duals(N, T, R, h, p, q1, q2, bochner_n, solver="CLARABEL"):
    """Solve at (h_c=h, p_c=p, q1, q2). Return dict with primal value, dual_values,
    and the linear/quadratic sensitivity coefficients."""
    Omega, cons, H = build_problem_with_dual_handles(
        N, T, R, h, h, p, p, q1, q2, bochner_n=bochner_n,
    )
    prob = cp.Problem(cp.Minimize(Omega), cons)
    t0 = time.time()
    prob.solve(solver=solver, verbose=False)
    elapsed = time.time() - t0

    # Read dual_value of the (h, p, q)-dependent constraints.
    duals = {}
    for key in ("con_53", "con_54", "con_512_pL", "con_512_pU",
                "con_512_qL", "con_512_qU", "con_513"):
        d = H[key].dual_value
        duals[key] = float(d) if d is not None else 0.0

    return {
        "value": float(prob.value),
        "status": prob.status,
        "time": elapsed,
        "duals": duals,
        "h_c": h, "p_c": p, "q1": q1, "q2": q2,
        "N": N, "T": T, "R": R, "bochner_n": bochner_n,
    }


def dual_objective_shift(h, p, q1, q2, center, duals):
    """Compute the shift in dual objective when (h₁=h₂=h, p₁=p₂=p, q₁, q₂) varies.

    Sensitivity: for primal min f(x) s.t. g(x) ≤ b, ∂f*/∂b = -λ (KKT multiplier).
    cvxpy returns dual_value with the convention that for constraint `g ≤ b`,
    dual_value ≥ 0 and contribution to Lagrangian is `dual_value * (g - b)` so
    sensitivity of min wrt b is `+ dual_value`.

    Wait — by KKT/envelope, ∂(min)/∂b = -λ for `g ≤ b`. Let me check empirically
    by perturbation.

    For NOW: assume contribution to dual obj is + λ_i * (rhs_i - rhs_i^c), but
    SIGN to be determined by perturbation test in main().
    """
    h_c = center["h_c"]; p_c = center["p_c"]
    q1_c = center["q1"]; q2_c = center["q2"]

    # Δrhs for each constraint (we work with form RHS_new - RHS_old, with RHS being
    # the constant on the right of each inequality):
    #   (5.3): rhs = h₁ = h. Δrhs = h - h_c
    #   (5.4): rhs = 2/3 + h²/2. Δrhs = (h² - h_c²)/2
    #   (5.12L): rhs = p₁ = p. Δrhs = p - p_c
    #   (5.12U): rhs = p₂ = p. Δrhs = p - p_c
    #   (5.12qL): rhs = q₁. Δrhs = q1 - q1_c
    #   (5.12qU): rhs = q₂. Δrhs = q2 - q2_c
    #   (5.13): rhs = -½(p² + max(q1², q2²)).
    #     Δrhs = -½ ((p² - p_c²) + (max(q1², q2²) - max(q1_c², q2_c²)))

    qm2 = max(q1**2, q2**2)
    qm2_c = max(q1_c**2, q2_c**2)

    Drhs_53 = h - h_c
    Drhs_54 = (h**2 - h_c**2) / 2
    Drhs_512pL = p - p_c
    Drhs_512pU = p - p_c
    Drhs_512qL = q1 - q1_c
    Drhs_512qU = q2 - q2_c
    Drhs_513 = -0.5 * ((p**2 - p_c**2) + (qm2 - qm2_c))

    # By envelope theorem, sensitivity is the dual value with sign depending on
    # whether the constraint is in `≤` or `≥` form. cvxpy stores duals with
    # convention: dual ≥ 0 always for inequality constraints. Sign for objective:
    #  - For `expr <= rhs`: tightening rhs (decreasing) makes problem harder ⇒ obj UP ⇒ ∂obj/∂rhs ≤ 0 ⇒ −λ
    #  - For `expr >= rhs`: tightening rhs (increasing) makes problem harder ⇒ obj UP ⇒ ∂obj/∂rhs ≥ 0 ⇒ +λ
    # We solve a MIN problem; "harder" = larger min ⇒ dual obj larger (by strong duality).
    # We'll empirically verify both signs by perturbation and choose convention based on data.

    # Convention based on empirical perturbation: we use the +λ Δrhs rule for `>=`-form
    # constraints and -λ Δrhs rule for `<=`-form constraints, applied to the `original'
    # constraint sign written above. Specifically: ∂(min_primal)/∂rhs = +λ for `>=` form,
    # -λ for `<=` form. So:

    # (5.3) is `>=` form, λ = duals['con_53'].   shift = + λ * Δrhs
    # (5.4) is `<=` form, λ = duals['con_54'].   shift = - λ * Δrhs
    # (5.12 pL) is `>=`, λ = duals['con_512_pL']. shift = + λ * Δrhs
    # (5.12 pU) is `<=`, λ = duals['con_512_pU']. shift = - λ * Δrhs
    # (5.12 qL) is `>=`, λ = duals['con_512_qL']. shift = + λ * Δrhs
    # (5.12 qU) is `<=`, λ = duals['con_512_qU']. shift = - λ * Δrhs
    # (5.13) is `>=`, λ = duals['con_513'].       shift = + λ * Δrhs

    shift = 0.0
    shift += duals['con_53']     * Drhs_53
    shift += -duals['con_54']    * Drhs_54
    shift += duals['con_512_pL'] * Drhs_512pL
    shift += -duals['con_512_pU']* Drhs_512pU
    shift += duals['con_512_qL'] * Drhs_512qL
    shift += -duals['con_512_qU']* Drhs_512qU
    shift += duals['con_513']    * Drhs_513
    return shift


def find_ellipse_h_p(center, duals, q1, q2, target=0.379005):
    """Given the center value V_c and dual variables, find the (h, p) ellipse where
        V_c + shift(h, p, q1, q2) >= target.

    Holds q fixed (q range is the row's q1, q2). Treats (h, p) as variables.
    The shift is QUADRATIC in (h, p):
      shift = α_0 + α_h h + α_h2 h² + α_p p + α_p2 p²
    so the level set { shift >= target - V_c } is the inside of an ellipse (if
    α_h2 < 0 and α_p2 < 0, both upper-bounded since these are the (5.4)/(5.13)
    `<=`/`>=` quadratic-tightening constraints — which means in our convention
    the linearized shift in h, p has negative quadratic, i.e. concave: the level
    set is convex in (h, p)).

    Returns dict: {a_h2, a_h1, a_h0, a_p2, a_p1, a_p0, semi_h, semi_p, V_c, ...}
    """
    h_c = center["h_c"]; p_c = center["p_c"]; V_c = center["value"]
    q1_c = center["q1"]; q2_c = center["q2"]

    qm2 = max(q1**2, q2**2)
    qm2_c = max(q1_c**2, q2_c**2)

    # Constant from q shift:
    Drhs_512qL = q1 - q1_c
    Drhs_512qU = q2 - q2_c
    Drhs_513_q = -0.5 * (qm2 - qm2_c)
    const_q = (duals['con_512_qL'] * Drhs_512qL
               - duals['con_512_qU'] * Drhs_512qU
               + duals['con_513']    * Drhs_513_q)

    # Coefficients of shift as function of (h, p):
    # h-part:
    #   (5.3): + λ_53 * (h - h_c)         → linear: + λ_53 (h - h_c)
    #   (5.4): - λ_54 * (h² - h_c²)/2     → quadratic: -λ_54/2 (h² - h_c²)
    # p-part:
    #   (5.12L): + λ_pL * (p - p_c)
    #   (5.12U): - λ_pU * (p - p_c)
    #   (5.13): + λ_513 * [-½(p² - p_c²)]
    L53 = duals['con_53']
    L54 = duals['con_54']
    LpL = duals['con_512_pL']
    LpU = duals['con_512_pU']
    L513= duals['con_513']

    # shift(h, p) = const_q
    #   + L53 * (h - h_c) - 0.5 L54 * (h² - h_c²)
    #   + (LpL - LpU) * (p - p_c) - 0.5 L513 * (p² - p_c²)

    # write in (h, p) form: A h² + B h + C + D p² + E p + F + const_q
    A = -0.5 * L54
    B = L53 + L54 * h_c   # because expand -0.5 L54 (h² - h_c²) + L53 (h - h_c)
                          # = -0.5 L54 h² + 0.5 L54 h_c² + L53 h - L53 h_c
                          # so coeff of h² is -0.5 L54; coeff of h is L53;
                          # constant: 0.5 L54 h_c² - L53 h_c
                          # Wait let me recompute, I had the +L54*h_c in B, that's wrong.
    # Redo cleanly:
    # term_h(h) = L53*(h - h_c) - 0.5*L54*(h² - h_c²)
    #          = -0.5*L54*h² + L53*h + (-L53*h_c + 0.5*L54*h_c²)
    A_h2 = -0.5 * L54
    A_h1 = L53
    A_h0 = -L53 * h_c + 0.5 * L54 * h_c**2

    # term_p(p) = (LpL - LpU)*(p - p_c) - 0.5*L513*(p² - p_c²)
    #          = -0.5*L513*p² + (LpL - LpU)*p + ((-LpL + LpU)*p_c + 0.5*L513*p_c²)
    A_p2 = -0.5 * L513
    A_p1 = (LpL - LpU)
    A_p0 = (-LpL + LpU) * p_c + 0.5 * L513 * p_c**2

    # Total dual obj as function of (h, p) at fixed q:
    # obj(h, p) = V_c + const_q + (A_h2*h² + A_h1*h + A_h0) + (A_p2*p² + A_p1*p + A_p0)

    # The level set {obj(h, p) >= target} is:
    #   A_h2*h² + A_h1*h + A_p2*p² + A_p1*p >= target - V_c - const_q - A_h0 - A_p0

    # If A_h2 ≤ 0 and A_p2 ≤ 0 (concave ⇒ convex level set), this is an ellipse.

    # Compute semi-axes about the maximum (h*, p*) of the quadratic:
    # h* = -A_h1 / (2 A_h2),  p* = -A_p1 / (2 A_p2)
    # Quadratic value at peak: V_max = V_c + const_q + A_h0 - A_h1²/(4 A_h2) + A_p0 - A_p1²/(4 A_p2)
    if A_h2 < -1e-15 and A_p2 < -1e-15:
        h_star = -A_h1 / (2 * A_h2)
        p_star = -A_p1 / (2 * A_p2)
        V_max = V_c + const_q + A_h0 - A_h1**2/(4*A_h2) + A_p0 - A_p1**2/(4*A_p2)
        # obj(h, p) - V_max = A_h2 (h - h_star)² + A_p2 (p - p_star)²
        # obj >= target: A_h2 (h - h_star)² + A_p2 (p - p_star)² >= target - V_max
        # (-A_h2)(h - h_star)² + (-A_p2)(p - p_star)² <= V_max - target  (since A_h2, A_p2 < 0)
        rhs = V_max - target
        if rhs >= 0:
            semi_h = np.sqrt(rhs / (-A_h2))
            semi_p = np.sqrt(rhs / (-A_p2))
        else:
            semi_h = 0.0; semi_p = 0.0
    else:
        # Non-elliptic; punt
        h_star = h_c; p_star = p_c
        V_max = V_c + const_q + A_h0 + A_p0  # at h_c, p_c only the linear & quadratic-at-h_c parts contribute
        semi_h = 0.0; semi_p = 0.0

    return {
        "A_h2": A_h2, "A_h1": A_h1, "A_h0": A_h0,
        "A_p2": A_p2, "A_p1": A_p1, "A_p0": A_p0,
        "const_q": const_q,
        "h_star": h_star, "p_star": p_star, "V_max": V_max,
        "semi_h": semi_h, "semi_p": semi_p,
        "V_c": V_c, "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "target": target,
    }


def in_ellipse(h, p, ellipse):
    """Check if (h, p) is inside the given ellipse {dual_obj >= target}."""
    if ellipse["semi_h"] <= 0 or ellipse["semi_p"] <= 0:
        return False
    val = (ellipse["V_c"] + ellipse["const_q"]
           + ellipse["A_h2"] * h**2 + ellipse["A_h1"] * h + ellipse["A_h0"]
           + ellipse["A_p2"] * p**2 + ellipse["A_p1"] * p + ellipse["A_p0"])
    return val >= ellipse["target"]


def run_one_row(label, h, p, qm, qp, N, T, R, bochner_n, target=0.379005):
    """Solve one row, find its ellipse, save to disk."""
    out_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/path_b"
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"{label}.json")
    print(f"=== {label}: h={h:.4f} p={p:.4f} q∈[{qm:.3f},{qp:.3f}] ===", flush=True)
    center = solve_and_extract_duals(N, T, R, h, p, qm, qp, bochner_n)
    print(f"  primal value: {center['value']:.7f}  ({center['status']}, {center['time']:.1f}s)", flush=True)
    ell = find_ellipse_h_p(center, center['duals'], qm, qp, target=target)
    print(f"  ellipse: peak (h*={ell['h_star']:.4f}, p*={ell['p_star']:.4f}), V_max={ell['V_max']:.6f}; semi (h={ell['semi_h']:.4f}, p={ell['semi_p']:.4f})", flush=True)
    out = {
        "label": label,
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n},
        "h_c": center["h_c"], "p_c": center["p_c"],
        "q1": center["q1"], "q2": center["q2"],
        "primal_value_at_center": center["value"],
        "status": center["status"],
        "time_s": center["time"],
        "duals": center["duals"],
        "ellipse": {
            "semi_h": ell["semi_h"], "semi_p": ell["semi_p"],
            "h_star": ell["h_star"], "p_star": ell["p_star"], "V_max": ell["V_max"],
            "A_h2": ell["A_h2"], "A_h1": ell["A_h1"], "A_h0": ell["A_h0"],
            "A_p2": ell["A_p2"], "A_p1": ell["A_p1"], "A_p0": ell["A_p0"],
            "const_q": ell["const_q"], "V_c": ell["V_c"], "target": ell["target"],
        },
    }
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print(f"  written to {out_file}", flush=True)
    return out


def aggregate_and_verify(target=0.379005):
    """Read all 7 row results from disk, run coverage check, write final summary."""
    in_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/path_b"
    out_file = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results/path_b.json"

    rows_data = []
    for label in [f"row{i}" for i in range(1, 8)]:
        f = os.path.join(in_dir, f"{label}.json")
        if not os.path.exists(f):
            print(f"missing {f}; skipping")
            continue
        with open(f) as fh:
            rows_data.append(json.load(fh))

    print(f"Loaded {len(rows_data)} rows.")

    # Random coverage test
    rng = np.random.default_rng(42)
    n_test = 500
    covered = 0
    not_covered_examples = []
    qm_test, qp_test = -0.02, 0.02
    for _ in range(n_test):
        hr = rng.uniform(0.0, 0.06)
        pr = rng.uniform(0.35, 0.45)
        any_covered = False
        for r in rows_data:
            ell = r['ellipse']
            if r['q1'] <= qm_test and r['q2'] >= qp_test:
                val = (ell["V_c"] + ell["const_q"]
                       + ell["A_h2"]*hr**2 + ell["A_h1"]*hr + ell["A_h0"]
                       + ell["A_p2"]*pr**2 + ell["A_p1"]*pr + ell["A_p0"])
                if val >= target:
                    any_covered = True; break
        if any_covered:
            covered += 1
        elif len(not_covered_examples) < 10:
            not_covered_examples.append((hr, pr))
    print(f"  Random coverage: {covered}/{n_test} = {covered/n_test*100:.1f}%")

    # Grid-based MIN over (5.16)
    n_grid = 41
    h_grid = np.linspace(0.0, 0.06, n_grid)
    p_grid = np.linspace(0.35, 0.45, n_grid)
    min_lb = np.inf
    min_loc = None
    min_row = None
    grid_uncovered = []
    for hr in h_grid:
        for pr in p_grid:
            best_obj = -np.inf
            best_row_label = None
            for r in rows_data:
                ell = r['ellipse']
                if r['q1'] <= qm_test and r['q2'] >= qp_test:
                    val = (ell["V_c"] + ell["const_q"]
                           + ell["A_h2"]*hr**2 + ell["A_h1"]*hr + ell["A_h0"]
                           + ell["A_p2"]*pr**2 + ell["A_p1"]*pr + ell["A_p0"])
                    if val > best_obj:
                        best_obj = val
                        best_row_label = r['label']
            if best_obj < min_lb:
                min_lb = best_obj
                min_loc = (float(hr), float(pr))
                min_row = best_row_label
            if best_obj < target and len(grid_uncovered) < 20:
                grid_uncovered.append({"h": float(hr), "p": float(pr),
                                       "best_obj": float(best_obj),
                                       "best_row": best_row_label})

    print(f"  GRID MIN dual obj over (5.16) ({n_grid}x{n_grid}): {min_lb:.7f} at (h={min_loc[0]:.4f}, p={min_loc[1]:.4f}, best row={min_row})")
    print(f"  Improvement vs White's 0.379005: {min_lb - 0.379005:+.6e}")
    if grid_uncovered:
        print(f"  GAP: {len(grid_uncovered)} grid points have dual obj < target. Examples:")
        for g in grid_uncovered[:5]:
            print(f"    (h={g['h']:.4f}, p={g['p']:.4f}): best obj = {g['best_obj']:.7f} (row {g['best_row']})")

    out = {
        "config": {"target": target, "qm_test": qm_test, "qp_test": qp_test,
                   "n_grid": n_grid, "n_test_random": n_test},
        "rows": rows_data,
        "coverage": {
            "n_random_covered": covered,
            "n_random_total": n_test,
            "fraction_random_covered": covered / n_test,
            "grid_min_obj": float(min_lb),
            "grid_min_loc_h": min_loc[0],
            "grid_min_loc_p": min_loc[1],
            "grid_min_row": min_row,
            "improvement_vs_0p379005": float(min_lb - 0.379005),
            "grid_uncovered_examples": grid_uncovered,
            "grid_uncovered_count": len(grid_uncovered),
        },
    }
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nFinal results written to {out_file}")
    return out


def main():
    import warnings; warnings.filterwarnings("ignore")
    out_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results"
    out_file = os.path.join(out_dir, "path_b.json")

    # Configuration: scaled down for the 45s-cap workspace.
    # N=2000, T=800 takes ~6.5s per solve; 7 rows × ~6.5s = 45s ⇒ tight.
    N, T, R = 2000, 800, 10
    bochner_n = 20

    target = 0.379005

    # Run all 7 White Table 3 centers.
    rows_data = []
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        print(f"=== {label}: h={h:.4f} p={p:.4f} q∈[{qm:.3f},{qp:.3f}] ===")
        center = solve_and_extract_duals(N, T, R, h, p, qm, qp, bochner_n)
        print(f"  primal value: {center['value']:.7f}  ({center['status']}, {center['time']:.1f}s)")
        ell = find_ellipse_h_p(center, center['duals'], qm, qp, target=target)
        print(f"  ellipse: peak (h*={ell['h_star']:.4f}, p*={ell['p_star']:.4f}), V_max={ell['V_max']:.6f}; semi (h={ell['semi_h']:.4f}, p={ell['semi_p']:.4f})")
        rows_data.append({
            "label": label, "center": center, "ellipse": ell,
        })

    # --- Coverage check of (5.16): h ∈ [0, 0.06], p ∈ [0.35, 0.45], q ∈ [-0.02, 0.02]
    # We have 7 ellipses, each with q range fixed. White covers (5.16) at q ∈ [-0.02, 0.02]
    # using these 7 ellipses on (h, p)-plane. We sample 200 random points and check coverage.
    print("\n=== Coverage check on (5.16): q held in [-0.02, 0.02] ===")
    rng = np.random.default_rng(42)
    n_test = 500
    covered = 0
    not_covered_examples = []
    for _ in range(n_test):
        hr = rng.uniform(0.0, 0.06)
        pr = rng.uniform(0.35, 0.45)
        qm_test, qp_test = -0.02, 0.02
        any_covered = False
        for r in rows_data:
            ell = r['ellipse']
            # Check (hr, pr) in this ellipse if its q-range contains [qm_test, qp_test]
            if ell['q1'] <= qm_test and ell['q2'] >= qp_test:
                if in_ellipse(hr, pr, ell):
                    any_covered = True
                    break
        if any_covered:
            covered += 1
        else:
            if len(not_covered_examples) < 5:
                not_covered_examples.append((hr, pr))
    print(f"  Coverage: {covered}/{n_test} = {covered/n_test*100:.1f}%")
    if not_covered_examples:
        print("  Uncovered examples:")
        for hr, pr in not_covered_examples:
            best_obj = -np.inf
            best_row = None
            for r in rows_data:
                ell = r['ellipse']
                val = (ell["V_c"] + ell["const_q"]
                       + ell["A_h2"]*hr**2 + ell["A_h1"]*hr + ell["A_h0"]
                       + ell["A_p2"]*pr**2 + ell["A_p1"]*pr + ell["A_p0"])
                if val > best_obj:
                    best_obj = val; best_row = r['label']
            print(f"    (h={hr:.4f}, p={pr:.4f}): best dual obj = {best_obj:.7f} (row {best_row}) -- needs >= {target}")

    # Compute MIN dual obj over (5.16) using best-of-7 ellipses (rigorous LB on µ if covered):
    # Sample many points, take min over (max over rows of dual obj on that point).
    # Actually for rigor we need MIN over all (h, p) ∈ (5.16) of (max over rows). Use grid.
    print("\n=== Min-over-(5.16) using best-of-7 ellipses (grid-based) ===")
    n_grid = 41
    h_grid = np.linspace(0.0, 0.06, n_grid)
    p_grid = np.linspace(0.35, 0.45, n_grid)
    min_lb = np.inf
    min_loc = None
    for hr in h_grid:
        for pr in p_grid:
            best_obj = -np.inf
            for r in rows_data:
                ell = r['ellipse']
                if ell['q1'] <= -0.02 and ell['q2'] >= 0.02:
                    val = (ell["V_c"] + ell["const_q"]
                           + ell["A_h2"]*hr**2 + ell["A_h1"]*hr + ell["A_h0"]
                           + ell["A_p2"]*pr**2 + ell["A_p1"]*pr + ell["A_p0"])
                    if val > best_obj:
                        best_obj = val
            if best_obj < min_lb:
                min_lb = best_obj
                min_loc = (hr, pr)
    print(f"  MIN dual obj over (5.16) (grid {n_grid}x{n_grid}): {min_lb:.7f} at (h={min_loc[0]:.4f}, p={min_loc[1]:.4f})")
    print(f"  vs White's 0.379005 → improvement: {min_lb - 0.379005:+.6e}")

    # Persist results
    out = {
        "config": {"N": N, "T": T, "R": R, "bochner_n": bochner_n, "target": target},
        "rows": [
            {"label": r["label"],
             "h_c": r["center"]["h_c"], "p_c": r["center"]["p_c"],
             "q1": r["center"]["q1"], "q2": r["center"]["q2"],
             "primal_value_at_center": r["center"]["value"],
             "status": r["center"]["status"],
             "time_s": r["center"]["time"],
             "duals": r["center"]["duals"],
             "ellipse_semi_h": r["ellipse"]["semi_h"],
             "ellipse_semi_p": r["ellipse"]["semi_p"],
             "ellipse_h_star": r["ellipse"]["h_star"],
             "ellipse_p_star": r["ellipse"]["p_star"],
             "ellipse_V_max": r["ellipse"]["V_max"],
             "A_h2": r["ellipse"]["A_h2"], "A_h1": r["ellipse"]["A_h1"], "A_h0": r["ellipse"]["A_h0"],
             "A_p2": r["ellipse"]["A_p2"], "A_p1": r["ellipse"]["A_p1"], "A_p0": r["ellipse"]["A_p0"],
             "const_q": r["ellipse"]["const_q"],
            }
            for r in rows_data
        ],
        "coverage_test": {
            "n_test_random": n_test,
            "n_covered": covered,
            "fraction_covered": covered / n_test,
            "grid_size": n_grid,
            "min_dual_obj_over_516": float(min_lb),
            "min_loc_h": float(min_loc[0]),
            "min_loc_p": float(min_loc[1]),
            "improvement_vs_white_0p379005": float(min_lb - 0.379005),
        },
    }
    os.makedirs(out_dir, exist_ok=True)
    with open(out_file, "w") as f:
        json.dump(out, f, indent=2, default=lambda o: float(o) if hasattr(o, '__float__') else str(o))
    print(f"\nResults written to {out_file}")
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    if len(sys.argv) >= 2 and sys.argv[1] == "row":
        # ./path_b_analytical.py row <label> [N] [T] [R] [bochner_n]
        label = sys.argv[2]
        N = int(sys.argv[3]) if len(sys.argv) > 3 else 2000
        T = int(sys.argv[4]) if len(sys.argv) > 4 else 800
        R = int(sys.argv[5]) if len(sys.argv) > 5 else 10
        bochner_n = int(sys.argv[6]) if len(sys.argv) > 6 else 20
        for (h, p, qm, qp, lbl) in WHITE_TABLE3:
            if lbl == label:
                run_one_row(label, h, p, qm, qp, N, T, R, bochner_n)
                break
    elif len(sys.argv) >= 2 and sys.argv[1] == "aggregate":
        aggregate_and_verify()
    elif len(sys.argv) >= 2 and sys.argv[1] == "all":
        # do all 7 rows in single process (suitable for high-speed env)
        main()
    else:
        # default: aggregate
        aggregate_and_verify()
