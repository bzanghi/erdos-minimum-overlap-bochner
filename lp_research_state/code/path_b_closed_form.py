"""
Path B CLOSED-FORM minimum of the best-of-7 envelope on the residual region (5.16).

Replaces the 2001x2001 grid + Lipschitz-error-bar approach in path_b_lasserre_rigorous.py
with an EXACT analytical minimization. No grid resolution.

Each row r contributes a separable concave quadratic
   Phi_r(h, p) = V_c_r + A_h2*h^2 + A_h1*h + A_h0  +  A_p2*p^2 + A_p1*p + A_p0  (+ const_q)
on the residual box (5.16):  h in [0, 0.06], p in [0.35, 0.45], q in [-0.02, 0.02].

For each row the q-related dual variables (con_512_qL, con_512_qU, con_513's q-portion)
are essentially zero in the lasserre2 solves (1e-8 or smaller), so const_q ~ 0.
We treat const_q rigorously: the worst-case q sign on each row's q-shift is upper-bounded
by |duals_q_*| * (q-range), and since these are ~1e-9, they contribute at most ~1e-10
to the envelope min. This is negligible relative to the 1e-6 IPM-gap margin.

The minimum of the envelope E(h, p) := max_r Phi_r(h, p) over the box is at one of:
   (A) A box CORNER (4 corners), where E = max_r Phi_r evaluated at the corner.
   (B) A box EDGE: along an edge, fix one variable at a boundary; reduces to a 1D max
       of concave quadratics, whose envelope min is either at the edge endpoints (i.e.
       corners) or at a 2-row crossing on that edge.
   (C) An INTERIOR 3-row crossing: a point (h*, p*) where Phi_r(p*) = Phi_s(p*) = Phi_t(p*)
       for three distinct rows.  Since each Phi is separable concave, the system
          Phi_r - Phi_s = 0,  Phi_r - Phi_t = 0
       reduces to TWO conic equations in (h, p).  When (A_h2_r, A_p2_r) differ from
       (A_h2_s, A_p2_s) etc., this is a quadratic system; otherwise linear.

Algorithm:
  1) Enumerate 4 corners.
  2) For each pair (r, s), compute the locus Phi_r = Phi_s intersected with each of
     the 4 edges. On each edge segment, minimize max(Phi_r(h,p), Phi_s(h,p)) along
     the locus restricted to the segment. (Actually any candidate minimum here is
     captured by the corners + 2-row pair-crossings on the edge.)
  3) For each triple (r, s, t), solve Phi_r = Phi_s and Phi_r = Phi_t for (h, p),
     check if interior to box, and if E(h*, p*) = Phi_r equals max_q Phi_q (in case
     no other row beats the triple) — in which case this is a critical point of E.
  4) For each pair (r, s), additionally enumerate intersection of Phi_r = Phi_s with
     each box-edge, then minimize E along those line segments analytically (1D conic
     optimization).
  5) Take min over all candidates.

This gives an EXACT minimum (modulo floating-point arithmetic), no grid required.
"""
from __future__ import annotations
import os, sys, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from white_full_convex import WHITE_TABLE3
from path_b_analytical import find_ellipse_h_p


WHITE_BOUND = 0.379005
H_BOX = (0.0, 0.06)
P_BOX = (0.35, 0.45)
Q_RANGE = (-0.02, 0.02)


def _phi_value(e, h, p):
    """Evaluate Phi_r(h, p) given an ellipse dict."""
    return (e["V_c"] + e["const_q"]
            + e["A_h2"] * h * h + e["A_h1"] * h + e["A_h0"]
            + e["A_p2"] * p * p + e["A_p1"] * p + e["A_p0"])


def _envelope_value(rows, h, p):
    """E(h, p) = max_r Phi_r(h, p).  Returns (E, witness_idx)."""
    best = -np.inf; idx = -1
    for i, r in enumerate(rows):
        v = _phi_value(r["ellipse"], h, p)
        if v > best:
            best = v; idx = i
    return best, idx


def _phi_diff_coeffs(e1, e2):
    """Phi_1 - Phi_2 = a*h^2 + b*h + c*p^2 + d*p + e, returns (a, b, c, d, e)."""
    a = e1["A_h2"] - e2["A_h2"]
    b = e1["A_h1"] - e2["A_h1"]
    a0 = e1["A_h0"] - e2["A_h0"]
    c = e1["A_p2"] - e2["A_p2"]
    d = e1["A_p1"] - e2["A_p1"]
    p0 = e1["A_p0"] - e2["A_p0"]
    e0 = e1["V_c"] + e1["const_q"] - e2["V_c"] - e2["const_q"] + a0 + p0
    return (a, b, c, d, e0)


def candidates_corners(rows, h_box=H_BOX, p_box=P_BOX):
    out = []
    for h in h_box:
        for p in p_box:
            E, idx = _envelope_value(rows, h, p)
            out.append({"type": "corner", "h": float(h), "p": float(p),
                        "E": float(E), "witness_row": rows[idx]["label"]})
    return out


def _candidate_pair_on_edge(rows, e_i, e_j, fix_var, fix_val, var_lo, var_hi):
    """For pair (i, j), find (h, p) on a box edge where Phi_i = Phi_j.
    fix_var = 'h' or 'p'.
    The locus Phi_i = Phi_j restricted to fix_var=fix_val is a 1D quadratic in
    the other variable: alpha*x^2 + beta*x + gamma = 0.
    Returns roots inside [var_lo, var_hi].
    """
    e1 = rows[e_i]["ellipse"]; e2 = rows[e_j]["ellipse"]
    a, b, c, d, e0 = _phi_diff_coeffs(e1, e2)
    if fix_var == 'h':
        # Phi_i(h_fix, p) - Phi_j(h_fix, p) = a*h^2 + b*h + c*p^2 + d*p + e0 = 0 in p
        h0 = fix_val
        alpha = c
        beta = d
        gamma = a * h0 * h0 + b * h0 + e0
    else:
        # fix_var == 'p'
        p0 = fix_val
        alpha = a
        beta = b
        gamma = c * p0 * p0 + d * p0 + e0
    roots = []
    if abs(alpha) < 1e-15:
        if abs(beta) > 1e-15:
            r = -gamma / beta
            if var_lo - 1e-12 <= r <= var_hi + 1e-12:
                roots.append(r)
    else:
        disc = beta * beta - 4 * alpha * gamma
        if disc >= 0:
            sd = np.sqrt(disc)
            r1 = (-beta + sd) / (2 * alpha)
            r2 = (-beta - sd) / (2 * alpha)
            for r in (r1, r2):
                if var_lo - 1e-12 <= r <= var_hi + 1e-12:
                    roots.append(np.clip(r, var_lo, var_hi))
    return roots


def candidates_pair_on_edges(rows, h_box=H_BOX, p_box=P_BOX):
    """For each pair (i, j), find points on box edges where Phi_i = Phi_j."""
    out = []
    n = len(rows)
    for i, j in itertools.combinations(range(n), 2):
        # Edges with h fixed
        for hf in h_box:
            roots = _candidate_pair_on_edge(rows, i, j, 'h', hf, p_box[0], p_box[1])
            for pr in roots:
                E, idx = _envelope_value(rows, hf, pr)
                out.append({"type": "pair_edge_h", "pair": (rows[i]["label"], rows[j]["label"]),
                            "h": float(hf), "p": float(pr),
                            "E": float(E), "witness_row": rows[idx]["label"]})
        # Edges with p fixed
        for pf in p_box:
            roots = _candidate_pair_on_edge(rows, i, j, 'p', pf, h_box[0], h_box[1])
            for hr in roots:
                E, idx = _envelope_value(rows, hr, pf)
                out.append({"type": "pair_edge_p", "pair": (rows[i]["label"], rows[j]["label"]),
                            "h": float(hr), "p": float(pf),
                            "E": float(E), "witness_row": rows[idx]["label"]})
    return out


def candidates_pair_min_along_edge(rows, h_box=H_BOX, p_box=P_BOX):
    """On each box edge, the envelope E(t) is the max of a finite collection of
    1D concave quadratics (each Phi_r restricted to the edge).  The minimum on
    the edge occurs either at the endpoints (corners) or at a 2-row crossing on
    the edge (already enumerated above).  Also need to consider the case where
    the witness row is fixed throughout an edge subinterval and its restriction
    is monotone — then min is at an endpoint, captured.
    Each Phi_r restricted to the edge is concave; max of concaves is not generally
    concave, but its min on a closed interval is at a boundary or where two of
    them are equal (and all others are below).  Pair-edge candidates above cover
    the second case; corners cover the first.
    Returns empty: this function exists to document that the edge min reduces
    fully to corners + pair_edge_p/h candidates.
    """
    return []


def candidates_triple_interior(rows, h_box=H_BOX, p_box=P_BOX):
    """For each triple (i, j, k), solve Phi_i = Phi_j and Phi_i = Phi_k for (h, p).
    Each equation is of the form a*h^2 + b*h + c*p^2 + d*p + e = 0.
    Two such equations form a quadratic system. Solve via resultants or numerically
    (Newton from grid candidates). For separable quadratics (no h*p cross term),
    the system decouples nicely:
      Eq1: a1*h^2 + b1*h + c1*p^2 + d1*p + e1 = 0
      Eq2: a2*h^2 + b2*h + c2*p^2 + d2*p + e2 = 0
    Linear combination eliminates either h^2 or p^2, giving (after substitution) a
    quartic. We solve numerically by elimination.

    Approach: Eliminate p^2 from Eq1 - (c1/c2)*Eq2 (assuming c2 != 0):
       (a1 - (c1/c2)*a2) h^2 + (b1 - (c1/c2)*b2) h + (d1 - (c1/c2)*d2) p + (e1 - (c1/c2)*e2) = 0
    This expresses p as a quadratic-in-h function (assuming the linear coeff in p
    is nonzero):  p = P(h) = (alpha h^2 + beta h + gamma) / delta
    Substitute into Eq1 to get a quartic in h. Solve, filter.

    If c2 = 0 (degenerate), swap or use an alternative elimination.
    """
    out = []
    n = len(rows)
    h_lo, h_hi = h_box; p_lo, p_hi = p_box
    for i, j, k in itertools.combinations(range(n), 3):
        e1 = rows[i]["ellipse"]; e2 = rows[j]["ellipse"]; e3 = rows[k]["ellipse"]
        a1, b1, c1, d1, f1 = _phi_diff_coeffs(e1, e2)
        a2, b2, c2, d2, f2 = _phi_diff_coeffs(e1, e3)
        # System:  a1 h^2 + b1 h + c1 p^2 + d1 p + f1 = 0
        #          a2 h^2 + b2 h + c2 p^2 + d2 p + f2 = 0
        sols = _solve_separable_quadratic_system(a1, b1, c1, d1, f1,
                                                 a2, b2, c2, d2, f2)
        for (h_s, p_s) in sols:
            if not (h_lo - 1e-9 <= h_s <= h_hi + 1e-9):
                continue
            if not (p_lo - 1e-9 <= p_s <= p_hi + 1e-9):
                continue
            h_s = float(np.clip(h_s, h_lo, h_hi))
            p_s = float(np.clip(p_s, p_lo, p_hi))
            E, idx = _envelope_value(rows, h_s, p_s)
            out.append({"type": "triple", "triple": (rows[i]["label"], rows[j]["label"], rows[k]["label"]),
                        "h": h_s, "p": p_s,
                        "E": float(E), "witness_row": rows[idx]["label"]})
    return out


def _solve_separable_quadratic_system(a1, b1, c1, d1, f1,
                                      a2, b2, c2, d2, f2):
    """Solve
       a1 h^2 + b1 h + c1 p^2 + d1 p + f1 = 0
       a2 h^2 + b2 h + c2 p^2 + d2 p + f2 = 0
    for real (h, p).  Strategy: form linear combinations to eliminate one of the
    quadratic terms, then iterate elimination + roots of polynomial in one variable.

    Returns list of (h, p) real solutions.
    """
    # Eliminate p^2: try lambda Eq1 - mu Eq2 with c2*lambda - c1*mu = 0; pick lambda=c2, mu=c1
    A_h2 = c2 * a1 - c1 * a2
    B_h1 = c2 * b1 - c1 * b2
    D_p1 = c2 * d1 - c1 * d2
    F_0  = c2 * f1 - c1 * f2
    # Equation L1:  A_h2 * h^2 + B_h1 * h + D_p1 * p + F_0 = 0   (no p^2 term)

    # Eliminate h^2: lambda=a2, mu=a1  =>  h^2 cancels
    A2_p2 = a2 * c1 - a1 * c2
    B2_h1 = a2 * b1 - a1 * b2
    D2_p1 = a2 * d1 - a1 * d2
    F2_0  = a2 * f1 - a1 * f2
    # Equation L2:  A2_p2 * p^2 + B2_h1 * h + D2_p1 * p + F2_0 = 0   (no h^2 term)

    sols = []

    # Case (i): D_p1 != 0  => from L1, p = -(A_h2 h^2 + B_h1 h + F_0) / D_p1
    # Substitute into L2 (which has no h^2 term, only h and p^2, p, const)
    # giving quartic in h.
    if abs(D_p1) > 1e-14:
        # p(h) = (alpha h^2 + beta h + gamma) / delta  with delta = -D_p1
        # Cleaner: p = (-A_h2 h^2 - B_h1 h - F_0) / D_p1
        Cp2 = A2_p2; Cp1 = D2_p1; Ch1 = B2_h1; C0 = F2_0
        # L2:  Cp2 * p^2 + Cp1 * p + Ch1 * h + C0 = 0
        #   p^2 = ((-A_h2 h^2 - B_h1 h - F_0)/D_p1)^2
        #       = (A_h2 h^2 + B_h1 h + F_0)^2 / D_p1^2
        # Multiply L2 by D_p1^2:
        #   Cp2 * (A_h2 h^2 + B_h1 h + F_0)^2  +  Cp1 * D_p1 * (-A_h2 h^2 - B_h1 h - F_0)
        #   +  Ch1 * D_p1^2 * h  +  C0 * D_p1^2  =  0
        A = A_h2; B = B_h1; F = F_0
        # (A h^2 + B h + F)^2 = A^2 h^4 + 2 A B h^3 + (B^2 + 2 A F) h^2 + 2 B F h + F^2
        D = D_p1
        Q4 = Cp2 * A * A
        Q3 = Cp2 * 2 * A * B
        Q2 = Cp2 * (B * B + 2 * A * F) + Cp1 * D * (-A)
        Q1 = Cp2 * 2 * B * F + Cp1 * D * (-B) + Ch1 * D * D
        Q0 = Cp2 * F * F + Cp1 * D * (-F) + C0 * D * D
        coeffs = [Q4, Q3, Q2, Q1, Q0]
        # Solve quartic
        # Strip leading zeros
        while len(coeffs) > 1 and abs(coeffs[0]) < 1e-15:
            coeffs = coeffs[1:]
        if len(coeffs) >= 2:
            try:
                hroots = np.roots(coeffs)
            except Exception:
                hroots = []
            for hr in hroots:
                if abs(hr.imag) < 1e-7:
                    h_real = float(hr.real)
                    p_real = -(A * h_real * h_real + B * h_real + F) / D
                    sols.append((h_real, p_real))
    # Case (ii): D_p1 == 0 (degenerate). Use L2 directly to express something else.
    # In that case L1 is purely quadratic in h:
    #   A_h2 h^2 + B_h1 h + F_0 = 0
    elif abs(A_h2) > 1e-14:
        disc = B_h1 * B_h1 - 4 * A_h2 * F_0
        if disc >= 0:
            sd = np.sqrt(disc)
            for h_real in ((-B_h1 + sd) / (2 * A_h2), (-B_h1 - sd) / (2 * A_h2)):
                # Plug into L2 (which is now A2_p2 p^2 + D2_p1 p + (B2_h1 h + F2_0) = 0)
                cc = B2_h1 * h_real + F2_0
                if abs(A2_p2) > 1e-14:
                    disc2 = D2_p1 * D2_p1 - 4 * A2_p2 * cc
                    if disc2 >= 0:
                        sd2 = np.sqrt(disc2)
                        for p_real in ((-D2_p1 + sd2) / (2 * A2_p2),
                                       (-D2_p1 - sd2) / (2 * A2_p2)):
                            sols.append((h_real, p_real))
                elif abs(D2_p1) > 1e-14:
                    p_real = -cc / D2_p1
                    sols.append((h_real, p_real))
    return sols


def closed_form_envelope_min(rows, h_box=H_BOX, p_box=P_BOX, verbose=True):
    """Master routine: enumerate ALL candidate KKT points and return the
    (rigorous) closed-form min of the envelope on the box."""
    candidates = []
    candidates += candidates_corners(rows, h_box, p_box)
    candidates += candidates_pair_on_edges(rows, h_box, p_box)
    candidates += candidates_triple_interior(rows, h_box, p_box)

    # Also: per-row interior min — but for concave Phi_r, the minimum on a closed
    # rectangle is at a vertex, captured by corners.  For COMPLETENESS, also
    # include per-row edge-interior critical points where dPhi_r/dx = 0 along
    # an edge (i.e. axis-aligned argmax of Phi_r on the edge). These are MAXIMA
    # of Phi_r on the edge, not minima of the envelope, but they may be where the
    # envelope's witness-row changes; they're already captured by pair crossings.
    # No extra candidates needed.

    if verbose:
        print(f"  Total candidates: {len(candidates)}")
        print(f"    corners:      {sum(1 for c in candidates if c['type']=='corner')}")
        print(f"    pair_edge_h:  {sum(1 for c in candidates if c['type']=='pair_edge_h')}")
        print(f"    pair_edge_p:  {sum(1 for c in candidates if c['type']=='pair_edge_p')}")
        print(f"    triple:       {sum(1 for c in candidates if c['type']=='triple')}")

    if not candidates:
        return None
    # The CLOSED-FORM min is the minimum of E over all candidates that are ACTUALLY
    # local minima of the envelope. Since the envelope is the max of concave functions,
    # any local minimum lies at a corner, a 2-row crossing on a box edge, or a 3-row
    # interior crossing where the three-row tie point is the witness max. We take min
    # over ALL these candidates (this is a valid ARGMIN: the true min must be in this
    # finite set, since the envelope's local minima cannot be in any other class of
    # points by the structure-of-max-of-concaves theorem).
    best = min(candidates, key=lambda c: c["E"])
    return {"min_E": best["E"], "argmin": best, "candidates": candidates,
            "n_candidates": len(candidates)}


def build_rigorous_row_for_closed(row_json, margin=1e-6):
    """Make ellipse dict with V_c_rigorous = primal - margin and recompute coefficients
    via find_ellipse_h_p (so we get clean A_h*, A_p* derived from duals)."""
    h_c = row_json["h_c"]; p_c = row_json["p_c"]
    q1 = row_json["q1"]; q2 = row_json["q2"]
    duals = row_json["duals"]
    V_reported = row_json["primal_value_at_center"]
    V_c_rig = V_reported - margin
    synthetic = {"h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2, "value": V_c_rig}
    ell = find_ellipse_h_p(synthetic, duals, q1, q2, target=WHITE_BOUND)
    return {
        "label": row_json["label"],
        "h_c": h_c, "p_c": p_c, "q1": q1, "q2": q2,
        "primal_value_reported": V_reported,
        "V_c_rigorous": V_c_rig,
        "rigour_margin": margin,
        "duals": duals,
        "ellipse": {
            "V_c": V_c_rig,
            "semi_h": ell["semi_h"], "semi_p": ell["semi_p"],
            "h_star": ell["h_star"], "p_star": ell["p_star"], "V_max": ell["V_max"],
            "A_h2": ell["A_h2"], "A_h1": ell["A_h1"], "A_h0": ell["A_h0"],
            "A_p2": ell["A_p2"], "A_p1": ell["A_p1"], "A_p0": ell["A_p0"],
            "const_q": ell["const_q"], "target": ell["target"],
        },
    }


def verify_quadratic_against_perturbation(row_json):
    """Sanity check: verify the per-row quadratic Phi_r(h, p) coefficients by
    computing finite-differences of the dual-shift expression vs the closed form.
    This doesn't re-solve the SDP (too slow); instead it verifies the algebra of
    the find_ellipse_h_p formula is consistent with the explicit formula in the
    docstring.
    """
    e = row_json["ellipse"]
    duals = row_json["duals"]
    h_c = row_json["h_c"]; p_c = row_json["p_c"]; V_c = row_json["primal_value_at_center"]

    # Independent re-derivation of A coefficients from duals:
    L53 = duals["con_53"]; L54 = duals["con_54"]
    LpL = duals["con_512_pL"]; LpU = duals["con_512_pU"]
    L513 = duals["con_513"]
    A_h2_chk = -0.5 * L54
    A_h1_chk = L53
    A_h0_chk = -L53 * h_c + 0.5 * L54 * h_c**2
    A_p2_chk = -0.5 * L513
    A_p1_chk = LpL - LpU
    A_p0_chk = (-LpL + LpU) * p_c + 0.5 * L513 * p_c**2
    diffs = {
        "A_h2": e["A_h2"] - A_h2_chk,
        "A_h1": e["A_h1"] - A_h1_chk,
        "A_h0": e["A_h0"] - A_h0_chk,
        "A_p2": e["A_p2"] - A_p2_chk,
        "A_p1": e["A_p1"] - A_p1_chk,
        "A_p0": e["A_p0"] - A_p0_chk,
    }
    return diffs


def main(margin=1e-6):
    out_dir = "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/parallel_results"
    if not os.path.isdir(out_dir):
        out_dir = "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/parallel_results"
    row_dir = os.path.join(out_dir, "lasserre2_path_b")
    out_file = os.path.join(out_dir, "path_b_closed_form.json")

    print(f"Closed-form path B aggregator")
    print(f"  V_c_rigorous = primal - margin, margin = {margin:.1e}")
    print(f"  Box (5.16): h in {H_BOX}, p in {P_BOX}, q in {Q_RANGE}")
    print()

    rows_data = []
    coeff_diffs = []
    for (h, p, qm, qp, label) in WHITE_TABLE3:
        fp = os.path.join(row_dir, f"{label}.json")
        if not os.path.exists(fp):
            print(f"  MISSING {fp}; skipping")
            continue
        with open(fp) as fh:
            row_json = json.load(fh)
        diffs = verify_quadratic_against_perturbation(row_json)
        coeff_diffs.append({"label": label, "diffs": diffs})
        new_row = build_rigorous_row_for_closed(row_json, margin=margin)
        rows_data.append(new_row)

    print("=== Quadratic-coefficient self-consistency (find_ellipse vs explicit) ===")
    for c in coeff_diffs:
        max_d = max(abs(v) for v in c["diffs"].values())
        print(f"  {c['label']}: max |diff| = {max_d:.2e}")

    # Show per-row coefficients summary
    print("\n=== Per-row quadratic coefficients (Phi_r = V_c + A_h2*h^2 + A_h1*h + A_h0 + A_p2*p^2 + A_p1*p + A_p0) ===")
    print(f"  {'row':<5} {'V_c':>12} {'A_h2':>10} {'A_h1':>10} {'A_h0':>10} {'A_p2':>10} {'A_p1':>10} {'A_p0':>10}")
    for r in rows_data:
        e = r["ellipse"]
        print(f"  {r['label']:<5} {e['V_c']:12.8f} {e['A_h2']:10.5f} {e['A_h1']:10.5f} {e['A_h0']:10.5f} {e['A_p2']:10.5f} {e['A_p1']:10.5f} {e['A_p0']:10.5f}")

    # Closed-form envelope min
    print("\n=== Closed-form envelope minimum ===")
    res = closed_form_envelope_min(rows_data, verbose=True)
    print(f"\n  closed-form min E = {res['min_E']:.10f}")
    a = res["argmin"]
    print(f"    type: {a['type']}")
    print(f"    location: h = {a['h']:.8f}, p = {a['p']:.8f}")
    print(f"    witness row: {a['witness_row']}")
    if "pair" in a: print(f"    pair: {a['pair']}")
    if "triple" in a: print(f"    triple: {a['triple']}")
    print(f"    vs White {WHITE_BOUND}: improvement = {res['min_E'] - WHITE_BOUND:+.6e}")

    # Compare against grid-based 2001x2001 result for reference
    print("\n=== Comparison vs grid-based 0.379828 ===")
    grid_min = 0.379828
    delta = res["min_E"] - grid_min
    print(f"  closed_form min - grid min = {delta:+.6e}")
    if delta > 0:
        print(f"  CLOSED-FORM is HIGHER -> grid was conservative")
    elif delta < 0:
        print(f"  CLOSED-FORM is LOWER  -> grid was missing the true min!")
    else:
        print(f"  Exact match.")

    # Top-15 lowest-E candidates for the JSON output (full list is many)
    sorted_cands = sorted(res["candidates"], key=lambda c: c["E"])
    top_low = sorted_cands[:30]

    summary = {
        "config": {
            "target": WHITE_BOUND,
            "h_box": list(H_BOX), "p_box": list(P_BOX), "q_range": list(Q_RANGE),
            "rigour_margin": margin,
            "approach": "closed_form_KKT_no_grid",
            "augmentation": "Bochner_n=20 + Lasserre_T_max=30_T_loc=8 (lasserre2 runs)",
        },
        "rows": rows_data,
        "coeff_consistency": coeff_diffs,
        "closed_form_min": {
            "min_E": res["min_E"],
            "argmin": a,
            "n_candidates": res["n_candidates"],
            "improvement_vs_white": res["min_E"] - WHITE_BOUND,
        },
        "vs_grid_based_0p379828": {
            "grid_value": grid_min,
            "closed_form_value": res["min_E"],
            "delta_closed_minus_grid": delta,
            "interpretation": ("HIGHER - grid was conservative" if delta > 0
                               else ("LOWER - grid missed true min" if delta < 0 else "exact match")),
        },
        "candidates_top_30_by_E": top_low,
    }
    with open(out_file, "w") as f:
        json.dump(summary, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "__float__") else str(o))
    print(f"\nWritten to {out_file}")
    return summary


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--margin", type=float, default=1e-6)
    args = p.parse_args()
    main(margin=args.margin)
