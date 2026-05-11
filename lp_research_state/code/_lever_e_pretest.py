"""Lever E pretest: three cheap tests to decide whether to invest weeks in an
M-side SDP encoding (minimizing Together's M-functional instead of/alongside
White's Ω).

Test 1 — Compute Together's M-functional on our LP-optimal f̃.
Test 2 — Run existing mside_bochner_n SDP at row 4.
Test 3 — Report what findings.md / existing code already establish.

Outputs JSON at lp_research_state/data/lever_e_pretest_results.json and is
intended to be summarized in LEVER_E_PRETEST.md at repo root.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Set up imports.
CODE_DIR = Path(__file__).resolve().parent
REPO_ROOT = CODE_DIR.parent.parent
DATA_DIR = REPO_ROOT / "lp_research_state" / "data"
sys.path.insert(0, str(CODE_DIR))


def cell_step_M(values_on_pos_half: np.ndarray, L: float) -> float:
    """Together's M-functional on a step function on [0,2] with cells of width L.

    M(h) = L * max_k Σ_i h_i (1 - h_{i+k}), zero-extension outside [0,2].
    """
    h = np.asarray(values_on_pos_half, dtype=np.float64)
    corr = np.correlate(h, 1.0 - h, mode="full")
    return float(L * corr.max())


def test_1_M_of_f_tilde() -> dict:
    """Apply Together's M to our LP-optimal f̃ under multiple interpretations.

    f̃ lives on [-2, 2] (10000-point grid). It violates [0,1]:
        min ≈ -0.97, max ≈ 1.55.
    We report:
      - M_strict_pos: M of f̃|_{[0,2]} (positive-half slice; no clipping).
      - M_strict_fold: M of |f̃(x)| + |f̃(-x)| ? Actually we use the
        even-symmetrized positive half g(x) = (f̃(x) + f̃(-x))/2 on [0,2].
      - M_strict_abs: M of |f̃|_{[0,2]} (defensive — at least nonneg).
      - M_clipped: M of clip(f̃, 0, 1) restricted to [0,2].
      - For each: also report ∫ f̃-slice (Together expects ∫ h = 1).
    """
    npz = np.load(DATA_DIR / "row4_f_tilde.npz")
    x = npz["x"]  # (10000,) on [-2, 2]
    f = npz["f_tilde"]  # (10000,)
    n = len(x)
    assert n == 10000 and abs(x[0] - (-2.0)) < 1e-12 and abs(x[-1] - 2.0) < 1e-9

    # Build a uniform cell-step interpretation: midpoint-based.
    # Cells are [-2 + j*Δx, -2 + (j+1)*Δx) for Δx = 4/n. The saved samples
    # likely lie at cell midpoints (j+0.5)*Δx - 2. Use them as cell values.
    dx = 4.0 / n  # cell width across [-2, 2]
    L_pos = dx     # cell width on positive half is identical
    assert abs(L_pos * (n // 2) - 2.0) < 1e-12, (
        "expected exactly n/2 cells on [0, 2]"
    )

    # Identify the positive-half cells (last n/2 entries map to x >= 0).
    half = n // 2
    f_pos = f[half:]  # values on [0, 2]
    f_neg = f[:half]  # values on [-2, 0]

    # 1a: directly use f̃ on [0, 2], no clipping.
    M_pos = cell_step_M(f_pos, L_pos)
    int_pos = float(L_pos * f_pos.sum())

    # 1b: even-fold. g(x) = (f̃(x) + f̃(-x))/2 on [0, 2]. f_neg in left-to-right
    # order corresponds to x = -2 + (j+0.5)*dx for j=0..half-1; its reverse
    # corresponds to x = -2 + (n - j - 0.5)*dx, i.e. positive-side cells.
    # So the partner of f_pos[i] is f_neg[half - 1 - i].
    g_fold = (f_pos + f_neg[::-1]) / 2.0
    M_fold = cell_step_M(g_fold, L_pos)
    int_fold = float(L_pos * g_fold.sum())

    # 1c: |f̃| on [0, 2] — non-physical but ensures nonneg.
    M_abs = cell_step_M(np.abs(f_pos), L_pos)
    int_abs = float(L_pos * np.abs(f_pos).sum())

    # 1d: clip(f̃, 0, 1) on [0, 2].
    f_clip_pos = np.clip(f_pos, 0.0, 1.0)
    M_clip_pos = cell_step_M(f_clip_pos, L_pos)
    int_clip_pos = float(L_pos * f_clip_pos.sum())

    # 1e: even-fold then clip.
    g_fold_clip = np.clip(g_fold, 0.0, 1.0)
    M_clip_fold = cell_step_M(g_fold_clip, L_pos)
    int_clip_fold = float(L_pos * g_fold_clip.sum())

    # Together's claimed bound.
    M_together = 0.380871

    return {
        "n_grid": n,
        "dx_cell_width": dx,
        "f_tilde_min": float(f.min()),
        "f_tilde_max": float(f.max()),
        "f_tilde_int_full": float(dx * f.sum()),
        "M_together_claimed": M_together,
        # Interpretation 1a: positive-half slice, no transform
        "M_strict_pos_slice": M_pos,
        "int_pos_slice": int_pos,
        # 1b: even-fold
        "M_strict_even_fold": M_fold,
        "int_even_fold": int_fold,
        # 1c: |f| on positive half
        "M_strict_abs_pos": M_abs,
        "int_abs_pos": int_abs,
        # 1d: clip on positive half
        "M_clipped_pos": M_clip_pos,
        "int_clipped_pos": int_clip_pos,
        # 1e: even-fold then clip
        "M_clipped_even_fold": M_clip_fold,
        "int_clipped_even_fold": int_clip_fold,
    }


def test_2_mside_bochner(n_M_list=(8,), include_schur=False, include_lasserre=False) -> dict:
    """Run row 4 build_problem with mside_bochner_n set, extract rigorous LB."""
    import cvxpy as cp
    from white_full_convex import build_problem
    from dual_extractor import solve_with_dual_extraction

    results = []

    # Baseline row 4 — same as CLAUDE.md snippet, bochner_n=20 only.
    common_args = dict(
        N=10000, T=4000, R=10, h1=0.004, h2=0.004,
        p1=0.3875, p2=0.3875, q1=-0.02, q2=0.02,
        bochner_n=20,
    )

    print("=" * 70)
    print("Test 2a: baseline (bochner_n=20, no M-side) — sanity check")
    print("=" * 70)
    t0 = time.time()
    Omega, w, v, c, d, eps, dlt, cons = build_problem(**common_args)
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    elapsed = time.time() - t0
    print(f"  status={res['status']}, reported={res['reported_value']}, "
          f"rigorous_LB={res['rigorous_dual_LB']}, elapsed={elapsed:.1f}s")
    results.append({
        "config": "bochner_n=20 only (baseline)",
        "mside_bochner_n": 0,
        "elapsed_s": elapsed,
        "status": res["status"],
        "reported_value": res["reported_value"],
        "rigorous_dual_LB": res["rigorous_dual_LB"],
        "dual_residual_at_LB": res.get("dual_residual_at_LB"),
    })

    # SOC variants.
    for n_M in n_M_list:
        print("=" * 70)
        print(f"Test 2b: mside_bochner_n={n_M} (SOC) + bochner_n=20")
        print("=" * 70)
        t0 = time.time()
        Omega, w, v, c, d, eps, dlt, cons = build_problem(
            **common_args, mside_bochner_n=n_M
        )
        prob = cp.Problem(cp.Minimize(Omega), cons)
        try:
            res = solve_with_dual_extraction(prob)
            elapsed = time.time() - t0
            print(f"  status={res['status']}, "
                  f"reported={res['reported_value']}, "
                  f"rigorous_LB={res['rigorous_dual_LB']}, "
                  f"elapsed={elapsed:.1f}s")
            results.append({
                "config": f"mside_bochner_n={n_M} (SOC) + bochner_n=20",
                "mside_bochner_n": n_M,
                "elapsed_s": elapsed,
                "status": res["status"],
                "reported_value": res["reported_value"],
                "rigorous_dual_LB": res["rigorous_dual_LB"],
                "dual_residual_at_LB": res.get("dual_residual_at_LB"),
            })
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  FAILED: {e}")
            results.append({
                "config": f"mside_bochner_n={n_M} (SOC) + bochner_n=20",
                "mside_bochner_n": n_M,
                "elapsed_s": elapsed,
                "error": str(e),
            })

    # Schur variants (if requested).
    if include_schur:
        for n_M in n_M_list:
            print("=" * 70)
            print(f"Test 2c: mside_bochner_schur_n={n_M} + bochner_n=20")
            print("=" * 70)
            t0 = time.time()
            Omega, w, v, c, d, eps, dlt, cons = build_problem(
                **common_args, mside_bochner_schur_n=n_M
            )
            prob = cp.Problem(cp.Minimize(Omega), cons)
            try:
                res = solve_with_dual_extraction(prob)
                elapsed = time.time() - t0
                print(f"  status={res['status']}, "
                      f"reported={res['reported_value']}, "
                      f"rigorous_LB={res['rigorous_dual_LB']}, "
                      f"elapsed={elapsed:.1f}s")
                results.append({
                    "config": f"mside_bochner_schur_n={n_M} + bochner_n=20",
                    "mside_bochner_schur_n": n_M,
                    "elapsed_s": elapsed,
                    "status": res["status"],
                    "reported_value": res["reported_value"],
                    "rigorous_dual_LB": res["rigorous_dual_LB"],
                    "dual_residual_at_LB": res.get("dual_residual_at_LB"),
                })
            except Exception as e:
                elapsed = time.time() - t0
                print(f"  FAILED: {e}")
                results.append({
                    "config": f"mside_bochner_schur_n={n_M} + bochner_n=20",
                    "mside_bochner_schur_n": n_M,
                    "elapsed_s": elapsed,
                    "error": str(e),
                })

    return {"row4_runs": results}


def main():
    out = {
        "test_1_M_of_f_tilde": None,
        "test_2_mside_bochner": None,
        "test_3_prior_findings": {
            "summary": (
                "findings.md (2026-05-10) documents that the existing SOC-relaxed "
                "M-side Bochner (mside_bochner_n) at n_M=5 produced Δ=+1.4e-9 vs "
                "baseline and at n_M=10 produced Δ=+1.65e-8 — empirically dead, "
                "labelled 'CANCEL' for higher n_M. Mechanism: the SOC slack "
                "U_m ≥ |f̂(m)|² absorbs the M-side PSD content without "
                "constraining (c,d). mside_bochner_schur defines the SAME convex "
                "set (per its docstring), so we expect identical behaviour. "
                "An exact-bilinear variant via Lasserre-lifted moments "
                "(mside_via_lasserre.py) exists but requires lasserre_T_max > 0; "
                "Lasserre level-2 is itself documented as 'non-rigorous and "
                "withdrawn' in CLAUDE.md and its localizing tail bound kills the "
                "gain quantitatively (communications/lasserre_tail_bound.md)."
            ),
            "prior_run_n5": {
                "config": "mside_bochner_n=5 (SOC), row4 N=2000",
                "Omega": 0.37627652427597486,
                "delta_vs_baseline": 1.4e-9,
                "status": "optimal",
            },
            "prior_run_n10": {
                "config": "mside_bochner_n=10 (SOC), row4 N=2000",
                "Omega": 0.3762765394,
                "delta_vs_baseline": 1.65e-8,
                "status": "optimal_inaccurate",
            },
            "phase_5_baseline": 0.380128,
            "white_baseline": 0.379005,
            "together_UB": 0.380871,
        },
    }

    print("\n###### TEST 1 ######\n")
    try:
        t1 = test_1_M_of_f_tilde()
        out["test_1_M_of_f_tilde"] = t1
        print(json.dumps(t1, indent=2))
    except Exception as e:
        import traceback
        traceback.print_exc()
        out["test_1_M_of_f_tilde"] = {"error": str(e)}

    print("\n###### TEST 2 ######\n")
    try:
        # Try n_M=8 only first (per anti-pattern note: 5-min cap)
        t2 = test_2_mside_bochner(n_M_list=(8,), include_schur=False)
        out["test_2_mside_bochner"] = t2
    except Exception as e:
        import traceback
        traceback.print_exc()
        out["test_2_mside_bochner"] = {"error": str(e)}

    out_path = DATA_DIR / "lever_e_pretest_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nSaved results to {out_path}")
    return out


if __name__ == "__main__":
    main()
