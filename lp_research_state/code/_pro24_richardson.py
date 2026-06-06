"""PRO-24: Richardson extrapolation on Ω*(N) trajectory.

Fit the model
    Ω*(N) = μ + c1 / N^α + c2 / N^(2α)
to the row-4 binding-center trajectory at fixed bn ∈ {20, 30, 40}, and the
Phase-5 cover trajectory.

Outputs:
- Best-fit μ (the N → ∞ ceiling)
- Discrimination between α=1 (pure power law) and α<1 (log corrections)
- Comparison to current rigorous LB headline 0.3803027228
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / "lp_research_state" / "experiments_done.json"
PARA = ROOT / "lp_research_state" / "parallel_results"
OUT = ROOT / "lp_research_state" / "data" / "pro24_richardson.json"


def load_row4_trajectory(bn: int) -> list[dict]:
    """Return all row-4 results at the given bochner_n."""
    with open(EXP) as f:
        results = json.load(f)["results"]
    out = []
    for r in results:
        if r.get("row") != "row4":
            continue
        if r.get("bochner_n") != bn:
            continue
        if r.get("value") is None:
            continue
        out.append({
            "N": int(r["N"]),
            "T": int(r["T"]),
            "value": float(r["value"]),
            "rigorous": r.get("rigorous_dual_LB"),
            "kind": r.get("kind"),
        })
    # dedupe by N — keep dual-extracted (lp_run_bochner_dual) if available
    by_N = {}
    for r in out:
        N = r["N"]
        if N not in by_N or "dual" in r["kind"]:
            by_N[N] = r
    out = sorted(by_N.values(), key=lambda x: x["N"])
    return out


def load_phase5_trajectory() -> list[dict]:
    """Return Phase-5 rigorous_LB at each N."""
    import re
    pts = []
    for p in sorted(PARA.glob("phase5_N*.json")):
        with open(p) as f:
            d = json.load(f)
        # Extract N from content if present, else parse from filename
        N = d.get("N")
        if N is None:
            # Parse filenames like "phase5_N15000.json" or "phase5_N20K_bn40.json"
            m = re.search(r"phase5_N(\d+)(K)?", p.name)
            if m:
                base = int(m.group(1))
                if m.group(2) == "K":
                    base *= 1000
                N = base
        rig = d.get("rigorous_LB")
        grid = d.get("grid_min")
        bn = d.get("bochner_n", 30)
        if N and rig:
            pts.append({"N": N, "rigorous_LB": rig, "grid_min": grid, "bn": bn,
                        "source": p.name})
    pts = sorted(pts, key=lambda x: x["N"])
    return pts


def model_alpha(N, mu, c1, c2, alpha):
    """Ω*(N) = mu + c1/N^alpha + c2/N^(2 alpha)."""
    return mu + c1 / N**alpha + c2 / N**(2 * alpha)


def model_alpha1(N, mu, c1, c2):
    """Ω*(N) = mu + c1/N + c2/N^2 (pure power-law, α=1)."""
    return mu + c1 / N + c2 / N**2


def model_log_correction(N, mu, c1, c2):
    """Ω*(N) = mu + c1 * log(N)/N + c2/N (one log term)."""
    return mu + c1 * np.log(N) / N + c2 / N


def fit_all(Ns: np.ndarray, ys: np.ndarray, label: str) -> dict:
    """Fit three models and report extrapolated μ for each."""
    out = {"label": label, "N_data": Ns.tolist(), "y_data": ys.tolist()}

    # Model 1: free alpha
    try:
        p_free, cov_free = curve_fit(
            model_alpha, Ns, ys,
            p0=[ys.max(), -1.0, 0.1, 1.0],
            maxfev=10000,
            bounds=([0.0, -np.inf, -np.inf, 0.5], [1.0, np.inf, np.inf, 2.0]),
        )
        err_free = np.sqrt(np.diag(cov_free))
        resid = ys - model_alpha(Ns, *p_free)
        out["free_alpha"] = {
            "mu": float(p_free[0]), "c1": float(p_free[1]),
            "c2": float(p_free[2]), "alpha": float(p_free[3]),
            "mu_err": float(err_free[0]), "alpha_err": float(err_free[3]),
            "residual_rms": float(np.sqrt(np.mean(resid**2))),
            "residuals": resid.tolist(),
        }
    except Exception as e:
        out["free_alpha"] = {"error": str(e)}

    # Model 2: alpha = 1 (constrained)
    try:
        p_a1, cov_a1 = curve_fit(
            model_alpha1, Ns, ys,
            p0=[ys.max(), -1.0, 0.1],
            maxfev=10000,
        )
        err_a1 = np.sqrt(np.diag(cov_a1))
        resid = ys - model_alpha1(Ns, *p_a1)
        out["alpha_1"] = {
            "mu": float(p_a1[0]), "c1": float(p_a1[1]), "c2": float(p_a1[2]),
            "mu_err": float(err_a1[0]),
            "residual_rms": float(np.sqrt(np.mean(resid**2))),
            "residuals": resid.tolist(),
        }
    except Exception as e:
        out["alpha_1"] = {"error": str(e)}

    # Model 3: log correction
    try:
        p_log, cov_log = curve_fit(
            model_log_correction, Ns, ys,
            p0=[ys.max(), -1.0, 0.1],
            maxfev=10000,
        )
        err_log = np.sqrt(np.diag(cov_log))
        resid = ys - model_log_correction(Ns, *p_log)
        out["log_correction"] = {
            "mu": float(p_log[0]), "c1": float(p_log[1]), "c2": float(p_log[2]),
            "mu_err": float(err_log[0]),
            "residual_rms": float(np.sqrt(np.mean(resid**2))),
            "residuals": resid.tolist(),
        }
    except Exception as e:
        out["log_correction"] = {"error": str(e)}

    return out


def main():
    print("=" * 78)
    print("PRO-24: Richardson extrapolation on Ω*(N)")
    print("=" * 78)

    all_fits = {}

    # ===== Row 4 trajectory at bn=20 =====
    traj20 = load_row4_trajectory(bn=20)
    print(f"\n--- Row 4, bn=20 ({len(traj20)} points) ---")
    for r in traj20:
        rig = f" (rig {r['rigorous']:.10f})" if r["rigorous"] else ""
        print(f"  N={r['N']:>5}  Ω* = {r['value']:.10f}{rig}")
    if len(traj20) >= 3:
        Ns = np.array([r["N"] for r in traj20], dtype=float)
        ys = np.array([r["rigorous"] if r["rigorous"] else r["value"]
                       for r in traj20], dtype=float)
        all_fits["row4_bn20"] = fit_all(Ns, ys, "row4_bn20")
        report(all_fits["row4_bn20"])

    # ===== Row 4 trajectory at bn=30 =====
    traj30 = load_row4_trajectory(bn=30)
    print(f"\n--- Row 4, bn=30 ({len(traj30)} points) ---")
    for r in traj30:
        rig = f" (rig {r['rigorous']:.10f})" if r["rigorous"] else ""
        print(f"  N={r['N']:>5}  Ω* = {r['value']:.10f}{rig}")
    if len(traj30) >= 3:
        Ns = np.array([r["N"] for r in traj30], dtype=float)
        ys = np.array([r["rigorous"] if r["rigorous"] else r["value"]
                       for r in traj30], dtype=float)
        all_fits["row4_bn30"] = fit_all(Ns, ys, "row4_bn30")
        report(all_fits["row4_bn30"])

    # ===== Phase 5 cover trajectory =====
    phase5 = load_phase5_trajectory()
    print(f"\n--- Phase 5 cover ({len(phase5)} points) ---")
    for r in phase5:
        print(f"  N={r['N']:>5}  rigorous_LB = {r['rigorous_LB']:.10f}  bn={r['bn']}  ({r['source']})")
    if len(phase5) >= 3:
        Ns = np.array([r["N"] for r in phase5], dtype=float)
        ys = np.array([r["rigorous_LB"] for r in phase5], dtype=float)
        all_fits["phase5"] = fit_all(Ns, ys, "phase5")
        report(all_fits["phase5"])

    # ===== Combined row-4 trajectory (bn=20 and bn=30 merged via bn-shift) =====
    # Empirically (findings.md), the bn=20→30 vertical shift is ~+2.55e-4 at large N.
    # Pool both bn series into a single (N, Ω_normalized) cloud using this shift.
    print("\n--- Combined row-4 (bn=20 + bn=30 shifted) ---")
    bn_shift = 2.55e-4  # empirical Δ from bn=20 to bn=30
    pool = []
    for r in traj30:
        pool.append({"N": r["N"], "value": r["rigorous"] if r["rigorous"] else r["value"]})
    for r in traj20:
        # Shift bn=20 values up by bn_shift to align with bn=30 series
        v = r["rigorous"] if r["rigorous"] else r["value"]
        pool.append({"N": r["N"], "value": v + bn_shift})
    pool = sorted(pool, key=lambda x: x["N"])
    for r in pool:
        print(f"  N={r['N']:>5}  Ω* (bn30-equiv) = {r['value']:.10f}")
    if len(pool) >= 4:
        Ns = np.array([r["N"] for r in pool], dtype=float)
        ys = np.array([r["value"] for r in pool], dtype=float)
        all_fits["row4_pooled"] = fit_all(Ns, ys, "row4_pooled")
        report(all_fits["row4_pooled"])

    # Save
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(all_fits, f, indent=2)
    print(f"\nFits saved to: {OUT}")


def report(fit_result: dict) -> None:
    """Print a Richardson fit's summary."""
    print()
    for model in ["free_alpha", "alpha_1", "log_correction"]:
        f = fit_result.get(model)
        if not f or "error" in f:
            continue
        mu = f["mu"]
        mu_err = f.get("mu_err", float("nan"))
        rms = f.get("residual_rms", float("nan"))
        extra = ""
        if model == "free_alpha":
            extra = f"  α = {f['alpha']:.4f} ± {f.get('alpha_err', 0):.4f}"
        print(f"  [{model:<16}] μ = {mu:.10f} ± {mu_err:.2e}  RMS = {rms:.2e}{extra}")


if __name__ == "__main__":
    main()
