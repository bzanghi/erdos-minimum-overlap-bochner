"""
Single-experiment driver invoked by the cron Claude session.

Reads experiments_queue.json, picks the highest-priority pending experiment,
runs it, appends to experiments_done.json, and prints a concise summary.

Usage:
    python3 cron_runner.py
    python3 cron_runner.py --dry-run     # show next experiment without running

The `experiments_queue.json` schema is documented in findings.md. The cron
Claude is responsible for reading this driver's output, updating findings.md
with any new insight, and (optionally) appending fresh experiments to the
queue.
"""
from __future__ import annotations
import json, os, sys, time, traceback, importlib.util
from pathlib import Path

STATE_DIR = Path(os.environ.get(
    "LP_STATE_DIR",
    "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state"
))
CODE_DIR = Path(os.environ.get(
    "LP_CODE_DIR",
    "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code"
))

# Future cron sessions may see the persistent code dir at a different mount
# point. Try a few candidates and use whichever exists.
def _safe_exists(p: str) -> bool:
    try:
        return Path(p).exists()
    except (PermissionError, OSError):
        return False

import glob as _glob
_dynamic_code = _glob.glob("/sessions/*/mnt/Erdos/lp_research_state/code")
_dynamic_state = _glob.glob("/sessions/*/mnt/Erdos/lp_research_state")

for candidate in _dynamic_code + [
    "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state/code",
    "/sessions/Erdos/lp_research_state/code",
    "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state/code",
]:
    if _safe_exists(candidate):
        CODE_DIR = Path(candidate)
        break

for candidate in _dynamic_state + [
    "/sessions/keen-magical-meitner/mnt/Erdos/lp_research_state",
    "/sessions/Erdos/lp_research_state",
    "/Users/benzanghi/Documents/Claude/Projects/Erdos/lp_research_state",
]:
    if _safe_exists(candidate):
        STATE_DIR = Path(candidate)
        break

QUEUE = STATE_DIR / "experiments_queue.json"
DONE = STATE_DIR / "experiments_done.json"
LOG = STATE_DIR / "cron_log.txt"

POINTS_BY_ROW = {
    "row1": (0.015, 0.381,  -0.02, 0.02),
    "row2": (0.015, 0.385,  -0.02, 0.02),
    "row3": (0.020, 0.375,  -0.02, 0.02),
    "row4": (0.004, 0.3875, -0.02, 0.02),
    "row5": (0.000, 0.4,    -0.02, 0.02),
    "row6": (0.000, 0.381,  -0.02, 0.02),
    "row7": (0.030, 0.375,  -0.02, 0.02),
}


def log(msg: str):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    try:
        with open(LOG, "a") as f:
            f.write(line + "\n")
    except Exception:
        pass


def load_json(p: Path, default):
    if p.exists():
        return json.loads(p.read_text())
    return default


def save_json(p: Path, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2))


def import_white_module():
    spec = importlib.util.spec_from_file_location(
        "white_full_convex", CODE_DIR / "white_full_convex.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load white_full_convex from {CODE_DIR}")
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(CODE_DIR))
    spec.loader.exec_module(mod)
    return mod


def already_done(done, exp):
    """Return True if an experiment with same params is already in the done log."""
    for r in done.get("results", []):
        if (
            r.get("kind") == exp.get("kind")
            and r.get("N") == exp["params"].get("N")
            and r.get("T") == exp["params"].get("T")
            and r.get("R") == exp["params"].get("R")
            and r.get("row") in exp["params"].get("rows", [])
            and r.get("alpha") == exp.get("alpha")
            and r.get("use_T5p") == exp.get("use_T5p")
        ):
            return True
    return False


def run_simple_t5p(N, T, R, rows, use_T5p_values, white):
    """Run baseline and +T5p on the listed rows. Return list of result dicts."""
    out = []
    for row in rows:
        h, p, qm, qp = POINTS_BY_ROW[row]
        for use_T5p in use_T5p_values:
            t0 = time.time()
            try:
                res = white.solve_full_program(
                    N, T, R, h, h, p, p, qm, qp, use_T5p=use_T5p
                )
                out.append({
                    "kind": "lp_run",
                    "N": N, "T": T, "R": R,
                    "row": row, "h": h, "p": p, "q_range": [qm, qp],
                    "use_T5p": use_T5p,
                    "value": float(res["value"]) if res["value"] is not None else None,
                    "status": res["status"],
                    "time": time.time() - t0,
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
            except Exception as e:
                out.append({
                    "kind": "lp_run",
                    "N": N, "T": T, "R": R,
                    "row": row, "use_T5p": use_T5p,
                    "error": f"{type(e).__name__}: {e}",
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
    return out


def run_bochner(N, T, R, rows, bochner_ns, white, also_T5p=False):
    """Run Bochner-augmented LP for a list of rows and Bochner levels."""
    out = []
    for row in rows:
        h, p, qm, qp = POINTS_BY_ROW[row]
        for n_b in bochner_ns:
            t0 = time.time()
            try:
                res = white.solve_full_program(
                    N, T, R, h, h, p, p, qm, qp,
                    use_T5p=also_T5p, bochner_n=n_b,
                )
                out.append({
                    "kind": "lp_run_bochner",
                    "N": N, "T": T, "R": R,
                    "row": row, "h": h, "p": p, "q_range": [qm, qp],
                    "use_T5p": also_T5p,
                    "bochner_n": n_b,
                    "value": float(res["value"]) if res["value"] is not None else None,
                    "status": res["status"],
                    "time": time.time() - t0,
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
            except Exception as e:
                out.append({
                    "kind": "lp_run_bochner",
                    "N": N, "T": T, "R": R,
                    "row": row, "bochner_n": n_b, "use_T5p": also_T5p,
                    "error": f"{type(e).__name__}: {e}",
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
    return out


def run_bochner_dual(N, T, R, rows, bochner_ns, white, also_T5p=False):
    """Run Bochner-augmented LP with CLARABEL verbose dual-extraction.

    Persists rigorous_dual_LB = reported_value - last_gap (the precise rigorous
    LB on the LP optimum, ≈ 6+ orders of magnitude tighter than the 1e-4
    inaccurate-status safety convention).
    """
    import cvxpy as cp
    # Find dual_extractor.py — it lives in the same code dir as white_full_convex.
    spec = importlib.util.spec_from_file_location(
        "dual_extractor", CODE_DIR / "dual_extractor.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load dual_extractor from {CODE_DIR}")
    de = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(de)

    out = []
    for row in rows:
        h, p, qm, qp = POINTS_BY_ROW[row]
        for n_b in bochner_ns:
            t0 = time.time()
            try:
                Omega, w_v, v_v, c_v, d_v, eps_v, dlt_v, cons = white.build_problem(
                    N, T, R, h, h, p, p, qm, qp,
                    use_T5p=also_T5p, bochner_n=n_b,
                )
                prob = cp.Problem(cp.Minimize(Omega), cons)
                res = de.solve_with_dual_extraction(prob)
                last_gap = None
                if res["raw_iterations"]:
                    last_gap = res["raw_iterations"][-1]["gap"]
                rigorous_LB_vmg = None
                if res["reported_value"] is not None and last_gap is not None:
                    rigorous_LB_vmg = res["reported_value"] - last_gap
                out.append({
                    "kind": "lp_run_bochner_dual",
                    "N": N, "T": T, "R": R,
                    "row": row, "h": h, "p": p, "q_range": [qm, qp],
                    "use_T5p": also_T5p,
                    "bochner_n": n_b,
                    "value": res["reported_value"],
                    "status": res["status"],
                    "rigorous_dual_LB": rigorous_LB_vmg,
                    "rigorous_dual_LB_low_precision": res["rigorous_dual_LB"],
                    "last_gap": last_gap,
                    "dual_residual_at_LB": res["dual_residual_at_LB"],
                    "best_iter": res["best_iter"],
                    "n_iters_total": res["n_iters_total"],
                    "time": time.time() - t0,
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
            except Exception as e:
                out.append({
                    "kind": "lp_run_bochner_dual",
                    "N": N, "T": T, "R": R,
                    "row": row, "bochner_n": n_b, "use_T5p": also_T5p,
                    "error": f"{type(e).__name__}: {e}",
                    "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                })
    return out


def run_mside_bochner(N, T, R, rows, mside_bochner_ns, white,
                      bochner_ns=(0,), also_T5p=False):
    """Run M-side (SOC-relaxed) Bochner-augmented LP for a list of rows and
    M-side levels. Optionally combine with f-side Bochner via `bochner_ns`.

    Persists kind='lp_run_mside_bochner' result rows; one per
    (row, mside_bochner_n, bochner_n) triple.
    """
    out = []
    for row in rows:
        h, p, qm, qp = POINTS_BY_ROW[row]
        for n_M in mside_bochner_ns:
            for n_b in bochner_ns:
                t0 = time.time()
                try:
                    res = white.solve_full_program(
                        N, T, R, h, h, p, p, qm, qp,
                        use_T5p=also_T5p,
                        bochner_n=n_b,
                        mside_bochner_n=n_M,
                    )
                    out.append({
                        "kind": "lp_run_mside_bochner",
                        "N": N, "T": T, "R": R,
                        "row": row, "h": h, "p": p, "q_range": [qm, qp],
                        "use_T5p": also_T5p,
                        "bochner_n": n_b,
                        "mside_bochner_n": n_M,
                        "value": float(res["value"]) if res["value"] is not None else None,
                        "status": res["status"],
                        "time": time.time() - t0,
                        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                    })
                except Exception as e:
                    out.append({
                        "kind": "lp_run_mside_bochner",
                        "N": N, "T": T, "R": R,
                        "row": row, "bochner_n": n_b,
                        "mside_bochner_n": n_M, "use_T5p": also_T5p,
                        "error": f"{type(e).__name__}: {e}",
                        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
                    })
    return out


def run_alpha_sweep(N, T, R, rows, alphas, white):
    """For each alpha, build the LP with a custom phi_alpha constraint and solve.

    NOT YET IMPLEMENTED in white_full_convex.py: we'd need to extend it to
    accept arbitrary alpha. For now, this returns a placeholder noting the
    cron task should extend the solver."""
    return [{
        "kind": "alpha_sweep_skipped",
        "reason": "white_full_convex.py does not yet accept arbitrary alpha; extend it to add the convex constraint c.T M_alpha [c;d] <= 0.5.",
        "params": {"N": N, "T": T, "R": R, "rows": rows, "alphas": alphas},
        "ts": time.strftime("%Y-%m-%d %H:%M:%S"),
    }]


def main():
    args = sys.argv[1:]
    dry_run = "--dry-run" in args

    queue_obj = load_json(QUEUE, {"next": []})
    done = load_json(DONE, {"results": []})

    pending = sorted(queue_obj.get("next", []), key=lambda e: e.get("priority", 99))
    if not pending:
        log("Queue empty. Nothing to do.")
        return 0

    # Pick first whose work is not redundant.
    def is_done(cand):
        if cand["kind"] in ("lp_run_bochner", "lp_run_bochner_sweep"):
            ns = cand["params"].get("bochner_ns", [cand["params"].get("bochner_n", 0)])
            for row in cand["params"].get("rows", []):
                for n_b in ns:
                    hit = False
                    for r in done["results"]:
                        if (
                            r.get("N") == cand["params"]["N"]
                            and r.get("T") == cand["params"]["T"]
                            and r.get("R") == cand["params"]["R"]
                            and r.get("row") == row
                            and r.get("bochner_n") == n_b
                        ):
                            hit = True; break
                    if not hit: return False
            return True
        elif cand["kind"] == "lp_run_bochner_dual":
            # Distinct from lp_run_bochner: requires a result row that has the
            # rigorous_dual_LB field (not just `value`).
            ns = cand["params"].get("bochner_ns", [cand["params"].get("bochner_n", 0)])
            for row in cand["params"].get("rows", []):
                for n_b in ns:
                    hit = False
                    for r in done["results"]:
                        if (
                            r.get("kind") == "lp_run_bochner_dual"
                            and r.get("N") == cand["params"]["N"]
                            and r.get("T") == cand["params"]["T"]
                            and r.get("R") == cand["params"]["R"]
                            and r.get("row") == row
                            and r.get("bochner_n") == n_b
                            and r.get("rigorous_dual_LB") is not None
                        ):
                            hit = True; break
                    if not hit: return False
            return True
        elif cand["kind"] == "lp_run_mside_bochner":
            ns_M = cand["params"].get("mside_bochner_ns",
                                      [cand["params"].get("mside_bochner_n", 0)])
            ns_b = cand["params"].get("bochner_ns",
                                      [cand["params"].get("bochner_n", 0)])
            for row in cand["params"].get("rows", []):
                for n_M in ns_M:
                    for n_b in ns_b:
                        hit = False
                        for r in done["results"]:
                            if (
                                r.get("kind") == "lp_run_mside_bochner"
                                and r.get("N") == cand["params"]["N"]
                                and r.get("T") == cand["params"]["T"]
                                and r.get("R") == cand["params"]["R"]
                                and r.get("row") == row
                                and r.get("mside_bochner_n") == n_M
                                and r.get("bochner_n", 0) == n_b
                            ):
                                hit = True; break
                        if not hit: return False
            return True
        elif cand["kind"] == "lp_run":
            for row in cand["params"].get("rows", []):
                for use_T5p in (False, True):
                    hit = False
                    for r in done["results"]:
                        if (
                            r.get("N") == cand["params"]["N"]
                            and r.get("T") == cand["params"]["T"]
                            and r.get("R") == cand["params"]["R"]
                            and r.get("row") == row
                            and r.get("use_T5p") == use_T5p
                            and r.get("bochner_n", 0) == 0
                        ):
                            hit = True; break
                    if not hit: return False
            return True
        return False  # unknown kinds always 'not done' so they're skipped via the kind check below

    exp = None
    for cand in pending:
        if cand["kind"] in ("infra", "alpha_sweep"):
            continue  # cron-Claude must handle these manually
        if not is_done(cand):
            exp = cand
            break

    if exp is None:
        log("All queued experiments appear already done — nothing to run.")
        return 0

    log(f"Running experiment {exp['id']}: {exp['description']}")

    if dry_run:
        log("(dry-run; not solving)")
        return 0

    if exp["kind"] not in ("diagnose", "scale_run", "lp_run", "lp_run_bochner",
                           "lp_run_bochner_sweep", "lp_run_bochner_dual",
                           "lp_run_mside_bochner"):
        log(f"Skipping — kind '{exp['kind']}' needs cron-Claude code extension.")
        return 0

    white = import_white_module()
    rows = exp["params"].get("rows", list(POINTS_BY_ROW.keys()))

    if exp["kind"] == "lp_run_mside_bochner":
        ns_M = exp["params"].get("mside_bochner_ns",
                                 [exp["params"].get("mside_bochner_n", 0)])
        ns_b = exp["params"].get("bochner_ns",
                                 [exp["params"].get("bochner_n", 0)])
        also_T5p = exp["params"].get("use_T5p", False)
        new_results = run_mside_bochner(
            exp["params"]["N"], exp["params"]["T"], exp["params"]["R"],
            rows, ns_M, white, bochner_ns=ns_b, also_T5p=also_T5p,
        )
    elif exp["kind"] == "lp_run_bochner_dual":
        ns = exp["params"].get("bochner_ns", [exp["params"].get("bochner_n", 0)])
        also_T5p = exp["params"].get("use_T5p", False)
        new_results = run_bochner_dual(
            exp["params"]["N"], exp["params"]["T"], exp["params"]["R"],
            rows, ns, white, also_T5p=also_T5p,
        )
    elif exp["kind"] in ("lp_run_bochner", "lp_run_bochner_sweep"):
        ns = exp["params"].get("bochner_ns", [exp["params"].get("bochner_n", 0)])
        also_T5p = exp["params"].get("use_T5p", False)
        new_results = run_bochner(
            exp["params"]["N"], exp["params"]["T"], exp["params"]["R"],
            rows, ns, white, also_T5p=also_T5p,
        )
    else:
        use_T5p_vals = (False, True)
        if exp.get("use_T5p_only"):
            use_T5p_vals = (True,)
        new_results = run_simple_t5p(
            exp["params"]["N"], exp["params"]["T"], exp["params"]["R"],
            rows, use_T5p_vals, white,
        )

    done["results"].extend(new_results)
    save_json(DONE, done)

    # Print summary for the cron Claude to read.
    log("--- new results ---")
    for r in new_results:
        if "error" in r:
            log(f"  ERROR  N={r['N']} T={r['T']} {r['row']} T5p={r['use_T5p']}: {r['error']}")
        elif r.get("kind") == "lp_run_mside_bochner":
            log(f"  N={r['N']} T={r['T']} R={r['R']} {r['row']} "
                f"f-side n={r.get('bochner_n', 0)} M-side n={r.get('mside_bochner_n')}: "
                f"Ω*={r.get('value', 'NA'):.10g} ({r.get('status', 'NA')}, {r.get('time', 0):.1f}s)")
        elif r.get("kind") == "lp_run_bochner_dual":
            rlb = r.get("rigorous_dual_LB")
            lg = r.get("last_gap")
            rlb_s = f"{rlb:.10g}" if rlb is not None else "None"
            lg_s = f"{lg:.2e}" if lg is not None else "None"
            log(f"  N={r['N']} T={r['T']} R={r['R']} {r['row']} n_b={r.get('bochner_n')}: "
                f"Ω*={r['value']:.10g} rigorousLB={rlb_s} last_gap={lg_s} "
                f"({r['status']}, {r['time']:.1f}s)")
        else:
            log(f"  N={r['N']} T={r['T']} R={r['R']} {r['row']} T5p={r['use_T5p']}: "
                f"Ω*={r['value']:.7f} ({r['status']}, {r['time']:.1f}s)")
    log(f"Total experiments done so far: {len(done['results'])}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        log(f"FATAL {type(e).__name__}: {e}")
        log(traceback.format_exc())
        sys.exit(1)
