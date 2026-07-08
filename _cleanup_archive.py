"""One-off repo cleanup (PRO-15): move clearly-stale files to _archive/ dirs.
SAFE: reference-checks every code/JSON candidate against files that REMAIN; skips
anything imported or path-referenced. .md memos have no code refs (safe to move).
Reports moved/skipped; caller then smoke-tests imports."""
import os, re, shutil, glob
ROOT = "/Users/benzanghi/Documents/Claude/Projects/Erdos"
CODE = f"{ROOT}/lp_research_state/code"
PR   = f"{ROOT}/lp_research_state/parallel_results"
for d in (f"{CODE}/_archive", f"{PR}/_archive", f"{ROOT}/docs/archive"):
    os.makedirs(d, exist_ok=True)

# ---- KEEP at repo root; every other top-level *.md -> docs/archive ----
KEEP_MD = {"README.md", "REPRODUCE.md", "REPO_AUDIT.md", "CLAUDE.md",
           "erdos_lower_bound_research_note.md", "SUBMISSION_CHECKLIST.md",
           "PROGRESS_AND_SIGNIFICANCE.md"}

# ---- stale scratch _*.py to archive (import-checked before moving) ----
ARCHIVE_PY_PREFIX = ("_lever_", "_redteam_", "_pro4_", "_pro26_", "_pro28_", "_pro29_",
    "_lasserre_scan", "_compare_lasserre", "_run_lasserre3_test", "_run_row1",
    "_run_row4_schur", "_run_one_rigorous", "_path_b_quick", "_check_status",
    "_quick_compare", "_sanity_closed_form", "_sat_Mn", "_brute_force_Mn_extended",
    "_lifted_density_compare", "_together_projection_independent")

# ---- stale JSONs to archive (path-ref-checked) ----
ARCHIVE_JSON = {f"row{i}_bochner.json" for i in range(1, 8)} | {
    "row1_test_10k_2500.json", "row1_test_10k_3000.json", "row1_test_10k_4000.json",
    "row5_bochner_N1000.json", "row5_bochner_N1500.json", "row5_bochner_n10.json",
    "row5_bochner_n15.json", "lasserre2_path_b_rigorous.json", "lasserre2_path_b_summary.json",
    "lasserre3.json", "mside_lasserre.json", "mside_schur_results.json", "path_b.json",
    "path_b_N2000.json", "path_b_closed_form.json", "path_b_rigorous.json",
    "phase4b_plus_T5p.json", "phase5_T5p_full.json", "symmetric_conditional.json",
    "symmetric_high_n.json", "sdpa_gmp.json", "push_high_n_summary.json",
    "cde_phase1_row8_h0_p394.json", "_r7_feas_map.json", "_r7_feas_map2.json",
    "_r7_floor_pruned.json", "_r9_probe_cache.json"}

def remaining_py_text(exclude):
    txt = []
    for f in glob.glob(f"{CODE}/*.py"):
        if os.path.basename(f) in exclude:
            continue
        try: txt.append(open(f, encoding="utf-8", errors="ignore").read())
        except Exception: pass
    return "\n".join(txt)

moved, skipped = [], []

# --- .md (safe) ---
for f in glob.glob(f"{ROOT}/*.md"):
    b = os.path.basename(f)
    if b not in KEEP_MD:
        shutil.move(f, f"{ROOT}/docs/archive/{b}"); moved.append(f"md/{b}")

# --- _*.py (import-checked) ---
py_cands = [f for f in glob.glob(f"{CODE}/_*.py")
            if os.path.basename(f).startswith(ARCHIVE_PY_PREFIX)]
cand_names = {os.path.basename(f) for f in py_cands}
rem = remaining_py_text(exclude=cand_names)
for f in py_cands:
    name = os.path.basename(f)[:-3]  # strip .py
    if re.search(rf"\b(import|from)\s+{re.escape(name)}\b", rem):
        skipped.append(f"py/{name} (imported)"); continue
    shutil.move(f, f"{CODE}/_archive/{name}.py"); moved.append(f"py/{name}")

# --- JSON (path-ref-checked) ---
rem_all = remaining_py_text(exclude=set())
for b in sorted(ARCHIVE_JSON):
    p = f"{PR}/{b}"
    if not os.path.exists(p):
        continue
    if b in rem_all:
        skipped.append(f"json/{b} (referenced)"); continue
    shutil.move(p, f"{PR}/_archive/{b}"); moved.append(f"json/{b}")

print(f"MOVED {len(moved)}:")
for m in moved: print("  +", m)
print(f"\nSKIPPED {len(skipped)} (kept — referenced):")
for s in skipped: print("  -", s)
print(f"\nremaining root *.md: {len(glob.glob(f'{ROOT}/*.md'))}  "
      f"code _*.py: {len(glob.glob(f'{CODE}/_*.py'))}  "
      f"parallel_results *.json: {len(glob.glob(f'{PR}/*.json'))}")
