# Repo audit & archive plan (PRO-15)

Snapshot 2026-06-01. Repo has accumulated **54 top-level `.md`**, **51 library
modules**, **92 `_*.py` scratch scripts**, **69 `parallel_results/*.json`**.
This catalog separates load-bearing/canonical files from archivable detritus.
**Nothing is moved by this document** — it is the plan; archiving is deferred
to explicit confirmation (several files are path-coupled; see Protected).

## ⛔ Protected — imported as machinery, do NOT archive

These underscore files are `import`ed by other code (verified via grep) and/or
are the certification-of-record for the verified full-space bound:

```
_fullspace_eval.py        (imported by _fs_recompute, _fs_recon, _fs_standard_floor, _fs_certify_R16)
_fs_recompute.py          (imported by _fs_standard_floor, _fs_certify_R16; globs fullspace_promote_R*.json)
_fs_recon.py              _verify_cover_dualext.py   _fullspace_stage2_solve.py   _r9_verify.py
_fs_certify_R16.py        _fs_standard_floor.py      _harden_loadbearing.py (active, PRO-44)
```

Canonical library stack (keep): `white_full_convex.py`, `white_full_convex_exact.py`,
`path_b_analytical.py`, `path_b_rigorous.py`, `path_b_independent.py`,
`poly_moment.py`, `bochner.py`, `bochner_independent.py`, `dual_extractor.py`,
`path_b_with_polymoment.py`, `iterate_centers_pm.py`, `cde_evaluate.py`,
`together_loader.py`, `sdpa_gmp_wrapper.py`, `pslq_hunt.py`, `arxiv_search.py`,
`wolfram_alpha.py`.

## `parallel_results/` — canonical vs archivable

**KEEP (load-bearing / referenced by code or the verified result):**
`phase5_N20K_bn40_dualext.json` (the 12 core anchors), `cde_phase5_corrected_tail.json`
(corrected-tail headline), `fullspace_promote_*.json`, `fullspace_promote_final.json`,
`fullspace_recon.json`, `fullspace_stage2_centers.json`, `fullspace_stage2_halo_centers.json`,
`verify_region_R*.json`, `verify_scope.json`, `harden_loadbearing_N24K.json` (incoming).

**KEEP (canonical history):** `cde_phase{1,2,3,4a,4b,5}*.json`, `cde_iter*.json`,
`cde_rigorous.json`, `phase5_N{15000,20000}.json`, `phase5_N20K_bn40.json`,
`fullspace_stage1{,_gate}.json`, `fullspace_rigor.json`.

**ARCHIVE candidates** (early/superseded snapshots, not referenced):
`row{1..7}_bochner*.json`, `row1_test_10k_*.json`, `row5_bochner_*.json`,
`path_b{,_N2000,_closed_form}.json`, `lasserre*.json` / `mside_*.json` (withdrawn line),
`_r7_feas_map*.json`, `_r9_probe_cache.json`, `_r7_floor_pruned.json`,
`phase5_T5p_full.json`, `phase4b_plus_T5p.json`, `symmetric_*.json`, `sdpa_gmp.json`,
`push_high_n_summary.json`. → ~24 files to `parallel_results/_archive/`.

## `_*.py` scratch — group actions (92 files)

| group | examples | action |
|-------|----------|--------|
| full-space machinery (this result) | `_fs_*`, `_fullspace_*`, `_verify_region`/`_eval_r*`/`_verify_r*`/`_promote_R*`, `_r9_*` | **KEEP** (provenance of PRO-38; some imported). Re-evaluate after preprint. |
| White correction | `_white_corr_*` | KEEP (provenance of PRO-43) |
| withdrawn/old levers | `_lever_*`, `_redteam_*`, `_pro4_*`, `_pro26_*`, `_pro28_*`, `_pro29_*`, `_lasserre_scan`, `_compare_lasserre`, `_run_lasserre3_test` | **ARCHIVE** → `code/_archive/` (~30 files) |
| one-shot diagnostics | `_check_status`, `_quick_compare`, `_diag_r17`, `_map_r17_dip`, `_sanity_closed_form`, `_sat_Mn`, `_run_row1`, `_run_row4_schur`, `_run_one_rigorous`, `_path_b_quick`, `_brute_force_Mn_extended`, `_lifted_density_compare`, `_together_projection_independent` | **ARCHIVE** → `code/_archive/` (~13 files) |
| keep-for-now | `_pro14_verifier`, `_pro24_richardson`, `_verify_mu`, `_verify_perturbation`, `_erd9_verify_step_e`, `_erd10_*` | KEEP (referenced by findings / reusable verifiers) |

## Top-level `.md` — consolidate (54 files)

KEEP at root: `README.md`, `REPRODUCE.md`, `CLAUDE.md`, `erdos_lower_bound_research_note.md`,
`FULLSPACE_VERIFICATION.md`-class memos, `SUBMISSION_CHECKLIST.md`, `PROGRESS_AND_SIGNIFICANCE.md`.
ARCHIVE → `docs/archive/`: the `LEVER_*` series (~16), `OUT_OF_BOX_*` (5),
`SESSION_*` (4), `PRO*_*.md` per-issue memos (~8) once their results are folded
into `findings.md` / the preprint. → ~37 files.

## Safe archive execution (deferred — run on confirmation)

```bash
cd lp_research_state
mkdir -p parallel_results/_archive code/_archive
# move ONLY the ARCHIVE-candidate lists above (never the Protected/KEEP sets)
# then: grep -rl "<moved-name>" . to confirm no live reference broke
```

Recommend doing this **after** the preprint is drafted (so memos are still
handy) and the PRO-44 hardening lands. Net effect when executed: `_*.py` count
~92 → ~35 (−62%), `parallel_results` split into canonical + `_archive/`,
top-level `.md` ~54 → ~17.
