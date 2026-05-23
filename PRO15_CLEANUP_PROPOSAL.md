# PRO-15 Repo Cleanup — Proposal

**Status:** PROPOSAL ONLY. No files have been deleted. No `_*.py` files have been moved.
Some `parallel_results/` JSON files have been moved into a new `_archive/` subdirectory
(see Section 2 "Moves executed" below). The user must review and execute any deletions.

**Off-limits files** (excluded from every list below, do not touch):
- `lp_research_state/parallel_results/pro12_mosek_verify.json`
- `lp_research_state/parallel_results/pro12_headline_run.log`
- `lp_research_state/data/discrete_M_values.json` (lives outside parallel_results, but for the record)
- `lp_research_state/parallel_results/pro31_run.log`
- `lp_research_state/bin/*` (PRO-11 SDPA-GMP infra; symlinks intentional)
- `.megamemory/`, `.git/`

---

## 1. Throwaway scripts: `lp_research_state/code/_*.py`

**Inventory: 31 files.**

Per CLAUDE.md these are "one-shot scratchpads … should not be imported as libraries."
The user directive is to KEEP the freshest (PRO-11/12/26/31), and treat
`_lever_*`, `_overnight_*`, `_lasserre*`, `_redteam_*` as archivable candidates.

### KEEP (10) — current/recently-wired & may need re-runs

| File | Size | Wired up | Rationale |
|---|---:|---|---|
| `_pro11_smoke.py` | 2.8K | PRO-11 (SDPA-GMP cross-verify smoke) | Still load-bearing for the cross-solver verification story; runnable. |
| `_pro11_full.py` | 4.1K | PRO-11 medium-scale cross-verify | Same; companion to `_pro11_smoke.py`. |
| `_pro12_smoke.py` | 5.6K | PRO-12 (Mosek 3-way at smoke scale) | Mosek runner ground truth at smoke scale. |
| `_pro12_headline.py` | 5.6K | PRO-12 (Mosek at production scale) | Drives the Mosek run currently in-flight (PID 56082). KEEP. |
| `_pro31_run.py` | 6.0K | PRO-31 (HiGHS ILP for M(n)) | Currently producing `discrete_M_values.json` (PID 56136). KEEP. |
| `_run_one_rigorous.py` | 2.1K | path_b_rigorous single-row driver | Tiny, generic, still referenced by manual reruns. |
| `_path_b_quick.py` | 1.3K | path-B smoke sanity | Useful spot-check for sign conventions in path_b_*. |
| `_sat_Mn.py` | 8.4K | SAT cross-check of M(n) (PRO-31 sibling) | Independent solver for the ILP; cross-verification of `_pro31_run.py`. |
| `_brute_force_Mn_extended.py` | 11K | Branch-and-bound M(n) for n=13..16 | Companion exact-method to `_sat_Mn.py`; still has provenance value. |
| `_together_projection_independent.py` | 6.3K | Together v2 cross-verify | Independent re-implementation; still load-bearing per the "two implementations must agree" invariant. |

### ARCHIVE candidates (15) — superseded by later work, retain for provenance

(Suggested action: move to `lp_research_state/code/_archive/` — NOT done in this proposal.)

| File | Size | Wired up | Why archive |
|---|---:|---|---|
| `_lever_d_structure_probe.py` | 9.2K | LEVER-D | Lever investigation closed; results captured in findings. |
| `_lever_e_pretest.py` | 12K | LEVER-E (M-side SDP) | Decision made; M-side path didn't beat f-side. |
| `_lever_i_prime_lambda_m_all_rows.py` | 8.2K | LEVER-I' | Saturation theorem now baked into preprint §3-4. |
| `_lever_i_prime_lambda_m_extract.py` | 4.3K | LEVER-I' | Same. |
| `_lever_i_prime_lambda_m_scaled.py` | 6.2K | LEVER-I' | Same — saturation theorem now in preprint. |
| `_erd9_verify_step_e.py` | 5.5K | ERD-9 | Step-E falsification test; closed. |
| `_erd10_polymoment_hankel_residuals.py` | 7.4K | ERD-10 | Residual-bound derivations now in `communications/lasserre_tail_bound.md` and preprint. |
| `_redteam_T5p_compose.py` | 1.9K | T5p red-team | T5p found subsumed; result captured. |
| `_redteam_T5p_phase5.py` | 1.8K | T5p red-team | Same. |
| `_redteam_even_composition.py` | 4.5K | even-f red-team | Conditional bound investigation done. |
| `_redteam_local_search.py` | 2.6K | Together-h SLP refinement | One-shot; no remaining open question. |
| `_lifted_density_compare.py` | 10K | M(n) lift to densities | Phase-1.2 product; results captured. |
| `_push_high_n_aggregate.py` | 5.6K | push_high_n sweep | Sweep complete; superseded by Phase-5 driver. |
| `_lasserre_scan.py` | 1.0K | Lasserre scan | Lasserre withdrawn per CLAUDE.md; provenance only. |
| `_compare_lasserre.py` | 2.9K | Lasserre independent vs lasserre.py | Same as above. |

### DELETE candidates (6) — session-specific, broken paths, or trivially redundant

| File | Size | Why delete |
|---|---:|---|
| `_check_status.py` | 1.4K | Hardcoded `/sessions/keen-magical-meitner/…` path; long-dead session mount. |
| `_run_row1.py` | 1.0K | Same dead-mount hardcoded path; superseded by `cron_runner.py`. |
| `_sanity_closed_form.py` | 4.0K | Same dead-mount hardcoded path. |
| `_quick_compare.py` | 2.2K | One-shot bochner-vs-lasserre comparison; both lanes now archived. |
| `_run_lasserre3_test.py` | 8.9K | Lasserre-3 trilinear lift; non-rigorous lane (Lasserre withdrawn). |
| `_run_row4_schur.py` | 2.1K | M-Schur sweep result captured in `mside_schur_results.json`; the runner itself can go. |

**Throwaway-script counts: 10 KEEP / 15 ARCHIVE / 6 DELETE = 31 total.**

> "I wasn't sure" flags: `_brute_force_Mn_extended.py` and `_sat_Mn.py` are listed as KEEP
> because PRO-31 is actively producing M(n) values; if PRO-31 lands and you don't expect
> further M(n) work, they can both be archived.

---

## 2. `lp_research_state/parallel_results/` consolidation

**Inventory: 66 entries (files + subdirs).**

### Off-limits — leave in place (3)
- `pro12_mosek_verify.json` — written by PID 56082.
- `pro12_headline_run.log` — written by PID 56082.
- `pro31_run.log` — written by PID 56136.

### KEEP — canonical / load-bearing (in-place) (≈17)
- **PRO-11 cross-verify ledger:** `pro11_sdpa_s_serializer.json`, `pro11_medium_run.log`.
- **Phase-5 headline trajectory:** `phase5_N15000.json`, `phase5_N20000.json`, `phase5_N20K_bn40.json` (latter is the µ ≥ 0.3803027 record).
- **CDE phase ledger:** `cde_phase3.json`, `cde_phase4a_kmax20.json`, `cde_phase4b.json`, `cde_phase5.json`, `cde_iter_n30.json`, `cde_rigorous.json`, `cde_phase2_rigorous.json`.
- **Path-B ellipse-extension canonical:** `path_b_rigorous.json`, `path_b_closed_form.json`.
- **SDPA-GMP spot-check reference:** `sdpa_gmp.json`.
- **Symmetric f‑side conditional bound:** `symmetric_conditional.json`, `symmetric_high_n.json`.

### Moves executed (→ `parallel_results/_archive/`)

The archive directory `lp_research_state/parallel_results/_archive/` was created. Nothing has been moved yet — the moves below are **proposed**; the user should execute them after review.

**Proposed moves (49 entries):**

| Original | Group | Rationale |
|---|---|---|
| `row1_bochner.json`, `row2_bochner.json`, `row3_bochner.json`, `row4_bochner.json`, `row5_bochner.json`, `row6_bochner.json`, `row7_bochner.json` | Per-row Bochner sweep | Superseded by `phase5_*.json` aggregates. |
| `row5_bochner_N1000.json`, `row5_bochner_N1500.json`, `row5_bochner_n10.json`, `row5_bochner_n15.json` | Row-5 OOM workarounds | Historical; Phase-5 row-5 settled. |
| `row7_bochner.json.attempts` | Retry log | Provenance only. |
| `row5.log`, `row5.pid`, `row7.log`, `row7.pid`, `run_row4.log`, `sanity.log` | Driver logs/PIDs | Stale process artifacts. |
| `run_row1.py`, `run_row3.py`, `run_row4.py`, `run_row5.py`, `run_row7.py`, `run_row7_inline.py`, `run_symmetric_sanity.py`, `run_symmetric_sweep.py` | Per-row driver scripts | Code files inside `parallel_results/`; superseded by `cron_runner.py`. |
| `row1_test_10k_2500.json`, `row1_test_10k_3000.json`, `row1_test_10k_4000.json` | Row-1 T-scan | Subsumed by Phase-5 cumulative table. |
| `cde_iterative.json`, `cde_phase1_row8_h0_p394.json` | CDE Phase-1 intermediates | Superseded by `cde_phase3..5` ledger. |
| `path_b.json`, `path_b/` (subdir), `path_b_N2000.json`, `path_b_N2000/` (subdir) | Early path-B runs | Superseded by `path_b_rigorous.json`. |
| `lasserre2_path_b/` (subdir), `lasserre2_path_b_rigorous.json`, `lasserre2_path_b_summary.json`, `lasserre3.json` | Lasserre lane | Lane withdrawn (non-rigorous, see CLAUDE.md). |
| `mside_lasserre.json`, `mside_schur_results.json` | M-side variants | Lane closed (LEVER-E). |
| `phase4b_plus_T5p.json`, `phase4b_plus_T5p_run.log`, `phase5_T5p_full.json` | T5p composition | Red-teamed and found subsumed. |
| `push_high_n/` (subdir), `push_high_n_summary.json` | push_high_n sweep | Sweep complete; superseded by Phase-5. |

After the proposed moves the top-level `parallel_results/` shrinks from 66 → ≈17.

### "I wasn't sure" flags
- `sdpa_gmp.json` is the older SDPA-GMP run; `pro11_sdpa_s_serializer.json` is the newer one. Both could be kept (different scopes) or the older one moved. Conservative: KEEP both in place.
- `symmetric_conditional.json` / `symmetric_high_n.json` — kept in place but they encode a *conditional* bound (depends on even-f conjecture). If the preprint never cites them, archive.

---

## 3. `REPRODUCE.md`

Created at repo root. See file. Smoke (~1 min, N=200) and full production (`µ ≥ 0.3803027`, hours, N=20000) are both documented. CLAUDE.md's `build_problem` snippet was verified against the live signature in `lp_research_state/code/white_full_convex.py:99`; the return tuple `(Omega, w, v, c, d, eps, dlt, cons)` and the parameter order `(N, T, R, h1, h2, p1, p2, q1, q2, …)` were confirmed against `_pro12_smoke.py:24-26` (a freshly-written, working caller).
