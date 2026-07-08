# Reproducing the results

This walks a fresh reader from clone → environment → the headline numbers.
Verification convention: **independent re-implementations agreeing to 10+
digits** and `rigorous_dual_LB = value − last_gap` (CLARABEL dual extraction),
not unit tests.

## 0. Environment

```bash
cd /path/to/Erdos
python3 -m venv .venv && source .venv/bin/activate
pip install cvxpy clarabel numpy scipy mpmath sympy
# arbitrary-precision cross-checks (optional): lp_research_state/bin/sdpa_gmp
```

All commands below run from `lp_research_state/code` with the repo `.venv`
python; scripts do `sys.path.insert(0, '.')`.

## 1. Core single-row smoke (~30 s) — sanity that the program builds & solves

```bash
cd lp_research_state/code
python3 -c "
import sys; sys.path.insert(0,'.')
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
import cvxpy as cp
Omega,w,v,c,d,eps,dlt,cons = build_problem(10000,4000,10, 0.004,0.004, 0.3875,0.3875, -0.02,0.02, bochner_n=20)
print(solve_with_dual_extraction(cp.Problem(cp.Minimize(Omega),cons))['rigorous_dual_LB'])  # >= 0.379653
"
```

## 2. Core-region headline µ ≥ 0.3802973 (corrected-tail) / 0.380284 (conservative)

Full Phase-5 cover at production config (12 centers; ~30 min, ~8 GB RAM):

```bash
python3 path_b_with_polymoment.py --N 20000 --T 4000 --R 10 \
        --bochner_n 40 --pm_k_max 20 --hankel_n 6
```

Persisted result: `../parallel_results/cde_phase5_corrected_tail.json`
(rigorous LB `0.3802973`, binding witness `cde_n30_iter3` at h≈0.0039, p≈0.3922).
The conservative `primal − 1e-5` anchor convention gives `0.380284`. The poly-moment
tail bound must be the **fixed** version (analytic remainder; `poly_moment.py:even_moment_tail_bound`)
— this is load-bearing for rigor (see `findings.md`, 2026-05-22 entry).

## 3. Full-space promotion µ ≥ 0.380284 (PRO-38, verified)

The augmented dual cover clears the core value over all 18 of White's Table-2
outside regions. Pure evaluation of saved duals — **no SDP solves**.

```bash
# per-region floors (12 core anchors only — shows the gate)
python3 _fullspace_eval.py            # -> ../parallel_results/fullspace_stage1.json

# recompute with the FULL union of 121 centers (core + stage2 + halo + promotion)
python3 _fs_recompute.py              # -> ../parallel_results/fullspace_promote_final.json
# expect: independently_certified_floor = 0.3802838, binding region = core

# DECISIVE rigor check of the tightest region R16 (adaptive subdivision)
python3 _fs_certify_R16.py
# expect: raw cover never < 0.3804026 anywhere; tight-box cover_min_lb >= 0.3802838
```

**Why adaptive subdivision is required:** the plain single-grid `cover_min_over_box`
gives ~0.289 on White's wide boxes — an `eps_grid = L_max·half_diag` artifact
(`L_max≈7.7`). Subdividing drives `eps_grid → 0` and recovers the true cover
infimum (`grid_min`), which is ≥ target everywhere. Verification of record:
`../FULLSPACE_VERIFICATION.md`. **Caveat:** load-bearing on the fresh R16/R17
poly-moment centers and the poly-moment cuts; thin margins (R16 +1.2e-4).
Hardening (N≥24000 re-solve + Farkas certs) tracked in Linear PRO-44.

## 4. White's 8→4 correction (PRO-43) — neutral; bound unchanged

Constraints 5.6/5.7 RHS coefficient is now the parameter `mside_sin_coeff`
(default `4.0` = White-corrected; pass `8.0` to reproduce the old behavior):

```bash
python3 _white_corr_row4.py    # compares coeff 8 vs 4 at row 4; delta ~ -3.6e-8 (neutral)
```

Analysis: `../WHITE_EMAIL_CORRECTION.md`. Verdict: neutral at the binding center;
the old `8` was conservative (valid-but-looser), never an overclaim.

## 5. Upper bound (Together, for the gap)

`µ ≤ 0.380871`, verified independently to `0.3808703106…`
(`together_loader.py` / `together_diagnostic.py`). Open gap
`[0.380284, 0.380871]`, width ≈ 5.87 × 10⁻⁴.

## Cross-checks

- Bochner constraint: `bochner.py` vs `bochner_independent.py` (bit-for-bit).
- Ellipse extension: `path_b_analytical.py` / `path_b_rigorous.py` / `path_b_independent.py` (agree 10+ digits).
- Precision spot-check: `sdpa_gmp_wrapper.py` (SDPA-GMP, feasibility error ~10⁻⁷⁵).
