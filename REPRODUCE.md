# REPRODUCE — Erdős minimum overlap µ ≥ 0.3803027

This document walks a cold reader from a clean clone to **µ ≥ 0.3803027**,
the headline rigorous lower bound established by this repo (PRO-21 Phase 5,
N = 20 000, bochner_n = 40).

> **Critical invariant.** A single-point SDP solve at a row center is **NOT**
> a rigorous bound on µ. Per `CLAUDE.md`, the unconditional bound requires
> White's §5.1 ellipse-extension argument — implemented in
> `lp_research_state/code/path_b_*.py`. The fast smoke recipe in §4 only
> demonstrates that the SDP encoding is correctly wired up. The full
> production recipe in §5 is what actually reproduces µ ≥ 0.3803027.

---

## 1. Prerequisites

- Python ≥ 3.10
- Linux or macOS; ≥ 16 GB RAM recommended for the full run (peak ~8 GB).
- No GPU required. CLARABEL is the default solver and ships with cvxpy.
- (Optional, for cross-verification): Mosek with a valid license; SDPA-GMP
  built from source. Neither is needed for the headline number.

## 2. Setup (≈ 3 minutes)

```bash
git clone <repo-url> erdos
cd erdos
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install cvxpy numpy scipy
```

That is the entire dependency surface for the headline reproduction. If you
also want the optional cross-checks add `pip install mosek pysat highspy`.

Quick sanity check that the import surface works:

```bash
.venv/bin/python -c "
import sys; sys.path.insert(0, 'lp_research_state/code')
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction
print('OK — encoding + dual extractor importable')
"
```

## 3. What's in the box

The core SDP-construction module is `lp_research_state/code/white_full_convex.py`.
Its `build_problem(...)` constructs White's §5 convex program plus all the
augmentations (Bochner-PSD, polynomial-moment, etc.) toggled via keyword args.

- `dual_extractor.solve_with_dual_extraction(prob)` — runs the SDP, parses
  CLARABEL's iteration log, and returns the **rigorous** dual lower bound as
  `result["rigorous_dual_LB"] = reported_value − last_gap`. This recovers
  ~+1 × 10⁻⁴ over the naive `prob.value` since CLARABEL's
  `optimal_inaccurate` flag is a labeling artifact (true gaps are ~10⁻⁷).
- `path_b_rigorous.py`, `path_b_with_polymoment.py`, `path_b_independent.py`
  — three implementations of the ellipse-extension argument. They agree to
  10+ significant digits.
- `iterate_centers_pm.py` — the Phase-5 driver. Iteratively adds dual-feasibility
  centers at the current binding point of the envelope until the rigorous
  lower bound stops moving.

## 4. Fast smoke (≈ 1 minute) — does the encoding work?

This solves **row 4** (the binding row at White's Table-3 ellipse centers)
at small scale and verifies the rigorous-dual-extraction pattern returns a
sensible value. It is **NOT** a bound on µ — it is a wiring check.

```bash
.venv/bin/python -c "
import sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, 'lp_research_state/code')

import cvxpy as cp
from white_full_convex import build_problem
from dual_extractor import solve_with_dual_extraction

# Row 4 (binding): (h, p, q) = (0.004, 0.3875, ±0.02), bochner_n=10
Omega, w, v, c, d, eps, dlt, cons = build_problem(
    1000, 400, 10,           # N, T, R
    0.004, 0.004,            # h1, h2
    0.3875, 0.3875,          # p1, p2
    -0.02, 0.02,             # q1, q2
    bochner_n=10,
)
prob = cp.Problem(cp.Minimize(Omega), cons)
res = solve_with_dual_extraction(prob)
print(f'rigorous_dual_LB = {res[\"rigorous_dual_LB\"]:.10f}')
print(f'reported_value   = {res[\"reported_value\"]:.10f}')
print(f'last_gap         = {res[\"last_gap\"]:.2e}')
print(f'status           = {res[\"status\"]}   time = {res[\"time\"]:.2f}s')
"
```

Expected: `rigorous_dual_LB ≈ 0.378–0.379`, status `optimal_inaccurate`,
`last_gap < 1e-6`, total wall time < 60 s on a laptop.

For a stronger single-point sanity (matches the value `≥ 0.379653` quoted in
`CLAUDE.md`), bump to `N=10000, T=4000, bochner_n=20`. Wall time ≈ 8 min,
peak memory ≈ 2 GB. Still a single-point bound — not µ.

## 5. Full production reproduction — µ ≥ 0.3803027 (hours)

The headline number comes from `path_b_with_polymoment.py`'s Phase-5 cover
iteration at `N=20000, T=4000, bochner_n=40, pm_k_max=20, hankel_n=6`. It

1. solves the augmented SDP at each of the 12 cover centers (7 White Table-3
   centers + 5 CDE-discovered centers from prior phases),
2. assembles the per-center duals into the ellipse-extension envelope,
3. takes the max over centers (rigorous, since each ellipse is independently
   dual-feasible), subtracts a uniform 1 × 10⁻⁶ margin and a 2.17 × 10⁻⁶
   Lipschitz envelope correction, and
4. prints the resulting unconditional lower bound on µ.

**Resource budget.** Wall time ≈ 60–90 min on a recent x86 laptop with
≥ 16 GB RAM. Peak RAM ≈ 8 GB. **Row 5 is memory-heavier than the rest** —
if you OOM on row 5 specifically, the workaround is to step `bochner_n` down
to 15 on that row only (see CLAUDE.md note).

```bash
cd lp_research_state/code
python3 path_b_with_polymoment.py \
    --N 20000 --T 4000 --R 10 \
    --bochner_n 40 --pm_k_max 20 --hankel_n 6
```

Output appends to `lp_research_state/parallel_results/phase5_N20K_bn40.json`.
The last printed line includes the rigorous lower bound; the JSON dump
includes the per-center dual data and the binding witness.

Expected outcome:
- `grid_min ≈ 0.3803049`
- `eps_grid ≈ 2.17 × 10⁻⁶` (Lipschitz envelope correction)
- post-margin LB = **0.3803027** (binding witness: `cde_n30_iter3` at
  `(h, p) ≈ (0.00385, 0.39222)`)

That is the headline. The committed reference JSON
`phase5_N20K_bn40.json` lets you diff your run line-by-line.

## 6. Cross-verification (optional, recommended)

The project's epistemic policy is that any rigorous claim must be confirmed
by an independent re-implementation agreeing to ≥ 10 digits.

- **SDPA-GMP at small N.** `lp_research_state/code/_pro11_full.py` runs the
  same problem through SDPA-GMP (multi-precision, independent solver) at
  smoke and medium scale. Output: `parallel_results/pro11_sdpa_s_serializer.json`.
- **Mosek at production scale.** `lp_research_state/code/_pro12_headline.py`
  re-solves row 4 at `(N=10000, T=4000, bn=20)` via Mosek, comparing
  primal/dual objectives. Output: `parallel_results/pro12_mosek_verify.json`.
- **Path-B independent.** `path_b_independent.py` is an independently-written
  ellipse-extension; `_run_one_rigorous.py 4` runs it on row 4. Agreement
  with `path_b_rigorous.py` must be ≥ 10 digits.

## 7. What you can NOT do with this recipe

- **Push past 0.3803027 with the same technique stack.** Per
  `erdos_lower_bound_research_note.md` and `lp_research_state/findings.md`, the
  Bochner-PSD + polynomial-moment + Hankel-PSD + ellipse-extension stack
  saturates near this value at currently-tractable SDP scale. The framework
  ceiling is ≈ 0.380558 (see `PRO6_COMPLEMENTARITY_PROOF.md`).
- **Quote any Lasserre-augmented bound as rigorous.** Lasserre is documented
  as withdrawn (`communications/lasserre_tail_bound.md`).
- **Quote a single-row solve as a bound on µ.** The ellipse-extension step
  in `path_b_*.py` is non-optional; per-row solves only bound the LP optimum
  at one parameter point.

## 8. Reference numbers

| Quantity | Value | Source |
|---|---|---|
| White (2023) lower bound | µ ≥ 0.379005 | arXiv:2201.05704 |
| This repo, prior headline (PRO-3, N=20K, bn=30) | µ ≥ 0.3802994 | `phase5_N20000.json` |
| **This repo, current headline (PRO-21, N=20K, bn=40)** | **µ ≥ 0.3803027** | `phase5_N20K_bn40.json` |
| Together (2026) upper bound | µ ≤ 0.380871 | github.com/togethercomputer/erdos-minimum-overlap |
| Open gap | 5.68 × 10⁻⁴ | |

## 9. Where to look next

- `lp_research_state/findings.md` — rolling lab notebook; the leading paragraph is always the
  freshest result.
- `communications/preprint_draft.tex` — the v2 preprint, including the
  saturation theorem (the methodological core).
- `CLAUDE.md` — orientation for new contributors / agents.
- `PRO15_CLEANUP_PROPOSAL.md` — provenance map for the auxiliary scripts and
  per-row result files.
