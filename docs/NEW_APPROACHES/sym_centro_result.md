# Approach ③ / PRO-49 — REAL Z/2 centrosymmetric Bochner block split

**Date:** 2026-06-06
**Lane:** the in-cvxpy realizable solver-cost speedup (THREE_APPROACHES.md §3).
**Honest framing:** a **solver-cost multiplier / enabler**, NOT a value-changer. It does
**not** pass the framework ceiling C_∞ ≈ 0.380558. Its purpose is to make the
high-`bochner_n` regime (LEVER-F3) and Approach ① cheaper on the same hardware.

## Headline

The complex-Hermitian route (`bochner_hermitian.py`) gave **ZERO** win — cvxpy's mandatory
`Complex2Real` lowers `hermitian=True` PSD straight back to the identical `2(n+1)×2(n+1)` real
embedding before any solver sees it (confirmed last iteration, `sym_reduction_result.{json,md}`).

The **real Z/2 centrosymmetric split DOES materialize the win**: CLARABEL actually receives a
PSD cone of side **`n+1` instead of `2(n+1)`** — verified directly from
`prob.get_problem_data(CLARABEL)`. The Hermitian-Toeplitz symmetry is even stronger than the
generic "two half-size blocks": the two centrosymmetric blocks **collapse to the same matrix**, so
a **single** `(n+1)×(n+1)` real PSD constraint replaces the full real form.

## The math (validated this session)

Bochner real form `RF = [[A, −B], [B, A]]`, `A = Re Mₙ` (symmetric), `B = Im Mₙ` (antisymmetric);
`RF ⪰ 0 ⟺ Mₙ = A + iB ⪰ 0`. `RF` is centrosymmetric; the orthogonal involution
`Q = (1/√2)[[I, J], [J, −I]]` (J = anti-diagonal reversal, `Q = Qᵀ = Q⁻¹`) gives
`Q·RF·Qᵀ = diag(B₁, B₂)`. For the Hermitian-Toeplitz structure `J·B = −B·J` (verified, **0.0**),
which forces

> **B₁ = B₂ = Bk := A + J·B**  (real, symmetric), and `RF ⪰ 0  ⟺  Bk ⪰ 0`,

with the `RF` spectrum equal to the `Bk` spectrum **doubled** (the U(1)/SO(2) multiplicity-2 fact).

Closed-form single block (verified bit-for-bit, **0.0** error vs numeric `A + J@B`):
```
Bk[j,k] = A_part + B_part,   0 ≤ j,k ≤ n
  A_part = 1/2                  if j==k
         = sign·(1/2)·c_{|j−k|} otherwise
  B_part = Im Mₙ[n−j, k];  with rk = (n−j)−k:
         = 0                    if rk==0
         = −sign·(1/2)·d_{|rk|} if rk>0
         = +sign·(1/2)·d_{|rk|} if rk<0
```

**Validation (2000 random trials, orders n=1..8, both signs):**

| check | result |
|---|---|
| `max |Bk − Bkᵀ|` | **0.0** (Bk symmetric → valid real PSDTriangle cone) |
| `max |J·B + B·J|` | **0.0** (JB = −BJ ⇒ blocks collapse to one) |
| `max |eig(RF) − eig(Bk)×2|` | **3.6e-15** (RF spectrum = Bk spectrum doubled) |
| `max |eig(A+iB) − eig(Bk)|` | **2.0e-15** (Bk is the complex Hermitian, real-similar) |
| closed-form Bk vs numeric `A+J@B` | **0.0** |
| PSD-consistency failures | **0 / 2000** |

## (A) Cone inventory — CLARABEL receives the half-size cone (the test the complex route FAILED)

From `prob.get_problem_data(cp.CLARABEL)` at row4, `bochner_n=20`, N=1500/T=600/R=10:

| form | PSD cone sides | A shape | A nnz |
|---|---|---|---|
| **REAL** `[[Re,−Im],[Im,Re]]` | **`[42, 42]`** (= 2(n+1)) | (12860, 5483) | 240168 |
| **CENTRO** `Bk ⪰ 0` | **`[21, 21]`** (= n+1) | (11516, 5483) | 239348 |

Side length `42 → 21` per block — a genuine **2× side reduction** (≈ **8× PSD-cone flops**,
`2·(½)³` per block), plus fewer rows/nnz. Contrast the complex-Hermitian route, which CLARABEL
received as `[42, 42]` **plus** hundreds of extra `zero`-cone tie equalities (net heavier).

## (B) Exact equivalence at the binding row4 — EXACT to machine precision

Spectra / residuals (NOT objective digits — CLARABEL solve-to-solve noise caps objective agreement
at ~7–9 digits per the prior finding). N=2000/T=800/R=10, `bochner_n=20`:

- Bochner-block spectra at **both** optima, **both** signs: `|RF spectrum − Bk spectrum doubled| =`
  **2.1e-15** (machine precision).
- The active `f≥0` block sits at the PSD boundary in BOTH forms (`min_eig ≈ 1e-9`), the same
  binding constraint.
- Objective real-vs-centro diff = **7.5e-9** — at CLARABEL's solve-to-solve noise floor, confirming
  same program (not a relaxation). Optimal `(c,d)` agree to **3.8e-6**.

**Exactness is proven IPM-independently** (eigenvalue spectra at the optimum), exactly as the prior
finding prescribed.

## (C) Win measurement — wall-time / peak-RSS factor and max bochner_n

Row4, fixed N=5000/T=2000/R=10; one child process per solve, peak RSS via
`resource.getrusage` (macOS = bytes). Full data in `sym_centro_result.json` (`win_measurements`).

| bochner_n | PSD side real → centro | wall real | wall centro | **time ×** | peak real | peak centro | **mem ×** | obj agree |
|---|---|---|---|---|---|---|---|---|
| 20 | 42 → **21** | 29.9 s | 2.8 s | 10.6× | 822 MB | 644 MB | 1.28× | 2.1e-8 |
| 30 | 62 → **31** | 39.7 s | 23.9 s | 1.66× | 1389 MB | 736 MB | 1.89× | 3.1e-8 |
| 40 | 82 → **41** | 64.8 s | 30.2 s | 2.15× | **2781 MB** | **859 MB** | **3.24×** | 3.9e-8 |

**Reading (honest):**
- **PSD side is exactly halved at every bn** — the structural win, deterministic.
- **Memory factor grows with bn** (1.28× → 1.89× → **3.24×**) as the PSD blocks come to dominate
  RSS. This is the win that matters: at bn=40 the real form needs **2.78 GB** vs centro's **0.86 GB**.
- **Time factor is noisier** (1.66×–10.6×): the bn=20 real solve happened to take an unusually large
  IPM iteration count (CLARABEL run-to-run variance), inflating that ratio. The robust large-bn
  point is **~2.15×** at bn=40 — consistent with the "2× side → ≳2× wall, ~8× PSD-cone flops"
  expectation. (Wall is dominated by per-iteration PSD linear algebra **and** iteration count; only
  the former is halved deterministically.)
- **Objective agreement 2–4e-8 across all bn** = CLARABEL solve-to-solve noise floor (the prior
  finding), NOT an encoding error. Rigorous dual LBs (dual_extractor) at bn=40 are **identical**
  (real 0.37965, centro 0.37965).

**Max bochner_n at a fixed memory budget.** The two forms produce the same-shape program apart from
the PSD cone, so the *only* thing that changes the memory wall is the PSD side. With the measured
~3.2× peak-RSS factor at bn=40, centro fits **roughly 1.7–2× the bochner_n** of the real form in the
same RAM. The real form's peak is 2.78 GB already at bn=40; centro at bn=40 is 0.86 GB. Centro-only
high-bn probe (to avoid OOM-ing the co-running L2 solve with the heavy real form):

| centro bochner_n | PSD side | peak RSS | wall | row4 value | rigorous dual LB |
|---|---|---|---|---|---|
| 60 | **61** | 1237 MB | 36.0 s | 0.3797342 | 0.37973 |
| 80 | **81** | 2729 MB | 63.6 s | 0.3797683 | 0.37977 |

**Equal-budget headline:** `real bn=40` (PSD `[82,82]`) needs **2781 MB**; `centro bn=80`
(PSD `[81,81]`) fits in **2729 MB** — *the same RAM at double the Bochner order*. The
side-length halving converts directly into ~2× the reachable `bochner_n` at fixed memory, and the
row4 diagnostic keeps climbing toward C_∞ with the extra levels
(bn=40 → 0.3796535, bn=60 → 0.3797342, bn=80 → 0.3797683), **without passing C_∞ ≈ 0.380558**.

## (D) Bound at the highest feasible bochner_n

Row4 rigorous dual LB via `dual_extractor` (matched N=5000/T=2000/R=10):

| bochner_n | row4 value (real) | row4 value (centro) | rigorous dual LB |
|---|---|---|---|
| 20 | 0.3792916 | 0.3792916 | — |
| 30 | 0.3795473 | 0.3795473 | — |
| 40 | 0.3796535 | 0.3796536 | **0.37965** (identical both forms) |

The centro highbn probe values are in the table above. As flagged, the row4 bound rises toward but
does **not** pass C_∞ ≈ 0.380558 — a cost lever / enabler, **not** a value lever. (These N=5000
values are below the production headline because production uses N=20000; the point here is the
real-vs-centro *equivalence* and the *cost*, not a new record.)

**CRITICAL caveat (repo policy):** a single-row-center SDP solve is **NOT** a rigorous bound on µ.
Converting per-row dual objectives into an unconditional µ bound requires the `path_b_*` ellipse
extension over the residual `(h,p,q)` region. These numbers are MIN-over-row-center *diagnostics*,
not a bound on µ.

## Honest conclusion

- **Exact-equivalence:** EXACT to machine precision. `RF ⪰ 0 ⟺ Bk ⪰ 0` with `RF` spectrum =
  `Bk` spectrum doubled (2.1e-15 at the row4 optimum, both signs; 3.6e-15 over 2000 random trials).
  Proven via eigenvalue spectra / residuals, NOT objective digits (CLARABEL nondeterminism caps
  objective agreement at ~7–9 digits — observed 2–4e-8).
- **The win materializes (unlike the complex route):** CLARABEL receives PSD side `n+1`, not
  `2(n+1)` — verified from `get_problem_data`. The complex-Hermitian `hermitian=True` route gave
  ZERO because cvxpy's `Complex2Real` re-expands it to the full `2(n+1)` embedding; the
  centrosymmetric block is built directly as a real `cp.bmat`, so nothing re-expands it.
- **Memory factor:** grows with bn, **3.24× at bn=40** (2.78 GB → 0.86 GB). At the real form's
  bn=40 RAM budget, centro reaches **bn=80** (same ~2.7 GB).
- **Time factor:** ~2.15× at bn=40 (noisier at small bn due to IPM iteration-count variance; the
  per-iteration PSD linear algebra is what halves deterministically).
- **Bound gain:** none — same optimum (rigorous LBs identical). Does **NOT** pass C_∞ ≈ 0.380558.
  This is a **cost multiplier / enabler** for Approach ① (the bispectrum/3-point lift, THREE_APPROACHES §1),
  exactly as flagged — not a value-changer.

So the centrosymmetric split is the **realizable in-cvxpy speedup** the complex-Hermitian negative
result pointed to: a deterministic 2× PSD-side reduction → ~2–3× memory and ~2× wall at high bn →
~2× the reachable `bochner_n`. The `bochner_hermitian.py` exact-equivalence encoding is reused with
confidence; the two predicted "half-size blocks" turn out to be identical (Toeplitz `JB = −BJ`), so a
single block does the job.

## Files (all NON-DESTRUCTIVE — `white_full_convex.py`/`bochner.py`/`_jansson_verify.py` untouched)

- `lp_research_state/code/bochner_centro.py` — the single-block centrosymmetric Bochner constraint
  `add_bochner_centro_constraint(cons, c, d, n, sign)` + numeric cross-check helper
  `make_centro_block(...)`. Self-test cross-validates against both `bochner_independent.make_real_form`
  and `bochner_hermitian.make_hermitian_matrix` (max err 3.6e-15 / 1.3e-15; 0/100 PSD-consistency
  failures; cvxpy feasibility smoke test passes).
- `lp_research_state/code/white_full_convex_centro.py` — `build_problem_centro(...)`, a thin wrapper
  delegating to the verified `build_problem(..., bochner_n=0)` and bolting on the centrosymmetric
  block (rest of the SDP byte-for-byte the verified program).
- `lp_research_state/code/_centro_check.py` — cone-inventory + exact-equivalence scratchpad.
- `lp_research_state/code/_centro_win.py` — real-vs-centro wall-time / peak-RSS / max-bochner_n
  sweep (writes the JSON `win_measurements`).
- `lp_research_state/code/_centro_highbn.py` — centro-only high-bn probe (writes `highbn_measurements`).
- `docs/NEW_APPROACHES/sym_centro_result.json` — all measurements (config, win sweep, highbn probe).
