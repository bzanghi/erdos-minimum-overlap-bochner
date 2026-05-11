# Fourier and SDP convention in `white_full_convex.py`

Reference document for the Together-as-primal SDP diagnostic (Task 2 of
the 12-task plan; depends on Task 1 from commit `688eb76`).

All line references are to `lp_research_state/code/white_full_convex.py`
as of commit `688eb76b48cc00216061f18da6fa4873fedbbc4c` (HEAD at the
time of writing).

---

## 1. Domain of `f`

`f` is defined on the domain **`[-2, 2]`** (length **`4`**, i.e. the
fundamental period of the Fourier basis `cos(πmx/2), sin(πmx/2)`).

The discretization (lines 122–134) introduces cell width

    L = 2.0 / N                          (line 122)
    j = 1, 2, …, N                        (line 123, `np.arange(1, N+1)`)

and two cvxpy variables of length `N`:

    w = cp.Variable(N)                   (line 133)
    v = cp.Variable(N)                   (line 134)

The cell layout (already documented in `together_loader.py:38–41`, which
itself cites `white_full_convex.py:99–152`) is:

| Variable        | Cell (positive index `j ∈ {1..N}`) |
|-----------------|------------------------------------|
| `w[j-1]`        | `[(j-1)L,  jL]`   ⊂  `[0, 2]`      |
| `v[j-1]`        | `[-jL, -(j-1)L]`  ⊂  `[-2, 0]`     |

So `w` covers the positive half `[0, 2]` and `v` covers the negative
half `[-2, 0]`, with `2N` cells of width `L = 2/N` tiling `[-2, 2]`.
The total domain length is `N · L + N · L = 4` ✓.

The mass constraint at line 142 is

    L * cp.sum(w + v) == 1               (line 142)

i.e. the trapezoidal-style sum `L · Σ(w_j + v_j) = ∫_{-2}^{2} f = 1`,
so `f` is an L¹-normalized density on `[-2, 2]`.

The pointwise constraint at line 141 is

    w >= 0, v >= 0, w <= Omega, v <= Omega, Omega <= 1

so `0 ≤ f ≤ Ω ≤ 1` cell-by-cell.

---

## 2. Variables

| Name       | Type                  | Length | Meaning                                                                  |
|------------|-----------------------|--------|--------------------------------------------------------------------------|
| `Omega`    | `cp.Variable()`       | 1      | SDP objective; pointwise upper bound on `f`; equal to `2·M̂(0)` (§5).   |
| `w[j-1]`   | `cp.Variable(N)`      | `N`    | Cell-average value of `f` on `[(j-1)L, jL]`. Cited: lines 133, 141–142. |
| `v[j-1]`   | `cp.Variable(N)`      | `N`    | Cell-average value of `f` on `[-jL, -(j-1)L]`. Same citation.            |
| `c[k]`     | `cp.Variable(T)`      | `T`    | Cosine Fourier modes of `f`. `f̂(k) = (c[k-1] − i·d[k-1])/2` for `k ≥ 1` (line 230). `c[0]` is parameter-bounded (see below). |
| `d[k]`    | `cp.Variable(T)`      | `T`    | Sine (imaginary) Fourier modes of `f`. Same citation as `c`.             |
| `eps[r-1]` | `cp.Variable(R)`      | `R`    | Tail-truncation slack for the cosine relation at odd `m = 2r-1` (lines 164–168, bound at line 195). |
| `dlt[r-1]` | `cp.Variable(R)`      | `R`    | Tail-truncation slack for the sine relation at odd `m = 2r-1` (lines 169–172, bound at line 196). |

Special parameter-bounds on the lowest Fourier coefficients (line 201):

    c[0] ∈ [p1, p2],   d[0] ∈ [q1, q2]

These bounds are *parameters of the problem* — they pin `(c_1, d_1)` to
a small rectangle in `(p, q)`-space, which is how White (2023) sweeps
the 7-row residual region. They are NOT free Fourier coefficients in
the usual sense; the rest of `c[1..T-1]` and `d[1..T-1]` are free except
for the box bound `|c|, |d| ≤ 2/π` (line 199) and the energy bound
`Σ(c² + d²) ≤ 1/2` (line 200).

A few additional shape constraints worth flagging:

- First moment of `f`: `L² · Σ_j (j·w_j − (j-1)·v_j) ≥ h1`  (line 151).
- Second moment: `L³ · Σ_j (j-1)² (w_j + v_j) ≤ 2/3 + h2²/2` (line 152).
- A 2-frequency cosine inequality: line 205 ((5.13) in White).

---

## 3. Inverse formula (`f` from its Fourier coefficients)

The Fourier basis is `cos(π m x / 2)` and `sin(π m x / 2)` for integer
`m`, period `4`, on the domain `[-2, 2]` (see lines 22, 26, 47, 51 and
the `np.cos(np.pi * m * x / 2)` / `np.sin(np.pi * m * x / 2)` evaluations
throughout `cos_cell_bounds_*` / `sin_cell_bounds_*`).

The mapping from Fourier modes to `f` consistent with line 230's
`f̂(0) = 1/2`, `f̂(k) = (c[k-1] − i·d[k-1])/2` is

    f(x) = 1/2  +  Σ_{m=1}^{∞} [ c[m-1] · cos(π m x / 2)
                                  + d[m-1] · sin(π m x / 2) ]   for x ∈ [-2, 2].

Reasoning:
- `f̂(0) = (1/(period)) ∫ f = ∫_{-2}^{2} f / 4 = 1/4`. But line 230
  STATES `f̂(0) = 1/2`. The reconciliation is that line 230's `f̂(0)`
  uses White's UNNORMALIZED convention where `f̂(m) := (1/2) ∫_{-2}^{2}
  f(x) e^{-i π m x / 2} dx` (the `1/2` prefactor corresponds to
  half-period normalization). Under that convention, `f̂(0) = (1/2)·∫f =
  1/2` matches line 142's `∫f = 1`.
- The two real Fourier modes `(c_m, d_m)` reconstruct `f` via the real
  expansion `1/2 + Σ_m (c_m cos(πmx/2) + d_m sin(πmx/2))`. The `1/2`
  constant term is what the Fourier-0 mode contributes back to `f` under
  the convention above (`f̂(0) e^{0} = 1/2`).
- `f̂(k) = (c[k-1] − i·d[k-1])/2` (line 230): the factor of `1/2` in
  this formula is the same half-period-normalization factor — i.e.
  White packages the real cos/sin coefficients without the `1/2`, then
  the complex coefficient absorbs it.

So the EXACT inverse formula used in the code is

    f(x) = 1/2 + Σ_{m=1}^{T} [ c[m-1] cos(π m x / 2) + d[m-1] sin(π m x / 2) ]
                                                                            (+ tail)

with the tail `Σ_{m > T}` absorbed (for the odd-`m` relations) into the
slack variables `eps` and `dlt` and bounded by `tail_bound_eps`,
`tail_bound_delta` (lines 91–96, applied at lines 195–196).

The cell-Fourier relations (lines 156–197) tie `(w, v)` to `(c, d, eps,
dlt)` through cell integrals of `cos(πmx/2)` and `sin(πmx/2)`, bounded
via `cos_cell_bounds_exact` / `sin_cell_bounds_exact` (lines 22–69) or
their Lipschitz analogues. The detailed encoding splits even and odd
`m`:

- Even `m = 2·half` (lines 158–160): the cell relation directly equates
  `a_m`, `b_m` to `c[half-1]/2`, `d[half-1]/2`.
- Odd `m` (lines 161–172): introduces tail-truncation slacks `eps`,
  `dlt` plus the `odd_coeff_factors(m, T)` weighting (defined lines
  84–88) that captures how odd-`m` cell integrals couple to the FULL
  sequence `c[k]`, `d[k]`.

---

## 4. Constraints (load-bearing subset)

- Pointwise / mass (lines 141–142): `0 ≤ w, v ≤ Ω ≤ 1`, `L·Σ(w+v) = 1`.
- Parameter bounds (line 201): `c[0] ∈ [p1, p2]`, `d[0] ∈ [q1, q2]`.
- Box / energy on Fourier modes (lines 199–200): `|c|, |d| ≤ 2/π`,
  `Σ c² + Σ d² ≤ 1/2`.
- Moment constraints (lines 151–152): first/second moments of `f` tied
  to parameters `h1, h2`.
- Tail slacks (lines 195–196): `|eps[r-1]| ≤ tail_bound_eps(2r-1, T)`,
  `|dlt[r-1]| ≤ tail_bound_delta(2r-1, T)`.
- Cell ↔ Fourier inequality for each `m ∈ {1..2R}` (lines 176–190):
  combines `cos_cell_bounds` / `sin_cell_bounds` with `a_m`, `b_m` to
  enforce White's `(5.10)`–`(5.13)`-style overlap inequalities.
- Bochner-PSD (lines 233–258, gated by `bochner_n > 0`): builds the
  Hermitian Toeplitz `M_n(f) = [f̂(j-k)]_{j,k=0..n}` with `f̂(0) = 1/2`,
  `f̂(k) = (c[k-1] − i·d[k-1])/2`, and imposes both `M_n(f) ⪰ 0` and
  `M_n(1−f) ⪰ 0` (sign flip on off-diagonals). Real-form embedding
  `[[Re, -Im],[Im, Re]]` at line 257.

(The code also has optional tightenings `T3`, `T5`, `T5'`, several
M-side Bochner variants, and a Lasserre level-2 lift — see lines
207–337. These are gated by their respective flags and are not part of
the default convention.)

---

## 5. The objective `Ω` — what it actually represents

The CVXPY problem is `cp.Minimize(Omega)` (line 360). At the SDP level
`Omega` is a single scalar variable that enters the problem in three
roles:

1. **Pointwise upper bound on `f`** (line 141): `w ≤ Ω, v ≤ Ω ≤ 1`.
   So every cell value of `f` is at most `Ω`.
2. **Implicit upper bound on the autocorrelation peak** via the per-`m`
   cell ↔ Fourier inequalities (lines 176–190). These inequalities
   collectively encode that
   `sup_t M(t) ≤ Ω/2` with `M(t) = ∫_{-2}^{2} f(x) f(x+t) dx`. In
   particular `Ω/2 = M̂(0)` is stated explicitly in the comment at line
   261 ("`M̂(0) = Ω/2`").
3. **The variable being minimized**, so the SDP's optimal value
   `Ω*` is the SMALLEST scalar that simultaneously satisfies (1) the
   cell-wise upper bound and (2) the integrated-Fourier upper bounds on
   `sup_t M(t)`, given all the other constraints (parameter rectangle,
   Bochner-PSD, etc.).

So `Ω` is NOT literally `sup_t ∫ f f(·+t) dx` evaluated at a fixed `f`
in closed form. It is a SDP-feasible upper bound on that supremum,
tightened by every constraint in the system. At a candidate primal
point `(w, v, c, d, eps, dlt)`, the smallest feasible `Ω` is

    Ω_feas(w, v, c, d, eps, dlt)  =
        max( max_j w_j,
             max_j v_j,
             max over all per-m inequalities (lines 176–182, 189–190)
                    of the implied LB on Ω ).

The SDP then minimizes this over all feasible primal points. This is
the natural reading of "Ω = primal objective" in `white_full_convex.py`.

The headline rigorous bound (`Ω* ≥ 0.379544` from the Bochner-augmented
program) is the dual lower bound on this Ω-minimization — see
`dual_extractor.py` and the project root `CLAUDE.md`.

---

## 6. The even-`f` case (`assume_even=True`)

Lines 144–150 add the constraints

    d == 0,    dlt == 0,    v == w                                  (line 150)

when `assume_even=True` (default `False` at line 111). This enforces:

- `d[k] = 0` for all `k ∈ {0..T-1}` (sine Fourier coefficients vanish).
- All sine tail slacks vanish.
- `v_j = w_j` for all `j ∈ {1..N}` (mirror cells: the cell value on
  `[-jL, -(j-1)L]` equals the cell value on `[(j-1)L, jL]`).

Together these enforce `f(-x) = f(x)`. The resulting bound is
CONDITIONAL on the (open) conjecture that the µ-optimal `f` is even
(White 2023, §4 and §6). Without `assume_even`, White's program allows
any non-even `f` on `[-2, 2]`.

---

## 7. White's `Ω` vs Together's `M` — **load-bearing distinction**

This is the most consequential note in this document for downstream
diagnostic tasks (Tasks 6 and 12 especially). White's `Ω` and
Together's `M(h)` are DIFFERENT functionals — not just different
discretizations.

### 7.1 The two functionals

White's program minimizes `Ω`, which (as discussed in §5) is an SDP
upper bound on

    M_W(t)  :=  ∫_{-2}^{2} f(x) f(x+t) dx          (autocorrelation of f
                                                    on the half-period
                                                    domain [-2, 2])
    sup_t M_W(t)  ≤  Ω/2.

Together's program minimizes

    M(h)  :=  max over real shifts k of  ∫_0^2 h(x) (1 - h(x+k)) dx
                                              (zero extension off [0, 2];
                                               Haugland 2016 formulation)

as documented in `together_loader.py:21–35`.

The two suprema coincide WITH THE ERDŐS CONSTANT `µ` at their
respective optimizers. That is, `inf_f sup_t M_W(t) = µ = inf_h M(h)`.
But at a fixed admissible input they give DIFFERENT NUMBERS.

### 7.2 Numerics for Together's `h*`

From Task 1 (`together_loader.py:55–73`, `verify_white_embedding`):

| Input                                                              | White's `Ω` (autocorr) | Together's `M(h)` |
|--------------------------------------------------------------------|------------------------|-------------------|
| Together's `h*` embedded as `f(x) = h(\|x\|)/2` on `[-2, 2]` (even) | `0.387337`             | `0.380871`        |
| Together's `h*` embedded as `f(x) = h(x)/2` on `[0, 2]`, `0` else  | `0.774675`             | `0.380871`        |

The 0.38734 vs 0.38087 gap is NOT a bug. It is the gap between
`sup_t M_W(t)` and `sup_k M_T(k)` evaluated at this specific `h*`.

(The 0.774675 entry for the asymmetric embedding is roughly `2 ×
0.38734`, reflecting that the asymmetric embedding has only half the
mass on the positive side — its autocorrelation peak at `t = 0` is
`∫ f² = ∫_0^2 (h/2)² ≈ 0.387/2`-scale but the constraint `∫f = 1` is
violated; that asymmetric variant is NOT admissible as a White primal
unless the embedding is renormalized differently — see
`together_loader.py:97–114`.)

### 7.3 Implication for downstream tasks

- **"`Ω(f*) ≈ 0.380871`" is NOT the right success criterion for Task 6
  or any successor.** Together's published value `0.380871` is `M(h*)`,
  not `Ω(f*)`.
- The correct internal comparison: solve White's SDP at the same
  parameter setting to get `Ω*_SDP`. The current rigorous LB is
  `Ω*_SDP ≥ 0.379544` (Bochner-augmented program; see project root
  `CLAUDE.md` and `erdos_lower_bound_research_note.md`). Then
  `Ω(Together_h*)` vs `Ω*_SDP` quantifies how much SLACK Together's
  point has UNDER WHITE'S OBJECTIVE.
- Task 1 already showed `Ω(f_even from h*) = 0.387337 > 0.379544`. So
  **Together's `h*` (embedded into White's space, evenly) is a STRICTLY
  WORSE Ω-point than the SDP-optimal `f`**. The structural lever
  Together exploits to push their bound to `0.380871` is, evidently,
  on the M-side, not the Ω-side.
- This is itself a research finding: it suggests that constraining the
  SDP to look like Together's `f` would NOT tighten White's LB on
  `µ` via the current Ω-functional, since Together's `f` is sub-optimal
  for that functional. Any lever for improving White's bound by
  importing Together's structure must operate on the M-side (e.g.
  through the M-side Bochner constraints — `mside_bochner.py`,
  `mside_bochner_schur.py`, `mside_via_lasserre.py`) or on a different
  cost functional.

### 7.4 Open question for Task 12

Do we need to extend the SDP machinery to optimize over a different
functional (closer in spirit to Together's `M`)? Or is the diagnosis
simply "Together's `M`-minimizer is a worse-than-optimal Ω-point, and
the gap is real but inert under our current objective"? Task 6 must
report `Ω(f*_direct)` and `Ω(f*_even)` BOTH, alongside `Ω*_SDP`, and
Task 12 must interpret the comparison in light of §7.3 above.

---

## Appendix A: Conventions cheat-sheet

| Symbol         | Definition                                                              | Source line  |
|----------------|-------------------------------------------------------------------------|--------------|
| `f`            | unknown density on `[-2, 2]`, period 4 extension                        | (implicit)   |
| `L`            | cell width `2/N`                                                        | 122          |
| `N`            | number of cells per half-domain                                         | 99           |
| `T`            | number of Fourier modes carried explicitly (`c[0..T-1]`, `d[0..T-1]`)   | 99           |
| `R`            | number of `m`-inequalities (`m = 1..2R`), and tail-slack length         | 99           |
| `Ω`            | objective: SDP UB on `sup_t M_W(t)`, also pointwise UB on `f`            | 132, 141     |
| `f̂(0)`        | `1/2` (White's half-period normalization; ties to `∫f = 1`)             | 230          |
| `f̂(k)`        | `(c[k-1] − i·d[k-1]) / 2` for `k ≥ 1`                                   | 230          |
| Basis          | `cos(π m x / 2)`, `sin(π m x / 2)`, period 4                            | 22–69        |
| `assume_even`  | sets `d = 0`, `dlt = 0`, `v = w` (conditional on even-`f` conjecture)   | 111, 149–150 |
| `bochner_n`    | adds PSD on `M_n(f)` and `M_n(1-f)`; THE rigorous augmentation          | 108, 233–258 |
