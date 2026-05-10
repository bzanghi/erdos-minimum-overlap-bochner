# Together-as-Primal SDP Diagnostic Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Plug Together Computer's piecewise-constant minimizer `f*` (which achieves µ ≤ 0.380871) into our SDP encoding to read off per-constraint slacks, then write a structural diagnosis memo selecting the next attack vector (Lukács SOS / combinatorial M(n) push / structural restriction theorem).

**Architecture:** Three new files under `lp_research_state/code/`: a loader/verifier for Together's `f*` (`together_loader.py`), a diagnostic driver that projects `f*` into White's Fourier basis and evaluates every constraint at it (`together_diagnostic.py`), and an independent re-implementation of the Fourier projection for cross-verification (`_together_projection_independent.py`). Outputs are JSON/NPZ artifacts plus the structural memo `TOGETHER_DIAGNOSTIC.md` at the repo root.

**Tech stack:** Python 3 / numpy / scipy / cvxpy 1.8.2 / CLARABEL (already installed in `.venv/`). No new dependencies. Verification by independent re-implementation to ≥10 digits (project standard).

**Verification model:** This is research code producing numerical findings, not unit-tested business logic. The TDD discipline adapts to: every quantitative output is computed two ways, and a step fails if the two computations disagree past their declared tolerance. The "expected value" in each verify step is a concrete numerical bound the engineer can check.

---

### Task 1: Discover and parse Together's `f*`

**Files:**
- Create: `lp_research_state/code/together_loader.py`
- Create: `lp_research_state/data/together_f_star.json`

- [ ] **Step 1.1: Investigate Together's repo**

Run:
```bash
cd /tmp && git clone --depth 1 https://github.com/togethercomputer/erdos-minimum-overlap.git together_repo
ls together_repo/
find together_repo -type f \( -name "*.py" -o -name "*.json" -o -name "*.csv" -o -name "*.md" -o -name "*.txt" -o -name "*.npy" \) | head -40
```
Expected: A small repo. Look for any artifact containing 600 numbers (the step values) or a script that produces them. Read the top-level README to understand their formulation (interval, normalization, value convention).

Record findings inline in a comment at the top of `together_loader.py`: domain (likely [0, 1] or [-1/2, 1/2]), number of steps, file path, claimed value.

- [ ] **Step 1.2: Document Together's problem formulation vs. White's**

Read the relevant sections of Together's README/paper. In `together_loader.py`, write a module docstring stating:
- Together's domain for `f`
- Together's overlap functional (the quantity they minimize)
- Together's claimed value (expected: 0.380871)
- Any normalization difference vs. White (factor of 2, mirror, etc.)

If the domain or normalization differs, document the explicit transformation `(domain_together, f_together) → (domain_white, f_white)`.

- [ ] **Step 1.3: Write the loader**

Create `lp_research_state/code/together_loader.py`:
```python
"""
Load Together Computer's piecewise-constant minimizer f* for the Erdős
minimum overlap problem.

Source: github.com/togethercomputer/erdos-minimum-overlap
Together's claimed bound: µ ≤ 0.380871
Together's domain/normalization: [DOCUMENTED IN STEP 1.2]
White's domain/normalization: f : [0, 2] → [0, 1] symmetric around x=1
    with f̂(0) = 1/2, f̂(k) = (c_k - i d_k)/2 — see white_full_convex.py:230
Transformation Together → White: [DOCUMENTED IN STEP 1.2]
"""
import json
import numpy as np
from pathlib import Path

DATA_PATH = Path(__file__).parent.parent / "data" / "together_f_star.json"


def load_together_raw():
    """Load Together's f* in their native representation.

    Returns:
        breakpoints: np.ndarray of shape (n_steps + 1,) — endpoints of each piece
        values:      np.ndarray of shape (n_steps,)     — f-value on each piece
        domain:      tuple (a, b) — domain in Together's convention
        meta:        dict — source provenance
    """
    # [IMPLEMENTATION BASED ON STEP 1.1 FINDINGS]
    raise NotImplementedError("Fill in once Step 1.1 identifies the artifact format")


def to_white_convention(breakpoints, values, domain):
    """Transform Together's f* into White's domain/normalization.

    Returns:
        wb: np.ndarray — breakpoints in White's domain
        wv: np.ndarray — values in White's convention
    """
    # [IMPLEMENTATION BASED ON STEP 1.2 DERIVATION]
    raise NotImplementedError("Fill in once Step 1.2 documents the transformation")
```
The two `NotImplementedError`s will be replaced in Step 1.4 once the format and transformation are known.

- [ ] **Step 1.4: Fill in `load_together_raw` and `to_white_convention`**

Implement both functions based on Step 1.1's discovery and Step 1.2's derivation. Save the parsed artifact:
```python
def save_canonical(out_path: Path = DATA_PATH):
    bp, vals, dom, meta = load_together_raw()
    wb, wv = to_white_convention(bp, vals, dom)
    out = {
        "together": {"breakpoints": bp.tolist(), "values": vals.tolist(), "domain": list(dom), "meta": meta},
        "white": {"breakpoints": wb.tolist(), "values": wv.tolist()},
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
```

Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_loader import save_canonical; save_canonical()"
ls -la lp_research_state/data/together_f_star.json
```
Expected: file exists, several hundred to several thousand bytes.

- [ ] **Step 1.5: Verify Together's claimed value directly from `f*`**

Add to `together_loader.py`:
```python
def compute_overlap_from_f(breakpoints, values):
    """Compute the Erdős overlap M(f) directly from a step function f.

    Uses Together's formulation: M(f) = max over t of ∫ f(x) f(x + t) dx
    (or whatever formulation is documented in Step 1.2).

    Returns: scalar — Together's claimed value if f* is their certificate.
    """
    # [IMPLEMENTATION using the formulation from Step 1.2]
    ...

def verify_together_value():
    bp, vals, _, _ = load_together_raw()
    mu = compute_overlap_from_f(bp, vals)
    print(f"Recomputed Together's value: {mu:.10f}")
    print(f"Together's claimed value:    0.380871")
    print(f"Difference:                  {abs(mu - 0.380871):.2e}")
    return mu
```

Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_loader import verify_together_value; verify_together_value()"
```
Expected: `|recomputed − 0.380871| < 1e-5`. **If this fails, stop and report — we don't understand their formulation. Do not proceed.**

- [ ] **Step 1.6: Commit**

```bash
cd /Users/benzanghi/Documents/Claude/Projects/Erdos
git add lp_research_state/code/together_loader.py lp_research_state/data/together_f_star.json
git commit -m "Together diagnostic: load and verify f* (claimed mu <= 0.380871)"
```

---

### Task 2: Document White's Fourier convention precisely

**Files:**
- Create: `lp_research_state/code/_fourier_convention_notes.md`

- [ ] **Step 2.1: Extract the convention from `white_full_convex.py`**

Read [white_full_convex.py:99-260](lp_research_state/code/white_full_convex.py). Identify:
- The variables `c` and `d` (length T+1 each).
- The relation between `(c, d)` and the Fourier coefficients of `f`.
- The domain of `f` (which interval, what symmetry assumption).
- The objective `Ω` as a function of `(c, d)` and any auxiliary `(w, v)`.

Key line is white_full_convex.py:230:
```
# M_n(f) := [f̂(j-k)]_{j,k=0..n} is PSD.  With f̂(0)=1/2, f̂(k)=(c_k - id_k)/2.
```
So `c_k = 2 Re f̂(k)` and `d_k = -2 Im f̂(k)` for `k ≥ 1`. Confirm `c[0]` and `d[0]`'s meaning (they are parameters with separate bounds `[p1, p2]`, `[q1, q2]` per white_full_convex.py:201 — NOT the same as `c_k` for `k ≥ 1`. Document carefully).

- [ ] **Step 2.2: Write the convention reference**

Create `lp_research_state/code/_fourier_convention_notes.md`:
```markdown
# Fourier convention in white_full_convex.py

Per white_full_convex.py:201, 230, and the assume_even branch:

## Domain
f is defined on [DOMAIN] (extract from the code — the cell width L and N cells determine this).

## Variables
- c[0], d[0]: parameters of the family (NOT Fourier coefficients of f).
- c[1..T], d[1..T]: encode f̂(k) = (c[k] - i·d[k]) / 2 for k = 1..T.
- f̂(0) = 1/2 (normalized).
- f̂(-k) = conj(f̂(k)) since f is real.

## Inverse formula (Fourier series for f)
f(x) = 1 + Σ_{k=1}^T [c[k] cos(πkx/?) + d[k] sin(πkx/?)] + tail
       (extract exact frequency from code: search for "cos(np.pi * m * x" terms)

## Objective Ω
[extract the exact CVXPY expression from build_problem]

## Even-f case
assume_even=True forces d[k] = 0 ∀k and v_j = w_j, restricting to even f.
```
Fill in every `[BRACKETED]` placeholder with concrete content from the source. This document will be referenced by every subsequent task.

- [ ] **Step 2.3: Commit**

```bash
git add lp_research_state/code/_fourier_convention_notes.md
git commit -m "Together diagnostic: document White's Fourier convention"
```

---

### Task 3: Implement primary Fourier projection of `f*`

**Files:**
- Create: `lp_research_state/code/together_diagnostic.py` (skeleton + projection)

- [ ] **Step 3.1: Closed-form Fourier coefficients of a step function**

For a step function `f(x) = v_i` on `[b_i, b_{i+1}]`, the cosine coefficient at frequency `ω_k = πk/L` (with L from Step 2.1) is:
```
c_k = (2/L) Σ_i v_i · [sin(ω_k b_{i+1}) - sin(ω_k b_i)] / ω_k
```
and similarly for `d_k` with `-(cos(...) - cos(...))/ω_k`. Derive the exact formulas from Step 2.1's convention; **do not guess**.

Add to `together_diagnostic.py`:
```python
"""
Diagnostic: evaluate every constraint in our SDP at Together's f*.

Outputs:
  lp_research_state/data/together_f_star_fourier.npz  — (c, d) projection
  lp_research_state/data/together_diagnostic_results.json  — constraint slacks
  lp_research_state/data/together_gap_function.npz  — f̃(x) - f*(x) samples
"""
import json
import numpy as np
from pathlib import Path
from .together_loader import load_together_raw, to_white_convention

# See _fourier_convention_notes.md for the exact convention.
# L = ... (from Step 2.1)


def project_step_function(breakpoints: np.ndarray, values: np.ndarray, T: int) -> tuple[np.ndarray, np.ndarray]:
    """Project a step function (in White's domain) onto the first T cos/sin Fourier modes.

    Returns:
        c: np.ndarray of shape (T+1,) — White's c[0..T]; c[0] is the constant term.
        d: np.ndarray of shape (T+1,) — White's d[0..T]; d[0] is the constant term.
    """
    # [CLOSED-FORM IMPLEMENTATION per Step 2.1's convention]
    ...
```

- [ ] **Step 3.2: Unit-test the projection on a known case**

Add a test function in the same file:
```python
def _test_projection_on_constant():
    """A constant function f(x) = c should give c[0] = 2c, all other coefficients = 0."""
    bp = np.array([0.0, 2.0])
    vals = np.array([0.5])
    c, d = project_step_function(bp, vals, T=10)
    # Expected per White's convention: f̂(0) = c[0]/? = 0.5 ⇒ c[0] = ?
    # (Use the EXACT relation derived in Step 2.1.)
    assert abs(c[0] - EXPECTED_C0) < 1e-12, f"c[0] = {c[0]} != expected {EXPECTED_C0}"
    assert np.allclose(c[1:], 0, atol=1e-12), f"c[1:] should be 0, got max {np.max(np.abs(c[1:]))}"
    assert np.allclose(d, 0, atol=1e-12), "d should be all 0 for a constant"
    print("✓ projection on constant function")
```
Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import _test_projection_on_constant; _test_projection_on_constant()"
```
Expected: prints `✓ projection on constant function` and exits 0. If the assertion on `c[0]` fails, the convention in Step 2.1 was misread — fix and re-run.

- [ ] **Step 3.3: Project Together's `f*` and save**

Add:
```python
def project_together_f_star(T: int = 4000) -> tuple[np.ndarray, np.ndarray]:
    bp_t, vals_t, dom, _ = load_together_raw()
    wb, wv = to_white_convention(bp_t, vals_t, dom)
    c, d = project_step_function(wb, wv, T=T)
    np.savez(
        Path(__file__).parent.parent / "data" / "together_f_star_fourier.npz",
        c=c, d=d, T=T,
    )
    return c, d
```
Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import project_together_f_star; c, d = project_together_f_star(); print('c[0:5]', c[:5]); print('|c|_max', np.max(np.abs(c))); print('|d|_max', np.max(np.abs(d)))"
```
Expected: `c[0]` matches White's convention for f̂(0) = 1/2; `|c|`, `|d|` decay roughly like 1/k for a step function.

- [ ] **Step 3.4: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py lp_research_state/data/together_f_star_fourier.npz
git commit -m "Together diagnostic: project f* into White's Fourier basis (T=4000)"
```

---

### Task 4: Independent Fourier projection (cross-check)

**Files:**
- Create: `lp_research_state/code/_together_projection_independent.py`

- [ ] **Step 4.1: Write an independent implementation using numerical quadrature**

Use a different algorithm than Task 3's closed-form: high-order Gaussian quadrature on a dense grid, or scipy.integrate.quad per coefficient. This MUST not share code with Task 3 except for the convention constants from `_fourier_convention_notes.md`.

```python
"""
Independent re-implementation of step-function → Fourier projection.

Method: per-cell numerical quadrature of f(x) cos(ω_k x) and f(x) sin(ω_k x).
Different algorithm from together_diagnostic.project_step_function (which uses
closed-form sin/cos integrals). Must agree to ≥10 digits — project rule.
"""
import numpy as np
from scipy.integrate import quad
from lp_research_state.code.together_loader import load_together_raw, to_white_convention

# Convention constants — KEEP IN SYNC WITH _fourier_convention_notes.md
L_DOMAIN = ...  # extract from Step 2.1
OMEGA = lambda k: ...  # extract from Step 2.1


def project_step_function_quad(breakpoints, values, T):
    c = np.zeros(T + 1)
    d = np.zeros(T + 1)
    # Per-cell quadrature against cos(ω_k x), sin(ω_k x)
    ...
    return c, d
```

- [ ] **Step 4.2: Cross-verify**

Add to `_together_projection_independent.py`:
```python
def cross_verify():
    bp, vals, dom, _ = load_together_raw()
    wb, wv = to_white_convention(bp, vals, dom)
    from lp_research_state.code.together_diagnostic import project_step_function
    c_fast, d_fast = project_step_function(wb, wv, T=4000)
    c_slow, d_slow = project_step_function_quad(wb, wv, T=4000)
    max_diff = max(np.max(np.abs(c_fast - c_slow)), np.max(np.abs(d_fast - d_slow)))
    print(f"Max abs diff between closed-form and quadrature: {max_diff:.3e}")
    assert max_diff < 1e-10, f"Cross-verify FAILED at {max_diff:.3e}; should be < 1e-10"
    print("✓ Fourier projections agree to ≥10 digits")
```
Run:
```bash
.venv/bin/python -c "from lp_research_state.code._together_projection_independent import cross_verify; cross_verify()"
```
Expected: prints `Max abs diff ... 1e-12` or thereabouts and `✓` line. **If it disagrees past 1e-10, stop and investigate — one implementation has a bug.**

- [ ] **Step 4.3: Commit**

```bash
git add lp_research_state/code/_together_projection_independent.py
git commit -m "Together diagnostic: independent Fourier projection cross-check (10+ digits)"
```

---

### Task 5: Rigorous truncation tail bound

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 5.1: Derive the tail bound for step functions**

A step function with `M` jumps has Fourier coefficients with magnitudes `|f̂(k)| ≤ V/k` where `V = total variation = Σ_i |v_i - v_{i-1}|`. By Parseval, `‖f - f_T‖_{L²}² = Σ_{k > T} |f̂(k)|² ≤ V² Σ_{k > T} 1/k² ≤ V² / T`. Hence `‖f - f_T‖_{L²} ≤ V/√T`.

For `L¹`: `‖f - f_T‖_{L¹} ≤ √(|domain|) · ‖f - f_T‖_{L²}` by Cauchy-Schwarz.

Confirm this derivation matches Together's `f*` shape — record the actual `V` (sum of absolute jumps) in the code.

- [ ] **Step 5.2: Implement and record**

Add to `together_diagnostic.py`:
```python
def truncation_tail_bound(breakpoints, values, T):
    """Rigorous bound on ||f - f_T||_{L²} and ||f - f_T||_{L¹} for a step function.

    Returns: dict with keys 'V', 'L2_bound', 'L1_bound', 'T'.
    """
    V = float(np.sum(np.abs(np.diff(values))))  # total variation
    domain_len = float(breakpoints[-1] - breakpoints[0])
    L2_bound = V / np.sqrt(T)
    L1_bound = np.sqrt(domain_len) * L2_bound
    return {"V": V, "L2_bound": L2_bound, "L1_bound": L1_bound, "T": T}
```

Run:
```bash
.venv/bin/python -c "
from lp_research_state.code.together_loader import load_together_raw, to_white_convention
from lp_research_state.code.together_diagnostic import truncation_tail_bound
bp_t, vals_t, dom, _ = load_together_raw()
wb, wv = to_white_convention(bp_t, vals_t, dom)
print(truncation_tail_bound(wb, wv, T=4000))
"
```
Expected: a dict with `L2_bound` roughly `10^-3` to `10^-2` (a 600-step function on a unit interval has V ≈ O(1)). Record exact value in commit message.

- [ ] **Step 5.3: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py
git commit -m "Together diagnostic: rigorous L2/L1 truncation tail bounds for f*"
```

---

### Task 6: Evaluate objective `Ω(f*)` at row 4

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 6.1: Identify row 4's parameters**

From CLAUDE.md and findings.md: row 4 is `(h_c, p_c, q_c) = (0.004, 0.3875, ±0.02)`. Look up the exact `build_problem` call signature in `path_b_rigorous.py` to mirror the row-4 invocation.

- [ ] **Step 6.2: Compute `Ω(f*)` by substitution**

Add to `together_diagnostic.py`:
```python
import cvxpy as cp
from .white_full_convex import build_problem


def evaluate_omega_at_f_star(c_proj, d_proj, row="row4", bochner_n=30,
                              N=10000, T=4000, R=10):
    """Substitute the projected f* into White's SDP and read off Ω.

    We FIX (c, d) at the projection and solve the auxiliary minimization over
    (w, v, dlt, eps, ...) only — this gives the value Ω(f*) under our encoding.
    """
    centers = {"row4": (0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02), ...}  # extract from path_b_rigorous.py
    h1, h2, p1, p2, q1, q2 = centers[row]
    Omega, w, v, c_var, d_var, eps, dlt, cons = build_problem(
        N, T, R, h1, h2, p1, p2, q1, q2, bochner_n=bochner_n,
    )
    # Pin c_var and d_var to the projection
    pin = [c_var[k] == c_proj[k] for k in range(min(T+1, len(c_proj)))]
    pin += [d_var[k] == d_proj[k] for k in range(min(T+1, len(d_proj)))]
    prob = cp.Problem(cp.Minimize(Omega), cons + pin)
    prob.solve(solver=cp.CLARABEL, verbose=False)
    return {
        "status": prob.status,
        "Omega_at_f_star": float(Omega.value) if Omega.value is not None else None,
        "infeasible": "infeasible" in (prob.status or ""),
    }
```

- [ ] **Step 6.3: Run and verify against Together's claimed value**

Run:
```bash
.venv/bin/python -u -c "
import numpy as np
from lp_research_state.code.together_diagnostic import evaluate_omega_at_f_star
data = np.load('lp_research_state/data/together_f_star_fourier.npz')
c, d = data['c'], data['d']
print(evaluate_omega_at_f_star(c, d, row='row4', bochner_n=30))
"
```
Expected outcomes — ALL are informative:
- `status='optimal'`, `Omega ≈ 0.380871`: our encoding correctly evaluates `f*`. Gap is structural.
- `status='optimal'`, `Omega` materially lower (e.g., 0.379): our objective under-encodes the true overlap. Bug or relaxation issue — flag for memo.
- `status='infeasible'`: a constraint rejects `f*`. **This is a major finding.** Move to Task 7 to find WHICH constraint.

Record the actual outcome.

- [ ] **Step 6.4: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py
git commit -m "Together diagnostic: evaluate Omega(f*) at row 4 [outcome: <actual>]"
```

---

### Task 7: Evaluate Bochner-PSD constraint at `f*`

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 7.1: Construct M_n(f*) directly**

Add to `together_diagnostic.py`:
```python
def bochner_matrix_at_f_star(c, d, n):
    """Construct the (n+1) x (n+1) Bochner moment matrix M_n(f*) from White's
    convention f̂(0)=1/2, f̂(k)=(c[k] - i d[k])/2.

    M_n[j, k] = f̂(j - k). Use f̂(-k) = conj(f̂(k)).
    Returns: complex (n+1, n+1) Hermitian matrix.
    """
    f_hat = np.empty(n + 1, dtype=complex)
    f_hat[0] = 0.5
    for k in range(1, n + 1):
        f_hat[k] = (c[k] - 1j * d[k]) / 2
    M = np.empty((n + 1, n + 1), dtype=complex)
    for j in range(n + 1):
        for k in range(n + 1):
            diff = j - k
            M[j, k] = f_hat[diff] if diff >= 0 else np.conj(f_hat[-diff])
    return M


def bochner_diagnostic(c, d, n=30):
    M_f = bochner_matrix_at_f_star(c, d, n)
    # 1 - f*: f̂_{1-f}(0) = 1 - 1/2 = 1/2, f̂_{1-f}(k) = -f̂(k) for k >= 1
    c1 = c.copy(); c1[1:] *= -1
    d1 = d.copy(); d1[1:] *= -1
    M_1mf = bochner_matrix_at_f_star(c1, d1, n)
    return {
        "lambda_min_M_n(f)":   float(np.min(np.linalg.eigvalsh(M_f))),
        "lambda_min_M_n(1-f)": float(np.min(np.linalg.eigvalsh(M_1mf))),
        "trace_M_n(f)":   float(np.real(np.trace(M_f))),
        "trace_M_n(1-f)": float(np.real(np.trace(M_1mf))),
        "n": n,
    }
```

- [ ] **Step 7.2: Cross-check against the SDP's encoding**

In `together_diagnostic.py` add a verification that the matrix above matches the matrix `build_problem` constructs internally (read white_full_convex.py:230–290 to see how it builds the matrix; build the same one and compare element-wise).

```python
def _bochner_xcheck(c, d, n=8):
    """Cross-check our M_n construction against the SDP's encoded version."""
    M_ours = bochner_matrix_at_f_star(c, d, n)
    # Construct M the way white_full_convex builds its constraint matrix
    M_white = ...  # mirror white_full_convex.py lines 230-290
    diff = np.max(np.abs(M_ours - M_white))
    assert diff < 1e-12, f"Bochner matrix mismatch: {diff:.3e}"
    print("✓ Bochner construction matches white_full_convex")
```
Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_diagnostic import _bochner_xcheck
data = np.load('lp_research_state/data/together_f_star_fourier.npz')
_bochner_xcheck(data['c'], data['d'], n=8)
"
```
Expected: prints `✓`.

- [ ] **Step 7.3: Run the diagnostic on `f*`**

Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_diagnostic import bochner_diagnostic
data = np.load('lp_research_state/data/together_f_star_fourier.npz')
print(bochner_diagnostic(data['c'], data['d'], n=30))
"
```
Expected: both `lambda_min` values should be ≥ 0 if `f*` truly satisfies `0 ≤ f ≤ 1`. A small negative value (down to ~`-1e-3`) is plausible from Fourier truncation; record actual values for the memo. A *large* negative value indicates either a bug or that Together's `f*` has features our T=4000 truncation can't represent.

- [ ] **Step 7.4: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py
git commit -m "Together diagnostic: Bochner-PSD slacks at f* (n=30)"
```

---

### Task 8: Evaluate poly-moment constraints at `f*`

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 8.1: Compute `m_{2k}(f*)` directly from the step function**

For a step function `f` on `[0, L_dom]` and `m_{2k}(f) = ∫ x^{2k} f(x) dx` (the formulation from `poly_moment.py` — check the exact integral being moment-bounded), the per-cell integral is closed form: `∫_a^b x^{2k} v dx = v · (b^{2k+1} - a^{2k+1}) / (2k+1)`.

Add to `together_diagnostic.py`:
```python
def poly_moments_at_f_star(breakpoints, values, k_list=(2, 4, 6, 8, 10, 12, 14)):
    """Directly compute m_{2k}(f*) = ∫ x^{2k} f*(x) dx for each k in k_list."""
    results = {}
    for k in k_list:
        m = 0.0
        for i in range(len(values)):
            a, b = breakpoints[i], breakpoints[i + 1]
            m += values[i] * (b**(2*k + 1) - a**(2*k + 1)) / (2*k + 1)
        results[2*k] = float(m)
    return results
```

- [ ] **Step 8.2: Look up tail bounds and compute slacks**

Read `lp_research_state/code/poly_moment.py` to find the tail bound formula `tail_k` it uses (Phase 3 uses k_max=14). For each `k`:
```python
def poly_moment_diagnostic(breakpoints, values, c, d, k_max=14, T=4000):
    from .poly_moment import compute_tail_bound  # extract actual name from poly_moment.py
    moments_direct = poly_moments_at_f_star(breakpoints, values, k_list=range(2, k_max + 1, 2))
    # Also compute m_{2k} via the Fourier representation (the formulation used in the SDP constraint)
    moments_via_fourier = ...  # use the integrals α_j^(k), β_j^(k) from poly_moment.py
    results = {}
    for k in range(2, k_max + 1, 2):
        tail = compute_tail_bound(k, T)
        direct = moments_direct[k]
        via_fourier = moments_via_fourier[k]
        agreement = abs(direct - via_fourier)
        slack = via_fourier - (-tail)  # constraint is m_k >= -tail
        results[k] = {
            "m_direct": direct, "m_via_fourier": via_fourier,
            "agreement": agreement, "tail_bound": tail, "slack": slack,
        }
    return results
```

- [ ] **Step 8.3: Run and tabulate**

Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_loader import load_together_raw, to_white_convention
from lp_research_state.code.together_diagnostic import poly_moment_diagnostic
bp_t, vals_t, dom, _ = load_together_raw()
wb, wv = to_white_convention(bp_t, vals_t, dom)
data = np.load('lp_research_state/data/together_f_star_fourier.npz')
results = poly_moment_diagnostic(wb, wv, data['c'], data['d'], k_max=14)
for k, info in results.items():
    print(f'k={k}: m={info[\"m_via_fourier\"]:.6e}, slack={info[\"slack\"]:.3e}, agreement={info[\"agreement\"]:.2e}')
"
```
Expected:
- `m_direct ≈ m_via_fourier` to ≥6 digits (cross-verify).
- All `slack` values ≥ 0 if Together's `f*` is moment-feasible.
- The k with the smallest slack is the one closest to binding — record for the memo.

- [ ] **Step 8.4: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py
git commit -m "Together diagnostic: poly-moment slacks at f* (k=2..14)"
```

---

### Task 9: Evaluate Hankel-PSD constraint at `f*`

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 9.1: Read the Hankel constraint encoding**

Read `lp_research_state/code/hankel_probe.py` and the Hankel section of `white_full_convex.py` to understand:
- What moments populate the Hankel matrix
- The matrix size `n=6` used in Phase 4B
- The expected PSD condition

- [ ] **Step 9.2: Build H(f*) and compute λ_min**

Add:
```python
def hankel_diagnostic(breakpoints, values, c, d, n=6):
    # Hankel matrix from the moments / Fourier coefs per hankel_probe.py
    H = ...  # mirror hankel_probe.py's construction
    eigs = np.linalg.eigvalsh(H)
    return {"lambda_min_H": float(eigs.min()), "lambda_max_H": float(eigs.max()), "n": n}
```

Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_loader import load_together_raw, to_white_convention
from lp_research_state.code.together_diagnostic import hankel_diagnostic
bp_t, vals_t, dom, _ = load_together_raw()
wb, wv = to_white_convention(bp_t, vals_t, dom)
data = np.load('lp_research_state/data/together_f_star_fourier.npz')
print(hankel_diagnostic(wb, wv, data['c'], data['d'], n=6))
"
```
Expected: `lambda_min_H ≥ -1e-6` if feasible; large negative is a finding.

- [ ] **Step 9.3: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py
git commit -m "Together diagnostic: Hankel-PSD slack at f* (n=6)"
```

---

### Task 10: Solve SDP at row 4 + extract dual + LP-optimal `f̃`

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`

- [ ] **Step 10.1: Solve the Phase 5 SDP at row 4**

Mirror the invocation from CLAUDE.md's "Reproducing the headline result" but with Phase 5's full constraint stack (`bochner_n=30` + poly_moment + Hankel — see `path_b_with_polymoment.py` for the exact call):

```python
def solve_row4_phase5():
    from .white_full_convex import build_problem
    from .dual_extractor import solve_with_dual_extraction
    import cvxpy as cp
    Omega, w, v, c_var, d_var, eps, dlt, cons = build_problem(
        10000, 4000, 10, 0.004, 0.004, 0.3875, 0.3875, -0.02, 0.02,
        bochner_n=30,
    )
    # Add poly_moment + Hankel (mirror path_b_with_polymoment.py)
    # ...
    prob = cp.Problem(cp.Minimize(Omega), cons)
    res = solve_with_dual_extraction(prob)
    c_opt = np.array([v.value for v in c_var])
    d_opt = np.array([v.value for v in d_var])
    return {"c_opt": c_opt, "d_opt": d_opt, "Omega_LB": res['rigorous_dual_LB'],
            "duals": [con.dual_value for con in cons]}
```

Run (takes 70–100s):
```bash
.venv/bin/python -u -c "
from lp_research_state.code.together_diagnostic import solve_row4_phase5
import numpy as np
r = solve_row4_phase5()
np.savez('lp_research_state/data/row4_phase5_primal.npz', c=r['c_opt'], d=r['d_opt'])
print('Omega_LB:', r['Omega_LB'])
"
```
Expected: `Omega_LB ≈ 0.38010` to `0.38013` (Phase 5 territory).

- [ ] **Step 10.2: Recover the LP-optimal density `f̃(x)` on a dense grid**

```python
def recover_f_tilde(c, d, x_grid):
    """Reconstruct f̃(x) from Fourier coefs using White's convention."""
    # f(x) = sum representation per Step 2.1
    f = np.full_like(x_grid, 1.0)  # the f̂(0) = 1/2 contribution times 2 (check convention)
    L = ...  # from Step 2.1
    for k in range(1, len(c)):
        f += c[k] * np.cos(np.pi * k * x_grid / L) + d[k] * np.sin(np.pi * k * x_grid / L)
    return f
```

Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_diagnostic import recover_f_tilde
data = np.load('lp_research_state/data/row4_phase5_primal.npz')
x = np.linspace(0, 2, 10000)  # adjust to White's domain
f_tilde = recover_f_tilde(data['c'], data['d'], x)
print('min f̃:', f_tilde.min(), 'max f̃:', f_tilde.max())
np.savez('lp_research_state/data/row4_f_tilde.npz', x=x, f_tilde=f_tilde)
"
```
Expected per probe.py history: `min f̃` is *very* negative (~−3.78 reported) — confirming Gibbs.

- [ ] **Step 10.3: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py lp_research_state/data/row4_phase5_primal.npz lp_research_state/data/row4_f_tilde.npz
git commit -m "Together diagnostic: solve Phase 5 SDP at row 4, recover LP-optimal f̃"
```

---

### Task 11: Compute gap function and plot

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`
- Create: `lp_research_state/data/together_gap_function.npz`
- Create: `lp_research_state/data/together_gap_plot.png`

- [ ] **Step 11.1: Sample `f*` on the same x-grid**

Add:
```python
def sample_f_star_on_grid(breakpoints, values, x_grid):
    """Evaluate the step function on a dense grid."""
    return np.array([values[max(0, min(len(values) - 1, np.searchsorted(breakpoints, xi) - 1))] for xi in x_grid])
```

- [ ] **Step 11.2: Compute gap and characterize**

```python
def gap_function(x_grid, f_tilde, f_star):
    g = f_tilde - f_star
    # Decompose: low-frequency vs high-frequency content via FFT
    G = np.fft.rfft(g)
    n = len(G)
    energy = np.abs(G) ** 2
    low_band = energy[:n // 20].sum()
    high_band = energy[n // 20:].sum()
    total = energy.sum()
    return {
        "max_abs_gap": float(np.max(np.abs(g))),
        "L2_gap": float(np.sqrt(np.mean(g ** 2))),
        "low_band_frac": float(low_band / total),
        "high_band_frac": float(high_band / total),
        "g": g,
    }
```

Run:
```bash
.venv/bin/python -c "
import numpy as np
from lp_research_state.code.together_loader import load_together_raw, to_white_convention
from lp_research_state.code.together_diagnostic import sample_f_star_on_grid, gap_function
bp_t, vals_t, dom, _ = load_together_raw()
wb, wv = to_white_convention(bp_t, vals_t, dom)
ftil = np.load('lp_research_state/data/row4_f_tilde.npz')
f_star = sample_f_star_on_grid(wb, wv, ftil['x'])
g_info = gap_function(ftil['x'], ftil['f_tilde'], f_star)
print('max|g|:', g_info['max_abs_gap'], 'L2(g):', g_info['L2_gap'])
print('low-band fraction:', g_info['low_band_frac'], 'high-band:', g_info['high_band_frac'])
np.savez('lp_research_state/data/together_gap_function.npz', x=ftil['x'], g=g_info['g'], f_tilde=ftil['f_tilde'], f_star=f_star)
"
```
Expected:
- `high_band_frac > 0.5`: Gibbs-dominated → memo recommends **A** (Lukács SOS / alt basis).
- `low_band_frac > 0.5`: structurally localized → memo recommends **D** (restriction theorem).
- Roughly equal: write up as ambiguous, default to A with caveat.

- [ ] **Step 11.3: Plot**

```python
def plot_gap(npz_path="lp_research_state/data/together_gap_function.npz",
             out="lp_research_state/data/together_gap_plot.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = np.load(npz_path)
    fig, ax = plt.subplots(3, 1, figsize=(10, 8))
    ax[0].plot(d["x"], d["f_tilde"], label="f̃ (our SDP)"); ax[0].legend(); ax[0].set_title("LP-optimal f̃")
    ax[1].plot(d["x"], d["f_star"], label="f* (Together)"); ax[1].legend(); ax[1].set_title("Together's f*")
    ax[2].plot(d["x"], d["g"], label="g = f̃ - f*"); ax[2].legend(); ax[2].set_title("Gap function")
    plt.tight_layout()
    plt.savefig(out, dpi=120)
    print(f"Wrote {out}")
```

Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import plot_gap; plot_gap()"
ls -la lp_research_state/data/together_gap_plot.png
```
Expected: PNG file exists.

- [ ] **Step 11.4: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py lp_research_state/data/together_gap_function.npz lp_research_state/data/together_gap_plot.png
git commit -m "Together diagnostic: gap function f̃ - f*, plot, frequency-band decomp"
```

---

### Task 12: Aggregate results JSON + write structural memo

**Files:**
- Modify: `lp_research_state/code/together_diagnostic.py`
- Create: `lp_research_state/data/together_diagnostic_results.json`
- Create: `TOGETHER_DIAGNOSTIC.md` (repo root)

- [ ] **Step 12.1: Aggregate every quantitative result into one JSON**

```python
def aggregate_results():
    import json
    from lp_research_state.code.together_loader import load_together_raw, to_white_convention, compute_overlap_from_f
    bp_t, vals_t, dom, meta = load_together_raw()
    wb, wv = to_white_convention(bp_t, vals_t, dom)
    data = np.load('lp_research_state/data/together_f_star_fourier.npz')
    c, d = data['c'], data['d']

    results = {
        "together_value_claimed": 0.380871,
        "together_value_recomputed": float(compute_overlap_from_f(bp_t, vals_t)),
        "truncation_tail": truncation_tail_bound(wb, wv, T=4000),
        "omega_at_f_star_row4": evaluate_omega_at_f_star(c, d, "row4", bochner_n=30),
        "bochner": bochner_diagnostic(c, d, n=30),
        "poly_moment": poly_moment_diagnostic(wb, wv, c, d, k_max=14),
        "hankel": hankel_diagnostic(wb, wv, c, d, n=6),
        # gap function summary
    }
    Path("lp_research_state/data/together_diagnostic_results.json").write_text(json.dumps(results, indent=2, default=float))
    return results
```

Run:
```bash
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import aggregate_results; import json; print(json.dumps(aggregate_results(), indent=2, default=float))"
```

- [ ] **Step 12.2: Write the structural memo**

Create `TOGETHER_DIAGNOSTIC.md` at repo root with this exact structure:

```markdown
# Together-as-Primal SDP Diagnostic

**Date:** 2026-05-10
**Predecessor:** [CDE_PHASE1_RESULT.md](CDE_PHASE1_RESULT.md) (µ ≥ 0.3801279, saturated stack)
**Together (2026):** [github.com/togethercomputer/erdos-minimum-overlap](https://github.com/togethercomputer/erdos-minimum-overlap) (µ ≤ 0.380871)

## TL;DR

[ONE PARAGRAPH: feasibility verdict + which lever the diagnostic selects + why]

## Q1: Is Together's f* feasible in our SDP encoding?
[Per-constraint table from `together_diagnostic_results.json`. Slack = ✓/×. Numbers to 4 sig figs.]

## Q2: What is Ω(f*) in our encoding?
[Value, comparison to 0.380871, comparison to our LB 0.3801279.]

## Q3: Which constraint has the most slack (over-engineered)?
[The family with the largest positive slack. Candidate to thin out.]

## Q4: Which constraint is binding (the active lever)?
[The family closest to violation. This is what to strengthen.]

## Q5: Gap function structure
[low_band_frac, high_band_frac from Task 11. Inline ref to gap plot.]

## Q6: The call
[A, C, or D, with one paragraph of justification.]

## Reproducing
```
.venv/bin/python -c "from lp_research_state.code.together_diagnostic import aggregate_results; aggregate_results()"
```

## Raw data
- `lp_research_state/data/together_f_star.json` — parsed f*
- `lp_research_state/data/together_f_star_fourier.npz` — Fourier projection
- `lp_research_state/data/together_diagnostic_results.json` — all slacks
- `lp_research_state/data/together_gap_function.npz` — f̃, f*, g samples
- `lp_research_state/data/together_gap_plot.png` — visual
```

Fill in every bracketed section with the concrete numbers from Step 12.1. Do not leave placeholders.

- [ ] **Step 12.3: Commit**

```bash
git add lp_research_state/code/together_diagnostic.py lp_research_state/data/together_diagnostic_results.json TOGETHER_DIAGNOSTIC.md
git commit -m "Together diagnostic: structural memo + aggregated results — next lever: <A/C/D>"
```

---

## Self-review (post-write check)

- **Spec coverage:** Each of the 5 protocol steps in the spec maps to tasks in this plan:
  - Spec Step 1 (fetch & verify) → Task 1
  - Spec Step 2 (Fourier projection + tail bound) → Tasks 2, 3, 4, 5
  - Spec Step 3 (evaluate constraints) → Tasks 6, 7, 8, 9
  - Spec Step 4 (dual & gap function) → Tasks 10, 11
  - Spec Step 5 (memo) → Task 12
- **Stopping criteria** from the spec are referenced inline (Step 1.5, Step 6.3).
- **Verification discipline** (≥10 digits for projection cross-check, ≥6 digits for Ω): Steps 4.2, 6.3, 8.3 enforce.
- **Placeholders:** Tasks 1.3, 2.2, 3.1 contain `[BRACKETED]` markers — these are *expected*: they get filled in during execution because the values depend on Together's actual format and White's exact convention text, which the engineer reads at execution time. Each one has a concrete source pointer (Step 1.1 findings, Step 2.1 extraction). They are placeholders inside steps, not placeholders for steps.

## Stopping & escalation

If at any point in Tasks 1–12 a verification step disagrees past its declared tolerance, STOP and produce a short writeup describing the discrepancy. Do not paper over numerical disagreements; in this project they are bugs.
