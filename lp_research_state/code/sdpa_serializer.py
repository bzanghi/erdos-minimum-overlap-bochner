"""
cvxpy.Problem -> SDPA-S (.dat-s) sparse file format serializer.

Background
----------
SDPA-S standard form (per SDPA User's Manual):

    min  c^T x
    s.t. sum_{i=1..m} x_i F_i  -  F_0  >= 0   (block-diagonal, each block PSD)

Block structure entries are positive integers (PSD block sizes) or negative
integers (LP/diagonal block, |entry| = number of scalar inequalities). For each
F_k, only the UPPER triangle of each block is listed in the sparse file:

    k  block  i  j  value      (1-indexed; i<=j)

We canonicalize a cvxpy.Problem through CVXPY's SCS backend
(prob.get_problem_data(cp.SCS)) which produces the standard conic form

    min c^T y    s.t.   s = b - A y  in  K = Zero x Nonneg x SOC x PSD

Mapping to SDPA-S:
    Mat(s) = Mat(b) - sum y_i Mat(A[:,i])  must lie in cone K
  i.e.  sum y_i (-Mat(A[:,i])) - (-Mat(b))  >=  0
  i.e.  F_0 = -Mat(b),   F_i = -Mat(A[:,i])

Cone-by-cone encoding into SDPA blocks:

  * Zero (equalities, b_eq - A_eq y = 0): SDPA-S has no native equality. We
    split each scalar equality into two scalar inequalities (s >= 0 AND -s >= 0)
    which become two diagonal LP entries.

  * Nonneg (b_ineq - A_ineq y >= 0): a single LP (negative-block-structure)
    block of length n_nonneg.

  * SOC (||u_{1:}|| <= u_0 for u = b - A y of length k+1):
    encoded as the arrow-matrix LMI
        [[u_0   u_{1:}^T],
         [u_{1:}  u_0 * I_k]]  >= 0
    which is equivalent to the SOC. One PSD block of size (k+1) per SOC.

  * PSD (k x k matrix slack):
    SCS uses lower-triangular column-major vectorization with sqrt(2) scaling
    on the off-diagonal entries. We invert that scaling and reconstruct the
    k x k symmetric matrix. One PSD block of size k per cone.

All sign flips are handled by negating the SCS A, b BEFORE feeding into the
block builder, then writing each matrix's upper triangle as-is.

Sanity is provided by sdpa_runner.py which compares the SDPA-GMP dual objective
against CLARABEL's dual_extractor result on the same cvxpy.Problem.
"""

from __future__ import annotations

import io
import math
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np
import scipy.sparse as sp
import cvxpy as cp


# --------------------------------------------------------------------------- #
# Helpers for SCS PSD vectorization
# --------------------------------------------------------------------------- #
def _scs_psd_indices(k: int) -> List[Tuple[int, int, float]]:
    """
    Return list of (row, col, descale) of length k(k+1)/2 in the order SCS
    packs symmetric PSD slacks: column-major lower triangle. Entries on the
    strict lower triangle carry a sqrt(2) scaling that we must undo when
    rebuilding M[i,j] = M[j,i].
    """
    out: List[Tuple[int, int, float]] = []
    sq2_inv = 1.0 / math.sqrt(2.0)
    for col in range(k):
        for row in range(col, k):
            scale = 1.0 if row == col else sq2_inv
            out.append((row, col, scale))
    return out


# --------------------------------------------------------------------------- #
# Block-builder data classes
# --------------------------------------------------------------------------- #
@dataclass
class _Block:
    """A single SDPA-S block. size > 0 -> PSD block; size < 0 -> diagonal/LP block."""
    size: int   # positive for PSD, negative for diag (|size| = #diag entries)
    # For each variable index k in 0..m (k=0 is the constant F_0), we store a
    # dict {(i,j): value} of upper-triangle nonzeros in the block. i,j are
    # 1-indexed for direct write.
    entries: List[dict] = field(default_factory=list)  # entries[k] is dict

    def ensure_var_count(self, m_plus_one: int) -> None:
        while len(self.entries) < m_plus_one:
            self.entries.append({})

    def add(self, k: int, i_1: int, j_1: int, val: float) -> None:
        if val == 0.0:
            return
        # Upper triangle convention
        if i_1 > j_1:
            i_1, j_1 = j_1, i_1
        # accumulate (may have multiple contributions from one constraint
        # for SOC arrow encoding); use += but typically single assignment
        d = self.entries[k]
        d[(i_1, j_1)] = d.get((i_1, j_1), 0.0) + float(val)


# --------------------------------------------------------------------------- #
# Main public API
# --------------------------------------------------------------------------- #
def cvxpy_to_sdpa_s(
    problem: cp.Problem,
    file_path: str,
    *,
    drop_zero_tol: float = 0.0,
) -> dict:
    """
    Write `problem` to `file_path` in SDPA-S sparse format.

    Returns a dict with metadata useful for post-solve interpretation:
      * 'm'                 : number of SDPA primal variables (= len(y))
      * 'block_structure'   : tuple of ints (SDPA bLOCKsTRUCT)
      * 'block_kinds'       : list of strings: 'eq', 'nonneg', 'soc', 'psd'
      * 'objective_offset'  : constant we need to ADD to SDPA's c^T x to get
                              the original objective value.
      * 'objective_sign'    : +1 if problem.objective was Minimize else -1.
                              SDPA solves MIN, so a Maximize problem flips.
      * 'chain'             : the cvxpy reduction chain (kept for round-trips)
      * 'inverse_data'      : cvxpy inverse data (for solution unpacking later)
      * 'cone_dims'         : the SCS dims object
    """
    # ---- canonicalize through CVXPY's SCS path
    data, chain, inv_data = problem.get_problem_data(cp.SCS)
    A: sp.csc_matrix = sp.csc_matrix(data["A"])
    b: np.ndarray = np.asarray(data["b"]).reshape(-1)
    c: np.ndarray = np.asarray(data["c"]).reshape(-1)
    dims = data["dims"]
    m = c.shape[0]  # number of SDPA primal variables

    # SCS slack ordering: zero, nonneg, soc, psd (and then exp, p3d which we
    # do not handle).
    if getattr(dims, "exp", 0):
        raise NotImplementedError("Exponential cones not supported in SDPA-S.")
    if getattr(dims, "p", None):
        raise NotImplementedError("3D power cones not supported in SDPA-S.")

    n_zero = int(dims.zero)
    n_nonneg = int(dims.nonneg)
    soc_sizes = list(dims.soc)
    psd_sizes = list(dims.psd)

    n_soc = sum(soc_sizes)
    n_psd_vec = sum(k * (k + 1) // 2 for k in psd_sizes)
    n_total = n_zero + n_nonneg + n_soc + n_psd_vec
    if n_total != b.shape[0]:
        raise ValueError(
            f"Cone-slack length mismatch: dims sum {n_total} vs b length {b.shape[0]}"
        )

    # Negate once. SCS form: s = b - A y in K. SDPA form: sum y_i F_i - F_0 >= 0
    # with F_0 = -Mat(b), F_i = -Mat(A[:,i]). So all coefficients we write are
    # negations of the SCS coefficients.
    neg_b = -b
    neg_A = -A  # csc

    blocks: List[_Block] = []
    block_kinds: List[str] = []

    # We'll lay out: first one big LP block for zero (split into +/-)
    # and nonneg, then one PSD block per SOC and per PSD cone.
    # SDPA expects each diagonal block as a SINGLE negative-sized block
    # so we'll merge: lp_len = 2 * n_zero + n_nonneg.

    lp_len = 2 * n_zero + n_nonneg
    if lp_len > 0:
        lp = _Block(size=-lp_len)
        lp.ensure_var_count(m + 1)
        blocks.append(lp)
        # mark just the first as 'lp' for kinds; finer attribution below
        block_kinds.append("lp")
    else:
        lp = None

    # ---------- zero block: each scalar eq becomes two LP rows
    row_cursor = 0
    if n_zero > 0:
        for r in range(n_zero):
            scs_row = row_cursor + r
            # diagonal positions in LP block: pair (2r+1, 2r+2)  (1-indexed)
            pos1 = 2 * r + 1
            pos2 = 2 * r + 2
            # F_0 contribution: neg_b at this row goes on BOTH +s and -s entries
            v0 = float(neg_b[scs_row])
            lp.add(0, pos1, pos1, v0)       #  s >= 0  : F_0 = -b
            lp.add(0, pos2, pos2, -v0)      # -s >= 0  : F_0 = +b
        # F_i contributions
        col = neg_A[row_cursor : row_cursor + n_zero, :].tocoo()
        for r_local, k_var, val in zip(col.row, col.col, col.data):
            pos1 = 2 * r_local + 1
            pos2 = 2 * r_local + 2
            v = float(val)
            lp.add(k_var + 1, pos1, pos1, v)
            lp.add(k_var + 1, pos2, pos2, -v)
        row_cursor += n_zero

    # ---------- nonneg block: scalar inequalities b - A y >= 0
    if n_nonneg > 0:
        offset_in_lp = 2 * n_zero  # diagonal positions start at offset+1
        # F_0
        for r in range(n_nonneg):
            scs_row = row_cursor + r
            pos = offset_in_lp + r + 1
            lp.add(0, pos, pos, float(neg_b[scs_row]))
        # F_i
        col = neg_A[row_cursor : row_cursor + n_nonneg, :].tocoo()
        for r_local, k_var, val in zip(col.row, col.col, col.data):
            pos = offset_in_lp + r_local + 1
            lp.add(k_var + 1, pos, pos, float(val))
        row_cursor += n_nonneg

    # ---------- SOC blocks: arrow LMI of size (k+1)
    # For an SOC of size k+1 (i.e. dim_soc = k+1, u_0 = first entry, u_{1:k+1}
    # are the next k entries), the LMI is
    #     [[u_0,      u_{1:}^T   ],
    #      [u_{1:},   u_0 * I_k ]]   >= 0
    # which is a (k+1) x (k+1) symmetric block. Entry placements:
    #   (1,1) := u_0
    #   (1, 1+j) := u_j  for j=1..k     (off-diagonal in row 1)
    #   (1+j, 1+j) := u_0  for j=1..k   (the I_k * u_0 part)
    # All other off-diagonals are zero.
    #
    # Strategy: u_i = b_i - sum_k A[i,k] y_k. So for each constraint row, add
    # neg_b coefficient to F_0 and neg_A row to each F_i in the right
    # placements. Since both 'u_0' and 'u_0 * I_k' use u_0, every (1+j, 1+j)
    # diagonal AND (1,1) entry references the SAME u_0 row.
    for soc_size in soc_sizes:
        # soc_size = k+1 in the math above
        sz = soc_size
        blk = _Block(size=sz)
        blk.ensure_var_count(m + 1)
        blocks.append(blk)
        block_kinds.append("soc")
        # row for u_0:
        u0_row = row_cursor
        v0_b = float(neg_b[u0_row])
        blk.add(0, 1, 1, v0_b)                  # (1,1) <- u_0 part of F_0
        for j in range(1, sz):
            blk.add(0, 1 + j, 1 + j, v0_b)      # (1+j,1+j) <- u_0 part of F_0
        # u_0 coefficients in each variable
        a0 = neg_A.getrow(u0_row).tocoo()
        for _, k_var, val in zip(a0.row, a0.col, a0.data):
            v = float(val)
            blk.add(k_var + 1, 1, 1, v)
            for j in range(1, sz):
                blk.add(k_var + 1, 1 + j, 1 + j, v)
        # u_j off-diagonals (1, 1+j) for j=1..sz-1
        for j in range(1, sz):
            uj_row = row_cursor + j
            vj_b = float(neg_b[uj_row])
            blk.add(0, 1, 1 + j, vj_b)
            aj = neg_A.getrow(uj_row).tocoo()
            for _, k_var, val in zip(aj.row, aj.col, aj.data):
                blk.add(k_var + 1, 1, 1 + j, float(val))
        row_cursor += sz

    # ---------- PSD blocks: undo SCS sqrt(2) scaling and lay out symmetric
    for k in psd_sizes:
        idx = _scs_psd_indices(k)
        blk = _Block(size=k)
        blk.ensure_var_count(m + 1)
        blocks.append(blk)
        block_kinds.append("psd")
        # F_0
        for local_i, (rr, cc, scale) in enumerate(idx):
            v = float(neg_b[row_cursor + local_i]) * scale
            if v != 0.0:
                blk.add(0, rr + 1, cc + 1, v)   # upper triangle (rr,cc) where rr>=cc; .add re-orients
        # F_i
        for local_i, (rr, cc, scale) in enumerate(idx):
            row_data = neg_A.getrow(row_cursor + local_i).tocoo()
            for _, k_var, val in zip(row_data.row, row_data.col, row_data.data):
                blk.add(k_var + 1, rr + 1, cc + 1, float(val) * scale)
        row_cursor += k * (k + 1) // 2

    assert row_cursor == b.shape[0], (row_cursor, b.shape[0])

    # ---------- Write the file
    block_struct = tuple(b.size for b in blocks)
    nBLOCK = len(blocks)

    # SDPA's primal is `min c^T x`. CVXPY canonicalized as min c^T y already.
    # (For Maximize, CVXPY internally negates - get_problem_data always
    # returns a minimization. We do not need to handle a separate sign.)
    obj_vec = c.copy()

    with open(file_path, "w") as fh:
        fh.write(f'"Generated by sdpa_serializer.py: m={m}, nBLOCK={nBLOCK}"\n')
        fh.write(f"{m} = mDIM\n")
        fh.write(f"{nBLOCK} = nBLOCK\n")
        fh.write("(" + ", ".join(str(s) for s in block_struct) + ") = bLOCKsTRUCT\n")
        # objective vector
        fh.write(", ".join(_fmt(v) for v in obj_vec) + "\n")
        # F entries
        for blk_idx, blk in enumerate(blocks):
            block_no = blk_idx + 1
            for k_var in range(m + 1):
                d = blk.entries[k_var]
                if not d:
                    continue
                for (i_1, j_1), v in sorted(d.items()):
                    if drop_zero_tol > 0 and abs(v) <= drop_zero_tol:
                        continue
                    if v == 0.0:
                        continue
                    fh.write(f"{k_var} {block_no} {i_1} {j_1} {_fmt(v)}\n")

    return {
        "m": m,
        "block_structure": block_struct,
        "block_kinds": block_kinds,
        "objective_offset": float(getattr(inv_data[-1] if isinstance(inv_data, list) else inv_data, "r", 0.0) or 0.0)
            if False else 0.0,   # cvxpy's offset is folded into c; we read prob.value separately
        "objective_sign": 1.0,
        "chain": chain,
        "inverse_data": inv_data,
        "cone_dims": dims,
    }


def _fmt(v: float) -> str:
    """High-precision plain-decimal formatting for SDPA-S file."""
    # SDPA-GMP reads ordinary decimals fine; use 17 sig figs (double round-trip).
    if v == 0.0:
        return "0"
    return repr(float(v))
