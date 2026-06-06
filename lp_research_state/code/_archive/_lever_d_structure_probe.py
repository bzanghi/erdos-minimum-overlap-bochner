"""
Lever D pre-investigation: empirical structure of optimal partitions.

Tests the hypothesis "the optimal f for the continuous Erdős problem lies
near piecewise-constant densities with O(1) breakpoints" by inspecting:

1. Together's claimed near-optimal h* (600 cells on [0,2]).
2. Brute-force-optimal small-n integer partitions A* for n=2..12.

Produces:
  - lp_research_state/data/together_h_structure.png
  - lp_research_state/data/Mn_optimizers_structure.png
  - Console output with all numerical findings.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lp_research_state" / "code"))

from together_loader import load_together_raw


# ---------- Part 1: Together's h* ----------

def block_count(h: np.ndarray, tol: float) -> int:
    """Count maximal runs where adjacent values differ by <= tol."""
    if len(h) == 0:
        return 0
    jumps = np.abs(np.diff(h)) > tol
    return int(jumps.sum()) + 1


def jump_count(h: np.ndarray, thresh: float) -> int:
    return int((np.abs(np.diff(h)) > thresh).sum())


def analyze_together(h: np.ndarray) -> dict:
    n = len(h)
    rounded_1e3 = np.round(h, 3)
    rounded_1e2 = np.round(h, 2)
    rounded_5e2 = (np.round(h / 0.05) * 0.05)
    rounded_1e1 = np.round(h, 1)
    out = {
        "n": n,
        "sum_h": float(h.sum()),
        "min": float(h.min()),
        "max": float(h.max()),
        "mean": float(h.mean()),
        "nonzero_cells_1e-6": int((h > 1e-6).sum()),
        "nonzero_cells_1e-3": int((h > 1e-3).sum()),
        "near_one_cells_1e-3": int((h > 1 - 1e-3).sum()),
        "mid_cells": int(((h > 0.05) & (h < 0.95)).sum()),
        "distinct_values_1e-3": int(len(np.unique(rounded_1e3))),
        "distinct_values_1e-2": int(len(np.unique(rounded_1e2))),
        "distinct_values_5e-2": int(len(np.unique(np.round(rounded_5e2, 5)))),
        "distinct_values_1e-1": int(len(np.unique(rounded_1e1))),
        "blocks_1e-3": block_count(h, 1e-3),
        "blocks_1e-2": block_count(h, 1e-2),
        "blocks_5e-2": block_count(h, 5e-2),
        "blocks_1e-1": block_count(h, 1e-1),
        "blocks_2e-1": block_count(h, 2e-1),
        "blocks_5e-1": block_count(h, 5e-1),
        "jumps_>0.01": jump_count(h, 0.01),
        "jumps_>0.05": jump_count(h, 0.05),
        "jumps_>0.1": jump_count(h, 0.1),
        "jumps_>0.2": jump_count(h, 0.2),
        "jumps_>0.5": jump_count(h, 0.5),
        "total_variation": float(np.abs(np.diff(h)).sum()),
    }
    return out


def plot_together(h: np.ndarray, out_path: Path) -> None:
    n = len(h)
    x = np.arange(n)
    fig, axes = plt.subplots(3, 1, figsize=(10, 9))

    # Top: h vs i
    axes[0].plot(x, h, lw=0.7, color="C0")
    axes[0].fill_between(x, 0, h, color="C0", alpha=0.25)
    axes[0].set_ylabel("h_i")
    axes[0].set_xlabel("cell index i (0..599)")
    axes[0].set_title("Together's h* on [0,2] (600 cells)")
    axes[0].set_xlim(0, n)
    axes[0].set_ylim(-0.02, 1.05)
    axes[0].grid(alpha=0.3)

    # Middle: histogram
    axes[1].hist(h, bins=20, range=(0, 1), color="C1", edgecolor="black")
    axes[1].set_xlabel("value h_i")
    axes[1].set_ylabel("# cells")
    axes[1].set_title("Histogram of h_i values (20 bins on [0,1])")
    axes[1].grid(alpha=0.3)

    # Bottom: cumulative total variation
    cumvar = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(h)))])
    axes[2].plot(np.arange(len(cumvar)), cumvar, color="C2")
    axes[2].set_xlabel("cell index i")
    axes[2].set_ylabel("cumulative TV = sum_{j<=i} |h_j - h_{j-1}|")
    axes[2].set_title(f"Total variation accumulation (final TV = {cumvar[-1]:.3f})")
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------- Part 2: small-n optimizers ----------

def A_to_vec(A_star, n: int) -> np.ndarray:
    """A_star is a list of 1-indexed positions in {1, ..., 2n}.
       Returns a 0/1 vector of length 2n (A=1, B=0)."""
    v = np.zeros(2 * n, dtype=np.int8)
    for j in A_star:
        v[j - 1] = 1
    return v


def count_blocks_int(v: np.ndarray) -> int:
    if len(v) == 0:
        return 0
    return int((np.diff(v) != 0).sum()) + 1


def is_palindrome(v: np.ndarray) -> bool:
    """A* is palindromic in the symmetry that maps j -> 2n+1-j (swap A/B)."""
    return bool(np.all(v == v[::-1]))


def is_antipalindrome(v: np.ndarray) -> bool:
    """A* is anti-palindromic: v[j] = 1 - v[2n-1-j] (A and B swap under reflection)."""
    return bool(np.all(v == 1 - v[::-1]))


def analyze_small_n(table: list) -> list:
    rows = []
    for entry in table:
        n = entry["n"]
        A_star = entry["A_star"]
        v = A_to_vec(A_star, n)
        rows.append({
            "n": n,
            "M": entry["M"],
            "ratio": entry["ratio"],
            "blocks": count_blocks_int(v),
            "palindrome": is_palindrome(v),
            "antipalindrome": is_antipalindrome(v),
            "vec": v.tolist(),
            "A_star": A_star,
        })
    return rows


def plot_small_n(rows: list, out_path: Path) -> None:
    n_max = max(r["n"] for r in rows)
    max_len = 2 * n_max
    fig, (ax_heat, ax_blocks) = plt.subplots(2, 1, figsize=(10, 7),
                                              gridspec_kw={"height_ratios": [3, 2]})

    # Heatmap: rows are n, columns are positions (padded with NaN)
    grid = np.full((len(rows), max_len), np.nan)
    for i, r in enumerate(rows):
        v = np.array(r["vec"], dtype=float)
        grid[i, :len(v)] = v
    cmap = matplotlib.colors.ListedColormap(["#fdebd0", "#1f4e79"])  # B, A
    cmap.set_bad("white")
    ax_heat.imshow(grid, aspect="auto", cmap=cmap, interpolation="nearest")
    ax_heat.set_yticks(range(len(rows)))
    ax_heat.set_yticklabels([f"n={r['n']}" for r in rows])
    ax_heat.set_xlabel("position j in {1, ..., 2n}")
    ax_heat.set_title("Brute-force-optimal partition A* (dark = A, light = B); rows padded to 2n_max")
    ax_heat.grid(False)

    # Block counts vs n
    ns = [r["n"] for r in rows]
    blocks = [r["blocks"] for r in rows]
    ax_blocks.plot(ns, blocks, "o-", color="C3", label="num blocks in A*")
    # Reference curves
    ax_blocks.plot(ns, [np.log2(n + 1) for n in ns], "--", color="grey",
                   label="log2(n+1)")
    ax_blocks.plot(ns, [np.sqrt(n) for n in ns], ":", color="black",
                   label="sqrt(n)")
    ax_blocks.set_xlabel("n")
    ax_blocks.set_ylabel("block count")
    ax_blocks.set_title("Block count vs n (compared with log2, sqrt growth)")
    ax_blocks.legend()
    ax_blocks.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


# ---------- Main ----------

def main():
    data_dir = ROOT / "lp_research_state" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Part 1 — Together's h*
    print("=" * 70)
    print("Part 1 — Together's h* structure")
    print("=" * 70)
    bp, h, dom, meta = load_together_raw()
    stats = analyze_together(h)
    for k, v in stats.items():
        if isinstance(v, float):
            print(f"  {k:30s} = {v:.6f}")
        else:
            print(f"  {k:30s} = {v}")

    plot_together(h, data_dir / "together_h_structure.png")
    print(f"\nSaved {data_dir / 'together_h_structure.png'}")

    # Part 2 — small-n
    print("\n" + "=" * 70)
    print("Part 2 — small-n brute-force optimizers")
    print("=" * 70)
    Mn_json = ROOT / "min_overlap_session_2026-05-09" / "Mn_brute.json"
    with open(Mn_json) as f:
        Mn_table = json.load(f)
    rows = analyze_small_n(Mn_table)

    print(f"\n{'n':>3} {'M(n)':>5} {'M/n':>7} {'blocks':>7} {'pal':>5} {'antipal':>8}  vector")
    for r in rows:
        v_str = "".join(str(x) for x in r["vec"])
        print(f"{r['n']:>3} {r['M']:>5} {r['ratio']:>7.4f} {r['blocks']:>7} "
              f"{str(r['palindrome']):>5} {str(r['antipalindrome']):>8}  {v_str}")

    # Quick growth check
    ns = np.array([r["n"] for r in rows])
    bls = np.array([r["blocks"] for r in rows])
    print("\nBlock-count growth diagnostics:")
    print(f"  blocks/n ratio: {(bls / ns).tolist()}")
    print(f"  blocks - log2(n+1): "
          f"{[round(b - np.log2(n + 1), 2) for n, b in zip(ns, bls)]}")
    # Linear-fit blocks vs n
    a, b = np.polyfit(ns, bls, 1)
    print(f"  linear fit:  blocks ≈ {a:.3f} * n + {b:.3f}")
    a_s, b_s = np.polyfit(np.sqrt(ns), bls, 1)
    print(f"  sqrt fit:    blocks ≈ {a_s:.3f} * sqrt(n) + {b_s:.3f}")
    a_l, b_l = np.polyfit(np.log(ns), bls, 1)
    print(f"  log fit:     blocks ≈ {a_l:.3f} * log(n) + {b_l:.3f}")

    plot_small_n(rows, data_dir / "Mn_optimizers_structure.png")
    print(f"\nSaved {data_dir / 'Mn_optimizers_structure.png'}")

    # Dump JSON summary too
    summary = {
        "together": stats,
        "small_n": rows,
    }
    (data_dir / "lever_d_structure_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(f"Saved {data_dir / 'lever_d_structure_summary.json'}")


if __name__ == "__main__":
    main()
