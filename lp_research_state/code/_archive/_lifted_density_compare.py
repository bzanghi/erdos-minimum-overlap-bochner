"""
Part 2/3 — Lift discrete M(n) optimizers to continuous densities and
compare to Together's h*.

The discrete A_n* ⊂ {1, ..., 2n} can be lifted to a density on [0, 1] via:

    f_n(x) := 1_{A_n*}(ceil(2n x))            (cell-indicator on 2n grid)

We then smooth f_n by convolution with a centered box of width 1/n.

Together's h* is on [0, 2] with 600 cells.  We compare on [0, 1] by
*folding* Together's h via:
    h_fold(x) := (h*(x) + h*(2 - x)) / 2,   x ∈ [0, 1]
Reason: A ↔ B complement symmetry in the discrete problem corresponds to
h ↔ 1 - h symmetry on [0, 2] (and to mirror symmetry around x=1 by
shift invariance of the autocorrelation argument).  We compute BOTH
"fold to [0,1]" (avg with mirror) AND simply restrict to [0, 1] for
comparison, since h* on [0,2] has asymmetry max |h_i - h_{n-1-i}| ≈ 0.53.

Note on normalization
---------------------
Discrete A_n* has |A| = n out of 2n positions, so f_n integrated over
[0, 1] equals... wait.  A ⊂ {1, .., 2n}, but we lift to [0, 1] using
the 2n-grid.  So ∫_0^1 f_n(x) dx = |A| / (2n) = 1/2.  This matches
Together's h_fold on [0, 1] which has integral 1/2 (since ∫_0^2 h = 1
and h_fold averages over a symmetric pair).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DATA = Path(__file__).parent.parent / "data"
MN_FILE = DATA / "Mn_optimizers_large.json"
TOGETHER_FILE = DATA / "together_f_star.json"


def load_Mn():
    with open(MN_FILE) as f:
        return json.load(f)


def load_h_star():
    with open(TOGETHER_FILE) as f:
        d = json.load(f)
    h = np.asarray(d["together"]["values"], dtype=np.float64)
    # bp on [0, 2] with 601 entries; cell i covers [2i/600, 2(i+1)/600]
    return h


def lift_to_density(A_star, n, x_grid):
    """f_n on [0, 1] as cell indicator on 2n-grid, sampled at x_grid."""
    # cell index for x in [0, 1]: cell j (1..2n) covers [(j-1)/(2n), j/(2n)]
    # x_grid in [0, 1]
    cell = np.minimum(np.floor(x_grid * 2 * n).astype(int), 2 * n - 1)
    indicator = np.zeros(2 * n, dtype=np.float64)
    for p in A_star:
        indicator[p - 1] = 1.0  # p is 1-indexed
    return indicator[cell]


def smooth_density(f, x_grid, width):
    """Box-convolve f on x_grid with a kernel of full-width `width`."""
    # Build kernel
    dx = x_grid[1] - x_grid[0]
    k_len = max(1, int(round(width / dx)))
    if k_len % 2 == 0:
        k_len += 1
    kernel = np.ones(k_len) / k_len
    # 'same' mode with reflection-style boundary
    fpad = np.concatenate([f[::-1], f, f[::-1]])
    sm = np.convolve(fpad, kernel, mode='same')
    return sm[len(f):2 * len(f)]


def h_fold(h, x_grid):
    """Sample h_fold(x) = (h(x) + h(2-x)) / 2 at x ∈ [0, 1]."""
    # h is on [0, 2] with 600 cells; cell i covers [2i/600, 2(i+1)/600]
    n_h = len(h)  # 600
    # value at x: cell index = floor(x * n_h / 2) for x in [0, 2]
    def sample(xs):
        ci = np.minimum(np.floor(xs * n_h / 2.0).astype(int), n_h - 1)
        ci = np.maximum(ci, 0)
        return h[ci]
    return (sample(x_grid) + sample(2.0 - x_grid)) / 2.0


def h_restrict_first_half(h, x_grid):
    """Just take h restricted to [0, 1]."""
    n_h = len(h)
    ci = np.minimum(np.floor(x_grid * n_h / 2.0).astype(int), n_h - 1)
    return h[ci]


def L1(f, g, x_grid):
    return float(np.trapz(np.abs(f - g), x_grid))


def L2(f, g, x_grid):
    return float(np.sqrt(np.trapz((f - g) ** 2, x_grid)))


def pearson(f, g):
    f = f - f.mean()
    g = g - g.mean()
    denom = np.sqrt((f * f).sum() * (g * g).sum())
    if denom < 1e-15:
        return 0.0
    return float((f * g).sum() / denom)


def ks_dist(f, g, x_grid):
    """Max distance between cumulative integrals."""
    Cf = np.concatenate([[0.0], np.cumsum((f[:-1] + f[1:]) / 2 * np.diff(x_grid))])
    Cg = np.concatenate([[0.0], np.cumsum((g[:-1] + g[1:]) / 2 * np.diff(x_grid))])
    return float(np.max(np.abs(Cf - Cg)))


def main():
    Mn_data = load_Mn()
    h = load_h_star()

    # Sample grid on [0, 1]
    NG = 1200
    x = np.linspace(0.0, 1.0, NG)

    h_fold_vals = h_fold(h, x)
    h_rest_vals = h_restrict_first_half(h, x)

    # Together's h represents a density on [0, 2] with integral 1.
    # h_fold on [0, 1] has integral 1/2 (matches f_n which integrates
    # to 1/2 on [0, 1]).  Verify:
    int_fold = np.trapz(h_fold_vals, x)
    int_rest = np.trapz(h_rest_vals, x)
    print(f"  ∫_0^1 h_fold = {int_fold:.6f}   (expect ~0.5)")
    print(f"  ∫_0^1 h_rest = {int_rest:.6f}   (no a priori expected value)")

    table_rows = []
    densities = {}
    for entry in Mn_data:
        n = entry["n"]
        A_star = entry["A_star"]
        f_A = lift_to_density(A_star, n, x)
        # The A↔B symmetry means 1_B = 1 - 1_A is equally optimal.
        # Try both orientations and report the better-matching one.
        f_A_sm = smooth_density(f_A, x, 1.0 / n)
        f_B_sm = 1.0 - f_A_sm  # smoothed complement

        # Cyclic-shift alignment: the discrete optimizer is only defined
        # up to cyclic shift on Z/2nZ in some formulations, but here A is
        # an asymmetric subset of {1..2n} so we should also consider all
        # cyclic shifts.  In addition to A↔B and reflection (which is
        # already canonicalized in the optimizer), the "lifting to [0,1]"
        # vs "Together's h* on [0,2]" has no canonical alignment, so we
        # search over translation + complement.

        def best_alignment(f_sm, target):
            """Search over cyclic shifts and complement; return best L1."""
            best = (np.inf, None, None)  # (L1, shift_idx, do_complement)
            n_grid = len(f_sm)
            for do_comp in (False, True):
                candidate = (1.0 - f_sm) if do_comp else f_sm
                for shift in range(0, n_grid, max(1, n_grid // 80)):
                    rolled = np.roll(candidate, shift)
                    d = L1(rolled, target, x)
                    if d < best[0]:
                        best = (d, shift, do_comp)
            return best

        L1_A_raw = L1(f_A_sm, h_fold_vals, x)
        L1_B_raw = L1(f_B_sm, h_fold_vals, x)
        L1_min_raw = min(L1_A_raw, L1_B_raw)

        # Best L1 over alignments
        align_fold = best_alignment(f_A_sm, h_fold_vals)
        align_rest = best_alignment(f_A_sm, h_rest_vals)

        # Apply best alignment to get matching f for plotting / correlation
        def apply_align(f_sm, align):
            _, shift, do_comp = align
            cand = (1.0 - f_sm) if do_comp else f_sm
            return np.roll(cand, shift)

        f_aligned_fold = apply_align(f_A_sm, align_fold)
        L1_fold = align_fold[0]
        L2_fold = L2(f_aligned_fold, h_fold_vals, x)
        corr_fold = pearson(f_aligned_fold, h_fold_vals)
        ks_fold = ks_dist(f_aligned_fold, h_fold_vals, x)

        f_aligned_rest = apply_align(f_A_sm, align_rest)
        L1_rest = align_rest[0]
        corr_rest = pearson(f_aligned_rest, h_rest_vals)

        densities[f"n_{n}"] = {
            "raw": f_A,
            "smoothed": f_A_sm,
            "aligned_to_fold": f_aligned_fold,
            "x": x,
        }
        row = {
            "n": n,
            "M": entry["M"],
            "blocks": entry["num_blocks"],
            "L1_fold_raw": L1_min_raw,
            "L1_fold": L1_fold,
            "L1_rest": L1_rest,
            "L2_fold": L2_fold,
            "corr_fold": corr_fold,
            "corr_rest": corr_rest,
            "ks_fold": ks_fold,
        }
        table_rows.append(row)
        print(f"n={n:2d}  M={entry['M']:2d}  blocks={entry['num_blocks']:2d}  "
              f"L1raw={L1_min_raw:.4f}  L1align={L1_fold:.4f}  "
              f"corr={corr_fold:+.3f}  KS={ks_fold:.4f}")

    # Save densities .npz
    npz_payload = {"x": x, "h_fold": h_fold_vals, "h_restrict": h_rest_vals}
    for k, v in densities.items():
        npz_payload[f"{k}_raw"] = v["raw"]
        npz_payload[f"{k}_smoothed"] = v["smoothed"]
    np.savez(DATA / "lifted_densities.npz", **npz_payload)
    print(f"\nSaved densities to {DATA / 'lifted_densities.npz'}")

    # Save table
    with open(DATA / "lifted_comparison_table.json", "w") as f:
        json.dump(table_rows, f, indent=2)

    # Plot: top = f_n smoothed for select n, with h_fold; bottom = L1 + corr
    fig, axes = plt.subplots(2, 1, figsize=(11, 9))
    ax_top, ax_bot = axes

    show_ns = [4, 8, 12, 14, 16]
    cmap = plt.get_cmap("viridis")
    for i, n in enumerate(show_ns):
        sm = densities[f"n_{n}"]["smoothed"]
        ax_top.plot(x, sm, color=cmap(i / max(1, len(show_ns) - 1)),
                    lw=1.4, alpha=0.85, label=f"f_{n} (smoothed, 1/n box)")
    ax_top.plot(x, h_fold_vals, color="crimson", lw=2.0,
                label="Together h_fold = (h*(x) + h*(2-x))/2")
    ax_top.plot(x, h_rest_vals, color="black", lw=1.0, ls="--", alpha=0.6,
                label="Together h*|[0,1]")
    ax_top.set_xlabel("x")
    ax_top.set_ylabel("density")
    ax_top.set_title("Lifted discrete optimizers f_n (smoothed) vs Together's h*")
    ax_top.legend(loc="upper right", fontsize=8)
    ax_top.grid(alpha=0.3)

    ns = np.array([r["n"] for r in table_rows])
    L1f = np.array([r["L1_fold"] for r in table_rows])
    L1r = np.array([r["L1_rest"] for r in table_rows])
    cf = np.array([r["corr_fold"] for r in table_rows])
    cr = np.array([r["corr_rest"] for r in table_rows])

    ax_bot.plot(ns, L1f, "o-", color="C0", label="L1(f_n, h_fold)")
    ax_bot.plot(ns, L1r, "s--", color="C0", alpha=0.5, label="L1(f_n, h*|[0,1])")
    ax_bot.set_xlabel("n")
    ax_bot.set_ylabel("L1 distance", color="C0")
    ax_bot.tick_params(axis="y", labelcolor="C0")
    ax_bot.grid(alpha=0.3)
    ax_bot.legend(loc="upper left", fontsize=8)

    ax_bot2 = ax_bot.twinx()
    ax_bot2.plot(ns, cf, "^-", color="C3", label="Pearson(f_n, h_fold)")
    ax_bot2.plot(ns, cr, "v--", color="C3", alpha=0.5,
                 label="Pearson(f_n, h*|[0,1])")
    ax_bot2.set_ylabel("Pearson correlation", color="C3")
    ax_bot2.tick_params(axis="y", labelcolor="C3")
    ax_bot2.legend(loc="upper right", fontsize=8)
    ax_bot.set_title("Convergence metrics vs n")

    fig.tight_layout()
    fig.savefig(DATA / "fn_vs_h_convergence.png", dpi=130)
    print(f"Saved plot to {DATA / 'fn_vs_h_convergence.png'}")

    return table_rows


if __name__ == "__main__":
    main()
