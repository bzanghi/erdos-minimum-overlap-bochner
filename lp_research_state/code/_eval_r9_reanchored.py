"""Reproduce R9's dedicated divide-and-conquer at the N=48000 anchors.

Same method as _eval_r9_combined.py: corehalo + R9's 3 fresh bn=40 strip
centers, cover_min_over_box on an 81x121x41 grid, p-range split into
LEFT[0,0.33] / STRIP[0.33,0.45] / RIGHT[0.45,1.0].  Only the core anchors
change.  All three sub-boxes are evaluated -- the stored result asserts LEFT and
RIGHT "already cleared by existing 23 corehalo centers", which is checked here
rather than assumed.
"""
import sys, json, time, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
CODE = Path("/Users/benzanghi/Documents/Claude/Projects/Erdos/.claude/worktrees/"
            "minimum-overlap-problem-5de0df/lp_research_state/code")
sys.path.insert(0, str(CODE))
PR = CODE.parent / "parallel_results"

import _fs_recon as FS
FS.DUALEXT = PR / "dualext_reanchored_N48000.json"      # <- the only change
from _fullspace_eval import cover_min_over_box

core = FS.load_core(); core = core[0] if isinstance(core, tuple) else core
halo = FS.load_halo()
corehalo = list(core) + list(halo)

d9 = json.load(open(PR / "fullspace_promote_R9.json"))
fresh = []
for c in d9.get("centers", []):
    c = dict(c)
    if c.get("dual_lb") is not None:
        c["primal"] = c["dual_lb"]      # producer convention: anchor = dual_lb - 1e-5
    fresh.append(c)
combined = corehalo + fresh
print(f"corehalo={len(corehalo)}  fresh={len(fresh)}  anchors from {FS.DUALEXT.name}\n")

TARGET = 0.3803954          # the N=48000 certified core floor
hr, qr = (0.0, 0.08), (0.025, 0.05)
worst = None
subboxes = {}
for name, pr_ in [("LEFT ", (0.0, 0.33)), ("STRIP", (0.33, 0.45)), ("RIGHT", (0.45, 1.0))]:
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        combined, "primal_m1e5", hr, pr_, qr, n_h=81, n_p=121, n_q=41)
    print(f"{name} p{pr_}: floor={lb:.7f} {'CLEARS' if lb>=TARGET else 'BELOW '} "
          f"grid_min={gmin:.7f} eps={eps:.2e} Lmax={Lm:.3f}")
    print(f"        worst@ h={pt[0]:.4f} p={pt[1]:.4f} q={pt[2]:.4f}  wit={wit}")
    subboxes[name.strip()] = {
        "p_range": list(pr_), "floor_single_grid": float(lb),
        "grid_min": float(gmin), "eps_grid": float(eps), "L_max": float(Lm),
        "worst_point": {"h": float(pt[0]), "p": float(pt[1]), "q": float(pt[2])},
        "witness": str(wit), "grid": [81, 121, 41]}
    worst = lb if worst is None else min(worst, lb)

print(f"\nR9 floor (min over the three sub-boxes) = {worst:.7f}")
print(f"stored (N=20000 anchors)                = 0.3803667")
print(f"target (N=48000 core floor)             = {TARGET}")
print("VERDICT:", "R9 clears -> core binds at 0.3803954" if worst >= TARGET
      else f"R9 still binds at {worst:.7f}")

# --- STRIP needs adaptive subdivision at the N=48000 anchors ------------------
# A single grid over the strip gives floor 0.3801263: grid_min is 0.3804300,
# comfortably above target, but eps = 3.04e-4 swamps it because the certified
# duals are steeper than the old convention's (Lmax 0.15 -> 0.39).  The same fix
# that recovered the core floor works here -- subdivide, recomputing Lmax per
# sub-box.  At base=41, depth=18 the strip clears in 27 leaves.
if __name__ == "__main__":
    from _eval_r6_box import adaptive_boxmin
    t0 = time.time()
    lb, gm, pt, wit, leaves = adaptive_boxmin(
        combined, (0.0, 0.08), (0.33, 0.45), (0.025, 0.05), TARGET,
        max_depth=18, base=(41,) * 3)
    dt = time.time() - t0
    print(f"\nSTRIP adaptive (base=41, depth=18): floor={lb:.7f} ceiling={gm:.7f} "
          f"leaves={leaves} {'CLEARS' if lb >= TARGET else 'BELOW'}")
    subboxes["STRIP"]["floor_adaptive"] = float(lb)
    subboxes["STRIP"]["adaptive"] = {
        "base": 41, "max_depth": 18, "leaves": int(leaves),
        "ceiling_grid_min": float(gm), "witness": str(wit),
        "worst_point": {"h": float(pt[0]), "p": float(pt[1]), "q": float(pt[2])},
        "secs": dt}

    # R9's floor is the min over the three p-sub-boxes.  LEFT and RIGHT clear on
    # the single grid; only STRIP needs the adaptive pass (a single grid there
    # fails on eps, not on the bound -- see the note above).
    region_floor = min(subboxes["LEFT"]["floor_single_grid"], lb,
                       subboxes["RIGHT"]["floor_single_grid"])
    print(f"R9 floor = min(LEFT, STRIP_adaptive, RIGHT) = {region_floor:.7f}")

    out = PR / "gate_region_R9_N48000.json"
    out.write_text(json.dumps({
        "region": 9,
        "anchors": FS.DUALEXT.name,
        "target": TARGET,
        "h_range": [0.0, 0.08], "p_range": [0.0, 1.0], "q_range": [0.025, 0.05],
        "method": ("Dedicated divide-and-conquer: p split LEFT[0,0.33] / "
                   "STRIP[0.33,0.45] / RIGHT[0.45,1.0] over 23 corehalo + 3 fresh "
                   "bn=40 R9 centers.  LEFT/RIGHT clear on a single 81x121x41 grid; "
                   "STRIP needs adaptive subdivision (base=41, depth=18) because "
                   "the certified duals are steeper (L_max 0.15 -> 0.39) so a "
                   "single grid's eps swamps an otherwise-fine grid_min."),
        "n_corehalo": len(corehalo), "n_fresh": len(fresh),
        "subboxes": subboxes,
        "region_floor": float(region_floor),
        "clears_target": bool(region_floor >= TARGET),
        "stored_previous_N20000_anchors": 0.38036671736066713,
    }, indent=2, default=float))
    print(f"saved -> {out}")
