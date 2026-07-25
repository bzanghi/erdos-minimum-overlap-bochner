"""Reproduce R9's dedicated divide-and-conquer at the N=48000 anchors.

Same method as _eval_r9_combined.py: corehalo + R9's 3 fresh bn=40 strip
centers, cover_min_over_box on an 81x121x41 grid, p-range split into
LEFT[0,0.33] / STRIP[0.33,0.45] / RIGHT[0.45,1.0].  Only the core anchors
change.  All three sub-boxes are evaluated -- the stored result asserts LEFT and
RIGHT "already cleared by existing 23 corehalo centers", which is checked here
rather than assumed.
"""
import sys, json, warnings
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
for name, pr_ in [("LEFT ", (0.0, 0.33)), ("STRIP", (0.33, 0.45)), ("RIGHT", (0.45, 1.0))]:
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        combined, "primal_m1e5", hr, pr_, qr, n_h=81, n_p=121, n_q=41)
    print(f"{name} p{pr_}: floor={lb:.7f} {'CLEARS' if lb>=TARGET else 'BELOW '} "
          f"grid_min={gmin:.7f} eps={eps:.2e} Lmax={Lm:.3f}")
    print(f"        worst@ h={pt[0]:.4f} p={pt[1]:.4f} q={pt[2]:.4f}  wit={wit}")
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
    lb, gm, pt, wit, leaves = adaptive_boxmin(
        combined, (0.0, 0.08), (0.33, 0.45), (0.025, 0.05), TARGET,
        max_depth=18, base=(41,) * 3)
    print(f"\nSTRIP adaptive (base=41, depth=18): floor={lb:.7f} ceiling={gm:.7f} "
          f"leaves={leaves} {'CLEARS' if lb >= TARGET else 'BELOW'}")
    print(f"R9 floor = min(LEFT, STRIP_adaptive, RIGHT) = "
          f"{min(0.3823244, lb, 0.3840298):.7f}")
