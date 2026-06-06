"""Independent CONSERVATIVE recheck: the STANDARD single-grid cover_min_over_box
floor over all 18 White Table-2 regions + core, using the harvested union of
point-centers. This deliberately AVOIDS the adaptive-subdivision / infeasibility
arguments that the workflow's headline (0.3802838) relies on, to establish the
defensible lower number. No SDP solves."""
import sys; sys.path.insert(0, '.')
import numpy as np
import _fs_recompute as M
from _fullspace_eval import (cover_min_over_box, reproduce_core_headline,
                             WHITE_TABLE2, WHITE_OUTSIDE_FLOOR)

_hv = M.harvest_centers()
centers = _hv[0]
sources = _hv[1] if len(_hv) > 1 else {}
# Normalize: anchor_value() needs c["primal"]; for centers carrying only dual_lb,
# set primal = dual_lb + 1e-5 so anchor = primal-1e-5 = dual_lb (conservative).
kept = []
for c in centers:
    if c.get("primal") is None:
        if c.get("dual_lb") is not None:
            c = dict(c); c["primal"] = float(c["dual_lb"]) + 1e-5
        else:
            continue
    kept.append(c)
print(f"harvested {len(centers)} centers, usable {len(kept)}; sources={sources}", flush=True)

core = reproduce_core_headline(kept, "primal_m1e5")
print(f"[CORE] standard rigorous_LB = {core['rigorous_LB']:.7f}  binding@{core['binding_point']}", flush=True)

worst = ("core", core["rigorous_LB"], None)
for idx, (hr, pr, qr, wb) in enumerate(WHITE_TABLE2, start=1):
    nh = 81 if (hr[1]-hr[0]) > 0.05 else 41
    npp = 161 if (pr[1]-pr[0]) > 0.2 else 81
    nq = 81 if (qr[1]-qr[0]) > 0.1 else 41
    lb, pt, wit, gmin, eps, Lm = cover_min_over_box(
        kept, "primal_m1e5", hr, pr, qr, n_h=nh, n_p=npp, n_q=nq)
    tag = "  <-- below core" if lb < core["rigorous_LB"] else ""
    print(f"R{idx:2d} standard={lb:.6f} grid_min={gmin:.6f} eps_grid={eps:.2e} "
          f"L_max={Lm:.2f} worst@({pt[0]:.3f},{pt[1]:.3f},{pt[2]:.3f}){tag}", flush=True)
    if lb < worst[1]:
        worst = (f"R{idx}", lb, pt)

print(f"\n==> STANDARD-METHOD full-space floor (min over core + 18 regions) = "
      f"{worst[1]:.7f} at {worst[0]} {worst[2]}", flush=True)
print(f"    vs White 0.379005: {worst[1]-0.379005:+.6e}", flush=True)
print(f"    vs core headline 0.3802838: {worst[1]-0.3802838:+.6e}", flush=True)
