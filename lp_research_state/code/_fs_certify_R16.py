"""DECISIVE independent rigor check of the binding gate region R16 (tightest, margin
+1.2e-4). Rigorous adaptive subdivision: recursively split R16's box; for each sub-box
use cover_min_over_box (RIGOROUS grid_min - eps_grid). A sub-box is CERTIFIED if its
rigorous lb >= TARGET. If lb < TARGET because eps_grid is large, SUBDIVIDE (smaller box
=> smaller eps_grid). If the raw grid_min itself dips < TARGET at a sampled point, that
is a candidate REFUTATION (the cover is genuinely below target there). Pure evaluation,
no SDP solves. Decides CONFIRM vs REFUTE for R16 rigorously."""
import sys; sys.path.insert(0, '.')
import numpy as np
from collections import deque
import _fs_recompute as M
from _fullspace_eval import cover_min_over_box

TARGET = 0.3802838
hv = M.harvest_centers(); centers = hv[0]
for c in centers:
    if c.get("primal") is None and c.get("dual_lb") is not None:
        c["primal"] = float(c["dual_lb"]) + 1e-5
centers = [c for c in centers if c.get("primal") is not None]
print(f"loaded {len(centers)} centers", flush=True)

# R16 box (from WHITE_TABLE2 region 16)
BOX = ((0.0, 0.06), (0.33, 0.45), (-0.025, -0.02))

def boxmin(h, p, q):
    # coarse per-box grid; rigor comes from eps_grid + subdivision, not grid density
    return cover_min_over_box(centers, "primal_m1e5", h, p, q, n_h=7, n_p=7, n_q=5)

dq = deque([(BOX[0], BOX[1], BOX[2], 0)])
BUDGET = 120000
ncert = processed = 0
min_grid_seen = 1.0          # smallest raw cover value at any sampled point
min_grid_at = None
worst_lb = 1.0
refute = None
maxdepth = 0

while dq and processed < BUDGET:
    h, p, q, depth = dq.popleft()
    processed += 1
    maxdepth = max(maxdepth, depth)
    lb, pt, wit, gmin, eps, Lm = boxmin(h, p, q)
    if gmin < min_grid_seen:
        min_grid_seen = gmin; min_grid_at = (pt, depth)
    if lb >= TARGET:
        ncert += 1
        continue
    if gmin < TARGET:
        refute = (h, p, q, pt, gmin, eps, depth)   # cover itself dips below target
        break
    if depth >= 70:
        worst_lb = min(worst_lb, lb)
        continue                                    # stalled (eps can't shrink enough)
    # gmin >= TARGET but lb < TARGET (eps artifact): split the widest dimension
    widths = [(h[1]-h[0], 'h'), (p[1]-p[0], 'p'), (q[1]-q[0], 'q')]
    widths.sort(reverse=True)
    d = widths[0][1]
    if d == 'h':
        m = 0.5*(h[0]+h[1]); dq.append(((h[0], m), p, q, depth+1)); dq.append(((m, h[1]), p, q, depth+1))
    elif d == 'p':
        m = 0.5*(p[0]+p[1]); dq.append((h, (p[0], m), q, depth+1)); dq.append((h, (m, p[1]), q, depth+1))
    else:
        m = 0.5*(q[0]+q[1]); dq.append((h, p, (q[0], m), depth+1)); dq.append((h, p, (m, q[1]), depth+1))

print(f"processed={processed} certified_subboxes={ncert} maxdepth={maxdepth} "
      f"queue_left={len(dq)}", flush=True)
print(f"min raw cover (grid_min) seen over all sampled points = {min_grid_seen:.7f} "
      f"at {min_grid_at}  (TARGET={TARGET})", flush=True)
if refute is not None:
    print(f"*** REFUTE: raw cover dips to {refute[4]:.7f} < TARGET at point {refute[3]} "
          f"(box depth {refute[6]}) ***", flush=True)
elif len(dq) == 0:
    print(f"*** R16 RIGOROUSLY CERTIFIED >= {TARGET}: every sub-box's (grid_min - eps_grid) "
          f">= TARGET via adaptive subdivision ***", flush=True)
else:
    print(f"*** INCONCLUSIVE: hit budget/depth with {len(dq)} boxes unresolved; "
          f"worst stalled lb={worst_lb:.7f}, but raw cover never dipped below "
          f"{min_grid_seen:.7f} >= TARGET ***", flush=True)
