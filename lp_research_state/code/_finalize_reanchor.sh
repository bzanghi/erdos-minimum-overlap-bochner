#!/bin/bash
# Finalize the Jansson re-anchoring: turn 12 certified centers into the two
# publishable bounds.  Run from lp_research_state/code.
#
#   TIER 1  (fully certified by us, plus White's PUBLISHED Table 2)
#       White's Table 2 column is headed "Optimum lower bound" and he states
#       "The data of the first line of Table 2 shows that either mu >= 0.38 or
#       E(M*) <= 0.75" — so his printed 0.38 is a rigorous >= 0.380000, used as
#       such by him.  Our 12 Jansson-certified core anchors cover the residual
#       region (5.16) AND lift White's weakest strip (region 18, his 0.37925)
#       above 0.380000.  Nothing else of ours is load-bearing.
#
#   TIER 2  (our own full-space cover)
#       Replaces White's Table 2 with our 121-center cover.  Stronger, but the
#       ~109 non-core centers are still anchored at the uncertified
#       `primal - 1e-5` convention, so it is a working frontier, not a theorem.
#
# The gate regions' stored floors were produced by adaptive subdivision that
# STOPS once it clears its target, so they must be re-run at the new, higher
# target — otherwise they look like hard caps when they are not.
set -euo pipefail

PY=../../.venv/bin/python
SCRATCH="${SCRATCH:-/tmp}"
REANCHORED=../parallel_results/phase5_N20K_bn40_dualext_reanchored.json

echo "=== 1. evaluate the re-anchored core envelope ==="
$PY -u _jansson_reanchor.py --evaluate | tee "$SCRATCH/final_core.log"

echo
echo "=== 2. emit the drop-in re-anchored center file ==="
$PY -u _jansson_reanchor.py --emit-dualext

CORE=$($PY - <<'EOF'
import json, pathlib
d = json.loads(pathlib.Path("../parallel_results/jansson_core12_reanchored.json").read_text())
print(repr(d["evaluation"]["reanchored"]["rigorous_LB"]))
EOF
)
echo "new certified core floor = $CORE"

echo
echo "=== 3. TIER 1: full space with White's published Table 2 ==="
LP_DUALEXT=$REANCHORED $PY -u _fullspace_eval.py | tee "$SCRATCH/final_tier1.log" | tail -25

echo
echo "=== 4. TIER 2: re-run the tightest gate region (R6) at the new target ==="
LP_DUALEXT=$REANCHORED LP_TARGET=$CORE $PY -u _eval_r6_box.py \
  | tee "$SCRATCH/final_r6.log" | tail -6

echo
echo "=== 5. TIER 2: full-space recompute (region floors read from stored JSONs) ==="
LP_DUALEXT=$REANCHORED LP_TARGET=$CORE $PY -u _fs_recompute.py \
  | tee "$SCRATCH/final_tier2.log" | tail -20

echo
echo "logs in $SCRATCH/final_*.log"
