#!/usr/bin/env bash
# Run this on your Mac to commit fresh.
# The Cowork sandbox couldn't finalize commits due to mount permissions
# on .git/index.lock — see the README for context.
set -euo pipefail
cd "$(dirname "$0")"

# Clear any orphaned git state from sandbox attempts.
if [ -d .git ]; then
    rm -rf .git
fi

# Fresh init and commit.
git init -b main
git config user.name "Ben Zanghi"
git config user.email "ben@benzanghi.com"
git add -A
git commit -m "Bochner-PSD strengthening of White (2023): mu >= 0.379544

+5.4e-4 rigorous improvement on E.P. White's published lower bound for the
Erdos minimum overlap constant mu.

  Prior:  White (2023, Acta Arithmetica)         mu >= 0.379005
  This:   Bochner-PSD + ellipse extension        mu >= 0.379544 (+5.4e-4)

Method: adjoin Bochner moment-matrix PSD constraints (M_n(f) and M_n(1-f)
PSD as Hermitian Toeplitz matrices) to White's Section-5 convex program,
then apply White's own Section 5.1 / Appendix II ellipse-extension argument
with the augmented dual objective. The 7 ellipses around White's Table-3
centers, recomputed with our augmented duals, fully cover the residual
region (5.16); their intersected minimum is 0.3795475 (closed-form Path B
KKT), reduced to 0.379544 with conservative 1e-6 IPM-gap margin.

Verification: three independently-written code paths agree to 10+ digits on
per-row SDPs; one independent re-encoding of the Bochner constraint by a
separate agent (no code-sharing) shows bit-for-bit agreement; SDPA-GMP
spot-check confirms CLARABEL is rigorous to ~5e-9 at small N. Status
optimal_inaccurate is a CLARABEL labeling artifact (actual gaps ~1e-7
verified via verbose-output dual extraction).

Note on a withdrawn extension: an earlier Lasserre level-2 augmentation
heuristically lifted the bound to 0.379828 (+8.2e-4) but the implementation
truncates the moment expansion of (f^2)^(m) without a tail bound. Adversarial
review correctly flagged this as non-rigorous. The Lasserre direction
remains promising once a Fejer-Riesz / Parseval tail bound is added; that
is explicit future work documented in the research note.

Open gap: mu in [0.379544, 0.380871], width ~1.3e-3, with the upper
bound from Together Computer (March 2026) via sequential-LP refinement.

Co-authored-by: Claude (Anthropic) <noreply@anthropic.com>"

echo
echo "Local commit done. To push to GitHub:"
echo "  Option A (gh CLI):"
echo "    gh repo create erdos-minimum-overlap-bochner --public --source=. --remote=origin --push"
echo
echo "  Option B (manual):"
echo "    1. Create new repo at github.com/<your-username>/erdos-minimum-overlap-bochner"
echo "    2. git remote add origin git@github.com:<your-username>/erdos-minimum-overlap-bochner.git"
echo "    3. git push -u origin main"
echo
echo "After pushing, update the README badge URL and the preprint repo URL."
