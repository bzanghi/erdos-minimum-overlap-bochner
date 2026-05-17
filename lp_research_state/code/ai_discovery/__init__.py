"""AI-driven constraint discovery for the Erdős minimum-overlap SDP.

Goal: find a new convex constraint family that pushes the SDP lower bound
on μ by ≥ 1×10⁻⁴, beyond the 10/10 prior levers.

Submodules:
  fast_eval     -- small-N SDP evaluator (~0.1s/solve at N=200)
  dsl           -- constraint DSL: parametric families generating cvxpy expr
  search        -- MLX-accelerated search loop over family parameters
  proposer      -- (Phase 4 fallback) MLX neural net proposer
"""
