"""
Phase 1.2 — Cross-check the continuous reformulation against discrete computation.

Test case: h = 1_{[1/2, 3/2]} (the 'rectangle' density on [0,2] with ∫h=1).

Two candidate continuous formulas:
  (a) Brief-as-written:   J_a(t) = ∫_0^2 h(x)(1 − h(x+t)) dx, h extended by 0 outside.
  (b) Corrected:           J_b(t) = ∫_0^{2−t} h(x)(1 − h(x+t)) dx, NO extension.

Discrete: A = middle half of [1, 2n], B = outer halves.
For this partition, max_k M_k(A,B)/n is computed.
"""
import numpy as np
from itertools import product

def discrete_M_over_n(n):
    """Rectangle partition: A = (n/2, 3n/2]."""
    A = set(range(n//2 + 1, 3*n//2 + 1))
    B = set(range(1, 2*n + 1)) - A
    assert len(A) == n and len(B) == n
    counts = {}
    for a in A:
        for b in B:
            counts[a-b] = counts.get(a-b, 0) + 1
    M = max(counts.values())
    kstar = max(counts, key=counts.get)
    return M / n, kstar / n

def J_brief(t, dx=1e-4):
    """Formula (a): integrate over full [0,2], extend h by 0."""
    xs = np.arange(0, 2, dx)
    h_x = ((xs >= 0.5) & (xs <= 1.5)).astype(float)
    xs_shift = xs + t
    h_xt = ((xs_shift >= 0.5) & (xs_shift <= 1.5)).astype(float)  # =0 outside [0,2] AND outside [1/2,3/2]
    return float(np.sum(h_x * (1 - h_xt)) * dx)

def J_corrected(t, dx=1e-4):
    """Formula (b): integrate only on [0, 2−t]; no extension needed."""
    if t < 0 or t >= 2:
        return 0.0
    xs = np.arange(0, 2 - t, dx)
    h_x = ((xs >= 0.5) & (xs <= 1.5)).astype(float)
    xs_shift = xs + t
    h_xt = ((xs_shift >= 0.5) & (xs_shift <= 1.5)).astype(float)
    return float(np.sum(h_x * (1 - h_xt)) * dx)

print("Discrete (rectangle partition):")
for n in [50, 200, 1000]:
    r, ts = discrete_M_over_n(n)
    print(f"  n={n}: M(A,B)/n = {r:.5f}, achieved at t=k/n={ts:.4f}")

print("\nContinuous formula (a) [brief, h extended by 0]: sup over t∈[0,2]")
ts = np.linspace(0.0, 1.999, 401)
Ja = [J_brief(t) for t in ts]
print(f"  sup_t J_a(t) ≈ {max(Ja):.5f}  (achieved t≈{ts[int(np.argmax(Ja))]:.3f})")

print("\nContinuous formula (b) [corrected, no extension]: sup over t∈[0,2]")
Jb = [J_corrected(t) for t in ts]
print(f"  sup_t J_b(t) ≈ {max(Jb):.5f}  (achieved t≈{ts[int(np.argmax(Jb))]:.3f})")
