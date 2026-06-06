"""PSLQ closed-form hunting for the Erdős minimum-overlap constant μ.

Usage:
    python3 pslq_hunt.py [value] [precision_digits]

Examples (from the venv):
    python3 lp_research_state/code/pslq_hunt.py 0.3803027 50
    python3 lp_research_state/code/pslq_hunt.py 0.380871 50

The script tries to find an integer-coefficient relation
    a_0 * μ + a_1 * 1 + a_2 * π + a_3 * e + a_4 * log(2) + ... = 0
across a small basis of standard constants.

CAVEAT: PSLQ needs enough decimal precision in μ to constrain coefficients
uniquely. Rule of thumb: to detect a relation with |coeff_max| ≤ C, you
need ~log10(C * basis_size) digits of precision. Our current LB headline
has ~10 reliable digits (CLARABEL `optimal_inaccurate` floor); SDPA-GMP
should yield ~30+ digits, making PSLQ much more powerful.
"""
from __future__ import annotations

import sys
from mpmath import mp, mpf, pi, exp, log, sqrt, pslq, gamma, zeta


def standard_basis(target: mpf) -> list[tuple[str, mpf]]:
    """Return a list of (name, value) pairs for the PSLQ search."""
    return [
        ("μ",      target),
        ("1",      mpf(1)),
        ("π",      pi),
        ("π²",     pi ** 2),
        ("π⁻¹",    1 / pi),
        ("e",      exp(1)),
        ("log 2",  log(2)),
        ("log 3",  log(3)),
        ("log π",  log(pi)),
        ("log(π/2)", log(pi / 2)),
        ("√2",     sqrt(2)),
        ("√3",     sqrt(3)),
        ("√5",     sqrt(5)),
        ("ζ(2)",   zeta(2)),
        ("ζ(3)",   zeta(3)),
        ("Γ(1/4)", gamma(mpf(1) / 4)),
        ("Γ(1/3)", gamma(mpf(1) / 3)),
    ]


def hunt(target_value: float | str, dps: int = 50,
         tol_exp: int = -30, maxcoeff: int = 10 ** 8) -> None:
    mp.dps = dps
    target = mpf(target_value)
    print(f"Hunting closed form for: {target}  (precision: {dps} digits)")
    print(f"PSLQ tolerance: 10^{tol_exp}, maxcoeff: {maxcoeff}\n")

    basis = standard_basis(target)
    values = [v for _, v in basis]
    names = [n for n, _ in basis]
    print("Basis:")
    for n, v in basis:
        print(f"  {n:<10} = {mp.nstr(v, 15)}")
    print()

    # Full-basis search
    print("=== Full-basis PSLQ ===")
    result = pslq(values, tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
    if result is None:
        print("No relation found in full basis.")
    else:
        # Pretty-print the relation
        print("Found relation (coefficients):")
        for n, c in zip(names, result):
            if c != 0:
                sign = " " if c > 0 else "-"
                print(f"  {sign}{abs(c):>10d} * {n}")
        # Verify the relation residual
        s = sum(c * v for c, v in zip(result, values))
        print(f"\nResidual: {mp.nstr(s, 5)} (target ≤ 10^{tol_exp})")

    # Subsets: pairs and triples involving μ
    print("\n=== Pair search: μ + c * x = 0 ===")
    for n, v in basis[1:]:  # skip μ itself
        r = pslq([target, v], tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
        if r is not None:
            a, b = r
            print(f"  {a:>4d} * μ + {b:>4d} * {n} = 0  →  μ = {-b}/{a} * {n}")

    print("\n=== Triple search: μ + a·c1 + b·c2 = 0 ===")
    for i in range(1, len(basis)):
        for j in range(i + 1, len(basis)):
            r = pslq([target, basis[i][1], basis[j][1]],
                     tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
            if r is not None and abs(r[0]) <= 100:
                print(f"  {r[0]:>4d}·μ  +  {r[1]:>4d}·{basis[i][0]}  +  "
                      f"{r[2]:>4d}·{basis[j][0]}  =  0")


if __name__ == "__main__":
    val = sys.argv[1] if len(sys.argv) > 1 else "0.3803027"
    dps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    hunt(val, dps=dps)
