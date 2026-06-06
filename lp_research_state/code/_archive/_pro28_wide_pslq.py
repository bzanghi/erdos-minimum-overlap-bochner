"""PRO-28: Wide-basis PSLQ for μ closed form.

Extends pslq_hunt.py with a much wider basis:
- Hurwitz zeta values ζ(s, q) for various rational q
- Polylogarithms Li_s(z) at small algebraic z
- Dirichlet L-function values L(s, χ_d)
- Bessel function zeros j_{ν,k}
- Apéry-like and known transcendentals
- Algebraic constants we might have missed

Targets:
- μ_UB at 50 digits (computed from Together's h*)
- Optionally: μ_LB-related quantities
"""
from __future__ import annotations
from pathlib import Path

from mpmath import (
    mp, mpf, pi, e, log, sqrt, exp, pslq, gamma, zeta, catalan, glaisher,
    euler, besselj, sin, cos, polylog, mpc
)

# Set high working precision
mp.dps = 60


# 50-digit-precise UB anchor (from PRO-26 scaffolding)
UB_50 = mpf("0.38087031058621710878661081496601738896393463045218")


def build_wide_basis(target: mpf) -> list[tuple[str, mpf]]:
    """Construct a wide basis of named constants."""
    basis = [
        ("μ",       target),
        # Algebraic / classical
        ("1",       mpf(1)),
        ("π",       pi),
        ("π²",      pi ** 2),
        ("π³",      pi ** 3),
        ("π⁴",      pi ** 4),
        ("1/π",     1 / pi),
        ("1/π²",    1 / pi ** 2),
        ("e",       e),
        ("e²",      e ** 2),
        ("eπ",      e * pi),
        ("e/π",     e / pi),
        # Logs
        ("log 2",   log(2)),
        ("log 3",   log(3)),
        ("log 5",   log(5)),
        ("log π",   log(pi)),
        ("log²π",   log(pi) ** 2),
        ("log(2π)", log(2 * pi)),
        ("log(π/2)", log(pi / 2)),
        # Roots
        ("√2",      sqrt(2)),
        ("√3",      sqrt(3)),
        ("√5",      sqrt(5)),
        ("√π",      sqrt(pi)),
        ("∛2",      mpf(2) ** (mpf(1) / 3)),
        # Riemann zeta
        ("ζ(2)",    zeta(2)),
        ("ζ(3)",    zeta(3)),
        ("ζ(4)",    zeta(4)),
        ("ζ(5)",    zeta(5)),
        ("ζ(7)",    zeta(7)),
        # Hurwitz zeta — new basis additions
        ("ζ(2,1/2)", zeta(2, mpf(1)/2)),
        ("ζ(2,1/3)", zeta(2, mpf(1)/3)),
        ("ζ(2,1/4)", zeta(2, mpf(1)/4)),
        ("ζ(3,1/2)", zeta(3, mpf(1)/2)),
        ("ζ(3,1/4)", zeta(3, mpf(1)/4)),
        # Gamma values
        ("Γ(1/4)",  gamma(mpf(1) / 4)),
        ("Γ(1/3)",  gamma(mpf(1) / 3)),
        ("Γ(1/6)",  gamma(mpf(1) / 6)),
        ("Γ(2/3)",  gamma(mpf(2) / 3)),
        ("Γ(3/4)",  gamma(mpf(3) / 4)),
        # Polylogs at 1/2
        ("Li₂(1/2)", polylog(2, mpf(1) / 2)),
        ("Li₃(1/2)", polylog(3, mpf(1) / 2)),
        ("Li₄(1/2)", polylog(4, mpf(1) / 2)),
        ("Li₅(1/2)", polylog(5, mpf(1) / 2)),
        # Polylogs at -1
        ("Li₂(-1)", polylog(2, mpf(-1))),  # = -π²/12
        ("Li₃(-1)", polylog(3, mpf(-1))),
        ("Li₄(-1)", polylog(4, mpf(-1))),
        # Special constants
        ("γ",       euler),         # Euler-Mascheroni
        ("Catalan", catalan),
        ("Glaisher", glaisher),
        # Bessel zeros (first few of J_0, J_1)
        ("j_{0,1}/π", besseljzero(0, 1) / pi if False else mpf(0)),
        # NOTE: Bessel zeros are non-trivial to compute reliably in mpmath
        # without the besseljzero function; use placeholder
    ]
    # Filter out placeholders / None
    return [(n, v) for n, v in basis if v != 0 or n == "0"]


def besseljzero(n: int, k: int) -> mpf:
    """Compute the k-th positive zero of J_n."""
    try:
        from mpmath import besseljzero as _bjz
        return _bjz(n, k)
    except Exception:
        return mpf(0)


def build_full_basis(target: mpf) -> list[tuple[str, mpf]]:
    """Final wide basis with Bessel zeros included."""
    base = [
        ("μ",       target),
        ("1",       mpf(1)),
        ("π",       pi),
        ("π²",      pi ** 2),
        ("1/π",     1 / pi),
        ("e",       e),
        ("eπ",      e * pi),
        ("log 2",   log(2)),
        ("log 3",   log(3)),
        ("log π",   log(pi)),
        ("log²π",   log(pi) ** 2),
        ("log(2π)", log(2 * pi)),
        ("√2",      sqrt(2)),
        ("√3",      sqrt(3)),
        ("√5",      sqrt(5)),
        ("ζ(2)",    zeta(2)),
        ("ζ(3)",    zeta(3)),
        ("ζ(5)",    zeta(5)),
        ("ζ(2,1/4)", zeta(2, mpf(1)/4)),
        ("ζ(3,1/2)", zeta(3, mpf(1)/2)),
        ("Γ(1/4)",  gamma(mpf(1) / 4)),
        ("Γ(1/3)",  gamma(mpf(1) / 3)),
        ("Γ(1/6)",  gamma(mpf(1) / 6)),
        ("Li₂(1/2)", polylog(2, mpf(1) / 2)),
        ("Li₃(1/2)", polylog(3, mpf(1) / 2)),
        ("Li₄(1/2)", polylog(4, mpf(1) / 2)),
        ("γ",       euler),
        ("Catalan", catalan),
        ("Glaisher", glaisher),
    ]
    # Try to add Bessel zeros — only if available
    j01 = besseljzero(0, 1)
    j02 = besseljzero(0, 2)
    j11 = besseljzero(1, 1)
    if j01 != 0:
        base.extend([
            ("j_{0,1}", j01),
            ("j_{0,2}", j02),
            ("j_{1,1}", j11),
            ("j_{0,1}/π", j01 / pi),
        ])
    return base


def run_pslq(basis, tol_exp=-45, maxcoeff=10**8, verbose=True):
    """Run full-basis PSLQ. Return found relation or None."""
    names = [n for n, _ in basis]
    vals = [v for _, v in basis]
    if verbose:
        print(f"  Basis size: {len(basis)} constants")
        print(f"  tol: 10^{tol_exp}, maxcoeff: 10^{int(round(__import__('math').log10(maxcoeff)))}")
    res = pslq(vals, tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
    if res is None:
        return None
    # Pretty-print
    if verbose:
        print("  Relation found:")
        for n, c in zip(names, res):
            if c != 0:
                sign = " " if c >= 0 else "-"
                print(f"    {sign}{abs(c):>8d} · {n}")
        # Residual
        s = sum(c * v for c, v in zip(res, vals))
        print(f"  Residual: {mp.nstr(s, 5)}")
    return res


def pair_search(target, basis, tol_exp=-45, maxcoeff=10**10, verbose=True):
    """All pair searches a·μ + b·c = 0."""
    hits = []
    for n, v in basis[1:]:
        r = pslq([target, v], tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
        if r is not None:
            a, b = r
            # Verify accuracy
            if a != 0:
                pred = -mpf(b) / a
                err = float(abs(pred - target))
                if err < 10 ** (tol_exp + 5):
                    hits.append((n, a, b, pred, err))
                    if verbose:
                        print(f"    PAIR: {a} · μ + {b} · {n} = 0  → μ = {mp.nstr(pred, 20)}  err={err:.2e}")
    return hits


def triple_search(target, basis, tol_exp=-45, maxcoeff=10**5, verbose=True):
    """Triple searches a·μ + b·c1 + d·c2 = 0 with small coefficients."""
    hits = []
    for i in range(1, len(basis)):
        for j in range(i + 1, len(basis)):
            r = pslq([target, basis[i][1], basis[j][1]],
                     tol=mpf(10) ** tol_exp, maxcoeff=maxcoeff)
            if r is not None:
                a, b, d = r
                if abs(a) <= 100 and abs(b) <= 200 and abs(d) <= 200 and a != 0:
                    # Compute the prediction and check
                    pred = -(mpf(b) * basis[i][1] + mpf(d) * basis[j][1]) / a
                    err = float(abs(pred - target))
                    # Filter out the redundant ones (involving only known relations)
                    if abs(err) < 10 ** (tol_exp + 5):
                        # Skip relations not involving μ
                        if a != 0:
                            hits.append((basis[i][0], basis[j][0], a, b, d, pred, err))
    return hits


def main():
    print(f"=== PRO-28: Wide-basis PSLQ for μ closed form ===\n")
    print(f"Target: μ_UB at 50 digits = {UB_50}\n")

    basis = build_full_basis(UB_50)
    print(f"=== Full-basis search ({len(basis)} constants, maxcoeff 1e8) ===")
    res = run_pslq(basis, tol_exp=-45, maxcoeff=10**8)
    if res is None:
        print("  No relation found.\n")

    print(f"=== Same basis, maxcoeff 1e10 (harder relations) ===")
    res2 = run_pslq(basis, tol_exp=-45, maxcoeff=10**10)
    if res2 is None:
        print("  No relation found.\n")

    print(f"=== Pair searches (maxcoeff 1e10) ===")
    pair_hits = pair_search(UB_50, basis, tol_exp=-45, maxcoeff=10**10)
    print(f"  ({len(pair_hits)} pair relations found)\n")

    print(f"=== Triple searches (small coeffs ≤ 100) ===")
    triple_hits = triple_search(UB_50, basis, tol_exp=-45, maxcoeff=10**5)
    # Filter: only μ-involving and high-precision
    μ_hits = [h for h in triple_hits if abs(h[2]) > 0]
    print(f"  ({len(μ_hits)} triple relations involving μ found)\n")
    for n1, n2, a, b, d, pred, err in μ_hits[:15]:
        print(f"    {a}·μ + {b}·{n1} + {d}·{n2} = 0  pred={mp.nstr(pred, 18)}, err={err:.2e}")

    # Save the basis used for reproducibility
    import json
    out_path = Path(__file__).parent.parent / "data" / "pro28_wide_pslq.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_data = {
        "target_value": str(UB_50),
        "target_label": "μ_UB at 50 digits (from Together's h*)",
        "basis_constants": [{"name": n, "value": mp.nstr(v, 30)} for n, v in basis],
        "full_basis_result_maxcoeff_1e8": "no relation" if res is None else str(res),
        "full_basis_result_maxcoeff_1e10": "no relation" if res2 is None else str(res2),
        "pair_hits": [{"const": n, "a": int(a), "b": int(b), "pred": mp.nstr(p, 30), "err": e}
                      for n, a, b, p, e in pair_hits],
        "triple_hits_involving_mu": [
            {"const1": n1, "const2": n2, "a": int(a), "b": int(b), "d": int(d),
             "pred": mp.nstr(p, 30), "err": e}
            for n1, n2, a, b, d, p, e in μ_hits[:50]
        ],
    }
    with open(out_path, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"\nResults saved to: {out_path}")


if __name__ == "__main__":
    main()
