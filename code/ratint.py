"""Deterministic integration of rational functions for the g-2 project.

SymPy's integrate() on rational functions with symbolic parameters is
slow and occasionally silently wrong (see 10-qed-g-2-NLO.md); this module
builds antiderivatives from explicit textbook formulas and *verifies them
by differentiation*, raising if the check fails.

Supported term shapes after apart(): polynomial, num/(linear)^k,
(alpha x + beta)/(quadratic)^k for k = 1, 2 with the quadratic
irreducible (over the rationals in the remaining symbols).
"""
from sympy import (Add, S, Poly, fraction, together, factor_list, apart,
                   integrate, log, atan, sqrt, diff, simplify, limit)


def anti_rational(expr, x, verify="sym"):
    """Antiderivative in x of a rational expr; verified by diff().

    verify="sym": symbolic zero-check (default);
    verify="numeric": exact rational-point checks (for coefficient fields
    where simplify() is too slow, e.g. extra symbolic parameters)."""
    terms = Add.make_args(apart(together(expr), x, full=False))
    F = S(0)
    for tm in terms:
        num, den = fraction(together(tm))
        p = Poly(den, x)
        if p.degree() == 0:
            F += integrate(tm, x)  # polynomial in x: safe
            continue
        f = factor_list(den, x)
        cc = f[0]
        facs = [(b, e) for (b, e) in f[1] if Poly(b, x).degree() > 0]
        assert len(facs) == 1, f"unexpected denominator {den}"
        base, k = facs[0]
        for (b, e) in f[1]:
            if Poly(b, x).degree() == 0:
                cc *= b ** e
        b = Poly(base, x)
        num = num / cc
        if b.degree() == 1:
            a1, c1 = b.all_coeffs()
            nn = Poly(num, x)
            assert nn.degree() == 0, f"apart left degree>0 over linear: {tm}"
            if k == 1:
                F += num / a1 * log(base)
            else:
                F += num / a1 * base ** (1 - k) / (1 - k)
            continue
        assert b.degree() == 2, f"unsupported base {base}"
        a2, b2, c2 = b.all_coeffs()
        D = 4 * a2 * c2 - b2 ** 2  # > 0 iff irreducible
        nn = Poly(num, x)
        assert nn.degree() <= 1, tm
        alpha = nn.coeff_monomial(x) if nn.degree() == 1 else S(0)
        beta = nn.coeff_monomial(1)
        I1 = 2 / sqrt(D) * atan((2 * a2 * x + b2) / sqrt(D))
        if k == 1:
            F += (alpha / (2 * a2) * log(base)
                  + (beta - alpha * b2 / (2 * a2)) * I1)
        elif k == 2:
            I2 = (2 * a2 * x + b2) / (D * base) + 2 * a2 / D * I1
            J2 = -(b2 * x + 2 * c2) / (D * base) - b2 / D * I1
            F += alpha * J2 + beta * I2
        else:
            raise NotImplementedError(f"power {k} over quadratic")
    check = together(diff(F, x) - expr)
    if verify == "sym":
        assert simplify(check) == 0, "antiderivative check failed"
    else:
        import random
        from sympy import Rational as _R
        random.seed(3)
        free = sorted(check.free_symbols, key=str)
        for _ in range(3):
            pt = {sym: _R(random.randint(5, 89), 97) for sym in free}
            val = check.subs(pt, simultaneous=True)
            assert simplify(val) == 0, "antiderivative check failed (num)"
    return F


def defint(expr, x, lo, hi):
    """Definite integral of a rational function via anti_rational."""
    from sympy import oo
    F = anti_rational(expr, x)
    Fhi = limit(F, x, hi, '-') if hi is not None else F
    Flo = limit(F, x, lo, '+') if lo is not None else F
    return Fhi - Flo
