"""One-loop momentum-integration tools (explicit components, metric +---).

Workflow for a one-loop (sub)integral:

1. combine propagator denominators with Feynman parameters
   (`feynman_shift` finds the shift s and gives Delta),
2. substitute k = l + s in the numerator, `symmetrize` in l
   (odd powers drop, l^a l^b -> g^{ab} l^2/4, etc.),
3. replace powers of l^2 by `loop_integral(a, n, Delta)`
   = int d^4 l/(2pi)^4 (l^2)^a / (l^2 - Delta)^n.

Divergent integrals (logarithmic only) are returned with the UV symbol
LUV = log Lambda^2 (Pauli-Villars scale, scheme constant absorbed): the
F_2 projection must be free of LUV, which we assert in the scripts.
"""
from itertools import combinations
from sympy import symbols, Symbol, expand, solve, factorial, I, pi, log, S

METRIC = (1, -1, -1, -1)
LUV = Symbol("LUV")  # log(Lambda^2) + scheme constant


def comps(name):
    """Four contravariant components of a momentum."""
    return list(symbols(f"{name}0 {name}1 {name}2 {name}3"))


def dot(a, b):
    return sum(METRIC[m] * a[m] * b[m] for m in range(4))


def msub(expr, k, val):
    """Substitute momentum components k -> val (both 4-lists)."""
    return expr.subs(dict(zip(k, val)), simultaneous=True)


def feynman_shift_general(D, k):
    """D(k) quadratic in k with k^2 coefficient A (possibly parameter
    dependent). Returns (s, Delta, A) with D(k = l + s) = A (l^2 - Delta).
    """
    from sympy import cancel
    D = expand(D)
    A = D.coeff(k[0], 2)
    s, Delta = feynman_shift(expand(cancel(D / A)), k)
    return s, Delta, A


def feynman_shift(D, k):
    """D(k) quadratic in components of k with unit k^2 coefficient.

    Returns (s, Delta) with D(k = l + s) = l^2 - Delta.
    """
    D = expand(D)
    s = comps("_s")
    l = comps("_l")
    shifted = expand(msub(D, k, [l[m] + s[m] for m in range(4)]))
    lin = [shifted.coeff(l[m], 1) for m in range(4)]
    sol = solve(lin, s, dict=True)
    assert len(sol) == 1
    sval = [sol[0][s[m]] for m in range(4)]
    # check quadratic part is exactly l^2
    quad = shifted - expand(msub(D, k, sval))
    for a in range(4):
        for b in range(4):
            c = quad.coeff(l[a]).coeff(l[b])
    Delta = -expand(msub(D, k, sval))
    check = expand(msub(D, k, [l[m] + sval[m] for m in range(4)])
                   - (dot(l, l) - Delta))
    if check != 0:
        from sympy import simplify, cancel, together
        check = simplify(cancel(together(check)))
    assert check == 0, "denominator not reduced to l^2 - Delta"
    return sval, Delta


def _pairings(idx):
    """All perfect matchings of a list of indices."""
    if not idx:
        yield []
        return
    a = idx[0]
    for j in range(1, len(idx)):
        rest = idx[1:j] + idx[j + 1:]
        for rest_pairs in _pairings(rest):
            yield [(a, idx[j])] + rest_pairs


def symmetrize(expr, k, l2):
    """Angular average over the direction of loop momentum k.

    Monomials odd in any component vanish; even monomials of degree 2n are
    replaced using

      <l^{a1}..l^{a2n}> = (l^2)^n / (2^n (n+1)!) * sum_matchings prod g^{ab}.
    """
    expr = expand(expr)
    terms = expr.as_ordered_terms()
    out = S(0)
    for t in terms:
        c = t
        powers = []
        for m in range(4):
            d = c.as_poly(k[m]).degree() if c.has(k[m]) else 0
            if d:
                c = c / k[m] ** d
                powers.extend([m] * d)
        c = expand(c)
        n2 = len(powers)
        if n2 == 0:
            out += t
            continue
        if n2 % 2 == 1:
            continue
        n = n2 // 2
        # sum over matchings of product of metric entries (diagonal metric:
        # a pairing contributes 0 unless it pairs equal indices)
        tot = S(0)
        for pairs in _pairings(powers):
            contrib = S(1)
            for (a, b) in pairs:
                contrib *= METRIC[a] if a == b else 0
            tot += contrib
        out += c * tot * l2 ** n / (S(2) ** n * factorial(n + 1))
    return out


def loop_integral(a, n, Delta):
    """int d^4 l/(2pi)^4 (l^2)^a / (l^2 - Delta)^n, as multiple of one.

    Finite (n - a - 2 > 0):
        i (-1)^{n+a}/(16 pi^2) (a+1)! (n-a-3)! / (n-1)!  Delta^{a+2-n}
    Log divergent (n - a - 2 == 0), Pauli-Villars:
        i (-1)^{n+a}/(16 pi^2) (a+1)!/(n-1)! (LUV - log Delta)
    """
    pref = I * S(-1) ** (n + a) / (16 * pi**2)
    if n - a - 2 > 0:
        return (pref * factorial(a + 1) * factorial(n - a - 3)
                / factorial(n - 1) * Delta ** (a + 2 - n))
    if n - a - 2 == 0:
        return pref * factorial(a + 1) / factorial(n - 1) * (LUV - log(Delta))
    raise ValueError(f"power divergence: a={a}, n={n}")
