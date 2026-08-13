"""IIc piece (a): s-integration of gchi (from g2_iic_analytic.py) over
s in (0, 1-t), using complex-linear letters.

Input:  g2_iic_a_gchi.pkl  (3-dim integrand, verified)
Output: g2_iic_a_Ga.pkl    (2-dim integrand Ga(t, u, R), R = sqrt(tu(1-u)),
                            complex-log/dilog form, verified numerically)

Method: atan(g) = (i/2)[log(D - iN) - log(D + iN)] for g = N/D with D > 0;
log(quadratic) = log(leading) + log(s - r+) + log(s - r-) with the exact
complex roots (rational in R); partial fractions in s over Q(t,u,R,i);
base integrals: elementary + the dilog primitive

    int log(s-b)/(s-a) ds
      = log(s-b) log((s-a)/(b-a)) + Li2((b-s)/(b-a)),

verified by differentiation, with numeric continuity checks per term.
"""
import pickle
import time

from sympy import (symbols, Symbol, S, I, pi, log, atan, sqrt, expand,
                   together, simplify, cancel, fraction, apart, factor,
                   lambdify, Rational, Poly, Add, polylog, diff, im,
                   factor_list, degree, Pow)
from mpmath import mp, quad as mquad, mpf, mpc

mp.dps = 25
T0 = time.time()


def tick(m):
    print("[%6.1fs] %s" % (time.time() - T0, m), flush=True)


s, t, u = symbols("s t u", positive=True)
R = Symbol("R", positive=True)

gchi = pickle.load(open("g2_iic_a_gchi.pkl", "rb"))
sq_args = list({b for a in gchi.atoms(Pow) for b, e in [a.as_base_exp()]
                if e in (Rational(1, 2), Rational(-1, 2))})
assert len(sq_args) == 1
big = sq_args[0]
g2 = gchi.subs(sqrt(big), 2 * R * (s * (1 - u) + t))

# ---- letters ----
alpha = R / (u * (1 - u))          # Q1 roots: +- i alpha
beta = R / (1 - u)                 # Qu roots: -t +- i beta
Q1 = u * (u - 1) ** 2 * s**2 + t * (1 - u)
Qu = u * (u - 1) ** 2 * (s + t) ** 2 + t * u**2 * (1 - u)
assert cancel(expand(Q1 - u * (u - 1) ** 2 * (s - I * alpha)
              * (s + I * alpha)).subs(R**2, t * u * (1 - u))) == 0
assert cancel(expand(Qu - u * (u - 1) ** 2 * (s + t - I * beta)
              * (s + t + I * beta)).subs(R**2, t * u * (1 - u))) == 0
tick("root structure verified")

# ---- rewrite transcendentals ----
rew = {}
for a_ in g2.atoms(atan):
    g = together(a_.args[0])
    N, D = fraction(g)
    # atan(N/D) = (i/2)[log(D - iN) - log(D + iN)], D > 0 on domain
    rew[a_] = (I / 2) * (log(expand(D - I * N)) - log(expand(D + I * N)))
for L in g2.atoms(log):
    arg = expand(L.args[0])
    if arg == expand(Q1):
        rew[L] = (log(u * (u - 1) ** 2) + log(s - I * alpha)
                  + log(s + I * alpha))
    elif arg == expand(Qu):
        rew[L] = (log(u * (u - 1) ** 2) + log(s + t - I * beta)
                  + log(s + t + I * beta))
    else:
        raise AssertionError(f"unknown log arg {arg}")
g3 = g2.subs(rew, simultaneous=True)

# numeric spot check of the rewrite
pt = {s: mpf("0.3"), t: mpf("0.2"), u: mpf("0.4")}
Rv = (pt[t] * pt[u] * (1 - pt[u])) ** mpf("0.5")
f2 = lambdify((s, t, u, R), g2, "mpmath")
f3 = lambdify((s, t, u, R), g3, "mpmath")
v2 = f2(pt[s], pt[t], pt[u], Rv)
v3 = f3(pt[s], pt[t], pt[u], Rv)
assert abs(complex(v2) - complex(v3)) < 1e-20, (v2, v3)
tick("letter rewrite verified")

# ---- split into rational * log terms, partial fraction in s ----
g3 = expand(g3)
LOGS = sorted(g3.atoms(log), key=str)
terms = []          # (rational_coefficient, log_atom or 1)
rest = g3
for L in LOGS:
    c = rest.coeff(L)
    if c != 0:
        terms.append((together(c), L))
        rest = expand(rest - c * L)
assert not rest.has(log)
terms.append((together(rest), S(1)))
tick("split: %d log groups + rational" % (len(terms) - 1))


def sub_bounds(F, lo, hi):
    return F.subs(s, hi) - F.subs(s, lo)


def anti_pole(a, k):
    """antiderivative of (s-a)^-k"""
    if k == 1:
        return log(s - a)
    return (s - a) ** (1 - k) / (1 - k)


def anti_log_pole(b, a, k):
    """antiderivative of log(s-b) * (s-a)^-k, verified by diff."""
    if k == 1:
        F = log(s - b) * log((s - a) / (b - a)) + polylog(2, (b - s) / (b - a))
    else:
        F = (log(s - b) * (s - a) ** (1 - k) / (1 - k)
             - anti_rat_over(b, a, k) / (1 - k))
    check = diff(F, s) - log(s - b) * (s - a) ** (-k)
    check = check.replace(lambda e: e.func == polylog and e.args[0] == 1,
                          lambda e: -log(1 - e.args[1]))
    check = simplify(cancel(together(check)))
    if check != 0:
        # log-combination residue: verify numerically instead
        import mpmath as _mp
        for (tv_, uv_, sv_) in [(0.2, 0.4, 0.11), (0.55, 0.7, 0.31)]:
            Rv_ = (tv_ * uv_ * (1 - uv_)) ** 0.5
            zz = dict(zip((s, t, u, R), (sv_, tv_, uv_, Rv_)))
            cv = complex(check.subs(zz))
            assert abs(cv) < 1e-9, (b, a, k, cv)
        check = S(0)
    assert check == 0, (b, a, k, check)
    return F


def anti_rat_over(b, a, k):
    """antiderivative of (s-a)^{1-k}/(s-b), k >= 2 (partial fractions)."""
    e = (s - a) ** (1 - k) / (s - b)
    return anti_rational_in_s(e)


def anti_rational_in_s(expr):
    from ratint import anti_rational
    return anti_rational(expr, s)


Ga = S(0)
lo_, hi_ = S(0), 1 - t
for (coef, L) in terms:
    tt0 = time.time()
    if L == 1:
        # plain rational in s
        F = anti_rational_in_s(coef)
        Ga += sub_bounds(F, lo_, hi_)
        tick("rational part done")
        continue
    # L = log(c1 (s - b)): normalize; branch offset eta checked numerically
    pol = Poly(L.args[0], s)
    assert pol.degree() == 1, L
    c1_, c0_ = pol.all_coeffs()
    b_ = cancel(-c0_ / c1_)
    if c1_ == 1:
        const_log = S(0)
    else:
        import mpmath as _mp
        etas = set()
        for (tv_, uv_) in [(0.2, 0.4), (0.55, 0.7), (0.35, 0.15)]:
            Rv_ = (tv_ * uv_ * (1 - uv_)) ** 0.5
            for sv_ in (0.1 * (1 - tv_), 0.5 * (1 - tv_), 0.9 * (1 - tv_)):
                zz = dict(zip((s, t, u, R), (sv_, tv_, uv_, Rv_)))
                lhs = complex(L.args[0].subs(zz))
                c1n = complex(c1_.subs(zz))
                bn = complex(b_.subs(zz))
                eta_ = (_mp.log(lhs) - _mp.log(c1n)
                        - _mp.log(complex(sv_) - bn)) / (2j * _mp.pi)
                etas.add(int(round(float(eta_.real))))
                assert abs(eta_.imag) < 1e-9
        assert len(etas) == 1, (L, etas)
        const_log = log(c1_) + 2 * I * pi * etas.pop()
    ap = apart(coef, s)
    contrib = S(0)
    if const_log != 0:
        Fc_ = anti_rational_in_s(coef)
        contrib += const_log * sub_bounds(Fc_, lo_, hi_)
    for tm in Add.make_args(ap):
        numf, denf = fraction(together(tm))
        p = Poly(denf, s)
        if p.degree() == 0:
            # polynomial * log(s-b): by parts, elementary
            P = Poly(numf / denf, s).integrate(s).as_expr()
            F = P * log(s - b_) - anti_rational_in_s(cancel(P / (s - b_)))
            contrib += sub_bounds(F, lo_, hi_)
            continue
        fl = factor_list(denf, s)
        cc = fl[0]
        for (bb, ee) in fl[1]:
            if Poly(bb, s).degree() == 0:
                cc *= bb ** ee
        facs = [(bb, ee) for (bb, ee) in fl[1] if Poly(bb, s).degree() > 0]
        assert len(facs) == 1 and Poly(facs[0][0], s).degree() == 1, tm
        base, k = facs[0]
        c1, c0 = Poly(base, s).all_coeffs()
        a_ = -c0 / c1
        const = numf / cc / c1 ** k
        assert not const.has(s)
        F = const * anti_log_pole(b_, a_, k)
        contrib += sub_bounds(F, lo_, hi_)
    Ga += contrib
    tick("log group %s done (%.0fs)" % (str(L)[:40], time.time() - tt0))

pickle.dump(Ga, open("g2_iic_a_Ga.pkl", "wb"))
tick("Ga assembled, ops = %d" % Ga.count_ops())

# emit Ga as Fortran for the quad-precision integrator g2_iic_a_quad.f90:
# polylog -> cdli2, log -> logc (complex-promoting), cmplx(0,1) -> CI
from sympy import fcode
_code = fcode(Ga, assign_to="gaval", source_format="free", standard=2008,
              user_functions={"polylog": "cdli2"})
_code = _code.replace("cmplx(0,1)", "CI")
import re as _re
_code = _re.sub(r"(?<![a-z_0-9])log\(", "logc(", _code)
with open("g2_iic_a_ga.inc", "w") as fh:
    fh.write("! generated by g2_iic_a_sint.py -- do not edit\n")
    fh.write(_code + "\n")
tick("wrote g2_iic_a_ga.inc")

# ---- numeric verification: s-quad of gchi vs Ga at two (t,u) points ----
gchif = lambdify((s, t, u), gchi, "mpmath")
Gaf = lambdify((t, u, R), Ga, "mpmath")
for tv, uv in [(mpf("0.2"), mpf("0.4")), (mpf("0.55"), mpf("0.7"))]:
    Rv = (tv * uv * (1 - uv)) ** mpf("0.5")
    direct = mquad(lambda ss: gchif(ss, tv, uv), [0, 1 - tv])
    mine = Gaf(tv, uv, Rv)
    print("check (t=%s,u=%s): %s vs %s" % (tv, uv, direct, mine), flush=True)
    assert abs(complex(direct) - complex(mine)) < 1e-15
tick("s-integration verified numerically")
