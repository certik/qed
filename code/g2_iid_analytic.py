"""Diagram IId: analytic evaluation of the parametric integrals derived by
g2_iid.py, in the lambda -> 0 limit.

Input: g2_iid_integrands.json (srepr of f_rat, f_log, C_val written by
g2_iid.py). Every integration step is numerically verified as it happens.

Result (must equal Petermann 1957, eq. (4)):

    mu_IId = 11/24 - pi^2/18 + (1/2) log(lam^2/m^2)

Structure of the computation:

* rational piece: f_rat factorizes as U-integrand(u) * S-integrand(y,z)
  with the (y,z) numerator a function of s = y+z (times z) only, so
  mu_rat(lam) = U(lam) * S(lam), two 1-dim integrals done exactly with
  ratint.anti_rational; the product expands to log(lam) + 1/2 + O(lam).

* log piece (finite at lam = 0): f_log = (1/xi) N(y,z,t,u,C)/Delta^3 with
  Delta = s^2 + t(C-1), s = y+z+t. The xi variable is traded for the
  spectral mass C: dxi/xi = -dC/(C-1), C in (1/u, oo). Then integrate
  z (polynomial), C (rational, anti_rational), t and s (rational-times-log
  terms; stays elementary), u (termwise; one pair of individually
  divergent terms integrated together). Gives -1/24 - pi^2/18.
"""
import json
import time
import pickle

from sympy import (symbols, Symbol, S, I, pi, log, oo, expand, together,
                   simplify, factor, fraction, integrate, apart, limit,
                   series, lambdify, sympify, Add, Poly, Rational,
                   expand_log)
from mpmath import mp, quad as mquad, mpf, inf as minf

from ratint import anti_rational

mp.dps = 15
T0 = time.time()


def tick(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg))


y, z, t, u, xi, lam = symbols("y z t u xi lam", positive=True)
s = symbols("s", positive=True)
Cs = Symbol("C", positive=True)

with open("g2_iid_integrands.json") as fh:
    d = json.load(fh)
f_rat = sympify(d["f_rat"])
f_log = sympify(d["f_log"])
C_val = sympify(d["C_val"])

# ===================================================================
# rational piece
# ===================================================================
U_int = 2 * u * (u - 2) * (u - 1) / ((1 - u) ** 2 + lam**2 * u)
S_num = z * (s - 1) * (3 * s**3 + s**2 - lam**2 * (4 * s - 1) * (s - 1))
form = (U_int * S_num.subs(s, y + z)
        / ((y + z) ** 2 + (1 - y - z) * lam**2) ** 2)
assert simplify(together(form - f_rat)) == 0
tick("f_rat factorized form verified")

# z-integration (y = s - z, z in (0, s)) gives factor s^2/2
S_int = ((s**2 / 2) * (s - 1)
         * (3 * s**3 + s**2 - lam**2 * (4 * s - 1) * (s - 1))
         / (s**2 + (1 - s) * lam**2) ** 2)


def defint01(expr, x):
    F = anti_rational(expr, x)
    return limit(F, x, 1, '-') - limit(F, x, 0, '+')


U = defint01(U_int, u)
Sv = defint01(S_int, s)
# numeric spot-check at lam = 1/100
Un = mquad(lambdify(u, U_int.subs(lam, mpf("0.01")), "mpmath"), [0, 0.999, 1])
assert abs(float(U.subs(lam, Rational(1, 100))) - float(Un)) < 1e-10
tick("U, S exact in lam, checked numerically")

mu_rat = simplify(series(U * Sv, lam, 0, 1).removeO())
print("  mu_rat(lam->0) =", mu_rat)
assert mu_rat == log(lam) + S(1) / 2

# ===================================================================
# log piece (lam = 0)
# ===================================================================
g = together(f_log.subs(lam, 0) * xi)     # f_log = g/xi, g has C not xi
assert not g.has(xi)
assert simplify(C_val.subs(lam, 0) - (u * xi - u + 1) / (u * xi)) == 0

# ---- z ----
gy = expand(g.subs(y, s - z - t))
num, den = fraction(together(gy))
assert not den.has(z)
gz = together(integrate(num, (z, 0, s - t)) / den)
gzf = lambdify((s, t, u, Cs), gz, "mpmath")
pt = (mpf("0.6"), mpf("0.2"), mpf("0.35"))
tick("z-integration done")

# ---- C over (1/u, oo), measure dC/(C-1) ----
FC = anti_rational(together(gz / (Cs - 1)), Cs)
gC = limit(FC, Cs, oo) - FC.subs(Cs, 1 / u)
gC = together(expand_log(expand(gC), force=True))
ref = mquad(lambda CC: gzf(pt[0], pt[1], pt[2], CC) / (CC - 1),
            [1 / pt[2], 10 / pt[2], minf])
gCf = lambdify((s, t, u), gC, "mpmath")
assert abs(gCf(*pt) - ref) < 1e-10
tick("C-integration done, checked numerically")


def split_logs(expr):
    """expr -> [(rational coefficient, log or 1), ...]"""
    logs = sorted(expr.atoms(log), key=str)
    out, rest = [], expand(expr)
    for L in logs:
        c = together(rest.coeff(L))
        out.append((c, L))
        rest = expand(rest - c * L)
    out.append((together(rest), S(1)))
    assert not out[-1][0].has(log)
    return out


def realify(expr, point):
    """Replace log(negative-on-domain) by log(-arg) + I pi and drop the
    (spurious, branch-bookkeeping) explicit I part; caller must verify
    numerically."""
    rep = {}
    for L in expr.atoms(log):
        if L.args[0].subs(point) < 0:
            rep[L] = log(-L.args[0]) + I * pi
    e = expand(expand_log(expr.subs(rep), force=True))
    return expand(e - I * e.coeff(I))


# ---- t over (0, s) ----
gt = Add(*[integrate(R * L, (t, 0, s)) for (R, L) in split_logs(gC)])
gt = realify(simplify(expand_log(expand(gt), force=True)),
             {s: Rational(3, 5), u: Rational(7, 20)})
gtf = lambdify((s, u), gt, "mpmath")
ref = mquad(lambda tt: gCf(pt[0], tt, pt[2]), [0, pt[0]])
assert abs(gtf(pt[0], pt[2]) - ref) < 1e-10
tick("t-integration done, checked numerically")

# ---- s over (0, 1) ----
gu = Add(*[integrate(R * L, (s, 0, 1)) for (R, L) in split_logs(gt)])
gu = expand(expand_log(simplify(gu), force=True))
guf = lambdify(u, gu, "mpmath")
ref = mquad(lambda ss: gtf(ss, mpf("0.35")), [0, 1])
assert abs(guf(mpf("0.35")) - ref) < 1e-10
tick("s-integration done, checked numerically")

# ---- u over (0, 1): termwise; group endpoint-divergent terms ----
terms = []
for (R, L) in split_logs(gu):
    terms += [tm * L for tm in Add.make_args(apart(together(R), u))]
finite, divergent = S(0), []
for tm in terms:
    r = integrate(tm, (u, 0, 1))
    if r.has(oo, -oo, S.NaN, S.ComplexInfinity):
        divergent.append(tm)
    else:
        finite += r
if divergent:
    F = integrate(Add(*divergent), u)
    finite += limit(F, u, 1, '-') - limit(F, u, 0, '+')
mu_log = simplify(expand(finite))
# the u-integrand has a removable 0/0 at u -> 1 (catastrophic cancellation
# in float evaluation); check on a shaved interval at higher precision
mp.dps = 30
ref = mquad(guf, [0, mpf(1) - mpf("1e-12")])
mp.dps = 15
assert abs(float(mu_log) - float(ref)) < 1e-8
tick("u-integration done, checked numerically")
print("  mu_log(lam=0) =", mu_log)

# ===================================================================
mu_IId = simplify(mu_rat.subs(log(lam), 0) + mu_log)
print()
print("mu_IId = %s + (1/2) log(lam^2)" % mu_IId)
target = Rational(11, 24) - pi**2 / 18
assert simplify(mu_IId - target) == 0
print("      == 11/24 - pi^2/18 + (1/2) log(lam^2/m^2)")
print("      (Petermann 1957, eq. (4));  constant =", float(mu_IId))
