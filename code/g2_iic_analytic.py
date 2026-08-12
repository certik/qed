"""Diagram IIc: analytic evaluation of the parametric integrals derived by
g2_iic.py, in the lambda -> 0 limit.  ** WORK IN PROGRESS **

Input: g2_iic_integrands.json. Numerically (g2_iic.f90) the pieces behave
as lambda -> 0 like

    mu_a -> 0.76976...   (finite)
    mu_b -> -log(lam) - 5/4        (proved exactly below)
    mu_c -> -0.08401...  (finite)

so mu_a + mu_c must equal (Petermann eq. (3) minus mu_b's constant)

    A* + C* = -37/24 + pi^2/18 - zeta(3)/2 + (pi^2/3) log 2
            = 0.685897...      (Fortran: 0.7697 - 0.0840 = 0.6857 OK)

Status:
  * piece (b): DONE exactly here: it factorizes into two 1-dim integrals
    over sig = u+v and s = y+z; mu_b = -log(lam) - 5/4 + O(lam).
  * piece (a): reduced (z and chi = u+v integrated, both verified
    numerically) to a 3-dim integral of gchi(s, t, u) over
    s in (0, 1-t), t, u in (0,1), where gchi is rational plus
    rational * log(Q_i) plus rational * atan(Mobius(s)) with
    Q = u(t(chi-1)+s(u-1))^2 + t(1-u)chi^2 the (irreducible,
    positive-definite) outer Delta at the chi-endpoints, and the atan
    scale R = sqrt(t u (1-u)) (rationalizable by t = tau^2,
    u = 1/(1+eta^2)). The remaining s, t, u integration is TODO.
  * piece (c): TODO (same program, with the spectral variable
    C in (1/u-ish, oo) from the xi-representation).
"""
import json
import pickle
import time

from sympy import (symbols, Symbol, S, pi, log, atan, expand, together,
                   simplify, factor, fraction, integrate, limit, series,
                   lambdify, sympify, Poly, Rational, factor_list,
                   expand_log)
from mpmath import mp, quad as mquad, mpf

from ratint import anti_rational

mp.dps = 15
T0 = time.time()


def tick(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


y, z, t, u, v, xi, lam = symbols("y z t u v xi lam", positive=True)
s, chi, sig = symbols("s chi sig", positive=True)
CIN = Symbol("CIN", positive=True)

with open("g2_iic_integrands.json") as fh:
    d = json.load(fh)


def defint01(expr, x):
    F = anti_rational(expr, x)
    return limit(F, x, 1, '-') - limit(F, x, 0, '+')


# ===================================================================
# piece (b): factorized cross term  -deltaF1(0)-rational-part x LO
# ===================================================================
fb = sympify(d["fb"])
num, den = fraction(together(fb))
# fb = 2 (y+z)(y+z-1) (sig^2+2sig-2) / [2 D0(sig) Delta(s)], sums only
form = (2 * (y + z) * (y + z - 1) * ((u + v) ** 2 + 2 * (u + v) - 2)
        / (2 * ((u + v) ** 2 + (1 - u - v) * lam**2)
           * ((y + z) ** 2 + (1 - y - z) * lam**2)))
assert simplify(together(form - fb)) == 0
tick("fb factorized form verified")

# simplex measures give one factor of the sum variable each
Iu = defint01(sig * (sig**2 + 2 * sig - 2) / (sig**2 + (1 - sig) * lam**2),
              sig)
Iy = defint01(s * s * (s - 1) / (s**2 + (1 - s) * lam**2), s)
mu_b = Iu * Iy
# numeric check against the Fortran b-column at lam = 0.001 (5.6424758)
val = mu_b.subs(lam, Rational(1, 1000)).evalf(12)
assert abs(float(val) - 5.6424758) < 1e-6
mu_b0 = simplify(expand(series(mu_b, lam, 0, 1).removeO()))
print("  mu_b(lam->0) =", mu_b0)
assert mu_b0 == -log(lam) - Rational(5, 4)
tick("piece (b) done: mu_b = -log(lam) - 5/4")

# ===================================================================
# piece (a): reduce to a 3-dim integral (z, chi integrated exactly)
# ===================================================================
fa = sympify(d["fa"]).subs(CIN, sympify(d["CIN_val"])).subs(lam, 0)
num, den = fraction(together(fa))
Q = [b for b, e in factor_list(den)[1] if e == 2][0]
Qs = expand(Q.subs(y, s - z).subs(v, chi - u))
# Q is positive definite: Q = u (t(chi-1)+s(u-1))^2 + t(1-u) chi^2
Qform = u * (t * (chi - 1) + s * (u - 1)) ** 2 + t * (1 - u) * chi**2
assert expand(Qs - Qform) == 0
tick("outer Delta structure verified (sum of squares)")

num_s = expand(num.subs(y, s - z))
numz = integrate(Poly(num_s, z).as_expr(), (z, 0, s))
numz = expand(numz.subs(v, chi - u))
ga = numz / (-2 * (u - 1) ** 15 * Qs**2)

# numeric check of the z-integration + variable change
gaf = lambdify((s, chi, t, u), ga, "mpmath")
faf = lambdify((y, z, t, u, v), fa, "mpmath")
sv, tv, uv, vv = mpf("0.3"), mpf("0.2"), mpf("0.4"), mpf("0.25")
direct = mquad(lambda zz: faf(sv - zz, zz, tv, uv, vv), [0, sv])
assert abs(complex(direct).real - complex(gaf(sv, uv + vv, tv, uv)).real) \
    < 1e-8
tick("z-integration verified")

Fchi = anti_rational(ga, chi)
gchi = expand_log(Fchi.subs(chi, 1) - Fchi.subs(chi, u), force=True)
gchif = lambdify((s, t, u), gchi, "mpmath")
direct2 = mquad(lambda cc: gaf(sv, cc, tv, uv), [uv, 1])
assert abs(complex(direct2).real - complex(gchif(sv, tv, uv)).real) < 1e-10
with open("g2_iic_a_gchi.pkl", "wb") as fh:
    pickle.dump(gchi, fh)
tick("chi-integration verified; gchi cached (g2_iic_a_gchi.pkl)")

print("""
TODO: s in (0,1-t), then t, u (rationalize R = sqrt(t u(1-u)) via
t = tau^2, u = 1/(1+eta^2)); then piece (c).""")
