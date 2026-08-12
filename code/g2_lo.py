"""LO anomalous moment: F_2(0) with a photon of mass lambda (mass m = 1).

F_2(0; lambda) = (alpha/pi) * K(t),   t = lambda^2/m^2,
K(t) = int_0^1 z(1-z)^2 / ((1-z)^2 + z t) dz,   K(0) = 1/2  (Schwinger).

This validates the parametric-integral pipeline and provides the massive-
photon kernel used by the vacuum-polarization diagram IIe (which needs it
for spectral masses t >= 4 only).

Note: asking SymPy for K(t) with t held symbolic returns a wrong result
(branch issues in the parametric integral), so for t > 4 we rationalize
the denominator first:

    (1-z)^2 + z t = (z + y)(z y + 1)/y   with   t = (1+y)^2/y,  0 < y <= 1.
"""
from sympy import symbols, integrate, Rational, simplify, nsimplify, log, S
from mpmath import mp, quad

z, y = symbols("z y", positive=True)

# --- Schwinger value: t = 0 => integrand reduces to z ---
K0 = integrate(z, (z, 0, 1))
print("K(0) =", K0, " => a_e = alpha/(2 pi)")
assert K0 == Rational(1, 2)

# --- massive photon, t = (1+y)^2/y > 4 ---
Ky = integrate(z * (1 - z) ** 2 * y / ((z + y) * (z * y + 1)), (z, 0, 1))
Ky = simplify(Ky)
print("K(y) =", Ky)

# numeric cross-check at t = 8  (y = 3 - 2*sqrt(2))
mp.dps = 25
t_num = 8
y_num = mp.mpf(3) - 2 * mp.sqrt(2)
direct = quad(lambda zz: zz * (1 - zz) ** 2 / ((1 - zz) ** 2 + zz * t_num), [0, 1])
closed = mp.mpf(str(Ky.subs(y, nsimplify(y_num, [S(2) ** Rational(1, 2)])).evalf(25)))
print("K(t=8) direct  =", direct)
print("K(t=8) closed  =", closed)
assert abs(direct - closed) < mp.mpf("1e-20")
print("OK: closed form matches numeric integration")
