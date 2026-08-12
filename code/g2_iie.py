"""Diagram IIe: vacuum-polarization insertion into the LO vertex photon.

Method (dispersion relation): the renormalized photon self-energy gives
the spectral representation of the dressed photon propagator,

    1/k^2  ->  1/k^2 + int_4^oo dt rho(t) / (k^2 - t),      (m = 1)
    rho(t) = (alpha/3 pi) (1/t) (1 + 2/t) sqrt(1 - 4/t),

so IIe is the LO diagram computed with a massive photon (the kernel K(t)
of g2_lo.py), folded with rho(t):

    mu_IIe = int_4^oo dt (1/3t)(1+2/t) sqrt(1-4/t) K(t)    [units (alpha/pi)^2]

Everything rationalizes with t = (1+y)^2/y, y in (0,1]:

    sqrt(t(t-4)) = (1-y^2)/y,     |dt/dy| = (1-y^2)/y^2,
    rho(t) dt -> (y^2+4y+1)(1-y)^2 / (3 y (1+y)^4) dy,
    (1-z)^2 + z t = (z+y)(z y+1)/y.

The double integral over (z, y) is done with the y (rational) integration
first; the remaining z integral has integrands of the form
rational * log(z), giving dilogarithms at 1 (i.e. pi^2).

SymPy notes: definite integrate() with a free parameter in the integrand
is slow/fragile here, so we take indefinite antiderivatives and evaluate
the limits at the endpoints ourselves.

Target (Petermann, Helv. Phys. Acta 30 (1957) 407, eq. (5)):
    mu_IIe = 119/36 - pi^2/3 = 0.0156874...
"""
from sympy import (symbols, integrate, apart, together, expand, expand_log,
                   simplify, limit, log, pi, Rational)

z, y = symbols("z y", positive=True)

# rho(t) dt in the y variable, and the massive-photon kernel integrand
weight = (y**2 + 4 * y + 1) * (1 - y) ** 2 / (3 * y * (y + 1) ** 4)
Kz = z * (1 - z) ** 2 * y / ((z + y) * (z * y + 1))

# ---- inner integral over y (rational function of y) ----
F = integrate(apart(together(weight * Kz), y), y)
Iy = simplify(F.subs(y, 1) - limit(F, y, 0, "+"))
Iy = expand(expand_log(Iy, force=True))
print("after y-integration:")
print("  I(z) =", simplify(Iy))

# ---- outer integral over z:  I(z) = A(z) log(z) + B(z), A, B rational ----
A = simplify(Iy.coeff(log(z)))
B = simplify(Iy - A * log(z))
G = integrate(expand(apart(A, z) * log(z)), z) + integrate(apart(B, z), z)
mu_IIe = expand(simplify(limit(G, z, 1, "-") - limit(G, z, 0, "+")))

print("mu_IIe =", mu_IIe)
print("       =", mu_IIe.evalf(20))
assert simplify(mu_IIe - (Rational(119, 36) - pi**2 / 3)) == 0
print("OK: mu_IIe = 119/36 - pi^2/3  (Petermann 1957, eq. (5))")
