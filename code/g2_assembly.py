"""Assembly of the order-alpha^2 anomalous moment from the per-diagram
values established in this chapter (units of (alpha/pi)^2; lam = photon
mass, L = log(lam^2/m^2) is the IR regulator log).

Each mu below was computed diagram by diagram in code/g2_*.py / *.f90 and
agrees with Petermann, Helv. Phys. Acta 30 (1957) 407, eqs. (1)-(6).
This script proves that the IR logs cancel and the sum is exactly

    A2 = 197/144 + pi^2/12 + (3/4) zeta(3) - (1/2) pi^2 log 2.
"""
from sympy import (Rational, pi, log, zeta, Symbol, simplify, nsimplify,
                   N)

L = Symbol('L')          # log(lam^2/m^2)

mu_I = Rational(1, 6) + Rational(13, 36)*pi**2 + Rational(5, 4)*zeta(3) \
    - Rational(5, 6)*pi**2*log(2)
mu_IIa = Rational(11, 48) + pi**2/18
mu_IIc = -Rational(67, 24) + pi**2/18 - zeta(3)/2 + pi**2*log(2)/3 \
    - L/2
mu_IId = Rational(11, 24) - pi**2/18 + L/2
mu_IIe = Rational(119, 36) - pi**2/3

total = simplify(mu_I + mu_IIa + mu_IIc + mu_IId + mu_IIe)
A2 = Rational(197, 144) + pi**2/12 + Rational(3, 4)*zeta(3) \
    - pi**2*log(2)/2

print("IR log cancellation: coeff of L in sum =", total.coeff(L))
print("sum  =", total)
print("A2   =", A2)
print("sum - A2 =", simplify(total - A2))
assert simplify(total - A2) == 0
print()
print("numeric values (units (alpha/pi)^2):")
for name, mu in [("I", mu_I), ("IIa", mu_IIa),
                 ("IIc", mu_IIc.subs(L, 0)), ("IId", mu_IId.subs(L, 0)),
                 ("IIe", mu_IIe)]:
    print("  mu_%-4s = %s" % (name, N(mu, 15)))
print("  A2      = %s" % N(A2, 15))
