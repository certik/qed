"""Stage 0: full LO vertex derivation from the Feynman rules.

Validates the whole pipeline (explicit Dirac algebra + F1/F2 projection +
Feynman parametrization + loop integration) by re-deriving

    F_2(0) = alpha/(2 pi)              (Schwinger), and
    F_2(0; lambda) = (alpha/pi) K(lambda^2/m^2)   (massive photon kernel,
                                                   cf. g2_lo.py)

from

  ubar dGamma^mu u = -i e^2 int d^4k/(2pi)^4
      [ubar gamma^nu (kslash' + m) gamma^mu (kslash + m) gamma_nu u]
      / [((k-p)^2 - lam^2) (k'^2 - m^2) (k^2 - m^2)],   k' = k + q.

Checks performed:
  * F_1 diverges (has LUV), F_2 is UV finite (no LUV)  [Ward-consistent]
  * F_2(0) = alpha/(2 pi) exactly (lam = 0)
  * F_2(0; lam) matches K(t) numerically at t = lam^2 = 8
"""
from sympy import (symbols, expand, simplify, factor, limit, integrate, I,
                   pi, Rational, S, zeros, nsimplify, Symbol)
from mpmath import mp, quad

from dirac import GAMMA, ID4, METRIC, slash, breit_frame, extract_F1_F2
from loops import comps, dot, msub, feynman_shift, symmetrize, loop_integral, LUV

m, w, lam = symbols("m w lam", positive=True)
x, y, z = symbols("x y z", positive=True)
e2 = Symbol("e2", positive=True)   # e^2 = 4 pi alpha
l2 = Symbol("l2", positive=True)   # l^2 after angular average

p, pp, q = breit_frame(m, w)
k = comps("k")
kp = [k[mu] + q[mu] for mu in range(4)]

# ---- Feynman-parametrized denominator: x photon, y k'-line, z k-line ----
D = (x * (dot([k[mu] - p[mu] for mu in range(4)],
              [k[mu] - p[mu] for mu in range(4)]) - lam**2)
     + y * (dot(kp, kp) - m**2)
     + z * (dot(k, k) - m**2))
D = D.subs(x, 1 - y - z)
shift, Delta = feynman_shift(D, k)
Delta = simplify(Delta)
print("Delta =", factor(Delta.subs(w, 0)), "  (at q^2 = 0)")

# ---- numerator matrix N^mu, shifted and angular-averaged ----
lsh = [k[mu] + shift[mu] for mu in range(4)]  # k -> l + s, l renamed k


def vertex(mu):
    N = zeros(4, 4)
    for nu in range(4):
        N += (METRIC[nu] * GAMMA[nu] * (slash(kp) + m * ID4) * GAMMA[mu]
              * (slash(k) + m * ID4) * GAMMA[nu])
    N = N.subs(dict(zip(k, lsh)), simultaneous=True)
    N = N.applyfunc(lambda entry: symmetrize(entry, k, l2))
    # -i e^2 * (Feynman factor 2) * loop integrals over l
    J0 = loop_integral(0, 3, Delta)
    J1 = loop_integral(1, 3, Delta)
    return N.applyfunc(
        lambda entry: -I * e2 * 2
        * (entry.coeff(l2, 0) * J0 + entry.coeff(l2, 1) * J1))


F1, F2 = extract_F1_F2(vertex, m, w)

# ---- checks ----
assert F1.has(LUV), "F1 should be UV divergent"
assert not expand(F2).has(LUV), "F2 must be UV finite"
print("UV structure OK: LUV in F1 only")

F2 = simplify(expand(F2.subs(w, 0)))          # q^2 -> 0
F2 = F2.subs(m, 1)
print("F2 integrand (q^2=0, m=1):", factor(F2))

# lam = 0: integrate over the triangle y in (0, 1-z), z in (0,1)
F2_0 = F2.subs(lam, 0)
res = integrate(integrate(F2_0, (y, 0, 1 - z)), (z, 0, 1))
res = simplify(res * 16 * pi**2 / e2)  # in units alpha/pi = e2/(4 pi^2)
print("F2(0) =", res, "* e2/(16 pi^2)  =>  a_e = alpha/(2 pi):",
      res == 2)
assert res == 2

# massive photon: compare against K(t) at t = 8 numerically
mp.dps = 20
F2_lam = F2.subs(lam**2, 8)
val = integrate(integrate(F2_lam, (y, 0, 1 - z)), (z, 0, 1))
val = float(simplify(val * 16 * pi**2 / e2)) / 4   # (alpha/pi) K
K8 = quad(lambda zz: zz * (1 - zz) ** 2 / ((1 - zz) ** 2 + 8 * zz), [0, 1])
print("K(8) from trace pipeline =", val)
print("K(8) direct              =", K8)
assert abs(val - float(K8)) < 1e-15
print("OK: massive-photon kernel reproduced")
