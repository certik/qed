"""Diagram IId, part 1: the self-energy subgraph and its counterterms.

Everything the worked chapter states about Sigma, delta_m and delta_Z2 is
computed here, in closed form in the photon mass lambda, and cross-checked:

  1. the gamma-matrix contraction identities used in the covariant
     reduction, verified against the explicit 4x4 matrices of dirac.py;
  2. Sigma_loop(k) = A(k^2) + B(k^2) kslash with Pauli-Villars;
  3. the on-shell counterterms delta_m and delta_Z2 in closed form,
     exhibiting the infrared logarithm of delta_Z2;
  4. the same in dimensional regularization, showing that the two schemes
     differ by a term linear in kslash with k-independent coefficients,
     which the on-shell subtraction annihilates identically: Sigma_R is
     scheme independent while the counterterms are not;
  5. the prediction that the whole infrared logarithm of mu_IId is the
     delta_Z2 subtraction, -2 delta_Z2 F2_LO = +(1/2) log(lam^2/m^2).

Conventions: metric (+,-,-,-); fermion propagator i(kslash+m)/(k^2-m^2);
photon propagator -i g_{mu nu}/(k^2-lam^2); vertex -i e gamma^mu;
c = e^2/(16 pi^2) = alpha/(4 pi); LUV = log(Lambda^2) + scheme constant.
"""
from sympy import (symbols, Symbol, S, Rational, expand, simplify, cancel,
                   factor, log, sqrt, pi, integrate, limit, series, O,
                   zeros, eye, together, apart, nsimplify)
from mpmath import mp, quad, mpf

from dirac import GAMMA, ID4, METRIC, slash
from ratint import anti_rational

u, lam, k2, eps = symbols("u lam k2 epsilon", positive=True)
m = Symbol("m", positive=True)
c = Symbol("c", positive=True)          # c = e^2/(16 pi^2) = alpha/(4 pi)
LUV = Symbol("LUV")
kk = symbols("k0 k1 k2c k3")            # generic momentum components

SEP = "=" * 70


def head(txt):
    print("\n" + SEP + "\n" + txt + "\n" + SEP)


# ----------------------------------------------------------------- 1.
head("1. gamma contraction identities (checked on explicit 4x4 matrices)")


def contract(build):
    """sum_nu g_{nu nu} gamma^nu (...) gamma^nu, with (...) = build(nu)."""
    tot = zeros(4, 4)
    for nu in range(4):
        tot += METRIC[nu] * GAMMA[nu] * build() * GAMMA[nu]
    return tot


for alpha_i in range(4):
    lhs = contract(lambda: GAMMA[alpha_i])
    rhs = -2 * GAMMA[alpha_i]
    assert (lhs - rhs) == zeros(4, 4)
print("gamma^nu gamma^a gamma_nu = -2 gamma^a                    OK")

for alpha_i in range(4):
    for beta_i in range(4):
        lhs = contract(lambda: GAMMA[alpha_i] * GAMMA[beta_i])
        rhs = 4 * (METRIC[alpha_i] if alpha_i == beta_i else 0) * ID4
        assert (lhs - rhs) == zeros(4, 4)
print("gamma^nu gamma^a gamma^b gamma_nu = 4 g^{ab}              OK")

for alpha_i in range(4):
    for beta_i in range(4):
        for gam_i in range(4):
            lhs = contract(lambda: GAMMA[alpha_i] * GAMMA[beta_i]
                           * GAMMA[gam_i])
            rhs = -2 * GAMMA[gam_i] * GAMMA[beta_i] * GAMMA[alpha_i]
            assert (lhs - rhs) == zeros(4, 4)
print("gamma^nu gamma^a gamma^b gamma^c gamma_nu = -2 g^c g^b g^a  OK")

# ----------------------------------------------------------------- 2.
head("2. the self-energy Sigma_loop(k) = A + B kslash (Pauli-Villars)")

# Delta_in = a - b k^2 from the Feynman shift (derived in the chapter)
a_in = (1 - u) * m**2 + u * lam**2
b_in = u * (1 - u)
D = a_in - b_in * k2
print("a =", a_in)
print("b =", factor(b_in))
print("D_in(k^2) = a - b k^2")

A = 4 * m * c * (LUV - log(D))
B = -2 * u * c * (LUV - log(D))
print("A(k^2) =", A)
print("B(k^2) =", B)

# ----------------------------------------------------------------- 3.
head("3. on-shell counterterms")

D0 = simplify(D.subs(k2, m**2))
print("D_0 = D_in(m^2) =", factor(D0))
assert simplify(D0 - ((1 - u)**2 * m**2 + u * lam**2)) == 0

dm_u = simplify((A + m * B).subs(k2, m**2))
dZ2_u = simplify((2 * m * A.diff(k2) + B + 2 * m**2 * B.diff(k2))
                 .subs(k2, m**2))
print("delta_m integrand   =", simplify(dm_u / c))
print("delta_Z2 integrand  =", simplify(dZ2_u / c))

# integrate over u in (0,1), at m = 1, exactly in lam
dm_i = simplify((dm_u / c).subs(m, 1))
dZ2_i = simplify((dZ2_u / c).subs(m, 1))
D0_1 = D0.subs(m, 1)

# WARNING (verified below): sympy's integrate() of P(u) log D0 over (0,1)
# returns a WRONG antiderivative here -- for P = 2u it claims
# lam^2 + 2 log(lam) - 3, which is -7.595 at lam = 0.1 while the true
# value is -2.477.  We integrate by parts instead, which leaves a rational
# integral handled deterministically by ratint.anti_rational, and check
# every result numerically.


def int01_Plog(P):
    """int_0^1 P(u) log D0 du, by parts; D0 = (1-u)^2 + u lam^2, D0(0)=1."""
    Q = integrate(P, u)                       # Q' = P, Q(0) = 0
    Q = Q - Q.subs(u, 0)
    boundary = Q.subs(u, 1) * log(lam**2)     # D0(1) = lam^2
    R = anti_rational(cancel(Q * D0_1.diff(u) / D0_1), u, verify="numeric")
    inner = limit(R, u, 1, "-") - limit(R, u, 0, "+")
    val = simplify(boundary - inner)
    # numeric verification at two values of lam
    for Ls in ["0.1", "0.01"]:
        l0 = mpf(Ls)
        Pf = lambda uu: float(P.subs(u, uu))
        num = quad(lambda uu: Pf(uu) * mp.log((1 - uu)**2 + uu * l0**2),
                   [0, 1])
        sym = mpf(str(val.subs(lam, Rational(Ls)).evalf(20)))
        assert abs(num - sym) < mpf("1e-15"), (P, Ls, num, sym)
    return val


mp.dps = 25
print("\n(all log-integrals below verified numerically at lam = 0.1, 0.01)")

# delta_m / c = int (4-2u)(LUV - log D0) du = 3 LUV - int (4-2u) log D0 du
dm_log = int01_Plog(4 - 2 * u)
dm_val = simplify(3 * LUV - dm_log)
print("\nint_0^1 (4-2u) log D0 du =", dm_log)
print("  limit lam->0:", limit(dm_log, lam, 0, "+"))
print("delta_m / c =", dm_val)
print("delta_m / c at lam -> 0 =", 3 * LUV - limit(dm_log, lam, 0, "+"))

# delta_Z2 / c = I1 - LUV + I2,
#   I1 = int 4u(1-u)(2-u)/D0 du     (rational, IR divergent)
#   I2 = int 2u log D0 du           (finite)
I1_anti = anti_rational(4 * u * (1 - u) * (2 - u) / D0_1, u,
                        verify="numeric")
I1 = limit(I1_anti, u, 1, "-") - limit(I1_anti, u, 0, "+")
I2 = int01_Plog(2 * u)
dZ2_val = I1 - LUV + I2
print("\nI1 = int 4u(1-u)(2-u)/D0 du,  lim (I1 + 4 log lam) =",
      limit(I1 + 4 * log(lam), lam, 0, "+"))
print("I2 = int 2u log D0 du,        lim I2 =", limit(I2, lam, 0, "+"))
print("delta_Z2 / c = I1 - LUV + I2")
print("  at lam -> 0:  -LUV - 4 log(lam) +",
      limit(I1 + 4 * log(lam), lam, 0, "+") + limit(I2, lam, 0, "+"))

for Ls in ["0.1", "0.01"]:
    l0 = mpf(Ls)
    n1 = quad(lambda uu: 4 * uu * (1 - uu) * (2 - uu)
              / ((1 - uu)**2 + uu * l0**2), [0, 1])
    s1 = mpf(str(I1.subs(lam, Rational(Ls)).evalf(20)))
    assert abs(n1 - s1) < mpf("1e-15"), (Ls, n1, s1)
print("  I1 verified numerically at lam = 0.1, 0.01")

# ----------------------------------------------------------------- 4.
head("4. dimensional regularization: same Sigma_R, different counterterms")

# d = 4 - 2 eps:  gamma^nu (lslash + m) gamma_nu = (2-d) lslash + d m
kslash = Symbol("ks")            # stands for kslash (commutes here: A + B ks)
num_d = (2 - (4 - 2 * eps)) * u * kslash + (4 - 2 * eps) * m
print("d-dim numerator after the shift:", expand(num_d))
ehat = Symbol("ehat")            # 1/eps - gamma_E + log(4 pi)
Sig_d = expand(c * num_d * (ehat - log(D)))
Sig_d = expand(Sig_d.subs(eps * ehat, 1))     # eps * (1/eps) -> 1 as eps->0
Sig_d = expand(Sig_d.subs(eps, 0))
print("Sigma_dimreg / c =", Sig_d)

Sig_pv = expand(c * (4 * m - 2 * u * kslash) * (ehat - log(D)))
diff = expand(Sig_d - Sig_pv)
print("Sigma_dimreg - Sigma_PV (with LUV <-> ehat) =", diff, "* c")

# the difference is E = alpha + beta * kslash with k-independent alpha,beta;
# the on-shell subtraction annihilates it identically
E = diff
E_m = E.subs(kslash, m)
E_p = E.diff(kslash)
assert expand(E - E_m - (kslash - m) * E_p) == 0
print("E - E(m) - (kslash - m) E'(m) =",
      expand(E - E_m - (kslash - m) * E_p), " => Sigma_R identical")

# ----------------------------------------------------------------- 5.
head("5. the infrared logarithm of mu_IId is exactly the delta_Z2 term")

# delta_Z2 = c * (dZ2_val); c = alpha/(4 pi) = (alpha/pi)/4.
# LO:  F2 = (alpha/pi) * 1/2.  The delta_Z2 insertion collapses to
# delta_Z2 * (LO diagram), and IId doubles it (mirror), so IId receives
#     -2 * delta_Z2 * F2_LO
#   = -2 * [(alpha/pi) dZ2_val / 4] * [(alpha/pi)/2]
#   = -(alpha/pi)^2 * dZ2_val / 4,
# i.e. a contribution -dZ2_val/4 to mu_IId.
dZ2_lim = -LUV - 4 * log(lam) + (limit(I1 + 4 * log(lam), lam, 0, "+")
                                 + limit(I2, lam, 0, "+"))
print("delta_Z2 / c ->", dZ2_lim, "  as lam -> 0")
ir = expand(-dZ2_lim / 4)
print("contribution to mu_IId from the delta_Z2 subtraction:", ir)
ir_log = ir.coeff(log(lam))
print("  coefficient of log(lam) =", ir_log,
      "  i.e. +(1/2) log(lam^2/m^2)")
assert ir_log == 1, ir_log
print("\nThis is exactly Petermann's infrared term for IId.  The rest of")
print("the diagram, F2[IId; Sigma - delta_m], therefore has to be infrared")
print("finite (and carries the compensating -LUV/4).  Sigma_R as a whole")
print("is LUV-free, which is what makes mu_IId well defined.")
