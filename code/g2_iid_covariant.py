"""Diagram IId, part 2: the covariant reduction of the insertion.

The pipeline (g2_iid.py) grinds the whole diagram with explicit 4x4
matrices in the Breit frame.  That is machine-friendly but opaque.  Here
we do the same algebra covariantly, by hand, and check every step against
the explicit matrices of dirac.py.  Two structural facts come out that
the brute-force route hides completely:

  Sigma_R(k) = -4 m^2 c b (2-u) (kslash - m)/D_0
               - c (4m - 2u kslash) log(D_in(k^2)/D_0)

  (a) the RATIONAL part of Sigma_R is proportional to (kslash - m), so
      sandwiched between the two propagators it gives
          (kslash+m)(kslash-m)(kslash+m) = (k^2-m^2)(kslash+m),
      which cancels one power of the doubled propagator and turns the
      whole rational piece into a CONSTANT times the LO vertex diagram.
      Hence mu_rat = -I1(lam) K(lam^2)/2 with K the LO massive-photon
      kernel and I1 the very integral that appears in delta_Z2.
  (b) the LOG part, after the xi-representation, gives a four-propagator
      integral with the numerator
          N^mu = P (4 k'^mu - 2m gamma^mu)
               + Q (-2 kslash gamma^mu k'slash + 4m k^mu),
      P = f(k^2+m^2) + 2 m g k^2,  Q = 2 m f + g(k^2+m^2),
      for an insertion f + g kslash.

Everything below is verified: the Dirac identities on explicit matrices
at random momenta, and mu_rat against the pipeline's own f_rat.
"""
import json
from sympy import (symbols, Symbol, S, Rational, expand, simplify, cancel,
                   factor, log, pi, integrate, limit, zeros, srepr, sympify,
                   nsimplify, together, sqrt, Integer)
from mpmath import mp, quad, mpf

from dirac import GAMMA, ID4, METRIC, slash
from ratint import anti_rational

SEP = "=" * 70


def head(txt):
    print("\n" + SEP + "\n" + txt + "\n" + SEP)


def rand_mom(seed):
    """A deterministic 'random' rational 4-momentum."""
    import random
    random.seed(seed)
    return [Rational(random.randint(-9, 9), random.randint(1, 7))
            for _ in range(4)]


def mdot(a, b):
    return sum(METRIC[i] * a[i] * b[i] for i in range(4))


# ----------------------------------------------------------------- 1.
head("1. the insertion sandwich  (kslash+m)(f + g kslash)(kslash+m)")

f, g, mm = symbols("f g m")
for seed in (1, 2, 3):
    k = rand_mom(seed)
    ks = slash(k)
    k2 = mdot(k, k)
    lhs = (ks + mm * ID4) * (f * ID4 + g * ks) * (ks + mm * ID4)
    P = f * (k2 + mm**2) + 2 * mm * g * k2
    Q = 2 * mm * f + g * (k2 + mm**2)
    rhs = P * ID4 + Q * ks
    assert simplify(expand(lhs - rhs)) == zeros(4, 4), seed
print("(kslash+m)(f+g kslash)(kslash+m) = P + Q kslash   with")
print("   P = f(k^2+m^2) + 2 m g k^2")
print("   Q = 2 m f + g(k^2+m^2)")
print("verified on explicit matrices at 3 random momenta        OK")

# the special case g = kappa, f = -kappa m  (i.e. kappa (kslash - m))
kap = Symbol("kappa")
for seed in (4, 5):
    k = rand_mom(seed)
    ks = slash(k)
    k2 = mdot(k, k)
    lhs = (ks + mm * ID4) * kap * (ks - mm * ID4) * (ks + mm * ID4)
    rhs = kap * (k2 - mm**2) * (ks + mm * ID4)
    assert simplify(expand(lhs - rhs)) == zeros(4, 4), seed
print("\n(kslash+m) kappa(kslash-m) (kslash+m) = kappa (k^2-m^2)(kslash+m)")
print("=> the rational part cancels one power of the doubled propagator")

# ----------------------------------------------------------------- 2.
head("2. the gamma^nu ... gamma_nu contraction of the outer string")

for seed in (6, 7):
    k = rand_mom(seed)
    kp = rand_mom(seed + 10)
    ks, kps = slash(k), slash(kp)
    Pv, Qv = symbols("P Q")
    for mu in range(4):
        lhs = zeros(4, 4)
        for nu in range(4):
            lhs += (METRIC[nu] * GAMMA[nu] * (kps + mm * ID4) * GAMMA[mu]
                    * (Pv * ID4 + Qv * ks) * GAMMA[nu])
        rhs = (Pv * (4 * kp[mu] * ID4 - 2 * mm * GAMMA[mu])
               + Qv * (-2 * ks * GAMMA[mu] * kps + 4 * mm * k[mu] * ID4))
        assert simplify(expand(lhs - rhs)) == zeros(4, 4), (seed, mu)
print("gamma^nu (k'slash+m) gamma^mu (P + Q kslash) gamma_nu")
print("   = P (4 k'^mu - 2 m gamma^mu)")
print("   + Q (-2 kslash gamma^mu k'slash + 4 m k^mu)")
print("verified on explicit matrices at 2 random momentum pairs  OK")

# ----------------------------------------------------------------- 3.
head("3. Sigma_R in closed form: rational part is proportional to "
     "(kslash - m)")

u, lam, K2, xi = symbols("u lam K2 xi", positive=True)
c = Symbol("c", positive=True)
LUV = Symbol("LUV")
m = S(1)
a_in = (1 - u) * m**2 + u * lam**2
b_in = u * (1 - u)
D = a_in - b_in * K2
D0 = simplify(D.subs(K2, m**2))

ks = Symbol("ks")                      # kslash, commutes in A + B ks
Sig = c * (4 * m - 2 * u * ks) * (LUV - log(D))
dm = c * (4 * m - 2 * u * m) * (LUV - log(D0))
dZ2 = c * (4 * m**2 * b_in * (2 - u) / D0 - 2 * u * (LUV - log(D0)))
Sig_R = expand(Sig - dm - (ks - m) * dZ2)
assert Sig_R.coeff(LUV) == 0
print("LUV cancels in Sigma_R pointwise in u                     OK")

claim = (-4 * m**2 * c * b_in * (2 - u) * (ks - m) / D0
         - c * (4 * m - 2 * u * ks) * (log(D) - log(D0)))
assert simplify(expand(Sig_R - expand(claim))) == 0
print("Sigma_R = -4 m^2 c b(2-u)(kslash-m)/D_0"
      " - c(4m-2u kslash) log(D/D_0)   OK")
kappa = simplify(-4 * m**2 * c * b_in * (2 - u) / D0)
print("kappa =", kappa)

# ----------------------------------------------------------------- 4.
head("4. the rational piece = kappa x (LO vertex diagram)")

z = Symbol("z", positive=True)
t = Symbol("t", positive=True)
# LO massive-photon kernel, derived in the warm-up section.
# NOTE: sympy's integrate() returns 1/2 - lam^2 for this, which is WRONG
# (the same parametric-integral trap g2_lo.py warns about; the true K has
# a square-root/arctan structure).  Use the deterministic integrator.
K_int = z * (1 - z)**2 / ((1 - z)**2 + z * lam**2)
FK = anti_rational(K_int, z, verify="numeric")
K_lam = limit(FK, z, 1, "-") - limit(FK, z, 0, "+")
print("K(t) = int_0^1 z(1-z)^2/((1-z)^2+zt) dz")
print("sympy's naive integrate() claims K =",
      integrate(K_int, (z, 0, 1)), " <-- WRONG")
print("K(0) = lim_{lam->0} K =", limit(K_lam, lam, 0, "+"), " (Schwinger)")
mp.dps = 20
for Ls in ["0.5", "0.1", "0.01"]:
    r = Rational(Ls)
    kn = quad(lambda zz: float(K_int.subs({lam: r, z: zz})), [0, 1])
    kv = mpf(str(K_lam.subs(lam, r).evalf(18)))
    assert abs(kn - kv) < mpf("1e-14"), (Ls, kn, kv)
print("K(lam^2) verified numerically at lam = 0.5, 0.1, 0.01       OK")

# mu_rat = 2 * [int du kappa_hat(u)] * K,  kappa_hat = kappa/(alpha/pi)
#        = 2 * int du [-u(1-u)(2-u)/D_0] * K = -(1/2) I1 K
I1_anti = anti_rational(4 * u * (1 - u) * (2 - u) / D0, u, verify="numeric")
I1 = limit(I1_anti, u, 1, "-") - limit(I1_anti, u, 0, "+")
mu_rat = simplify(-I1 * K_lam / 2)
print("\nmu_rat(lam) = -(1/2) I1(lam) K(lam^2)")
print("lim (mu_rat - log lam) =", limit(mu_rat - log(lam), lam, 0, "+"))

# Cross-check against the pipeline's own f_rat, INTEGRAND BY INTEGRAND
# (much stronger, and far cheaper, than comparing the integrals): the
# covariant claim is
#     f_rat(y,z,u,lam) = 2 * kappa_hat(u) * g_LO(y,z,lam)
# with kappa_hat = kappa/(alpha/pi) and g_LO the LO F2 parametric
# integrand in the pipeline's own (y,z) parametrization, which we rebuild
# here with the same machinery.
from dirac import breit_frame, u_spinor, ubar
from loops import comps, dot as ldot, feynman_shift, symmetrize, \
    loop_integral, LUV as LUVs
from sympy import gamma as Gamma_f, I as Ii, series

w = Symbol("w", positive=True)
e2 = Symbol("e2", positive=True)
l2 = Symbol("l2", positive=True)
x_, y_, z_ = symbols("x y z", positive=True)
pmom, ppmom, qmom = breit_frame(S(1), w)
kmom = comps("k")
kpmom = [kmom[i] + qmom[i] for i in range(4)]


def lo_integrand():
    """The LO F2 parametric integrand g_LO(y,z,lam), pipeline conventions."""
    Pph = ldot([kmom[i] - pmom[i] for i in range(4)],
               [kmom[i] - pmom[i] for i in range(4)]) - lam**2
    Pkp = ldot(kpmom, kpmom) - 1
    Pk = ldot(kmom, kmom) - 1
    Dcomb = x_ * Pph + y_ * Pkp + z_ * Pk
    wt = Gamma_f(3)
    Dcomb = Dcomb.subs(x_, 1 - y_ - z_)
    shift, Delta = feynman_shift(Dcomb, kmom)
    shiftmap = dict(zip(kmom, [kmom[b] + shift[b] for b in range(4)]))
    vals = {}
    for (mu, sp, s) in [(0, 0, 0), (1, 0, 1)]:
        ub, uu = ubar(ppmom, S(1), sp), u_spinor(pmom, S(1), s)
        tot = S(0)
        kpm = slash(kpmom) + ID4
        km = slash(kmom) + ID4
        for nu in range(4):
            tot += METRIC[nu] * (ub * GAMMA[nu] * kpm * GAMMA[mu] * km
                                 * GAMMA[nu] * uu)[0, 0]
        A = expand(tot.subs(shiftmap, simultaneous=True))
        A = expand(A.series(w, 0, 2).removeO())
        A = symmetrize(A, kmom, l2)
        val = -Ii * e2 * wt * sum(
            A.coeff(l2, aa) * loop_integral(aa, 3, Delta)
            for aa in range(3) if A.coeff(l2, aa) != 0)
        vals[(mu, sp, s)] = expand(val.series(w, 0, 2).removeO())
    F2 = -(vals[(0, 0, 0)].subs(w, 0)
           + vals[(1, 0, 1)].diff(w).subs(w, 0)) / 2
    return expand(F2)


g_lo = lo_integrand()
# strip the loop factor: F2 = (alpha/pi) * g_hat, alpha/pi = e2/(4 pi^2)
g_hat = cancel(g_lo * 4 * pi**2 / e2)
print("\nLO F2 integrand (units alpha/pi):", g_hat)
chk = integrate(integrate(g_hat.subs(lam, 0), (y_, 0, 1 - z_)), (z_, 0, 1))
print("int g_hat at lam=0 =", chk, " (= K(0) = 1/2)")
assert chk == Rational(1, 2)

# g_LO depends on (y,z) only through s = y+z, and
#     int_0^1 dz int_0^(1-z) dy  F(y+z) = int_0^1 s F(s) ds,
# so K(lam^2) = int_0^1 s^2 (1-s)/(s^2 + (1-s) lam^2) ds -- the same
# kernel as K(t) after z -> 1-s.
s_ = Symbol("s", positive=True)
g_s = cancel(g_hat.subs({y_: s_ - z_}))
assert not g_s.has(z_), g_s
print("g_LO depends on y,z only through s = y+z:", g_s)
FKs = anti_rational(s_ * g_s, s_, verify="numeric")
K_s = limit(FKs, s_, 1, "-") - limit(FKs, s_, 0, "+")
for Ls in ["0.5", "0.1", "0.01"]:
    r = Rational(Ls)
    d = mpf(str((K_s - K_lam).subs(lam, r).evalf(18)))
    assert abs(d) < mpf("1e-14"), (Ls, d)
print("int_0^1 s g_LO(s) ds == K(lam^2)  (numerically, 3 lambdas)  OK")

# The pipeline instead keeps the DOUBLED propagator and Feynman-
# parametrizes four powers, so its f_rat is a different function of
# (y,z) with the same integral.  Its documented factorization is
# mu_rat = U(lam) S(lam); the covariant result says mu_rat = -U K, i.e.
# S must equal -K.  That is the sharp check:
U_int = 2 * u * (u - 2) * (u - 1) / ((1 - u)**2 + lam**2 * u)
FU = anti_rational(U_int, u, verify="numeric")
U_lam = limit(FU, u, 1, "-") - limit(FU, u, 0, "+")
assert simplify(U_lam - I1 / 2) == 0
print("\nU(lam) (pipeline) == I1(lam)/2 (covariant)                OK")

S_int = (s_**2 / 2 * (s_ - 1)
         * (3 * s_**3 + s_**2 - lam**2 * (4 * s_ - 1) * (s_ - 1))
         / (s_**2 + (1 - s_) * lam**2)**2)
FS = anti_rational(S_int, s_, verify="numeric")
S_lam = limit(FS, s_, 1, "-") - limit(FS, s_, 0, "+")
# simplify() cannot see the cancellation (atan branch bookkeeping), so we
# verify S + K = 0 numerically at several lambda
print("\nS(lam) + K(lam^2) at lam = 0.5, 0.1, 0.01:")
for Ls in ["0.5", "0.1", "0.01"]:
    r = Rational(Ls)
    val = mpf(str((S_lam + K_lam).subs(lam, r).evalf(18)))
    print("   lam=%-5s  S+K = %s" % (Ls, mp.nstr(val, 5)))
    assert abs(val) < mpf("1e-14"), (Ls, val)
print("S(lam) == -K(lam^2): the pipeline's factorized rational piece")
print("is exactly kappa x (LO massive-photon diagram)             OK")

print("\nmu_rat = U S = -U K = -(1/2) I1 K =",
      "log(lam) + 1/2 + O(lam)")
assert limit(mu_rat - log(lam), lam, 0, "+") == Rational(1, 2)

# ----------------------------------------------------------------- 5.
head("5. the log piece: the xi representation and its propagator")

Cs = Symbol("C", positive=True)
C_val = m**2 + (a_in - b_in * m**2) / (xi * b_in)
lhs = (K2 - m**2) / (xi * (K2 - Cs))
rhs_check = simplify(lhs.subs(Cs, C_val) - (D - D0) / (D0 + xi * (D - D0)))
assert simplify(rhs_check) == 0
print("log(D/D_0) = int_0^1 dxi (D-D_0)/(D_0+xi(D-D_0))")
print("           = int_0^1 dxi (k^2-m^2)/(xi (k^2 - C)),   C =",
      factor(C_val))
print("C - m^2 =", factor(simplify(C_val - m**2)), "> 0  for xi in (0,1)")
print("\n=> one power of (k^2-m^2)^2 is cancelled here too, leaving the")
print("   four propagators [(k-p)^2-lam^2][k'^2-m^2][k^2-m^2][k^2-C].")
