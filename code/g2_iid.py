"""Diagram IId: one-loop self-energy inserted on an internal electron line
of the LO vertex (plus mirror), renormalized on shell.

Steps (all mechanical, using dirac.py + loops.py; m = 1):

1. Inner loop (derived, not assumed):
   Sigma(k) = e^2/(16 pi^2) (4m - 2u kslash)(LUV - log D_in),
   D_in = a - b k^2,  a = (1-u) m^2 + u lam^2,  b = u(1-u),  u in (0,1).

2. On-shell subtraction Sigma_R = Sigma - dm - (kslash - m) dZ2 with
   Sigma = A(k^2) + B(k^2) kslash,  dm = [A + m B](m^2),
   dZ2 = [2m A' + B + 2 m^2 B'](m^2). LUV cancels exactly (asserted):

   Sigma_R = M_rat + M_log log(D_in(k^2)/D0),  D0 = D_in(m^2).

3. Insertion: i S(k) -> i (kslash+m) Sigma_R (kslash+m) / (k^2-m^2)^2 in
   the LO integrand. The log is made rational by
   log(X/Y) = int_0^1 dxi (X-Y)/(Y + xi(X-Y)), i.e.
   log(D_in(k^2)/D0) -> (k^2 - m^2)/(xi (k^2 - C)),
   C = m^2 + (a - b m^2)/(xi b), cancelling one power of (k^2 - m^2).

4. Outer loop: Feynman-combine, shift, angular average, integrate, project
   F2 at q^2 -> 0. In the Breit frame (dirac.breit_frame) with sandwiches
   A_mu = ubar(p',s') Gamma^mu u(p,s),

   F2 = -m (w A_0^{(0,0)} + m A_1^{(0,1)}) / (2 w (m^2 + w^2))
     -> -(A_0|_{w=0} + d A_1/d w|_{w=0}) / 2      (m = 1, w -> 0).

   This shortcut is validated on the LO diagram inside this script.

Twice the result (mirror) must reproduce (Petermann 1957, eq. (4))

   mu_IId = 11/24 - pi^2/18 + (1/2) log(lam^2/m^2)   [(alpha/pi)^2 units].

This script produces the parametric integrands and writes them as Fortran
into g2_iid_integrand.inc (used by g2_iid.f90 for numeric integration);
the analytic integration is a follow-up script.
"""
from sympy import (symbols, Symbol, expand, simplify, factor, cancel, log,
                   I, pi, S, zeros, fcode, gamma as Gamma_f)

from dirac import (GAMMA, ID4, METRIC, slash, breit_frame, u_spinor, ubar)
from loops import (comps, dot, feynman_shift, symmetrize, loop_integral, LUV)

m = S(1)
w, lam = symbols("w lam", positive=True)
u, xi = symbols("u xi", positive=True)
x, y, z, t = symbols("x y z t", positive=True)
e2 = Symbol("e2", positive=True)
l2 = Symbol("l2", positive=True)
K2 = Symbol("K2")

p, pp, q = breit_frame(m, w)
k = comps("k")
kp = [k[mu] + q[mu] for mu in range(4)]

SPINS = [(0, 0, 0), (1, 0, 1)]  # (mu, s', s) sandwiches used for F2


def sandwich(mu, sp, s, mid):
    """ubar(p',s') [sum_nu g_nn gamma^nu (k'+m) gamma^mu (k+m) mid (k+m)
    gamma^nu] u(p,s), progressively, as a scalar."""
    ub = ubar(pp, m, sp)
    uu = u_spinor(p, m, s)
    tot = S(0)
    kpm = slash(kp) + m * ID4
    km = slash(k) + m * ID4
    for nu in range(4):
        left = ub * GAMMA[nu] * kpm * GAMMA[mu] * km
        if mid is not None:
            tot += METRIC[nu] * (left * mid * km * GAMMA[nu] * uu)[0, 0]
        else:
            tot += METRIC[nu] * (left * GAMMA[nu] * uu)[0, 0]
    return tot


def project_F2_q0(pieces):
    """pieces: dict {(mu,sp,s): A value (function of w)}; returns F2 at
    w -> 0 via F2 = -(A_0|_0 + dA_1/dw|_0)/2."""
    A0 = pieces[(0, 0, 0)]
    A1 = pieces[(1, 0, 1)]
    A0w0 = A0.subs(w, 0)
    dA1 = A1.diff(w).subs(w, 0)
    return -(A0w0 + dA1) / 2


def assemble(mid, denom_weights, extra_pref, name):
    """Outer-loop assembly; returns the F2(q^2=0) parametric integrand."""
    Dcomb = S(0)
    npow = 0
    wt = S(1)
    for (P, par, nu) in denom_weights:
        Dcomb += par * P
        npow += nu
        wt *= par ** (nu - 1) / Gamma_f(nu)
    wt *= Gamma_f(npow)
    first = denom_weights[0][1]
    others = sum(par for (_, par, _) in denom_weights[1:])
    Dcomb = Dcomb.subs(first, 1 - others)
    wt = wt.subs(first, 1 - others)
    shift, Delta = feynman_shift(Dcomb, k)
    Delta = expand(Delta)
    shiftmap = dict(zip(k, [k[b] + shift[b] for b in range(4)]))

    vals = {}
    for (mu, sp, s) in SPINS:
        A = sandwich(mu, sp, s, mid)
        A = expand(A.subs(shiftmap, simultaneous=True))
        # only w^0 and w^1 needed; series-truncate early
        A = A + Symbol("_O") * 0
        A = expand(A.series(w, 0, 2).removeO())
        A = symmetrize(A, k, l2)
        val = -I * e2 * wt * extra_pref * sum(
            A.coeff(l2, a) * loop_integral(a, npow, Delta)
            for a in range(3) if A.coeff(l2, a) != 0)
        vals[(mu, sp, s)] = expand(val.series(w, 0, 2).removeO())
    F2 = project_F2_q0(vals)
    F2 = expand(F2)
    assert not F2.has(LUV), f"{name}: LUV must not appear in F2"
    print(f"{name}: F2 integrand assembled")
    return F2


# ------------------------------------------------------------ validation
print("validating projector shortcut on the LO diagram ...")
Pph = dot([k[a] - p[a] for a in range(4)],
          [k[a] - p[a] for a in range(4)]) - lam**2
Pkp = dot(kp, kp) - m**2
Pk = dot(k, k) - m**2
F2_lo = assemble(None, [(Pph, x, 1), (Pkp, y, 1), (Pk, z, 1)], S(1), "LO")
from sympy import integrate
res = integrate(integrate(F2_lo.subs(lam, 0), (y, 0, 1 - z)), (z, 0, 1))
res = simplify(res * 16 * pi**2 / e2)
print("  LO F2(0) =", res, "* e2/(16 pi^2)")
assert res == 2, "LO validation failed"
print("  LO validation OK: a_e = alpha/(2 pi)")

# ---------------------------------------------------------------- step 1
k2v = comps("k2v")
Din_D = (u * (dot([k2v[a] - k[a] for a in range(4)],
                  [k2v[a] - k[a] for a in range(4)]) - lam**2)
         + (1 - u) * (dot(k2v, k2v) - m**2))
sh_in, Delta_in = feynman_shift(Din_D, k2v)
Ninner = zeros(4, 4)
for a in range(4):
    Ninner += METRIC[a] * GAMMA[a] * (slash(k2v) + m * ID4) * GAMMA[a]
Ninner = Ninner.subs(dict(zip(k2v, [k2v[b] + sh_in[b] for b in range(4)])),
                     simultaneous=True)
Ninner = Ninner.applyfunc(lambda entry: symmetrize(entry, k2v, l2))
assert not any(Ninner[i, j].has(l2) for i in range(4) for j in range(4))
Sigma = (-I * e2 * loop_integral(0, 2, Delta_in)) * Ninner
a_in = expand(Delta_in.subs(dict(zip(k, [0, 0, 0, 0]))))
b_in = -(expand(Delta_in - a_in)).coeff(k[0] ** 2)
assert expand(Delta_in - (a_in - b_in * dot(k, k))) == 0
print("D_in = a - b k^2 with a =", factor(a_in), ", b =", factor(b_in))

# replace log(Delta_in(k)) (any equivalent form) by log(a - b K2)
SigK = Sigma.applyfunc(lambda entry: entry.subs(
    {lg: log(a_in - b_in * K2) for lg in entry.atoms(log)}))
A_f = expand(SigK.trace() / 4)
BslashK = (SigK - A_f * ID4).applyfunc(expand)
B_f = cancel(BslashK[0, 0] / k[0])
assert (BslashK - B_f * slash(k)).applyfunc(expand) == zeros(4, 4)
assert not A_f.has(*k) and not B_f.has(*k)
print("A =", A_f)
print("B =", B_f)

# ---------------------------------------------------------------- step 2
dm = (A_f + m * B_f).subs(K2, m**2)
dZ2 = expand(2 * m * A_f.diff(K2) + B_f + 2 * m**2 * B_f.diff(K2)).subs(
    K2, m**2)
LR = Symbol("LR")  # log(D_in(k^2)/D0)
D0 = expand(a_in - b_in * m**2)
subsLR = {log(a_in - b_in * K2): LR + log(D0)}
Sigma_R_K = (SigK - dm * ID4 - (slash(k) - m * ID4) * dZ2).applyfunc(
    lambda entry: expand(entry.subs(subsLR)))
assert not Sigma_R_K.has(LUV), "LUV must cancel in on-shell Sigma_R"
print("on-shell subtraction OK: LUV cancelled")

M_log = Sigma_R_K.applyfunc(lambda entry: entry.coeff(LR, 1))
M_rat = Sigma_R_K.applyfunc(lambda entry: entry.coeff(LR, 0))
assert ((M_rat + M_log * LR - Sigma_R_K).applyfunc(expand)
        == zeros(4, 4))
M_rat = M_rat.applyfunc(lambda entry: expand(entry.subs(K2, dot(k, k))))
M_log = M_log.applyfunc(lambda entry: expand(entry.subs(K2, dot(k, k))))

# ---------------------------------------------------------------- step 3+4
# C is kept opaque during assembly (its value in terms of u, xi below);
# this keeps the intermediate expressions small.
Cs = Symbol("C", positive=True)
C_val = m**2 + (a_in - b_in * m**2) / (xi * b_in)
PC = dot(k, k) - Cs

print("assembling rational piece ...")
F2r = assemble(M_rat, [(Pph, x, 1), (Pkp, y, 1), (Pk, z, 2)], S(1), "rat")
print("assembling log piece ...")
F2l = assemble(M_log, [(Pph, x, 1), (Pkp, y, 1), (Pk, z, 1), (PC, t, 1)],
               1 / xi, "log")

# ------------------------------------------------------------- output
# mu_IId = 2 * (int F2r + int F2l) / (alpha/pi)^2, (alpha/pi)^2 = e2^2/(16 pi^4)
import pickle
with open("g2_iid_cache.pkl", "wb") as fh:
    pickle.dump({"F2r": F2r, "F2l": F2l, "C_val": C_val}, fh)
print("cached raw F2r, F2l -> g2_iid_cache.pkl")

from sympy import cancel as _cancel, Rational
f_rat = _cancel(2 * F2r * 16 * pi**4 / e2**2)
# F2l is exactly homogeneous of degree 2 in e2 and -4 in pi (two loop
# factors e2/(16 pi^2)); verify at a random point, then strip them instead
# of running an expensive cancel() on the large expression.
pt = {y: Rational(1, 7), z: Rational(2, 7), t: Rational(1, 5),
      u: Rational(3, 8), xi: Rational(2, 9), lam: Rational(1, 10),
      Cs: Rational(50, 7)}
v1 = F2l.subs(pt).subs([(e2, 1), (pi, 1)])
assert F2l.subs(pt).subs(pi, 1).subs(e2, 2) / v1 == 4
assert F2l.subs(pt).subs(e2, 1).subs(pi, 2) / v1 == Rational(1, 16)
f_log = 2 * 16 * F2l.subs(e2, 1).subs(pi, 1)
assert f_rat.free_symbols <= {y, z, u, lam}
assert f_log.free_symbols <= {y, z, t, u, xi, lam, Cs}
print("\nmu_IId = int f_rat dy dz du  +  int f_log dy dz dt du dxi,")
print("with C =", factor(C_val))
print("f_rat =", f_rat)
print("f_log = (%d-term expression, see g2_iid_flog.inc)"
      % len(f_log.as_ordered_terms()))

with open("g2_iid_frat.inc", "w") as f:
    f.write("! generated by g2_iid.py -- do not edit\n")
    f.write("! f_rat(y,z,u,lam): region y+z<1, u in (0,1)\n")
    f.write(fcode(f_rat, assign_to="f_rat", source_format="free",
                  standard=2008) + "\n")
with open("g2_iid_flog.inc", "w") as f:
    f.write("! generated by g2_iid.py -- do not edit\n")
    f.write("! f_log(y,z,t,u,xi,lam): region y+z+t<1, u,xi in (0,1)\n")
    f.write(fcode(C_val, assign_to="C", source_format="free",
                  standard=2008) + "\n")
    f.write(fcode(f_log, assign_to="f_log", source_format="free",
                  standard=2008) + "\n")
import pickle as _p
with open("g2_iid_final.pkl", "wb") as fh:
    _p.dump({"f_rat": f_rat, "f_log": f_log, "C_val": C_val}, fh)
import json as _json
from sympy import srepr as _srepr
with open("g2_iid_integrands.json", "w") as fh:
    _json.dump({"f_rat": _srepr(f_rat), "f_log": _srepr(f_log),
                "C_val": _srepr(C_val)}, fh)
print("\nwrote g2_iid_frat.inc, g2_iid_flog.inc, g2_iid_integrands.json")
