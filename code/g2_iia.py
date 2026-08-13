"""Diagram IIa (ladder): one-loop vertex part inserted at the EXTERNAL
photon vertex of the LO vertex diagram, with the same pointwise
renormalization as IIc: the inner vertex part is subtracted by its
on-shell value L(u,v) gamma^mu (integral = deltaF1(0)).

  ubar gamma^nu S(p'-k) gamma^rho S(p'-k-k2) gamma^mu S(p-k-k2)
       gamma_rho S(p-k) gamma_nu u

Inner loop (k2; u for the (p'-k-k2) line, v for the (p-k-k2) line, photon
rest). Outer propagators identical to IIc's -> same three assemblies;
NO mirror factor (the ladder is its own mirror image).

Target (Petermann 1957, eq. (2)): mu_IIa = 11/48 + pi^2/18 = 0.777478.
"""
import json
import time

from sympy import (symbols, Symbol, S, I, pi, log, oo, expand, together,
                   simplify, factor, cancel, fraction, integrate, limit,
                   lambdify, zeros, fcode, srepr, Rational, Poly,
                   gamma as Gamma_f)

from dirac import GAMMA, ID4, METRIC, slash, breit_frame, u_spinor, ubar
from loops import (comps, dot, feynman_shift, feynman_shift_general,
                   symmetrize, loop_integral, LUV)

T0 = time.time()


def tick(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


m = S(1)
w, lam = symbols("w lam", positive=True)
u, v, xi = symbols("u v xi", positive=True)
x, y, z, t = symbols("x y z t", positive=True)
e2 = Symbol("e2", positive=True)
l2 = Symbol("l2", positive=True)
CIN = Symbol("CIN", positive=True)
CXI = Symbol("CXI", positive=True)

p, pp, q = breit_frame(m, w)
# Only O(w) survives the q^2 -> 0 projection, so truncate the kinematics:
# sqrt(1+w^2) -> 1 + w^2/2 (on-shellness violated only at O(w^4)).
# This removes all square roots of w and speeds the algebra up greatly.
_sq = (1 + w**2) ** Rational(1, 2)


def truncw(e):
    return expand(S(e).subs(_sq, 1 + w**2 / 2))


p = [truncw(c) for c in p]
pp = [truncw(c) for c in pp]


def spinor_series(vec):
    return vec.applyfunc(
        lambda c: c.subs(_sq, 1 + w**2 / 2).series(w, 0, 3).removeO())


k = comps("k")
k2 = comps("k2")

# =================================================================
# inner loop: Lambda^nu(k)
# =================================================================
pkk2 = [pp[a] - k[a] - k2[a] for a in range(4)]   # p' - k - k2
mkk2 = [p[a] - k[a] - k2[a] for a in range(4)]    # p  - k - k2

D_in = (u * (dot(pkk2, pkk2) - 1) + v * (dot(mkk2, mkk2) - 1)
        + (1 - u - v) * (dot(k2, k2) - lam**2))
sh_in, Delta_in = feynman_shift(D_in, k2)
Delta_in = expand(Delta_in)
# D_in = ahat + Lk - bhat k^2, ahat k-free, Lk linear in k
ahat = expand(Delta_in.subs(dict(zip(k, [0, 0, 0, 0]))))
bhat = -(expand(Delta_in - ahat)).coeff(k[0] ** 2)
Lk = expand(Delta_in - ahat + bhat * dot(k, k))
assert Poly(Lk, *k).total_degree() <= 1 if Lk != 0 else True
assert expand(Delta_in - (ahat + Lk - bhat * dot(k, k))) == 0
assert not ahat.has(*k) and not bhat.has(*k)
tick("inner shift done: bhat = %s, Lk = %s" % (factor(bhat), Lk))

shmap_in = dict(zip(k2, [k2[b] + sh_in[b] for b in range(4)]))
Nin = {}
c_l2 = {}
for nu in range(4):
    Nn = zeros(4, 4)
    for rho in range(4):
        Nn += (METRIC[rho] * GAMMA[rho] * (slash(pkk2) + ID4) * GAMMA[nu]
               * (slash(mkk2) + ID4) * GAMMA[rho])
    Nn = Nn.subs(shmap_in, simultaneous=True)
    Nn = Nn.applyfunc(lambda entry: symmetrize(entry, k2, l2))
    Nin[nu] = Nn.applyfunc(lambda entry: entry.coeff(l2, 0))
    c_l2[nu] = Nn.applyfunc(lambda entry: entry.coeff(l2, 1))
    # the l^2 coefficient must be proportional to gamma^nu
    i0, j0 = next((i_, j_) for i_ in range(4) for j_ in range(4)
                  if GAMMA[nu][i_, j_] != 0)
    ratio = simplify(c_l2[nu][i0, j0] / GAMMA[nu][i0, j0])
    assert (c_l2[nu] - ratio * GAMMA[nu]).applyfunc(simplify) == zeros(4, 4)
    if nu == 0:
        c_const = ratio
    else:
        assert simplify(ratio - c_const) == 0
tick("inner numerator matrices done; c = %s" % c_const)

# inner prefactor: -i e2 * (Feynman 2) * (i/(16 pi^2)) J
#   J(0,3): coefficient -1/(2 Delta); J(1,3): (LUV - log Delta)
PREF_IN = -I * e2 * 2
J03 = loop_integral(0, 3, Symbol("_D"))     # -I/(32 pi^2)/_D
J13 = loop_integral(1, 3, Symbol("_D"))     # I/(16 pi^2)(LUV - log _D)
c03 = J03 * Symbol("_D")                     # scalar coefficient of 1/D
c13LUV = expand(J13).coeff(LUV)              # I/(16 pi^2)

# subtraction constant L(u,v): sandwich at k = 0, w = 0
p0 = [S(1), 0, 0, 0]
ub0 = ubar(p0, 1, 0)
uu0 = u_spinor(p0, 1, 0)
sub0 = dict(zip(k, [0, 0, 0, 0])) | {w: 0}
N00 = Nin[0].subs(sub0, simultaneous=True)
La = ((ub0 * N00 * uu0)[0, 0] / (ub0 * GAMMA[0] * uu0)[0, 0])
La = simplify(La * PREF_IN * c03)            # coefficient of 1/D0
D0 = expand(ahat.subs(sub0, simultaneous=True))
tick("subtraction: La = %s, D0 = %s" % (factor(La), factor(D0)))
# ahat is w-free up to O(w^4) (p'^2 = 1 + w^4/4 after truncation), which
# is beyond the O(w) we need: truncate, so CIN/CXI values are w-free and
# the opaque symbols can be kept through the whole assembly.
resid = expand(ahat - D0.subs(lam, lam))  # ahat - ahat(w=0)
resid = expand(ahat - expand(ahat.subs(w, 0)))
from sympy import degree
# for IIa, ahat has genuine O(w^2) terms (q flows through the inner
# loop); they cannot affect the O(w) projection, so truncation is valid
assert resid == 0 or all(mm[0] >= 2 for mm in Poly(resid, w).monoms())
ahat = expand(ahat.subs(w, 0))
# c-part of L matches Lambda's LUV part by construction (same c_const)

# =================================================================
# outer loop assembly
# =================================================================
# outer photon carries momentum k itself in this routing
Pph = dot(k, k) - lam**2
Ppk = dot([pp[a] - k[a] for a in range(4)],
          [pp[a] - k[a] for a in range(4)]) - 1
Pmk = dot([p[a] - k[a] for a in range(4)],
          [p[a] - k[a] for a in range(4)]) - 1
# 1/D_in  = -(1/bhat)   / PCIN, PCIN = k^2 - Lk/bhat - CIN, CIN = ahat/bhat
# 1/D_xi  = -(1/xi bhat)/ PCXI, PCXI = k^2 - Lk/bhat - CXI,
#           CXI = (D0(1-xi) + xi ahat)/(xi bhat)
PCIN = dot(k, k) - Lk / bhat - CIN
PCXI = dot(k, k) - Lk / bhat - CXI

SPINS = [(0, 0, 0), (1, 0, 1)]


def sandwich_iic(mu_, sp, s_, mid_of_nu, extra_num):
    """ubar gamma^nu S(p'-k) [mid^mu] S(p-k) gamma_nu u * extra_num."""
    ub = spinor_series(ubar(breit_frame(m, w)[1], m, sp))
    uu = spinor_series(u_spinor(breit_frame(m, w)[0], m, s_))
    left_a = slash([pp[a] - k[a] for a in range(4)]) + ID4
    right_b = slash([p[a] - k[a] for a in range(4)]) + ID4
    tot = S(0)
    for nu in range(4):
        tot += METRIC[nu] * (ub * GAMMA[nu] * left_a * mid_of_nu[mu_]
                             * right_b * GAMMA[nu] * uu)[0, 0]
    return expand(tot * extra_num)


def assemble_iic(mid_of_nu, extra_num, denoms, scal_pref, name):
    """Outer loop with numerator sandwich; returns F2 integrand at q^2->0.

    denoms: [(P, param, power)]. Opaque symbols (CIN/CXI) may appear in
    the propagators and scal_pref; their (w-free) values are substituted
    by the caller afterwards. No expand() is applied to the rational
    loop-integral factors (it explodes); F2 is returned unexpanded.
    """
    Dcomb, npow, wt = S(0), 0, S(1)
    for (P, par, nu_) in denoms:
        Dcomb += par * P
        npow += nu_
        wt *= par ** (nu_ - 1) / Gamma_f(nu_)
    wt *= Gamma_f(npow)
    first = denoms[0][1]
    others = sum(par for (_, par, _) in denoms[1:])
    Dcomb = Dcomb.subs(first, 1 - others)
    wt = wt.subs(first, 1 - others)
    shift, Delta, Acoef = feynman_shift_general(Dcomb, k)
    assert not any(sh.has(CIN, CXI) for sh in shift)
    shmap = dict(zip(k, [k[b] + shift[b] for b in range(4)]))

    assert Acoef == 1
    # expand Delta and the J-factors to O(w) once (Delta is polynomial
    # in w); avoids per-term series of identical Delta functions
    Dw0 = Delta.subs(w, 0)
    Dw1 = Delta.diff(w).subs(w, 0)

    def Jser(a_):
        # returns (J|_{w=0}, dJ/dw|_{w=0}) for loop_integral(a_, npow, .)
        Jgen = loop_integral(a_, npow, Symbol("_DD"))
        J0 = Jgen.subs(Symbol("_DD"), Dw0)
        J1 = Jgen.diff(Symbol("_DD")).subs(Symbol("_DD"), Dw0) * Dw1
        return J0, J1

    vals = {}
    for (mu_, sp, s_) in SPINS:
        A = sandwich_iic(mu_, sp, s_, mid_of_nu, extra_num)
        A = expand(A.subs(shmap, simultaneous=True))
        A = expand(A.series(w, 0, 2).removeO())
        A = symmetrize(A, k, l2)
        val0, val1 = S(0), S(0)
        for a_ in range(4):
            ca = A.coeff(l2, a_)
            if ca == 0:
                continue
            ca0 = ca.subs(w, 0)
            ca1 = ca.diff(w).subs(w, 0)
            J0, J1 = Jser(a_)
            val0 += ca0 * J0
            val1 += ca1 * J0 + ca0 * J1
        pref = (-I * e2) * wt * scal_pref
        vals[(mu_, sp, s_)] = (val0 * pref, val1 * pref)
        tick(f"    {name}: sandwich mu={mu_} assembled")
    A0 = vals[(0, 0, 0)][0]
    dA1 = vals[(1, 0, 1)][1]
    F2 = -(A0 + dA1) / 2
    if F2.has(LUV):
        # exact zero test of the LUV coefficient at random rational points
        import random
        random.seed(7)
        dF = F2.diff(LUV)
        free = sorted(dF.free_symbols, key=str)
        iszero = True
        for _ in range(3):
            pt = {sym: Rational(random.randint(3, 97), 101) for sym in free}
            val0 = dF.subs(pt, simultaneous=True)
            if val0 != 0:
                iszero = False
                break
        print(f"  [{name}] LUV coefficient zero: {iszero}", flush=True)
        if iszero:
            F2 = F2.subs(LUV, 0)
    print(f"  [{name}] F2 integrand done, LUV present: {F2.has(LUV)}",
          flush=True)
    return F2


CIN_val = cancel(ahat / bhat)
CXI_val = cancel((D0 * (1 - xi) + xi * ahat) / (xi * bhat))

# ---- piece (a): P^nu/D_in -> extra propagator PCIN, pref -1/bhat
tick("assembling piece (a) ...")
mid_a = {nu: PREF_IN * c03 * Nin[nu] for nu in range(4)}
F2a = assemble_iic(mid_a, S(1),
                   [(Pph, x, 1), (Ppk, y, 1), (Pmk, z, 1), (PCIN, t, 1)],
                   -1 / bhat, "a")

# ---- piece (b): -La/D0 gamma^nu, LO-type loop
tick("assembling piece (b) ...")
mid_b = {nu: -La / D0 * GAMMA[nu] for nu in range(4)}
F2b = assemble_iic(mid_b, S(1),
                   [(Pph, x, 1), (Ppk, y, 1), (Pmk, z, 1)],
                   S(1), "b")

# ---- piece (c): -c gamma^nu log(D_in/D0) via xi-representation
tick("assembling piece (c) ...")
assert c13LUV != 0 and c_const != 0
mid_c = {nu: -PREF_IN * c13LUV * c_const * GAMMA[nu] for nu in range(4)}
extra_num_c = expand(Delta_in - D0)
F2c = assemble_iic(mid_c, extra_num_c,
                   [(Pph, x, 1), (Ppk, y, 1), (Pmk, z, 1), (PCXI, t, 1)],
                   -1 / (xi * bhat), "c")

# =================================================================
# output: mu_IIc = 2*(F2a + F2b + F2c) in (alpha/pi)^2 units
# =================================================================
def stabilize(expr):
    """cancel() every negative-power base and log argument into a single
    polynomial ratio; avoids catastrophic (u-1) cancellation when the
    opaque CIN/CXI values are substituted for numeric evaluation."""
    from sympy import Pow as _Pow
    expr = expr.replace(
        lambda e: e.is_Pow and e.exp.is_negative,
        lambda e: _Pow(cancel(together(e.base)), e.exp))
    expr = expr.replace(
        lambda e: isinstance(e, log),
        lambda e: log(cancel(together(e.args[0]))))
    return expr


out = {"CIN_val": CIN_val, "CXI_val": CXI_val}
for nm, F2p, syms, opq in [("fa", F2a, (y, z, t, u, v, lam), [(CIN, CIN_val)]),
                           ("fb", F2b, (y, z, u, v, lam), []),
                           ("fc", F2c, (y, z, t, u, v, xi, lam),
                            [(CXI, CXI_val)])]:
    fexpr = 1 * 16 * F2p.subs(e2, 1).subs(pi, 1)  # no mirror for IIa
    extra = fexpr.free_symbols - set(syms) - {o[0] for o in opq}
    assert not extra, (nm, extra)
    out[nm] = fexpr
    fstab = stabilize(fexpr.subs(dict(opq)))
    with open(f"g2_iia_{nm}.inc", "w") as fh:
        fh.write("! generated by g2_iic.py -- do not edit\n")
        fh.write(fcode(fstab, assign_to=nm, source_format="free",
                       standard=2008) + "\n")
    print(nm, "terms:", len(fexpr.as_ordered_terms()), flush=True)

with open("g2_iia_integrands.json", "w") as fh:
    json.dump({nm: srepr(expr) for nm, expr in out.items()}, fh)
tick("wrote g2_iic_fa/fb/fc.inc and g2_iia_integrands.json")
print("""
mu_IIa(lam) = int fa + int fb + int fc   (same regions as IIc)
target: 11/48 + pi^2/18 = 0.777478...  (UV and IR finite)""")
