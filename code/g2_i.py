"""Diagram I (crossed ladder): the only irreducible fourth-order vertex
diagram. No subdivergences, no renormalization: a single two-loop
assembly.

String (from ubar(p') to u(p); photon A (momentum a) contracts vertices
1 and 4 counted from the ubar side, photon B (momentum b) vertices 2 and
5; external photon q in the middle):

  ubar g^al S(p'-a) g^be S(p'-a-b) g^mu S(p-a-b) g_al S(p-b) g_be u
    / [(a^2-lam^2)(b^2-lam^2) (fermion propagators)]

Inner loop over a: four a-dependent denominators (box):
  a^2-lam^2, (p'-a)^2-1, (p'-a-b)^2-1, (p-a-b)^2-1
=> Feynman parameters (u, v, r) with weight Gamma(4)=6; the a-numerator
is cubic, so only J(0,4) ~ 1/Da^2 and J(1,4) ~ 1/Da appear: the inner
loop is UV finite (no logs, no LUV).

Da(b) = ahat + L(b) - bhat b^2; the two branches 1/Da^2 and 1/Da are
rewritten with the propagator PA = b^2 - L/bhat - CA (CA = ahat/bhat
kept opaque) and combined with the outer denominators
  b^2-lam^2 [x], (p-b)^2-1 [y], PA [t, power 1 or 2].
The b-numerator reaches degree 4, so J(2, npow) log terms can appear at
npow = 5; LUV must cancel in the F2 projection (rational-point test).

Target (Petermann 1957, eq. (1)):
  mu_I = 1/6 + (13/36) pi^2 + (5/4) zeta(3) - (5/6) pi^2 log 2
       = -0.467645...    (lambda-finite; no mirror factor)
"""
import json
import time

from sympy import (symbols, Symbol, S, I, pi, log, expand, together,
                   simplify, factor, cancel, fraction, lambdify, zeros,
                   fcode, srepr, Rational, Poly, gamma as Gamma_f)

from dirac import GAMMA, ID4, METRIC, slash, breit_frame, u_spinor, ubar
from loops import (comps, dot, feynman_shift, feynman_shift_general,
                   symmetrize, loop_integral, LUV)

T0 = time.time()


def tick(m):
    print("[%6.1fs] %s" % (time.time() - T0, m), flush=True)


m = S(1)
w, lam = symbols("w lam", positive=True)
u, v, r = symbols("u v r", positive=True)
x, y, t = symbols("x y t", positive=True)
e2 = Symbol("e2", positive=True)
la2 = Symbol("la2", positive=True)
lb2 = Symbol("lb2", positive=True)
CA = Symbol("CA", positive=True)

p, pp, q = breit_frame(m, w)
_sq = (1 + w**2) ** Rational(1, 2)
p = [expand(S(c).subs(_sq, 1 + w**2 / 2)) for c in p]
pp = [expand(S(c).subs(_sq, 1 + w**2 / 2)) for c in pp]


def spinor_series(vec):
    return vec.applyfunc(
        lambda c: c.subs(_sq, 1 + w**2 / 2).series(w, 0, 3).removeO())


a = comps("a")
b = comps("b")

# ---------------- inner combination over a ----------------
pa = [pp[i] - a[i] for i in range(4)]
pab = [pp[i] - a[i] - b[i] for i in range(4)]
mab = [p[i] - a[i] - b[i] for i in range(4)]
mb = [p[i] - b[i] for i in range(4)]

D_in = (x * (dot(a, a) - lam**2) + u * (dot(pa, pa) - 1)
        + v * (dot(pab, pab) - 1) + r * (dot(mab, mab) - 1))
D_in = D_in.subs(x, 1 - u - v - r)
sh_in, Delta_in = feynman_shift(D_in, a)
Delta_in = expand(Delta_in)
ahat = expand(Delta_in.subs(dict(zip(b, [0, 0, 0, 0]))))
bhat = -(expand(Delta_in - ahat)).coeff(b[0] ** 2)
Lb = expand(Delta_in - ahat + bhat * dot(b, b))
assert expand(Delta_in - (ahat + Lb - bhat * dot(b, b))) == 0
assert not ahat.has(*b) and not bhat.has(*b)
# w-truncation of ahat (O(w^2) irrelevant for the O(w) projection)
resid = expand(ahat - expand(ahat.subs(w, 0)))
assert resid == 0 or all(mm[0] >= 2 for mm in Poly(resid, w).monoms())
ahat = expand(ahat.subs(w, 0))
tick("inner shift: bhat = %s" % factor(bhat))
shmap_in = dict(zip(a, [a[i] + sh_in[i] for i in range(4)]))

PA = dot(b, b) - Lb / bhat - CA
CA_val = cancel(ahat / bhat)
Pph = dot(b, b) - lam**2
Pmb = dot(mb, mb) - 1

SPINS = [(0, 0, 0), (1, 0, 1)]
PREF = -I * e2 * (-I * e2) * 6   # two photon+vertex-pair factors, Gamma(4)
# bookkeeping: (-i e^2) per loop as in the LO/IIc pattern; inner Feynman
# weight Gamma(4) = 6; outer weights added per branch below.


def sandwich_inner(mu_, sp, s_):
    """Full string sandwich, inner-loop (a) done: returns dict
    {1: coeff of 1/Da, 2: coeff of 1/Da^2} as polynomials in b, w-linear."""
    ub = spinor_series(ubar(breit_frame(m, w)[1], m, sp))
    uu = spinor_series(u_spinor(breit_frame(m, w)[0], m, s_))
    tot = S(0)
    for al in range(4):
        for be in range(4):
            tot += (METRIC[al] * METRIC[be] *
                    (ub * GAMMA[al] * (slash(pa) + ID4) * GAMMA[be]
                     * (slash(pab) + ID4) * GAMMA[mu_]
                     * (slash(mab) + ID4) * GAMMA[al]
                     * (slash(mb) + ID4) * GAMMA[be] * uu)[0, 0])
    tot = expand(tot.subs(shmap_in, simultaneous=True))
    tot = expand(tot.series(w, 0, 2).removeO())
    tot = symmetrize(tot, a, la2)
    # attach inner J's: J(j,4,Da): j=0 -> 1/Da^2 * (i/96pi^2-type), j=1 -> 1/Da
    out = {1: S(0), 2: S(0)}
    for j in range(3):
        cj = tot.coeff(la2, j)
        if cj == 0:
            continue
        assert j <= 1, "unexpected la2^2 term (would be UV divergent)"
        Jg = loop_integral(j, 4, Symbol("_D"))
        # Jg = const * _D**(j-2): split constant and power
        constJ = Jg * Symbol("_D") ** (2 - j)
        assert not constJ.has(Symbol("_D"))
        out[2 - j] += cj * constJ
    return out


def assemble_outer(cb, kpa, name):
    """Outer loop for numerator cb (poly in b, w-linear) with denominators
    [Pph x, Pmb y, PA t^kpa]; returns the F2-projection pieces."""
    Dcomb = x * Pph + y * Pmb + t * PA
    npow = 2 + kpa
    wt = Gamma_f(npow) / Gamma_f(kpa) * t ** (kpa - 1)
    Dcomb = Dcomb.subs(x, 1 - y - t)
    wt = wt.subs(x, 1 - y - t)
    shift, Delta, Acoef = feynman_shift_general(Dcomb, b)
    assert Acoef == 1
    shmap = dict(zip(b, [b[i] + shift[i] for i in range(4)]))
    Dw0 = Delta.subs(w, 0)
    Dw1 = Delta.diff(w).subs(w, 0)

    A = expand(cb.subs(shmap, simultaneous=True))
    A = expand(A.series(w, 0, 2).removeO())
    A = symmetrize(A, b, lb2)
    val0, val1 = S(0), S(0)
    for j in range(4):
        cj = A.coeff(lb2, j)
        if cj == 0:
            continue
        Jg = loop_integral(j, npow, Symbol("_DD"))
        J0 = Jg.subs(Symbol("_DD"), Dw0)
        J1 = Jg.diff(Symbol("_DD")).subs(Symbol("_DD"), Dw0) * Dw1
        val0 += cj.subs(w, 0) * J0
        val1 += cj.diff(w).subs(w, 0) * J0 + cj.subs(w, 0) * J1
    pref = PREF * (-1 / bhat) ** kpa * wt
    return val0 * pref, val1 * pref


vals = {}
for (mu_, sp, s_) in SPINS:
    tt0 = time.time()
    parts = sandwich_inner(mu_, sp, s_)
    tick("inner sandwich mu=%d done (%.0fs)" % (mu_, time.time() - tt0))
    v0tot, v1tot = S(0), S(0)
    for kpa in (1, 2):
        if parts[kpa] == 0:
            continue
        v0, v1 = assemble_outer(parts[kpa], kpa, f"mu{mu_}-PA^{kpa}")
        v0tot += v0
        v1tot += v1
    vals[(mu_, sp, s_)] = (v0tot, v1tot)
    tick("outer mu=%d done" % mu_)

A0 = vals[(0, 0, 0)][0]
dA1 = vals[(1, 0, 1)][1]
F2 = -(A0 + dA1) / 2
if F2.has(LUV):
    import random
    random.seed(11)
    dF = F2.diff(LUV)
    free = sorted(dF.free_symbols, key=str)
    iszero = True
    for _ in range(3):
        pt = {sym: Rational(random.randint(3, 97), 101) for sym in free}
        if dF.subs(pt, simultaneous=True) != 0:
            iszero = False
            break
    print("LUV coefficient zero:", iszero, flush=True)
    assert iszero
    F2 = F2.subs(LUV, 0)
tick("F2 projection done; LUV-free")

# units: (alpha/pi)^2; strip e2^2 pi^-4 (verified by homogeneity)
pt = {u: Rational(1, 7), v: Rational(1, 5), r: Rational(1, 6),
      y: Rational(2, 7), t: Rational(1, 4), lam: Rational(1, 10),
      CA: Rational(50, 7)}
v1_ = F2.subs(pt).subs([(e2, 1), (pi, 1)])
assert F2.subs(pt).subs(pi, 1).subs(e2, 2) / v1_ == 4
assert F2.subs(pt).subs(e2, 1).subs(pi, 2) / v1_ == Rational(1, 16)
fI = 16 * F2.subs(e2, 1).subs(pi, 1)
assert fI.free_symbols <= {u, v, r, y, t, lam, CA}, fI.free_symbols

with open("g2_i_integrands.json", "w") as fh:
    json.dump({"fI": srepr(fI), "CA_val": srepr(CA_val)}, fh)


def stabilize(expr):
    from sympy import Pow as _Pow
    expr = expr.replace(lambda e: e.is_Pow and e.exp.is_negative,
                        lambda e: _Pow(cancel(together(e.base)), e.exp))
    expr = expr.replace(lambda e: isinstance(e, log),
                        lambda e: log(cancel(together(e.args[0]))))
    return expr


fst = stabilize(fI.subs(CA, CA_val))
with open("g2_i_f.inc", "w") as fh:
    fh.write("! generated by g2_i.py -- do not edit\n")
    fh.write(fcode(fst, assign_to="fI", source_format="free",
                   standard=2008) + "\n")
tick("wrote g2_i_integrands.json, g2_i_f.inc")
print("""
mu_I(lam) = int fI du dv dr dy dt   (u+v+r<1, y+t<1)
target: 1/6 + 13/36 pi^2 + 5/4 z3 - 5/6 pi^2 log2 = -0.467645...""")
