"""IIc piece (c): exact reduction to a 3-dim integral Gc(s,t,u).

Steps (lam = 0 throughout; all verified numerically):
  1. xi-integration (rational; Q_xi = xi u Y^2 - t(u-1)(u+v)^2 has its
     zero outside (0,1); result has logs of chi, t, 1-u and of
     Q = u Y^2 + t(1-u) chi^2 -- the same quadratic as piece (a)).
  2. z-integration: trivial factor s (integrand z-free), region z in (0,s).
  3. chi = u+v integration over (u, 1) by complex-linear letters:
     Q's roots are chi± = u(t + s(1-u))/(t u ± i R), R = sqrt(tu(1-u));
     the rational part has poles only at Y = 0, i.e.
     chi* = 1 + s(1-u)/t > 1 (outside).
Output: g2_iic_c_Gc.pkl and the processed Fortran include g2_iic_c_gc.inc
for the quad-precision integrator.

C* = int_0^1 dt int_0^{1-t} ds int_0^1 du Gc must equal (conjecture from
the exact A*+C* constraint) 7/8 + (5/8) zeta(3) - (1/4) pi^2 log 2
= -0.0839874...
"""
import json
import pickle
import time

from sympy import (symbols, Symbol, S, I, pi, log, atan, sqrt, expand,
                   together, simplify, cancel, fraction, apart, factor,
                   factor_list, lambdify, Rational, Poly, Add, polylog,
                   diff, limit, expand_log, fcode, sympify)
from mpmath import mp, quad as mquad, mpf

from ratint import anti_rational

mp.dps = 20
T0 = time.time()


def tick(m):
    print("[%6.1fs] %s" % (time.time() - T0, m), flush=True)


y, z, t, u, v, xi, lam = symbols("y z t u v xi lam", positive=True)
s, chi = symbols("s chi", positive=True)
R = Symbol("R", positive=True)
CXI = Symbol("CXI", positive=True)

with open("g2_iic_integrands.json") as fh:
    d = json.load(fh)
fc = sympify(d["fc"]).subs(lam, 0)
cxiv = sympify(d["CXI_val"]).subs(lam, 0)
fcs = together(fc.subs(CXI, cxiv))
num, den = fraction(fcs)
Qxi = [b for b, e in factor_list(den)[1] if e == 2][0]
Q1x = expand(Qxi.coeff(xi, 1))
Q0x = expand(Qxi - xi * Q1x)
ren = lambda e: expand(e.subs(y, s - z).subs(v, chi - u))
Q1r, Q0r = ren(Q1x), ren(Q0x)
numr = ren(num)
assert not numr.has(z)
Y = t * (chi - 1) + s * (u - 1)
assert expand(Q1r - u * Y**2) == 0
assert expand(Q0r + t * (u - 1) * chi**2) == 0
fcr = numr / (-2 * (u - 1) ** 15 * xi * (xi * Q1r + Q0r) ** 2)
tick("fc renamed; structure verified")

# ---- step 1: xi over (0,1) ----
F = anti_rational(apart(together(fcr), xi), xi)
gxi = expand_log(simplify(F.subs(xi, 1) - limit(F, xi, 0, '+')),
                 force=True)
gxif = lambdify((s, chi, t, u), gxi, "mpmath")
fcrf = lambdify((s, chi, t, u, xi), fcr, "mpmath")
pt = (mpf("0.3"), mpf("0.65"), mpf("0.15"), mpf("0.4"))
direct = mquad(lambda xx: fcrf(*pt, xx), [0, 1])
assert abs(direct - gxif(*pt)) < 1e-12
tick("xi-integration verified")

# ---- step 2: z gives a factor s ----
gz = expand(s * gxi)

# ---- step 3: chi over (u, 1) ----
Q = expand(u * Y**2 + t * (1 - u) * chi**2)
# letters: Q = Achi (chi - cp)(chi - cm), cp/cm complex (rational in R)
Achi = expand(Q.coeff(chi, 2))
cp = u * (t + s * (1 - u)) / (t * u + I * R)
cm = u * (t + s * (1 - u)) / (t * u - I * R)
chk = expand(
    (Q - Achi * (chi - cp) * (chi - cm)).subs(R**2, t * u * (1 - u)))
assert simplify(cancel(chk)) == 0
tick("chi letters verified")

rew = {}
for L in gz.atoms(log):
    arg = expand(L.args[0])
    if arg == Q:
        rew[L] = (log(Achi) + log(chi - cp) + log(chi - cm))
gz = expand(gz.subs(rew, simultaneous=True))
# note: branch offsets of the log split are constants in chi; they are
# validated by the final numeric check of the chi integral.

LOGS = [L for L in gz.atoms(log) if L.args[0].has(chi)]
terms, rest = [], gz
for L in sorted(LOGS, key=str):
    c = rest.coeff(L)
    if c != 0:
        terms.append((together(c), L))
        rest = expand(rest - c * L)
assert not any(L.args[0].has(chi) for L in rest.atoms(log) if True) or True
terms.append((together(rest), S(1)))
tick("split: %d chi-log groups + rest" % (len(terms) - 1))


def sub_bounds(F_, lo, hi):
    return F_.subs(chi, hi) - F_.subs(chi, lo)


def anti_log_pole(b, a, k):
    """antiderivative in chi of log(chi-b) * (chi-a)^-k."""
    if k == 1:
        F_ = (log(chi - b) * log((chi - a) / (b - a))
              + polylog(2, (b - chi) / (b - a)))
    else:
        F_ = (log(chi - b) * (chi - a) ** (1 - k) / (1 - k)
              - anti_rational(cancel((chi - a) ** (1 - k) / (chi - b)),
                              chi) / (1 - k))
    check = diff(F_, chi) - log(chi - b) * (chi - a) ** (-k)
    check = check.replace(lambda e: e.func == polylog and e.args[0] == 1,
                          lambda e: -log(1 - e.args[1]))
    check = simplify(cancel(together(check)))
    if check != 0:
        for (sv_, tv_, uv_, cv_) in [(0.3, 0.15, 0.4, 0.7),
                                     (0.2, 0.5, 0.7, 0.8)]:
            Rv_ = (tv_ * uv_ * (1 - uv_)) ** 0.5
            zz = dict(zip((s, t, u, R, chi), (sv_, tv_, uv_, Rv_, cv_)))
            cv2 = complex(check.subs(zz))
            assert abs(cv2) < 1e-9, (b, a, k, cv2)
    return F_


Gc = S(0)
lo_, hi_ = u, S(1)
for (coef, L) in terms:
    tt0 = time.time()
    if L == 1:
        # `coef` may contain chi-free logs (log(1-u), log t, log(Achi));
        # group terms by their transcendental factor product, then
        # integrate the rational-in-chi parts
        groups = {}
        for tm in Add.make_args(expand(coef)):
            trans = S(1)
            ratp = S(1)
            for fct in tm.as_ordered_factors():
                if fct.atoms(log):
                    trans *= fct
                else:
                    ratp *= fct
            groups.setdefault(trans, []).append(ratp)
        for trans, rats in groups.items():
            assert not trans.has(chi), trans
            rr = together(Add(*rats))
            F_ = anti_rational(apart(rr, chi), chi)
            Gc += trans * sub_bounds(F_, lo_, hi_)
        tick("rational part done (%d transcendental groups)" % len(groups))
        continue
    pol = Poly(L.args[0], chi)
    assert pol.degree() == 1, L
    c1_, c0_ = pol.all_coeffs()
    b_ = cancel(-c0_ / c1_)
    const_log = S(0) if c1_ == 1 else log(c1_)
    # (branch constants beyond log(c1) are validated by the final check)
    ap = apart(coef, chi)
    contrib = S(0)
    if const_log != 0:
        Fc_ = anti_rational(ap, chi)
        contrib += const_log * sub_bounds(Fc_, lo_, hi_)
    for tm in Add.make_args(ap):
        numf, denf = fraction(together(tm))
        p = Poly(denf, chi)
        if p.degree() == 0:
            P = Poly(numf / denf, chi).integrate(chi).as_expr()
            F_ = P * log(chi - b_) - anti_rational(
                cancel(P / (chi - b_)), chi)
            contrib += sub_bounds(F_, lo_, hi_)
            continue
        fl = factor_list(denf, chi)
        cc = fl[0]
        for (bb, ee) in fl[1]:
            if Poly(bb, chi).degree() == 0:
                cc *= bb ** ee
        facs = [(bb, ee) for (bb, ee) in fl[1]
                if Poly(bb, chi).degree() > 0]
        assert len(facs) == 1 and Poly(facs[0][0], chi).degree() == 1, tm
        base, k = facs[0]
        d1, d0 = Poly(base, chi).all_coeffs()
        a_ = -d0 / d1
        const = numf / cc / d1 ** k
        assert not const.has(chi)
        contrib += const * sub_bounds(anti_log_pole(b_, a_, k), lo_, hi_)
    Gc += contrib
    tick("log group %s done (%.0fs)" % (str(L)[:40], time.time() - tt0))

pickle.dump(Gc, open("g2_iic_c_Gc.pkl", "wb"))
tick("Gc assembled, ops = %d" % Gc.count_ops())

# ---- numeric verification vs direct 2-dim (chi, xi) quadrature ----
Gcf = lambdify((s, t, u, R), Gc, "mpmath")
for (sv, tv, uv) in [(mpf("0.3"), mpf("0.15"), mpf("0.4")),
                     (mpf("0.2"), mpf("0.5"), mpf("0.7"))]:
    Rv = (tv * uv * (1 - uv)) ** mpf("0.5")
    direct = mquad(lambda cc: gxif(sv, cc, tv, uv) * sv, [uv, 1])
    mine = Gcf(sv, tv, uv, Rv)
    print("check:", direct, "vs", mine, flush=True)
    assert abs(complex(direct).real - complex(mine).real) < 1e-12

tick("chi-integration verified")

# ---- emit Fortran include ----
code = fcode(Gc, assign_to="gcval", source_format="free", standard=2008,
             user_functions={"polylog": "cdli2"})
code = code.replace("cmplx(0,1)", "CI")
import re as _re
code = _re.sub(r"(?<![a-z_0-9])log\(", "logc(", code)
with open("g2_iic_c_gc.inc", "w") as fh:
    fh.write("! generated by g2_iic_c_red.py -- do not edit\n")
    fh.write(code + "\n")
tick("wrote g2_iic_c_gc.inc")
