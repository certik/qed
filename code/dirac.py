"""Explicit Dirac algebra for the g-2 project.

Conventions: metric (+,-,-,-), Dirac representation. Momenta are lists of
4 *contravariant* components [p0, p1, p2, p3].

Kinematics used everywhere (Breit frame): the external photon carries
q = p' - p along z,

    p  = (E, 0, 0,  w),   p' = (E, 0, 0, -w),   E = sqrt(m^2 + w^2),
    q  = (0, 0, 0, -2w),  q^2 = -4 w^2,

so on-shellness is automatic and the q^2 -> 0 limit is w -> 0.

F_2 extraction: sandwich the vertex between explicit on-shell spinors and
match against the form-factor decomposition

    ubar(p') Gamma^mu u(p)
      = ubar(p') [ gamma^mu F_1 + i sigma^{mu nu} q_nu /(2m) F_2 ] u(p).

For mu = 1 (x) and spins (s', s) = (up, up) both structures are nonzero
and independent, so two spin/mu choices determine F_1 and F_2.
"""
from sympy import (Matrix, I, eye, zeros, sqrt, symbols, simplify, expand,
                   Rational)

METRIC = (1, -1, -1, -1)

g0 = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]])
g1 = Matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0], [-1, 0, 0, 0]])
g2 = Matrix([[0, 0, 0, -I], [0, 0, I, 0], [0, I, 0, 0], [-I, 0, 0, 0]])
g3 = Matrix([[0, 0, 1, 0], [0, 0, 0, -1], [-1, 0, 0, 0], [0, 1, 0, 0]])
GAMMA = [g0, g1, g2, g3]
ID4 = eye(4)


def dot(a, b):
    """Minkowski dot product of two 4-component (contravariant) momenta."""
    return sum(METRIC[m] * a[m] * b[m] for m in range(4))


def slash(p):
    """p-slash = gamma^mu p_mu = gamma^mu g_{mu nu} p^nu."""
    s = zeros(4, 4)
    for mu in range(4):
        s += METRIC[mu] * p[mu] * GAMMA[mu]
    return s


def sigma(mu, nu):
    """sigma^{mu nu} = (i/2)[gamma^mu, gamma^nu]."""
    return I * (GAMMA[mu] * GAMMA[nu] - GAMMA[nu] * GAMMA[mu]) / 2


def u_spinor(p, m, s):
    """On-shell u(p): (p-slash - m) u = 0, normalized ubar u = 2m.

    s = 0, 1: spin up/down two-spinor xi in the Dirac representation.
    p must satisfy p0 = E = sqrt(m^2 + |vec p|^2) (checked by caller).
    """
    E = p[0]
    xi = Matrix([1, 0]) if s == 0 else Matrix([0, 1])
    sig = [Matrix([[0, 1], [1, 0]]), Matrix([[0, -I], [I, 0]]),
           Matrix([[1, 0], [0, -1]])]
    sp = zeros(2, 2)
    for i in range(3):
        sp += p[i + 1] * sig[i]
    upper = sqrt(E + m) * xi
    lower = sp * xi / sqrt(E + m)
    return Matrix([upper[0], upper[1], lower[0], lower[1]])


def ubar(p, m, s):
    return u_spinor(p, m, s).H * g0


def breit_frame(m, w):
    """Return (p, pprime, q) in the Breit frame; q^2 = -4 w^2."""
    E = sqrt(m**2 + w**2)
    p = [E, 0, 0, w]
    pp = [E, 0, 0, -w]
    q = [0, 0, 0, -2 * w]
    return p, pp, q


def formfactor_basis(m, w, mu, sp, s):
    """(B1, B2) = ubar gamma^mu u, ubar i sigma^{mu nu} q_nu u / (2m)."""
    p, pp, q = breit_frame(m, w)
    ub, uu = ubar(pp, m, sp), u_spinor(p, m, s)
    B1 = (ub * GAMMA[mu] * uu)[0, 0]
    S = zeros(4, 4)
    for nu in range(4):
        S += sigma(mu, nu) * METRIC[nu] * q[nu]
    B2 = (ub * (I * S / (2 * m)) * uu)[0, 0]
    return simplify(B1), simplify(B2)


def extract_F1_F2(vertex_fn, m, w):
    """Given vertex_fn(mu) -> 4x4 matrix Gamma^mu (explicit components),
    return (F1, F2) by matching two independent spinor sandwiches.

    The mu = 1 spin-flip sandwich alone is degenerate (it measures the
    magnetic combination F1 + F2: B1 = B2 = -2w), so we pair it with the
    mu = 0 no-flip sandwich (B1 = 2m, B2 = -2w^2/m):

        (mu, s', s) = (0, 0, 0)  and  (1, 0, 1).
    """
    from sympy import linsolve, symbols as syms
    F1s, F2s = syms("F1s F2s")
    p, pp, q = breit_frame(m, w)
    eqs = []
    for (mu, sp, s) in [(0, 0, 0), (1, 0, 1)]:
        ub, uu = ubar(pp, m, sp), u_spinor(p, m, s)
        A = (ub * vertex_fn(mu) * uu)[0, 0]
        B1, B2 = formfactor_basis(m, w, mu, sp, s)
        eqs.append(expand(A - F1s * B1 - F2s * B2))
    sol = linsolve(eqs, [F1s, F2s])
    (f1, f2), = sol
    return simplify(f1), simplify(f2)
