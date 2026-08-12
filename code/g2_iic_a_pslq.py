"""IIc piece (a): high-precision 2-dim quadrature of Ga(t,u) and PSLQ
identification of A* = mu_a(lam=0).

Input: g2_iic_a_Ga.pkl (from g2_iic_a_sint.py).
Fortran reference: A* ~ 0.7697 (g2_iic.f90 a-column, lam -> 0).
"""
import pickle
import time
from multiprocessing import Pool

from sympy import symbols, Symbol, lambdify
from mpmath import (mp, quad as mquad, mpf, pslq, pi as mpi, zeta,
                    log as mlog, nstr)

DPS = 20
NP = 24
T0 = time.time()


def tick(m):
    print("[%7.1fs] %s" % (time.time() - T0, m), flush=True)


def make_f():
    t, u = symbols("t u", positive=True)
    R = Symbol("R", positive=True)
    Ga = pickle.load(open("g2_iic_a_Ga.pkl", "rb"))
    return lambdify((t, u, R), Ga, modules="mpmath", cse=True)


F = None


def init():
    global F
    mp.dps = DPS
    F = make_f()


def safe(tv, uv):
    EPS = mpf("1e-17")
    if tv < EPS or uv < EPS or 1 - uv < EPS or 1 - tv < EPS:
        return mp.zero
    try:
        return F(tv, uv, (tv * uv * (1 - uv)) ** mpf("0.5")).real
    except ZeroDivisionError:
        return mp.zero


def patch(k):
    mp.dps = DPS
    a, b = mpf(k) / NP, mpf(k + 1) / NP
    val = mquad(lambda tv: mquad(lambda uv: safe(tv, uv), [0, 1]), [a, b])
    return str(val)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] != "collect":
        # worker mode: compute one patch, write result
        k = int(sys.argv[1])
        init()
        val = patch(k)
        with open(f"g2_iic_a_patch_{k}.txt", "w") as fh:
            fh.write(val + "\n")
        tick(f"patch {k} done: {val[:30]}...")
        sys.exit(0)
    # collect mode
    mp.dps = DPS
    parts = [open(f"g2_iic_a_patch_{k}.txt").read().strip()
             for k in range(NP)]
    Astar = sum(mpf(p) for p in parts)
    tick("A* = %s" % nstr(Astar, DPS))
    with open("g2_iic_a_Astar.txt", "w") as fh:
        fh.write(str(Astar) + "\n")

    for name, basis in [
        ("[1, pi^2, z3, pi^2 log2]",
         [mpf(1), mpi**2, zeta(3), mpi**2 * mlog(2)]),
        ("extended",
         [mpf(1), mpi**2, zeta(3), mpi**2 * mlog(2), mlog(2)**2,
          mlog(2)**3, mpi**2 * mlog(2)**2, mlog(2), mpi**4]),
    ]:
        rel = pslq([Astar] + basis, maxcoeff=10**5, maxsteps=10**6)
        tick("PSLQ %s: %s" % (name, rel))
        if rel:
            break
