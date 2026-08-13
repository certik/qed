"""Systematic Wick-contraction enumeration for the QED vertex function.

The object expanded is the 3-point Green function (electron in an external
field, cf. 09-qed-g-2-LO.md and 10-qed-g-2-NLO-derivation.md)

    G^mu = < Omega | T  psi(p')  A^mu(-q)  psibar(-p) | Omega >,

whose order-e^n term (Gell-Mann-Low + Dyson series) is

    ((-i e)^n / n!) int d^4x_1 ... d^4x_n
    < 0 | T  psi(p') A^mu(-q) psibar(-p)
             (psibar A-slash psi)(x_1) ... (psibar A-slash psi)(x_n) | 0 >.

By Wick's theorem this vacuum expectation value is the sum over ALL
complete contractions: every psi paired with a psibar (fermion propagator
S), every A paired with another A (photon propagator D).  This script
enumerates every contraction exactly:

  * photon side: perfect matchings of {A_ext, A_1, ..., A_n},
  * fermion side: bijections  {psi_out, psi_1..psi_n} -> {psibar_in,
    psibar_1..psibar_n}  (psi_a paired with psibar_b gives S(x_a - x_b)),

computes the fermion anticommutation sign of each term (by bringing the
contracted pairs adjacent in the canonical operator order and counting
crossings), classifies the resulting topology, and tabulates the counts.

Checks made explicit in the output:

  * the count of every distinct topology is a multiple of n!
    (relabelings of the integration variables x_i), so after the 1/n! of
    the Dyson series every distinct diagram enters with integer weight --
    weight 1 per orientation/mirror copy, no fractional symmetry factors;
  * the sign of every term equals (-1)^(number of closed fermion loops)
    (asserted for all terms), the textbook loop-sign rule;
  * at n = 5 the connected, loop-free part is exactly the 15 open-line
    topologies (I, IIa, IIb/IIf, IIc x2, IId x2, external-leg-only), and
    the loop classes reproduce Karplus-Kroll's IIe, III, IV, V.

Run:  pixi run python code/wick_enum.py
"""
from itertools import permutations
from math import factorial


# ----------------------------------------------------------------- helpers
def pairings(items):
    """All perfect matchings of a list."""
    if not items:
        yield []
        return
    a, rest = items[0], items[1:]
    for j, b in enumerate(rest):
        for tail in pairings(rest[:j] + rest[j + 1:]):
            yield [(a, b)] + tail


def wick_sign(n, sigma):
    """Fermion sign of the contraction sigma (psi-slot -> psibar-slot).

    Canonical operator order of the fermion fields in the T-product:

        psi_out  psibar_in  (psibar_1 psi_1) (psibar_2 psi_2) ...

    A contracted pair encountered in the order (psi ... psibar) counts as
    +S; in the order (psibar ... psi) as -S (T-reordering of the pair).
    Bringing the pairs adjacent costs one factor (-1) per crossed
    still-uncontracted field.  Returns the total sign relative to the
    product of propagators S.
    """
    tokens = [("psi", "out"), ("psibar", "in")]
    for i in range(1, n + 1):
        tokens += [("psibar", i), ("psi", i)]
    pos_psi = {lab: k for k, (t, lab) in enumerate(tokens) if t == "psi"}
    pos_bar = {lab: k for k, (t, lab) in enumerate(tokens) if t == "psibar"}
    partner = {}
    for a, b in sigma.items():
        partner[pos_psi[a]] = pos_bar[b]
        partner[pos_bar[b]] = pos_psi[a]
    alive = [True] * len(tokens)
    sign = 1
    for L in range(len(tokens)):
        if not alive[L]:
            continue
        P = partner[L]
        between = sum(1 for k in range(L + 1, P) if alive[k])
        sign *= (-1) ** between
        if tokens[L][0] == "psibar":       # pair met in (psibar ... psi) order
            sign *= -1
        alive[L] = alive[P] = False
    return sign


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


# ----------------------------------------------------------- classification
def classify(n, phot, sigma):
    """Return (class name, number of closed fermion loops)."""
    # fermion chain from the outgoing electron down to the incoming one
    chain = []
    cur = sigma["out"]
    while cur != "in":
        chain.append(cur)
        cur = sigma[cur]
    online = set(chain)
    loops, seen = [], set(online)
    for i in range(1, n + 1):
        if i in seen:
            continue
        cyc, j = [i], sigma[i]
        seen.add(i)
        while j != i:
            cyc.append(j)
            seen.add(j)
            j = sigma[j]
        loops.append(cyc)
    L = len(loops)

    if any(len(c) == 1 for c in loops):
        return "tadpole loop  tr[gamma S(0)] = 0", L

    # connectivity (nodes: LINE, EXTA, vertices)
    uf = UF()
    uf.find("LINE"), uf.find("EXTA")
    for vtx in online:
        uf.union("LINE", vtx)
    for cyc in loops:
        for a2 in cyc:
            uf.union(cyc[0], a2)
    extv = None
    for a, b in phot:
        if a == "ext":
            extv, other = b, b
            uf.union("EXTA", b)
        elif b == "ext":
            extv = a
            uf.union("EXTA", a)
        else:
            uf.union(a, b)
    root = uf.find("LINE")
    disconnected = any(uf.find(v) != root for v in range(1, n + 1)) \
        or uf.find("EXTA") != root

    if not chain:
        return "forward line x (rest): no scattering / disconnected", L
    if disconnected:
        return "disconnected (vacuum bubble or detached blob)", L

    line = list(reversed(chain))          # incoming-side first
    posn = {v: k + 1 for k, v in enumerate(line)}
    ipos = posn.get(extv)                 # None if ext photon on a loop
    ppairs = []
    for a, b in phot:
        if "ext" in (a, b):
            continue
        ppairs.append((a, b))

    if L == 0:
        idx = [tuple(sorted((posn[a], posn[b]))) for (a, b) in ppairs]
        if not idx:
            return "tree vertex gamma^mu", L
        spans = [pr for pr in idx if pr[0] < ipos < pr[1]]
        if len(spans) == 2:
            (a1, b1), (a2, b2) = sorted(spans)
            if a1 < a2 < b1 < b2:
                return "I: crossed ladder", L
            return "IIa: ladder (vertex part at the external vertex)", L
        if len(spans) == 1:
            (a, b) = spans[0]
            if len(idx) == 1:
                return "LO one-loop vertex correction", L
            (c, d) = [pr for pr in idx if pr != (a, b)][0]
            if c > b or d < a:
                return "IIb/IIf: LO vertex x self-energy on an external leg", L
            if a < c and d < b:
                return "IId: self-energy on an internal line (+ mirror)", L
            return "IIc: corner, vertex part at an internal vertex (+ mirror)", L
        return "external-leg corrections only (tree vertex)", L

    sizes = sorted(len(c) for c in loops)
    if any(s % 2 == 1 for s in sizes):
        return "V: odd fermion loop (Furry pair, cancels)", L
    if sizes == [2]:
        lp = loops[0]
        if extv in lp:
            # photon from the loop to the line at line position j;
            # the line's own photon pair (c,d)
            j = None
            for a, b in ppairs:
                if a in lp and b in online:
                    j = posn[b]
                elif b in lp and a in online:
                    j = posn[a]
            rest = [tuple(sorted((posn[a], posn[b])))
                    for (a, b) in ppairs if not (a in lp or b in lp)]
            if not rest:
                return "III: VP on the external-potential line x tree vertex", L
            (c, d) = rest[0]
            if c < j < d:
                return "III: VP on the external-potential line x vertex corr.", L
            return "III: VP on the external-potential line x leg SE", L
        att = sorted(posn[a] if a in online else posn[b]
                     for (a, b) in ppairs if (a in lp) != (b in lp))
        if att[0] < ipos < att[1]:
            return "IIe: vacuum polarization in the internal photon", L
        return "external-leg SE with VP insertion", L
    if sizes == [4]:
        return "III: two-loop VP blob on the external-potential line", L
    if sizes == [2, 2]:
        return "IV: iterated one-loop VP chain on the external-potential line", L
    return f"unclassified loops {sizes}", L


# ------------------------------------------------------------------ driver
def enumerate_order(n):
    V = list(range(1, n + 1))
    table = {}
    total = 0
    for phot in pairings(["ext"] + V):
        for perm in permutations(["in"] + V):
            sigma = dict(zip(["out"] + V, perm))
            total += 1
            cls, L = classify(n, phot, sigma)
            s = wick_sign(n, sigma)
            assert s == (-1) ** L, (cls, sigma, s, L)
            cnt, sgns = table.get(cls, (0, set()))
            table[cls] = (cnt + 1, sgns | {s})
    nm = factorial(n)
    nmatch = 1
    for j in range(n, 0, -2):
        nmatch *= j
    assert total == nmatch * factorial(n + 1)
    print(f"order e^{n}:  {total} complete contractions "
          f"({nmatch} photon matchings x {n + 1}! fermion pairings)")
    print(f"  {'topology class':66s} {'count':>6s} {'/n!':>4s} sign")
    for cls in sorted(table, key=lambda c: (-table[c][0], c)):
        cnt, sgns = table[cls]
        # every single-topology class must come in multiples of n!
        # (buckets like "tadpole"/"disconnected" merge several structures)
        bucket = cnt % nm != 0
        assert not bucket or cls.startswith(("tadpole", "disconnected",
                                             "forward")), (cls, cnt)
        wtxt = f"{cnt // nm:4d}" if not bucket else "   -"
        stxt = "/".join(f"{s:+d}" for s in sorted(sgns, reverse=True))
        print(f"  {cls:66s} {cnt:6d} {wtxt} {stxt}")
    print(f"  {'total':66s} {total:6d} {total // nm:4d}")
    print("  all signs equal (-1)^(# closed fermion loops): OK")
    return table


if __name__ == "__main__":
    for n in (1, 3, 5):
        enumerate_order(n)
        print()
