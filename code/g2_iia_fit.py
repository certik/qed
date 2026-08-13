"""IIa piece (a): fit the lambda -> 0 limit A0 of the quad-precision
ladder produced by code/g2_iia_a3.f90 (values V(lam) = I(lam)
- (1/2)log lam, 30 digits printed, stable to ~18 digits between DE
levels 5 and 6).

The integrand is rational in lam (the photon mass enters as lam^2 times
rational factors, and the measure brings one odd power), so the tail
basis is integer powers of lam times powers of log lam.  Fits on nested
subsets show convergence of the constant term to the predicted

    A0 = 11/12 + pi^2/18 = 1.46497802228274...
"""
from mpmath import mp, mpf, log, pi, matrix, lu_solve, nstr

mp.dps = 30
data = [(mpf("1e-2"), mpf("1.397257953741492807873377568570")),
        (mpf("3e-3"), mpf("1.439499192379301986687186310182")),
        (mpf("1e-3"), mpf("1.454839988486778227163151036370")),
        (mpf("3e-4"), mpf("1.461381945469017179823820054682")),
        (mpf("1e-4"), mpf("1.463608481859397470290915827937")),
        (mpf("3e-5"), mpf("1.464510667192499352430857425917"))]
pred = mpf(11)/12 + pi**2/18


def basis(l):
    return [mpf(1), l, l*log(l), l*log(l)**2, l**2, l**2*log(l)]


def fit(pts, nb):
    M = matrix(len(pts), nb)
    rhs = matrix(len(pts), 1)
    for i, (l, v) in enumerate(pts):
        row = basis(l)
        for j in range(nb):
            M[i, j] = row[j]
        rhs[i] = v
    return lu_solve(M, rhs)[0]


print("A0 fits over {1, lam, lam log lam, lam log^2 lam, lam^2, "
      "lam^2 log lam}:")
for label, pts, nb in [("smallest 4 points, 4 terms", data[2:], 4),
                       ("smallest 5 points, 5 terms", data[1:], 5),
                       ("all 6 points,     6 terms", data, 6)]:
    c0 = fit(pts, nb)
    print("  %-27s A0 = %s   (A0 - pred = %s)"
          % (label, nstr(c0, 12), nstr(c0 - pred, 3)))
print("predicted A0 = 11/12 + pi^2/18 =", nstr(pred, 15))
