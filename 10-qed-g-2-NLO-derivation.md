# Magnetic moment of an electron (NLO): ab-initio derivation

This section derives, from the QED Lagrangian and with no steps assumed,
every integral that the NLO section (the next section, "Magnetic moment of
an electron (NLO)") evaluates. Nothing is *computed* here: the end
products are the exact loop/parameter integrals — the LO massive-photon
kernel, the five two-loop diagrams I, IIa, IIc, IId, IIe with their
subtractions, and the projector that turns each of them into a number —
together with a proof that their sum is the physical $A_2$ and an
accounting of every other term the expansion generates (they all cancel,
vanish, or renormalize parameters).

The logical chain, each link derived below:

1. QED Lagrangian with regulators (photon mass $\lambda$, Pauli–Villars)
   → free propagators and the interaction vertex.
2. The quantity: $a_e = F_2(0)$, with the form-factor decomposition and
   the explicit Breit-frame projector proved.
3. The 3-point Green function and the LSZ reduction *with interacting
   fields*: residue factors $Z_2$, amputation by full propagators,
   why external-leg diagrams (IIb/IIf) drop out.
4. Gell-Mann–Low formula (proved in the QFT section) + Dyson series +
   **Wick's theorem** (proved here, with the fermion sign rules)
   → the completely systematic algebraic expansion.
5. The expansion executed at orders $e^1$, $e^3$, $e^5$: every single
   contraction enumerated and classified (mechanically verified by
   `code/wick_enum.py`).
6. The loop-integration toolkit derived: Feynman parameters, the shift,
   Wick rotation, the master integral table of `code/loops.py`, angular
   averages, Pauli–Villars.
7. The one-loop building blocks derived: the massive-photon vertex
   (kernel $K$), the self-energy $\Sigma$, the vacuum polarization
   $\hat\Pi$ and its spectral function $\rho(t)$, the one-loop vertex
   parts $\Lambda^\nu$ with their Feynman-parameter denominators.
8. The five raw two-loop integrals written out explicitly.
9. Renormalization: the on-shell counterterms, the Ward identity
   $Z_1 = Z_2$, and the **assembly theorem**: the
   Karplus–Kroll-subtracted diagrams of the NLO section sum exactly to
   the ab-initio $a_e\big|_{\alpha^2}$.

A note on terminology, since the question arises naturally: the
"systematic algebraic expansion" has two distinct steps. Expanding the
time-ordered exponential $S = T\exp(i\int\mathcal{L}_1)$ in powers of $e$
is the **Dyson series** (an ordinary Taylor expansion under the
$T$-symbol). Converting each resulting vacuum expectation value of a
$T$-product of free fields into a sum of products of propagators is
**Wick's theorem**. Feynman diagrams are nothing but the bookkeeping of
the second step: one diagram per class of complete contractions. Both
steps are carried out explicitly below.

## 1. The theory

We work in the conventions of the QED and LO sections: metric
$(+,-,-,-)$, $\hbar = c = 1$, Dirac representation for $\gamma^\mu$,
electron charge $e$ (only $e^2 = 4\pi\alpha$ ever appears). The
Lagrangian density, including the two regulators used by
Karplus–Kroll/Petermann and by all our computations, is

$$\mathcal{L} = \bar\psi\left(i\gamma^\mu\partial_\mu - m_0\right)\psi
  - \frac14 F_{\mu\nu}F^{\mu\nu}
  + \frac{\lambda^2}{2}A_\mu A^\mu
  - \frac12\left(\partial_\mu A^\mu\right)^2
  - e\,\bar\psi\gamma^\mu\psi\, A_\mu.$$

* The last term is the interaction $\mathcal{L}_1 = -e\bar\psi\slashed{A}\psi$
  that follows from gauge coupling (QED section); in the Dyson series
  each factor $i\mathcal{L}_1$ produces the vertex factor
  $-ie\gamma^\mu$. Strictly the Lagrangian coupling is the *bare*
  charge $e_0$; the distinction ($e_0^2 = e^2\left[1-\Pi(0)\right]$,
  §3.3) matters only in the charge-renormalization bookkeeping of
  §3.3/§9, so we write $e$ throughout and reinstate $e_0$ exactly
  there.
* $\lambda$ is the **photon mass**, the infrared regulator. Together
  with the Feynman gauge-fixing term $-\frac12(\partial A)^2$ it gives
  the photon propagator its simple $-ig_{\mu\nu}$ numerator (derived
  below). This *defines* the IR regularization scheme of the 1950s
  papers: photon propagator $-ig_{\mu\nu}/(k^2-\lambda^2)$ everywhere,
  $\lambda \to 0$ at the end. (A genuine Proca photon would carry an
  extra $k_\mu k_\nu/\lambda^2$ term; by the collapse lemma of §2.4 such
  longitudinal pieces only produce $\gamma^\mu$-structures, i.e. they
  change $F_1$-type contact terms, and the $\lambda\to 0$ limit of the
  total $F_2$ is prescription independent. All per-diagram values in
  the NLO section refer to the $-ig_{\mu\nu}/(k^2-\lambda^2)$
  prescription.)
* $m_0 = m - \delta m$ is the bare mass. We expand around the
  *physical* mass $m$:

  $$\mathcal{L} \supset -m_0\bar\psi\psi = -m\bar\psi\psi
     + \delta m\,\bar\psi\psi,$$

  so propagators carry $m$ and the counterterm generates an extra
  two-point vertex $+i\,\delta m$ (one insertion is order $\alpha$,
  exactly what is needed at two loops). The on-shell condition fixing
  $\delta m$ is derived in §4.2.
* Ultraviolet regularization is **Pauli–Villars**, implemented at the
  level of the single logarithmically divergent master integral (§6.5):
  every log divergence is carried as the explicit symbol
  $L_{UV} = \log\Lambda^2 + \text{const}$, and every $F_2$ below must be
  (and is, assertedly) free of $L_{UV}$. For the fermion loop of the
  vacuum polarization the PV subtraction must be gauge invariant
  (two-mass PV); §7.3 explains where this matters and why the final
  answer is unambiguous.

### 1.1 Free propagators

The fermion propagator was derived in the QFT section:

$$\langle 0|T\,\psi(x)\bar\psi(y)|0\rangle = S(x-y)
  = \int\frac{\mathrm{d}^4p}{(2\pi)^4}\,\tilde S(p)\,e^{-ip(x-y)},
  \qquad
  \tilde S(p) = \frac{i\left(\slashed{p}+m\right)}{p^2-m^2+i\epsilon}.$$

For the photon we need the propagator that follows from *our* quadratic
Lagrangian (mass term plus Feynman gauge fixing). Integrating by parts
(dropping total derivatives, which do not change the action),

$$-\frac14F_{\mu\nu}F^{\mu\nu}
 = -\frac12\left(\partial_\mu A_\nu\partial^\mu A^\nu
              - \partial_\mu A_\nu\partial^\nu A^\mu\right)
 = \frac12 A_\nu\,\partial^2 A^\nu
   - \frac12 A_\nu\,\partial^\nu\partial_\mu A^\mu ,$$

$$-\frac12\left(\partial_\mu A^\mu\right)^2
 = \frac12 A_\nu\,\partial^\nu\partial_\mu A^\mu ,$$

so the gauge-fixing term cancels the second piece exactly and

$$\mathcal{L}_A = \frac12 A_\mu\left[g^{\mu\nu}\left(\partial^2
   + \lambda^2\right)\right]A_\nu .$$

Each of the four components of $A_\mu$ is a free scalar of mass
$\lambda$ (with a sign from the metric for the spatial ones), so the
propagator is the inverse of the operator in brackets with the causal
$i\epsilon$, i.e. in momentum space $-g^{\mu\nu}(k^2-\lambda^2)$
inverted times $i$:

$$\tilde D_{\mu\nu}(k) = \frac{-i g_{\mu\nu}}{k^2 - \lambda^2 + i\epsilon}.$$

This is the massless propagator of the QFT section with
$k^2 \to k^2-\lambda^2$, as promised.

## 2. The quantity: $a_e = F_2(0)$

### 2.1 Setup

As in the LO section, the electron scatters off a fixed classical
potential $A_\mu^{\mathrm{cl}}$; to first order in the potential and to
all orders in the radiative corrections the amplitude is

$$iM\,2\pi\,\delta(p'^0-p^0)
  = -ie\,\bar u(p')\,\Gamma^\mu(p',p)\,u(p)\,A^{\mathrm{cl}}_\mu(q),
  \qquad q = p'-p,$$

and the anomalous moment is $a_e = F_2(0)$ where $F_2$ is defined by the
decomposition of $\Gamma^\mu$ derived next. The connection between
$F_{1,2}(0)$ and the $g$-factor ($g = 2(F_1+F_2)$, Born-level matching
to $-\boldsymbol{\mu}\cdot\mathbf{B}$) was derived in full in the LO
section and is not repeated.

### 2.2 The form-factor decomposition, derived

**Claim.** For on-shell spinors ($\slashed{p}u(p) = m\,u(p)$,
$\bar u(p')\slashed{p}' = m\,\bar u(p')$, $p^2 = p'^2 = m^2$) and a
parity-conserving theory, the most general vertex satisfying current
conservation is

$$\bar u(p')\,\Gamma^\mu(p',p)\,u(p)
 = \bar u(p')\left[\gamma^\mu F_1(q^2)
 + \frac{i\sigma^{\mu\nu}q_\nu}{2m}F_2(q^2)\right]u(p),
 \qquad \sigma^{\mu\nu} = \frac{i}{2}[\gamma^\mu,\gamma^\nu].$$

*Step 1 — basis.* $\Gamma^\mu$ is a $4\times4$ matrix carrying one
Lorentz index, built from $\gamma^\mu$, $p$, $p'$, $g$, and (by parity
conservation, no $\gamma^5$ or $\epsilon$-tensor) nothing else. Any
product of $\slashed{p}$'s and $\slashed{p}'$'s sandwiched between the
on-shell spinors collapses: commute $\slashed{p}$ to the right
($\slashed{p}\to m$) and $\slashed{p}'$ to the left
($\slashed{p}'\to m$) using
$\slashed{a}\slashed{b} = 2a\cdot b - \slashed{b}\slashed{a}$; each
commutation only produces scalars $p\cdot p' = m^2 - q^2/2$. The
sandwiched vertex therefore reduces to

$$\bar u'\,\Gamma^\mu u = \bar u'\left[A\,\gamma^\mu + B\,P^\mu
  + C\,q^\mu\right]u,
  \qquad P^\mu \equiv (p+p')^\mu,$$

with scalar functions $A,B,C$ of $q^2$ (the only invariant). Structures
$\sigma^{\mu\nu}q_\nu$ and $\sigma^{\mu\nu}P_\nu$ are not extra: writing
$i\sigma^{\mu\nu}a_\nu = -\gamma^\mu\slashed{a}+a^\mu
= \slashed{a}\gamma^\mu - a^\mu$
(direct consequence of $\sigma^{\mu\nu} = i(\gamma^\mu\gamma^\nu -
g^{\mu\nu})$),

$$\bar u'\,i\sigma^{\mu\nu}q_\nu\,u
  = \bar u'\left[2m\gamma^\mu - P^\mu\right]u
  \qquad\text{(Gordon identity)},$$

$$\bar u'\,i\sigma^{\mu\nu}P_\nu\,u
  = \bar u'\left[P^\mu - 2p'^\mu\right]u = -\,q^\mu\,\bar u' u .$$

(The first line: $\bar u'(-\gamma^\mu\slashed{q}+q^\mu)u$ with
$\slashed{q}u = (\slashed{p}'-\slashed{p})u = \slashed{p}'u - mu$ and
$\gamma^\mu\slashed{p}' = 2p'^\mu - \slashed{p}'\gamma^\mu$, then
$\bar u'\slashed{p}' = m\bar u'$; both identities are verified
numerically below.)

*Step 2 — current conservation.* The Ward–Takahashi identity (derived
ab initio in §4.3) gives $q_\mu\,\bar u'\Gamma^\mu u = 0$ on shell.
Now $q\cdot P = p'^2 - p^2 = 0$ and
$\bar u'\slashed{q}u = (m-m)\bar u'u = 0$ hold identically, so
transversality forces $q^2\,C(q^2) = 0$, i.e. $C \equiv 0$.

*Step 3 — repackaging.* Use the Gordon identity to trade $P^\mu$ for
$\gamma^\mu$ and $i\sigma^{\mu\nu}q_\nu$:

$$A\gamma^\mu + BP^\mu
 = (A + 2mB)\,\gamma^\mu - B\,i\sigma^{\mu\nu}q_\nu
 \equiv \gamma^\mu F_1 + \frac{i\sigma^{\mu\nu}q_\nu}{2m}F_2 . \qquad\blacksquare$$

### 2.3 The explicit projector (as used by the pipeline)

All our computations extract $F_2$ by sandwiching the *matrix*
$\Gamma^\mu$ between explicit spinors in the Breit frame
(`code/dirac.py`),

$$p = (E,0,0,w),\quad p' = (E,0,0,-w),\quad q = (0,0,0,-2w),
\quad E = \sqrt{m^2+w^2},\quad q^2 = -4w^2 .$$

Dirac-representation spinors ($\xi_0 = \binom10$, $\xi_1 = \binom01$,
$\sigma^3\xi_s = \pm\xi_s$):

$$u^s(p) = \begin{pmatrix}\sqrt{E+m}\,\xi_s\\[2pt]
  \dfrac{\mathbf{p}\cdot\boldsymbol{\sigma}}{\sqrt{E+m}}\,\xi_s\end{pmatrix},
  \qquad \bar u = u^\dagger\gamma^0 .$$

Define the two sandwiches used everywhere ($A_\mu$ here is a number,
not a field):

$$A_0 \equiv \bar u^{0}(p')\,\Gamma^0\,u^{0}(p),\qquad
  A_1 \equiv \bar u^{0}(p')\,\Gamma^1\,u^{1}(p).$$

Evaluate the two basis structures on these sandwiches. With
$\mathbf{p}'\cdot\boldsymbol{\sigma} = -w\sigma^3$,
$\mathbf{p}\cdot\boldsymbol{\sigma} = w\sigma^3$:

* $(\mu,s',s)=(0,0,0)$:
  $\bar u'\gamma^0 u = u'^\dagger u
  = (E+m)\,\xi_0^\dagger\xi_0
  + \frac{(-w\sigma^3)(w\sigma^3)}{E+m}
  = (E+m) - \frac{w^2}{E+m} = 2m$.
  For the $\sigma$-structure only $\nu=3$ contributes
  ($q_3 = -q^3 = 2w$), and
  $\frac{i\sigma^{03}q_3}{2m} = \frac{i(i\gamma^0\gamma^3)(2w)}{2m}
  = -\frac{w}{m}\gamma^0\gamma^3$; using
  $u'^\dagger\gamma^3u = 2w$ (same block algebra) the value is
  $-\frac{2w^2}{m}$.
* $(\mu,s',s)=(1,0,1)$: with
  $\gamma^0\gamma^1 = \begin{pmatrix}0&\sigma^1\\ \sigma^1&0\end{pmatrix}$
  and $\sigma^1\xi_1 = \xi_0$:
  $\bar u'\gamma^1 u = u'^\dagger\gamma^0\gamma^1 u
  = -2w\,\xi_0^\dagger\sigma^1\xi_1 = -2w$; and
  $\frac{i\sigma^{13}q_3}{2m} = -\frac{w}{m}\gamma^1\gamma^3$ with
  $\gamma^1\gamma^3 = i\sigma^2\,\mathbb{1}_{2\times2\,\text{blocks}}$,
  $i\sigma^2\xi_1 = \xi_0$, gives
  $-\frac{w}{m}\left[(E+m)-\frac{w^2}{E+m}\right] = -2w$.

So on these two sandwiches

$$A_0 = 2m\,F_1 - \frac{2w^2}{m}F_2,\qquad
  A_1 = -2w\,(F_1+F_2),$$

— the $\mu=1$ flip sandwich alone measures only the *magnetic*
combination $F_1+F_2$, which is why it must be paired with the $\mu=0$
sandwich. Solving the $2\times2$ linear system:

$$\boxed{\;F_1 = \frac{m A_0 - w A_1}{2\left(m^2+w^2\right)},\qquad
  F_2 = -\,\frac{m\left(w A_0 + m A_1\right)}{2w\left(m^2+w^2\right)}\;}$$

exact at finite $w$. As $w\to0$ (i.e. $q^2\to0$), $A_1$ vanishes
linearly in $w$ (it is proportional to $w(F_1+F_2)$), so with $m=1$

$$F_2(0) = -\tfrac12\left(A_0\big|_{w=0} + \partial_w A_1\big|_{w=0}\right),$$

which is precisely the projection formula quoted and used in the NLO
section (`code/g2_iid.py` validates it on the LO diagram). Everything
above is verified symbolically:

    >>> import sys; sys.path.insert(0, 'code')
    >>> from sympy import symbols, simplify
    >>> from dirac import formfactor_basis
    >>> m, w = symbols("m w", positive=True)
    >>> [simplify(b) for b in formfactor_basis(m, w, 0, 0, 0)]
    [2*m, -2*w**2/m]
    >>> [simplify(b) for b in formfactor_basis(m, w, 1, 0, 1)]
    [-2*w, -2*w]

### 2.4 The collapse lemma (longitudinal photon pieces)

For later use (gauge/regulator independence, Ward manipulations): for
any internal-photon momentum $k$ attached to the electron line next to
an on-shell end,

$$\frac{1}{\slashed{p}-\slashed{k}-m}\,\slashed{k}\,u(p)
  = \frac{1}{\slashed{p}-\slashed{k}-m}
    \left[-(\slashed{p}-\slashed{k}-m) + (\slashed{p}-m)\right]u(p)
  = -\,u(p),$$

$$\bar u(p')\,\slashed{k}\,\frac{1}{\slashed{p}'-\slashed{k}-m}
  = -\,\bar u(p').$$

Iterating such telescoping steps, a $k_\mu k_\nu$ numerator on an
internal photon collapses the two adjacent propagator chains and leaves
only $\gamma^\mu$-type (i.e. $F_1$-type) structures. Both identities
are checked numerically in the pipeline conventions (they are two-line
consequences of $(\slashed{p}-m)u(p)=0$).

## 3. Green functions, LSZ with interacting fields, and amputation

### 3.1 The 3-point function

Exactly as in the LO section, the amplitude is obtained from

$$\tilde G^\mu(p',q,p) = \langle\Omega|\,T\,
   \tilde\psi(p')\,\tilde A^\mu(-q)\,\tilde{\bar\psi}(-p)\,|\Omega\rangle,$$

with the external legs stripped by the reduction rules of the QFT
section. Those rules were derived there with *free*-field residues; at
two loops we need the refinement that the interacting field creates the
physical one-electron state with amplitude

$$\langle\Omega|\psi(0)|p,s\rangle = \sqrt{Z_2}\;u^s(p).$$

$Z_2$ is defined by this equation. Repeating the free-field computation
of the propagator with this matrix element (insert a complete set of
states into $\langle T\psi\bar\psi\rangle$; the one-particle
contribution reproduces the free two-point integral times $Z_2$, with
the *physical* mass $m$ in the exponentials; multi-particle states start
at the cut $p^2 \ge (m+\lambda)^2$) gives the exact two-point function
near its pole:

$$\tilde S'(p) \equiv \int \mathrm{d}^4x\, e^{ipx}
  \langle\Omega|T\psi(x)\bar\psi(0)|\Omega\rangle
  \;\xrightarrow{\;\slashed{p}\to m\;}\;
  \frac{iZ_2\left(\slashed{p}+m\right)}{p^2-m^2+i\epsilon}
  + \text{regular}.$$

Consequently the correct reduction rule for each external fermion leg is
the QFT-section rule divided by $\sqrt{Z_2}$ (the free rule assumed unit
amplitude), with $m$ the physical mass.

### 3.2 Two-point functions from 1PI blocks

Let $-i\Sigma(\slashed{p})$ denote the sum of all one-particle-
irreducible (1PI) self-energy insertions (plus the $+i\delta m$
counterterm vertex, so $\Sigma_{\rm eff} = \Sigma - \delta m$; we write
$\Sigma$ for $\Sigma_{\rm eff}$ from here on). Summing the geometric
series of blobs strung on free propagators,

$$\tilde S' = \tilde S + \tilde S(-i\Sigma)\tilde S + \cdots
 = \frac{i}{\slashed{p}-m-\Sigma(\slashed{p})+i\epsilon}.$$

Matching to the pole form above:

$$\Sigma(\slashed{p}=m) = 0
 \;\;\Leftrightarrow\;\; \delta m = \Sigma_{\rm loop}(m)
 \qquad\text{(on-shell mass)},$$

$$Z_2 = \frac{1}{1-\Sigma'(m)} = 1 + \Sigma'(m) + O(\alpha^2),
 \qquad \delta Z_2 \equiv \Sigma'(m) = 
 \left.\frac{\partial\Sigma}{\partial\slashed{p}}\right|_{\slashed{p}=m}.$$

Identically for the photon: with $i\Pi^{\mu\nu}(q) =
(q^2g^{\mu\nu}-q^\mu q^\nu)\,i\Pi(q^2)$ the 1PI photon self-energy
(transversality is derived in §7.3), the dressed propagator sums to

$$\tilde D'_{\mu\nu}(q) = \frac{-ig_{\mu\nu}}{q^2\left(1-\Pi(q^2)\right)}
  + \left(q_\mu q_\nu\ \text{terms}\right)$$

(displayed at $\lambda=0$; the $q_\mu q_\nu$ terms never contribute by
§2.4 and the Ward identity).

### 3.3 Factorized form of $\tilde G^\mu$ and the S-matrix vertex

Every connected diagram of $\tilde G^\mu$ decomposes uniquely as
(blobs on the $p$ leg) ⊗ (blobs on the $p'$ leg) ⊗ (vacuum-polarization
chain on the photon leg) ⊗ (core vertex 1PI in all three legs). This is
not an assumption: at the orders we need it is *verified term by term*
by the complete enumeration of §5 — every contraction lands in exactly
one class (the odd-loop class V would be additional core diagrams, but
it cancels pairwise by Furry's theorem, §5). Therefore, exactly,

$$\tilde G^\mu(p',q,p) = (2\pi)^4\delta^{(4)}(p'-p-q)\;
  \tilde S'(p')\,\left[-ie_0\,\Gamma^\nu(p',p)\right]\,\tilde S'(p)\;
  \frac{{-ig_\nu{}^\mu}}{q^2\left(1-\Pi(q^2)\right)},$$

where $\Gamma^\nu$ is the **proper (1PI) vertex** — at tree level
$\gamma^\nu$; at one loop the LO diagram; at two loops the five diagrams
of §8.

Now apply the corrected reduction of §3.1. Each fermion leg contributes

$$\frac{1}{\sqrt{Z_2}}\;\bar u(p')\,\frac{1}{\tilde S(p')}\,\tilde S'(p')
 \;\xrightarrow{\;\text{on shell}\;}\; \frac{1}{\sqrt{Z_2}}\,
 \bar u(p')\,\frac{\slashed{p}'-m}{i}\,
 \frac{iZ_2(\slashed{p}'+m)}{p'^2-m^2} = \sqrt{Z_2}\;\bar u(p'),$$

(using $(\slashed{p}'-m)(\slashed{p}'+m) = p'^2-m^2$; the regular part
of $\tilde S'$ is killed by $(\slashed{p}'-m)u = 0$), and similarly
$\sqrt{Z_2}\,u(p)$ on the incoming side. The photon leg is the external
potential; its propagator and $1/(1-\Pi)$ dressing combine with the
charge as follows. On-shell **charge renormalization** defines

$$e^2 \equiv \frac{e_0^2}{1-\Pi(0)},\qquad
  \hat\Pi(q^2) \equiv \Pi(q^2) - \Pi(0)\ \ (\hat\Pi(0) = 0),$$

so that
$e_0^2/(1-\Pi(q^2)) = e^2\left[1+\hat\Pi(q^2)+O(\alpha^2)\right]$.
Collecting everything:

$$iM = -ie\,\bar u(p')\left[\,Z_2\,\Gamma^\mu(p',p)\,\right]u(p)\,
  A^{\mathrm{cl}}_\mu(q)\;\times\;\left[1+\hat\Pi(q^2)+\dots\right].$$

Three structural conclusions, each an exact statement about classes of
diagrams:

1. **External-leg self-energies (KK classes IIb, IIf and all their
   two-loop variants) never enter the vertex.** They are part of
   $\tilde S'(p^{(\prime)})$ and are consumed by the reduction, which
   leaves only the number $Z_2$ per leg — total factor $Z_2$.
2. **Vacuum-polarization chains on the external-potential line (KK
   classes III and IV)** are the factor $1/(1-\Pi(q^2))$. After charge
   renormalization they contribute $F_i^{\rm blob}(q^2)\,\hat\Pi(q^2)$;
   since $\hat\Pi(0) = 0$ and all $F_i^{\rm blob}$ are finite at
   $q^2 = 0$, they contribute **nothing** to $a_e = F_2(0)$: pure
   charge renormalization.
3. The physical form factors are those of $Z_2\,\Gamma^\mu$:

   $$a_e = F_2^{S}(0) = Z_2\;F_2(0),\qquad
     F_1^{S}(0) = Z_2\,F_1(0).$$

Expanding to order $\alpha^2$ with
$F_2(0) = F_2^{(2)}(0) + F_2^{(4)}(0) + \dots$ (superscripts count
powers of $e$) and $Z_2 = 1+\delta Z_2$:

$$\boxed{\;a_e\big|_{\alpha^2}
   = F_2^{(4)}(0) \;+\; \delta Z_2\,F_2^{(2)}(0),\qquad
   F_2^{(2)}(0) = \frac{\alpha}{2\pi}.\;}$$

This is the **master ab-initio formula**. By the Ward identity (§4.3)
$\delta Z_2 = -\delta F_1(0)$, so the second term is exactly the
"vertex renormalization cross term
$-\delta F_1(0)\cdot\frac{\alpha}{2\pi}$" of the NLO section.

## 4. The expansion machinery

### 4.1 Gell-Mann–Low and the Dyson series

From the QFT section (proved there),

$$\langle\Omega|T\,\phi_H\cdots|\Omega\rangle
 = \frac{\langle0|T\,\phi_I\cdots\,S\,|0\rangle}{\langle0|S|0\rangle},
 \qquad
 S = T\exp\left(i\int \mathrm{d}^4x\,\mathcal{L}_1(x)\right),$$

with here $\mathcal{L}_1 = -e\bar\psi\slashed{A}\psi
+ \delta m\,\bar\psi\psi$ (all fields interaction-picture free fields).
Taylor-expanding the exponential under the $T$ symbol (legitimate
because everything inside one $T$-product commutes by construction)
gives the order-$e^n$ term of the 3-point function:

$$\tilde G^\mu\Big|_{e^n} = \frac{(-ie)^n}{n!}
  \int \mathrm{d}^4x_1\cdots\mathrm{d}^4x_n\,
  \langle0|T\,\tilde\psi(p')\tilde A^\mu(-q)\tilde{\bar\psi}(-p)\,
  \left(\bar\psi\slashed{A}\psi\right)(x_1)\cdots
  \left(\bar\psi\slashed{A}\psi\right)(x_n)|0\rangle$$

plus, at each order in $\alpha$, the terms with $\delta m$ insertions.
Only odd $n$ survives (the photon fields must pair up: $n$ vertex
photons + 1 external photon must be even). The vacuum diagrams generated
simultaneously in $\langle0|S|0\rangle$ cancel: any complete contraction
factorizes uniquely into its externals-connected part (say $n_1$
vertices) and its vacuum part ($n_2 = n-n_1$ vertices); summing over
which labeled vertices go where turns $1/n!\binom{n}{n_1}$ into
$\frac{1}{n_1!}\frac{1}{n_2!}$, so the numerator factorizes as
(no-vacuum series) × $\langle0|S|0\rangle$ and the denominator removes
the second factor exactly.

### 4.2 Wick's theorem (with fermion signs)

**Definitions.** For a free field $\phi = \phi^+ + \phi^-$
($\phi^+$ annihilates, $\phi^-$ creates), the *normal product*
$N\{\cdots\}$ puts all $\phi^-$ left of all $\phi^+$, with a factor
$(-1)$ for each transposition of two fermionic fields performed in the
rearrangement. The *contraction* of two fields is the c-number

$$\underbrace{\phi(x)\phi(y)} \equiv T\{\phi(x)\phi(y)\} - N\{\phi(x)\phi(y)\}.$$

Because $[\phi^+(x),\phi^-(y)]_{\mp}$ is a c-number for free fields,
the contraction is a c-number, and taking $\langle0|\cdot|0\rangle$
(normal products of one or more fields have zero vacuum expectation
value) shows it equals the propagator:
$\underbrace{\psi_a(x)\bar\psi_b(y)} = S_{ab}(x-y)$,
$\underbrace{A_\mu(x)A_\nu(y)} = D_{\mu\nu}(x-y)$, and
$\underbrace{\psi\psi} = \underbrace{\bar\psi\bar\psi} = 0$.

**Theorem (Wick).**

$$T\{\phi_1\phi_2\cdots\phi_n\} = N\Big\{\phi_1\cdots\phi_n
 \;+\; \sum_{\text{all contractions}}(\pm)\,(\text{contracted pairs})
 \times(\text{uncontracted rest})\Big\},$$

where the sign of each term is the parity of the permutation that
brings every contracted pair into adjacent order (fermions
anticommute, bosons commute, in this counting).

*Proof (induction in $n$).* For $n=2$ this is the definition of the
contraction. Assume the statement for $n-1$ fields. Both sides change
by the same graded sign under any permutation of the fields (that is
how $T$, $N$, and the pair-adjacency sign are all defined), so we may
relabel such that $t_1$ is the *earliest* time; then
$T\{\phi_1\cdots\phi_n\} = T\{\phi_2\cdots\phi_n\}\,\phi_1$, and by the
induction hypothesis the right-hand side is [the signed contraction sum
of $\phi_2\cdots\phi_n$ inside normal products] $\times\,\phi_1$.
Split $\phi_1 = \phi_1^+ + \phi_1^-$. The annihilation part $\phi_1^+$
already stands in normal position at the far right. The creation part
$\phi_1^-$ must be moved to the far left; each transposition past the
annihilation part $\phi_k^+$ of an uncontracted field produces the
graded commutator $[\phi_k^+,\phi_1^-]_\mp$ — a c-number equal to
$\langle0|\phi_k\phi_1|0\rangle = \langle0|T\phi_k\phi_1|0\rangle$
(since $t_k \ge t_1$), i.e. exactly the contraction
$\underbrace{\phi_k\phi_1}$ — while the transposition signs accumulate
into precisely the pair-adjacency sign of that new contraction.
Collecting the terms with and without the new contraction reproduces
the claimed sum for $n$ fields. $\;\blacksquare$

Taking the vacuum expectation value kills every term with uncontracted
fields:

$$\langle0|T\{\phi_1\cdots\phi_n\}|0\rangle
 = \sum_{\text{complete contractions}} (\pm)\ \prod \text{propagators}.$$

**Sign rules** (consequences used constantly below):

* *Open fermion line.* In the canonical operator order
  $\psi_{\rm out}\,\bar\psi_{\rm in}\,
  (\bar\psi_1\psi_1)(\bar\psi_2\psi_2)\cdots$ every contraction whose
  fermion pairs form only an open line from $\bar\psi_{\rm in}$ to
  $\psi_{\rm out}$ has total sign $+1$.
* *Closed fermion loop.* Each closed loop contributes one factor
  $(-1)$ and a Dirac trace. Derivation on the shortest case, the
  two-vertex loop
  $\langle(\bar\psi\Gamma^{(1)}\psi)_x(\bar\psi\Gamma^{(2)}\psi)_y\rangle$:
  the only contraction is
  $\bar\psi_a(x)\Gamma^{(1)}_{ab}\psi_b(x)\,
  \bar\psi_c(y)\Gamma^{(2)}_{cd}\psi_d(y)$ with
  $\underbrace{\psi_b(x)\bar\psi_c(y)} = S_{bc}(x-y)$ and the pair
  $(\bar\psi_a,\ \psi_d)$ contracted *in reversed order*: bringing
  $\psi_d$ in front of $\bar\psi_a$ costs the T-reordering sign
  $\langle T\bar\psi_a\psi_d\rangle = -S_{da}(y-x)$, so the term equals
  $-\,S_{da}(y-x)\Gamma^{(1)}_{ab}S_{bc}(x-y)\Gamma^{(2)}_{cd}
  = -\,\mathrm{tr}\left[\Gamma^{(1)}S(x-y)\Gamma^{(2)}S(y-x)\right]$.
  The general loop works identically (one reversed pair per loop).

Both rules — and nothing else about signs — are verified for **every
single contraction** through two loops by `code/wick_enum.py` (§5),
which computes each term's permutation sign explicitly and asserts
$\text{sign} = (-1)^{\#\text{closed loops}}$.

### 4.3 The Ward–Takahashi identity, $Z_1 = Z_2$

The fermion equation of motion from §1,
$(i\slashed{\partial} - m_0 - e\slashed{A})\psi = 0$ and its conjugate,
give current conservation as an operator equation regardless of the
photon mass (the $\slashed{A}$ and mass terms cancel between the pair):

$$\partial_\mu j^\mu = 0,\qquad j^\mu = \bar\psi\gamma^\mu\psi .$$

Differentiate the time-ordered product (the $T$-symbol's
$\theta$-functions produce equal-time commutators):

$$\partial_\mu^z\,T\,j^\mu(z)\psi(x)\bar\psi(y)
 = T\,(\partial_\mu j^\mu)(z)\,\psi(x)\bar\psi(y)
 + \delta(z^0-x^0)\,T\left[j^0(z),\psi(x)\right]\bar\psi(y)
 + \delta(z^0-y^0)\,T\,\psi(x)\left[j^0(z),\bar\psi(y)\right].$$

The first term vanishes; the equal-time commutators follow from
$\{\psi_a(\mathbf{z}),\psi_b^\dagger(\mathbf{x})\} =
\delta_{ab}\delta^{(3)}(\mathbf{z}-\mathbf{x})$ via
$[AB,C] = A\{B,C\} - \{A,C\}B$:

$$\left[j^0(z),\psi(x)\right]_{z^0=x^0} = -\,\delta^{(3)}(\mathbf{z}-\mathbf{x})\,\psi(x),
\qquad
\left[j^0(z),\bar\psi(y)\right]_{z^0=y^0} = +\,\delta^{(3)}(\mathbf{z}-\mathbf{y})\,\bar\psi(y).$$

Fourier transforming (conventions of the QFT section) and amputating
the two fermion legs with full propagators turns this into the
momentum-space **Ward–Takahashi identity** for the proper vertex; the
overall normalization is fixed by the tree-level check
($\Gamma^\mu = \gamma^\mu$, $\tilde S'^{-1} = -i(\slashed{p}-m)$):

$$q_\mu\,\Gamma^\mu(p',p) = i\left[\tilde S'^{-1}(p')
  - \tilde S'^{-1}(p)\right]
  = \slashed{q} - \left[\Sigma(\slashed{p}') - \Sigma(\slashed{p})\right].$$

Two consequences used above and below:

* **On shell** both sides sandwiched with $\bar u(p')\dots u(p)$
  vanish: $q_\mu\bar u'\Gamma^\mu u = 0$ (used in §2.2).
* **$q\to0$:** write $\Sigma(\slashed{p}) = A(p^2)+B(p^2)\slashed{p}$
  and expand the identity to first order in $q$; sandwiching and using
  the $q=0$ Gordon identity $\bar u\, p^\mu\, u = m\,\bar u\gamma^\mu u$
  to convert the $P^\mu$-structures,

  $$F_1(0) = 1 - \left[B + 2mA' + 2m^2B'\right]_{p^2=m^2}
     = 1 - \Sigma'(m) = 1 - \delta Z_2 .$$

  Hence with $F_1(0) \equiv 1 + \delta F_1(0)$ and
  $Z_1^{-1} \equiv F_1(0)$:

  $$\boxed{\ \delta F_1(0) = -\,\delta Z_2,\qquad Z_1 = Z_2
    + O(\alpha^2).\ }$$

  (The identity holds to all orders; we only need it at order
  $\alpha$.) In particular $F_1^S(0) = Z_2F_1(0) = 1$ exactly: the
  charge is not renormalized by the vertex, and the LO section's
  renormalization condition "$F_1(0)=1$ to all orders" is automatic.

### 4.4 Momentum-space Feynman rules (assembled once)

Applying Wick's theorem to the Dyson term of §4.1 and Fourier
transforming every propagator, each $\int\mathrm{d}^4x_i$ produces
$(2\pi)^4\delta^{(4)}(\sum k_{\rm in})$ at vertex $i$; one delta
survives as overall momentum conservation and the remaining momentum
integrals are the loop integrals. The dictionary (all derived above):

| element | factor |
| --- | --- |
| vertex | $-ie\gamma^\mu$ |
| mass counterterm vertex | $+i\,\delta m$ |
| fermion propagator | $i(\slashed{k}+m)/(k^2-m^2+i\epsilon)$ |
| photon propagator | $-ig_{\mu\nu}/(k^2-\lambda^2+i\epsilon)$ |
| closed fermion loop | $(-1)\times\mathrm{tr}[\cdots]$, $\int\frac{\mathrm{d}^4l}{(2\pi)^4}$ |
| loop momentum | $\int\frac{\mathrm{d}^4k}{(2\pi)^4}$ |

**Combinatorial weight.** For a fixed topology with no vertex-exchange
automorphism, the $n!$ assignments of the labeled integration variables
$x_1\dots x_n$ to the $n$ vertex positions give identical integrals and
cancel the Dyson $1/n!$ exactly: every distinct diagram enters **with
weight one**. This is not left as folklore: `code/wick_enum.py` counts
the contractions of every topology at $n = 1,3,5$ and verifies that
each single-topology class appears exactly $n!$ times per orientation
or mirror copy.

## 5. The expansion, executed

`code/wick_enum.py` enumerates *every* complete contraction of

$$\langle0|T\,\tilde\psi(p')\,\tilde A^\mu(-q)\,\tilde{\bar\psi}(-p)\,
  \left(\bar\psi\slashed{A}\psi\right)(x_1)\cdots
  \left(\bar\psi\slashed{A}\psi\right)(x_n)|0\rangle,\qquad n = 1,3,5:$$

all perfect matchings of the $n{+}1$ photon fields times all $(n{+}1)!$
pairings of $\psi$'s with $\bar\psi$'s, then classifies the topology
and computes the fermion sign of each term. Output
(`pixi run python code/wick_enum.py`):

    order e^1:  2 complete contractions (1 photon matchings x 2! fermion pairings)
      topology class                                                      count  /n! sign
      tadpole loop  tr[gamma S(0)] = 0                                        1    1 -1
      tree vertex gamma^mu                                                    1    1 +1
      total                                                                   2    2
      all signs equal (-1)^(# closed fermion loops): OK

    order e^3:  72 complete contractions (3 photon matchings x 4! fermion pairings)
      topology class                                                      count  /n! sign
      tadpole loop  tr[gamma S(0)] = 0                                       39    - +1/-1
      external-leg corrections only (tree vertex)                            12    2 +1
      III: VP on the external-potential line x tree vertex                    6    1 -1
      LO one-loop vertex correction                                           6    1 +1
      forward line x (rest): no scattering / disconnected                     6    1 -1
      disconnected (vacuum bubble or detached blob)                           3    - -1
      total                                                                  72   12
      all signs equal (-1)^(# closed fermion loops): OK

    order e^5:  10800 complete contractions (15 photon matchings x 6! fermion pairings)
      topology class                                                      count  /n! sign
      tadpole loop  tr[gamma S(0)] = 0                                     6165    - +1/-1
      external-leg corrections only (tree vertex)                           840    7 +1
      forward line x (rest): no scattering / disconnected                   660    - +1/-1
      disconnected (vacuum bubble or detached blob)                         495    - +1/-1
      V: odd fermion loop (Furry pair, cancels)                             480    4 -1
      III: two-loop VP blob on the external-potential line                  360    3 -1
      III: VP on the external-potential line x leg SE                       240    2 -1
      IIb/IIf: LO vertex x self-energy on an external leg                   240    2 +1
      IIc: corner, vertex part at an internal vertex (+ mirror)             240    2 +1
      IId: self-energy on an internal line (+ mirror)                       240    2 +1
      external-leg SE with VP insertion                                     240    2 -1
      I: crossed ladder                                                     120    1 +1
      III: VP on the external-potential line x vertex corr.                 120    1 -1
      IIa: ladder (vertex part at the external vertex)                      120    1 +1
      IIe: vacuum polarization in the internal photon                       120    1 -1
      IV: iterated one-loop VP chain on the external-potential line         120    1 +1
      total                                                               10800   90
      all signs equal (-1)^(# closed fermion loops): OK

The classification logic is elementary once the fermion chain is laid
out: order the open-line vertices from the incoming electron as
$v_1\dots v_k$, note where the external photon attaches, and note which
pairs the two internal photons connect; closed fermion cycles are the
loops. Reading the $n=5$ table row by row gives the complete fate map of
the two-loop expansion:

* **Kept (the 1PI vertex cores):**
  * *I — crossed ladder* (weight 1): both photons span the external
    vertex and cross;
  * *IIa — ladder* (weight 1): both span, nested — a one-loop vertex
    subgraph at the external vertex;
  * *IIc — corner* (weight 2 = the two mirror reflections): outer
    photon spans the external vertex, inner photon spans one of its
    internal attachment points;
  * *IId* (weight 2 = mirrors): self-energy subgraph on one of the two
    internal fermion propagators of the LO diagram;
  * *IIe* (weight 1): fermion loop (with its $-1$ and trace) inserted
    in the internal photon.

  The mirror image of a diagram (reflect through the external vertex,
  i.e. read the string backwards and exchange $p\leftrightarrow p'$)
  has the same $F_2(0)$ — the Breit-frame projection of §2.3 is
  symmetric under $w\to-w$ combined with the relabeling — so the NLO
  section computes one orientation of IIc and IId and doubles it.
  The five *distinct* contributions are exactly Karplus–Kroll's.
* **External-leg classes** (IIb/IIf, the seven "external-leg-only"
  topologies, "leg SE with VP insertion"): all live inside
  $\tilde S'(p)$ or $\tilde S'(p')$ and are consumed by LSZ (§3.3),
  leaving the single factor $Z_2$. No integral to compute, ever.
* **External-potential classes III and IV** (VP blob, chain, or blob ×
  line correction on the $q$ line): the $1/(1-\Pi(q^2))$ dressing;
  after on-shell charge renormalization they contribute
  $\hat\Pi(0) = 0$ to $F_2(0)$ (§3.3). The weight 3 of the two-loop
  blob row is the three standard two-loop photon-self-energy graphs;
  none of this is ever integrated either.
* **Class V — odd fermion loops** (a fermion triangle attached to the
  potential and by two photons to the electron line, or attached once
  to the line and closed on itself by a photon with the potential on
  the line): vanish pairwise by **Furry's theorem**, proved as
  follows. The charge
  conjugation matrix $C = i\gamma^2\gamma^0$ satisfies
  $C^{-1}\gamma^\mu C = -(\gamma^\mu)^T$, hence
  $C^{-1}\tilde S(k)C = \tilde S(-k)^T$. Insert $CC^{-1}$ between all
  factors of an $\ell$-vertex loop trace and transpose:

  $$\mathrm{tr}\left[\gamma^{\mu_1}\tilde S(k_1)\cdots
    \gamma^{\mu_\ell}\tilde S(k_\ell)\right]
    = (-1)^\ell\,
    \mathrm{tr}\left[\gamma^{\mu_1}\tilde S(-k_\ell)\gamma^{\mu_\ell}
    \tilde S(-k_{\ell-1})\cdots\gamma^{\mu_2}\tilde S(-k_1)\right],$$

  which after renaming loop momenta is exactly the trace of the
  orientation-reversed loop. The enumeration confirms that odd loops
  come in orientation pairs (the count $480 = 2$ subclasses $\times$ an
  orientation pair $\times\,5!$); for odd $\ell$ each pair sums to
  zero. The $\ell=1$ tadpole is
  its own reversal and vanishes identically,
  $\mathrm{tr}[\gamma^\mu\tilde S(l)] \propto l^\mu \to 0$ by symmetric
  integration — consistent with the same $C$ argument.
* **Disconnected / forward classes:** cancelled by the
  $\langle0|S|0\rangle$ denominator (§4.1) or proportional to
  no-scattering $\delta$-functions; they never reach $\Gamma^\mu$.
* **$\delta m$ counterterm terms** (order $e^3\times\delta m$, not in
  the table): insertions on the two internal lines pair with IId
  (§8.4); insertions on external legs join the external-leg class and
  drop with it.

At order $e^3$ the same table shows the familiar one-loop result: the
proper vertex is the single LO diagram (weight 1) — everything else is
legs, potential dressing, tadpoles, or vacuum.

## 6. The loop-integration toolkit, derived

These are the five mechanical steps of `code/loops.py`, used by every
diagram below.

### 6.1 Feynman parameters

Base case, by explicit antiderivative:

$$\int_0^1\frac{\mathrm{d}x}{\left[xA+(1-x)B\right]^2}
 = \left[\frac{-1}{(A-B)\left[B+x(A-B)\right]}\right]_0^1
 = \frac{-1}{(A-B)}\left(\frac1A - \frac1B\right)
 = \frac{1}{AB}.$$

Induction step: assume the $(n{-}1)$-denominator formula; apply the base
case to combine the $(n{-}1)$-fold denominator with $A_n$, then rescale
the parameters — this yields

$$\frac{1}{A_1A_2\cdots A_n} = (n-1)!\int_0^1\mathrm{d}x_1\cdots\mathrm{d}x_n\,
  \frac{\delta\left(\textstyle\sum_i x_i-1\right)}{\left[\sum_i x_iA_i\right]^n}.$$

Repeated powers: differentiate w.r.t. $A_i$ ($\nu_i - 1$ times each) to
get the general weight used by `assemble()` in the scripts:

$$\frac{1}{\prod_i A_i^{\nu_i}}
 = \frac{\Gamma\!\left(\sum\nu_i\right)}{\prod_i\Gamma(\nu_i)}
   \int_0^1\prod_i\mathrm{d}x_i\,x_i^{\nu_i-1}\,
   \frac{\delta\left(\sum x_i-1\right)}{\left[\sum x_iA_i\right]^{\sum\nu_i}}.$$

### 6.2 Completing the square

With propagators $(k-a_i)^2 - m_i^2$ and weights $x_i$
($\sum x_i = 1$),

$$\sum_i x_i\left[(k-a_i)^2 - m_i^2\right] = (k-s)^2 - \Delta,
\qquad s = \sum_i x_i a_i,$$

$$\Delta = s^2 - \sum_i x_i a_i^2 + \sum_i x_i m_i^2 ,$$

by expanding both sides ($k^2$ coefficients match; linear terms fix
$s$; the rest is $\Delta$). The shift $k = \ell + s$ has unit Jacobian.
`feynman_shift` does exactly this and verifies the reduction to
$\ell^2 - \Delta$ symbolically.

### 6.3 Wick rotation

For $\Delta > 0$ the integrand poles in the $\ell^0$ plane sit at
$\ell^0 = \pm\left(\sqrt{\boldsymbol{\ell}^2+\Delta} - i\epsilon\right)$
— second and fourth quadrants. The closed contour (real axis, arc,
imaginary axis, arc) therefore encloses no pole, the arcs vanish for
convergent integrands, and

$$\int_{-\infty}^{\infty}\mathrm{d}\ell^0
 = \int_{-i\infty}^{i\infty}\mathrm{d}\ell^0
 = i\int_{-\infty}^{\infty}\mathrm{d}\ell_E^4,\qquad
 \ell^0 = i\ell_E^4,\quad \ell^2 = -\ell_E^2 .$$

### 6.4 The master table

Four-dimensional Euclidean spherical measure: from
$\pi^2 = \left(\int e^{-t^2}\mathrm{d}t\right)^4
= S_3\int_0^\infty r^3e^{-r^2}\mathrm{d}r = S_3/2$, the unit-sphere
surface is $S_3 = 2\pi^2$, so with $t = \ell_E^2$:

$$\int\mathrm{d}^4\ell_E\,f(\ell_E^2) = \pi^2\int_0^\infty t\,f(t)\,\mathrm{d}t .$$

For $n-a-2 > 0$ (convergent case),

$$\int_0^\infty\frac{t^{a+1}\,\mathrm{d}t}{(t+\Delta)^n}
 = \Delta^{a+2-n}\,B(a+2,\,n-a-2)
 = \Delta^{a+2-n}\,\frac{(a+1)!\,(n-a-3)!}{(n-1)!},$$

(the Beta integral, by the substitution $t = \Delta\,v/(1-v)$).
Assembling with the Wick-rotation factor $i\,(-1)^{a+n}$ from
$(\ell^2)^a = (-1)^a(\ell_E^2)^a$ and
$(\ell^2-\Delta)^n = (-1)^n(\ell_E^2+\Delta)^n$:

$$\int\frac{\mathrm{d}^4\ell}{(2\pi)^4}\,
  \frac{(\ell^2)^a}{(\ell^2-\Delta+i\epsilon)^n}
 = \frac{i\,(-1)^{a+n}}{16\pi^2}\,
   \frac{(a+1)!\,(n-a-3)!}{(n-1)!}\;\Delta^{a+2-n},$$

which is `loop_integral(a, n, Delta)` verbatim. The three cases used
below:

$$\int\frac{\mathrm{d}^4\ell}{(2\pi)^4}\frac{1}{(\ell^2-\Delta)^3}
  = \frac{-i}{32\pi^2\Delta},\qquad
  \int\frac{\mathrm{d}^4\ell}{(2\pi)^4}\frac{\ell^2}{(\ell^2-\Delta)^3}
  \ \text{and}\
  \int\frac{\mathrm{d}^4\ell}{(2\pi)^4}\frac{1}{(\ell^2-\Delta)^2}:
  \ \text{log divergent.}$$

### 6.5 Pauli–Villars for the logarithmic case

For $n-a-2 = 0$ subtract the identical integrand with
$\Delta \to \Delta_\Lambda$ ($\Delta$ shifted by the regulator mass
$\Lambda^2\to\infty$). The difference is convergent and elementary:
differentiating in $\Delta$ lands on the convergent case of §6.4,

$$\frac{\mathrm{d}}{\mathrm{d}\Delta}\int_0^\infty
 \frac{t^{a+1}\,\mathrm{d}t}{(t+\Delta)^{a+2}}
 = -(a+2)\int_0^\infty\frac{t^{a+1}\,\mathrm{d}t}{(t+\Delta)^{a+3}}
 = -(a+2)\,\frac{(a+1)!\;0!}{(a+2)!}\,\frac1\Delta = -\frac1\Delta ,$$

so the regularized radial integral is
$\int_\Delta^{\Delta_\Lambda}\mathrm{d}\Delta'/\Delta'
 = \log\left(\Delta_\Lambda/\Delta\right)$ and, restoring the
prefactors of §6.4,

$$\int\frac{\mathrm{d}^4\ell}{(2\pi)^4}\,
  \frac{(\ell^2)^a}{(\ell^2-\Delta)^n}\,\Bigg|_{\rm PV}
 = \frac{i\,(-1)^{a+n}}{16\pi^2}\,\frac{(a+1)!}{(n-1)!}
 \left(L_{UV} - \log\Delta\right),\qquad n-a-2 = 0,$$

with $L_{UV} \equiv \log\Lambda^2 + \text{const}$. The additive
constant (and the difference between subtracting at the integral level
and the historical propagator-level PV) is a scheme constant *inside*
$L_{UV}$; since $L_{UV}$ cancels identically in every renormalized
quantity below (asserted by the scripts), it never matters.

### 6.6 Angular averages

After the shift the denominator depends on $\ell^2$ only, so odd
powers of $\ell$ integrate to zero and even ones can be replaced by
their rotational averages (Euclidean rotation invariance after §6.3).
By Lorentz covariance $\langle\ell^\mu\ell^\nu\rangle = c\,g^{\mu\nu}\ell^2$;
tracing with $g_{\mu\nu}$ gives $4c = 1$:

$$\langle\ell^\mu\ell^\nu\rangle = \frac{g^{\mu\nu}}{4}\,\ell^2,\qquad
 \langle\ell^\mu\ell^\nu\ell^\rho\ell^\sigma\rangle
 = \frac{(\ell^2)^2}{24}\left(g^{\mu\nu}g^{\rho\sigma}
 + g^{\mu\rho}g^{\nu\sigma} + g^{\mu\sigma}g^{\nu\rho}\right),$$

the rank-4 coefficient fixed by the double trace ($16+4+4 = 24$). The
general even rank ($(\ell^2)^n/\left(2^n(n+1)!\right)\sum_{\rm
matchings}\prod g$) is what `symmetrize` implements; ranks 2 and 4 are
all the diagrams below need.

### 6.7 Dirac algebra contractions

From the Clifford algebra
$\{\gamma^\mu,\gamma^\nu\} = 2g^{\mu\nu}$ (each line uses the
previous):

$$\gamma_\nu\gamma^\nu = 4,\qquad
 \gamma_\nu\gamma^\alpha\gamma^\nu = -2\gamma^\alpha,\qquad
 \gamma_\nu\gamma^\alpha\gamma^\beta\gamma^\nu = 4g^{\alpha\beta},\qquad
 \gamma_\nu\gamma^\alpha\gamma^\beta\gamma^\gamma\gamma^\nu
 = -2\gamma^\gamma\gamma^\beta\gamma^\alpha,$$

$$\mathrm{tr}\,\mathbb{1} = 4,\qquad
 \mathrm{tr}\left[\gamma^\mu\gamma^\nu\right] = 4g^{\mu\nu},\qquad
 \mathrm{tr}\left[\gamma^\mu\gamma^\alpha\gamma^\nu\gamma^\beta\right]
 = 4\left(g^{\mu\alpha}g^{\nu\beta} - g^{\mu\nu}g^{\alpha\beta}
 + g^{\mu\beta}g^{\alpha\nu}\right),$$

odd traces vanishing. (The pipeline needs none of these as *rules* —
it multiplies explicit $4\times4$ matrices — but they are what make the
hand-derivations below short.)

## 7. One-loop building blocks

### 7.1 The LO vertex with a massive photon → the kernel $K(t)$

From the $n=3$ row of §5, the proper vertex at one loop is the single
diagram with the photon spanning the external vertex. Its Wick
contraction, Fourier transformed with the rules of §4.4 (routing as in
the LO section and `code/g2_lo_trace.py`: fermion momenta $k$, $k'=k+q$,
photon $k-p$):

$$\bar u(p')\,\delta\Gamma^\mu\,u(p) = -\,i e^2\int\frac{\mathrm{d}^4k}{(2\pi)^4}\;
 \frac{\bar u(p')\,\gamma^\nu\left(\slashed{k}'+m\right)\gamma^\mu
 \left(\slashed{k}+m\right)\gamma_\nu\,u(p)}
 {\left[(k-p)^2-\lambda^2\right]\left[k'^2-m^2\right]\left[k^2-m^2\right]},$$

the prefactor being $(-ie)^2\,i\,i\,(-i) = -ie^2$ from two vertices,
two fermion propagators and one photon propagator (the $g_{\nu\rho}$ of
the photon contracts the two vertex indices into
$\gamma^\nu\cdots\gamma_\nu$). The $i\epsilon$'s are kept implicit from
here on. Using §6.7 the numerator reduces exactly as in the LO section:

$$\gamma^\nu\left(\slashed{k}'+m\right)\gamma^\mu\left(\slashed{k}+m\right)\gamma_\nu
 = -2\slashed{k}\gamma^\mu\slashed{k}' + 4m\,(k+k')^\mu - 2m^2\gamma^\mu .$$

Feynman parameters (§6.1) with $x$ on the photon, $y$ on $k'$, $z$ on
$k$, $x = 1-y-z$; the shift (§6.2) with $a_x = p$, $a_y = -q$,
$a_z = 0$ gives $s = xp - yq$ and (using $p\cdot q = -q^2/2$,
$x+y+z=1$)

$$\Delta = (y+z)^2 m^2 \;-\; yz\,q^2 \;+\; (1-y-z)\,\lambda^2 ,$$

which at $q^2 = 0$ is exactly the `Delta` printed by
`pixi run lo-trace-sympy`. (The LO section's variables put the
parameter $z$ on the photon; the two labelings are related by
$z_{\rm LO} = 1-y-z$ and give identical integrals.) The numerator with
$k = \ell + s$ is expanded, averaged (§6.6), and integrated with the
master table (§6.4–6.5); the $F_2$ projection (§2.3) of the result —
performed symbolically by the pipeline, every intermediate checked —
gives at $q^2 = 0$, $m = 1$:

$$F_2\text{-integrand} = \frac{e^2}{4\pi^2}\,
  \frac{(y+z)(1-y-z)}{(y+z)^2 + \lambda^2(1-y-z)}$$

(the integrand printed by `pixi run lo-trace-sympy`),

while all of $L_{UV}$ lands in $F_1$ (asserted). Since the integrand
depends on $y,z$ only through $s = y+z$,
$\int_0^1\!\!\int_0^{1-z} \mathrm{d}y\,\mathrm{d}z\,f(y{+}z)
= \int_0^1 s f(s)\,\mathrm{d}s$, and with $s \to 1-z$:

$$F_2(0;\lambda) = \frac{\alpha}{\pi}\int_0^1
  \frac{s^2(1-s)}{s^2+\lambda^2(1-s)}\,\mathrm{d}s
 = \frac{\alpha}{\pi}\int_0^1
  \frac{z(1-z)^2}{(1-z)^2+z\,t}\,\mathrm{d}z
 = \frac{\alpha}{\pi}\,K(t),\qquad t = \frac{\lambda^2}{m^2},$$

**exactly the Stage-0 kernel of the NLO section**, with
$K(0) = \int_0^1 z\,\mathrm{d}z = \frac12$ reproducing Schwinger's
$a_e = \alpha/2\pi$ (the LO section's result).

The same computation defines the one-loop $F_1$: UV divergent
($L_{UV}$), IR divergent ($\log\lambda$), and its value at $q^2 = 0$,

$$\delta F_1(0;\lambda) \equiv F_1(0)\big|_{\alpha} ,$$

is the subtraction constant of diagrams IIa/IIc below. We never need
its value here, only its *integrand* (§7.4).

### 7.2 The electron self-energy $\Sigma(k)$

From the SE subgraph contraction (internal fermion momentum $l$,
photon $k-l$):

$$-\,i\Sigma_{\rm loop}(k) = -\,e^2\int\frac{\mathrm{d}^4l}{(2\pi)^4}\,
 \frac{\gamma^\nu\left(\slashed{l}+m\right)\gamma_\nu}
 {\left[(k-l)^2-\lambda^2\right]\left[l^2-m^2\right]}
 \quad\Longrightarrow\quad
 \Sigma_{\rm loop}(k) = -\,ie^2\int\frac{\mathrm{d}^4l}{(2\pi)^4}\,
 \frac{-2\slashed{l}+4m}
 {\left[(k-l)^2-\lambda^2\right]\left[l^2-m^2\right]},$$

the prefactor being $(-ie)^2\,i\,(-i) = -e^2$ (two vertices, one
fermion and one photon propagator) and the numerator reduced with
§6.7.

One Feynman parameter $u$ on the photon; $s = uk$;
$\Delta_{\rm in} = a - b\,k^2$ with

$$a = (1-u)\,m^2 + u\,\lambda^2,\qquad b = u(1-u)$$

(§6.2 with $a_{\rm ph} = k$, $m_{\rm ph} = \lambda$; $a_f = 0$,
$m_f = m$). The shifted numerator $-2(\slashed{\ell}+u\slashed{k})+4m$
loses its odd part, and the log-divergent master integral (§6.5) gives

$$\boxed{\;\Sigma_{\rm loop}(k) = \frac{e^2}{16\pi^2}\int_0^1\mathrm{d}u\,
 \left(4m - 2u\slashed{k}\right)
 \left(L_{UV} - \log\left(a - b\,k^2\right)\right)\;}$$

— the formula quoted in the NLO section and derived mechanically by
`pixi run iid-sympy`. Writing
$\Sigma = A(k^2) + B(k^2)\,\slashed{k}$, the on-shell constants of §3.2
are

$$\delta m = \left[A + mB\right]_{k^2=m^2},\qquad
 \delta Z_2 = \left[2mA' + B + 2m^2B'\right]_{k^2=m^2},$$

and the fully subtracted self-energy (double zero on shell)

$$\Sigma_R(k) = \Sigma_{\rm loop}(k) - \delta m
 - \left(\slashed{k}-m\right)\delta Z_2$$

is free of $L_{UV}$: the three $L_{UV}$ coefficients are
$(4m-2u\slashed{k})$, $-(4m-2um)$, and $+2u(\slashed{k}-m)$, which sum
to zero *pointwise in $u$* — the cancellation asserted by the script.
$\Sigma_R$ splits into a rational part and a
$\log\left(D_{\rm in}(k^2)/D_{\rm in}(m^2)\right)$ part; the log is
rationalized with the elementary identity (differentiate
$\log\left[Y+\xi(X-Y)\right]$ and integrate back)

$$\log\frac{X}{Y} = \int_0^1\mathrm{d}\xi\,\frac{X-Y}{Y+\xi(X-Y)},$$

which, applied to $X = a-bk^2$, $Y = a-bm^2$, produces the factor
$(k^2-m^2)$ that cancels one power of the doubled propagator and leaves
an ordinary propagator of squared mass
$C = m^2 + (a-bm^2)/(\xi b) > m^2$ — the NLO section's $\xi$-trick,
here derived.

### 7.3 The vacuum polarization $\hat\Pi$ and its spectral function

The fermion-loop contraction with its $(-1)$ and trace (§4.2):

$$i\Pi^{\mu\nu}(q) = (-1)\,(-ie)^2\,i^2\int\frac{\mathrm{d}^4l}{(2\pi)^4}\,
 \frac{\mathrm{tr}\left[\gamma^\mu\left(\slashed{l}+m\right)
 \gamma^\nu\left(\slashed{l}-\slashed{q}+m\right)\right]}
 {\left[l^2-m^2\right]\left[(l-q)^2-m^2\right]},$$

$$\mathrm{tr}\left[\cdots\right] = 4\left[l^\mu(l-q)^\nu + l^\nu(l-q)^\mu
 - g^{\mu\nu}\left(l\cdot(l-q) - m^2\right)\right] \qquad(\text{Sec. 6.7}).$$

**Transversality.** The Ward–Takahashi argument of §4.3 applied to
$\langle T j^\mu(z) j^\nu(y)\rangle$ (the equal-time commutator
$[j^0,j^\nu]$ vanishes up to Schwinger terms that a gauge-invariant
regulator cancels) gives $q_\mu\Pi^{\mu\nu} = 0$, i.e.

$$\Pi^{\mu\nu}(q) = \left(q^2g^{\mu\nu} - q^\mu q^\nu\right)\Pi(q^2).$$

Maintaining this in regularization is exactly why the fermion loop
needs *gauge-invariant* Pauli–Villars (two heavy masses $M_i$ with
coefficients satisfying $\sum c_i = 0$, $\sum c_iM_i^2 = 0$, which
remove the quadratically divergent non-transverse polynomial). Granted
transversality, $\Pi(q^2)$ can be read off from the $q^\mu q^\nu$
coefficient, which is only *log* divergent and hence unambiguous:
Feynman parameter $x$, shift $l = \ell + xq$,

$$l^\mu(l-q)^\nu + l^\nu(l-q)^\mu \to
 2\,\langle\ell^\mu\ell^\nu\rangle - 2x(1-x)\,q^\mu q^\nu,
 \qquad \Delta_x = m^2 - x(1-x)\,q^2,$$

so the $q^\mu q^\nu$ coefficient of $\Pi^{\mu\nu}$ is
$ie^2\cdot4\cdot\left(-2x(1-x)\right)$ integrated with the log master
integral, and matching to $-q^\mu q^\nu\,\Pi(q^2)$:

$$\Pi(q^2) = -\frac{2\alpha}{\pi}\int_0^1\mathrm{d}x\,x(1-x)
 \left(L_{UV} - \log\Delta_x\right),$$

$$\boxed{\;\hat\Pi(q^2) = \Pi(q^2) - \Pi(0)
 = \frac{2\alpha}{\pi}\int_0^1\mathrm{d}x\;x(1-x)\,
 \log\frac{m^2 - x(1-x)\,q^2}{m^2}\;}$$

— finite, with $\hat\Pi(0) = 0$: the once-subtracted (on-shell charge
renormalized) vacuum polarization of §3.3.

**Spectral (dispersion) form.** $\hat\Pi(q^2)/q^2$ is analytic in the
$q^2$ plane cut along $[4m^2,\infty)$ (the log argument
$m^2 - x(1-x)q^2$ can vanish only for $q^2 \ge 4m^2$ since
$x(1-x)\le\frac14$) and falls off like $\log(-q^2)/q^2$; the Cauchy
formula on a contour hugging the cut plus a large circle (whose
contribution vanishes) gives the unsubtracted dispersion relation

$$\frac{\hat\Pi(k^2)}{k^2}
 = \int_{4m^2}^{\infty}\mathrm{d}t\;\frac{\rho(t)}{k^2 - t},
 \qquad
 \rho(t) = -\frac{1}{\pi t}\,\mathrm{Im}\,\hat\Pi(t+i\epsilon).$$

The imaginary part is elementary: with $q^2 = t+i\epsilon$ the argument
$m^2 - x(1-x)t$ becomes negative for
$x\in(x_-,x_+)$, $x_\pm = \frac12(1\pm\beta)$,
$\beta = \sqrt{1-4m^2/t}$, where the log acquires $-i\pi$; hence

$$\mathrm{Im}\,\hat\Pi(t+i\epsilon)
 = -2\alpha\int_{x_-}^{x_+}x(1-x)\,\mathrm{d}x
 = -2\alpha\cdot\frac{\beta}{12}\left(3-\beta^2\right)$$

(substitute $x = \frac12(1+\beta v)$, $v\in(-1,1)$:
$\int = \frac{\beta}{8}\int_{-1}^1(1-\beta^2v^2)\mathrm{d}v
= \frac{\beta}{4}\left(1-\frac{\beta^2}{3}\right)$), and with
$3-\beta^2 = 2+4m^2/t$:

$$\boxed{\;\rho(t) = \frac{\alpha}{3\pi t}
 \left(1+\frac{2m^2}{t}\right)\sqrt{1-\frac{4m^2}{t}}\;}$$

— precisely the spectral density used for diagram IIe in the NLO
section, and the dressed-propagator statement quoted there:

$$\frac{1}{k^2} \;\to\; \frac{1}{k^2}
 + \int_{4m^2}^\infty\mathrm{d}t\,\frac{\rho(t)}{k^2-t}.$$

Both the derivative at zero and the full function check numerically
against the parametric form:

    >>> from mpmath import mp, quad, mpf, log, sqrt
    >>> mp.dps = 30
    >>> r = lambda t: (1 + 2/t)*sqrt(1 - 4/t)/(3*t)   # rho = (alpha/pi) r,  m = 1
    >>> Pihat = lambda q2: 2*quad(lambda x: x*(1-x)*log(1 - x*(1-x)*q2), [0, 1])
    >>> quad(lambda t: r(t)/t, [4, mp.inf])            # = -Pihat'(0) = 1/15
    mpf('0.0666666666666666666666666666666634')
    >>> mpf(1)/15
    mpf('0.0666666666666666666666666666666634')
    >>> Pihat(-5)/-5 - quad(lambda t: r(t)/(-5 - t), [4, mp.inf])
    mpf('0.0')

### 7.4 One-loop vertex parts with open legs: $\Lambda^\nu(k)$

Diagrams IIa and IIc contain a one-loop vertex subgraph whose photon
leg (momentum $k$) is internal and off shell. Two versions occur.

**IIc version** (one electron leg on shell at $p'$, one internal with
momentum $p'-k$; inner photon momentum $k_2$):

$$\Lambda^\nu(k) = -\,ie^2\int\frac{\mathrm{d}^4k_2}{(2\pi)^4}\,
 \frac{\gamma^\rho\left(\slashed{p}'-\slashed{k}_2+m\right)\gamma^\nu
 \left(\slashed{p}'-\slashed{k}-\slashed{k}_2+m\right)\gamma_\rho}
 {\left[k_2^2-\lambda^2\right]
  \left[(p'-k_2)^2-m^2\right]
  \left[(p'-k-k_2)^2-m^2\right]} .$$

Feynman parameters $u$ on the $(p'-k-k_2)$ line, $v$ on the $(p'-k_2)$
line, $1-u-v$ on the photon. The shift formula (§6.2) with
$a_u = p'-k$, $a_v = p'$, $a_{\rm ph} = 0$ gives
$s = (u+v)p' - uk$ and (using $p'^2 = m^2$; verified symbolically)

$$\Delta_{\rm in} = D_0 + L(k) - \hat b\,k^2,\qquad
 D_0 = (u+v)^2m^2 + (1-u-v)\lambda^2,$$

$$L(k) = 2u\,(1-u-v)\;p'\!\cdot k,\qquad \hat b = u(1-u).$$

$D_0$ is the LO $\Delta$ (at $q^2=0$) with $z = u+v$ — the origin of
the "$D_0$" in the NLO section. After the shift and angular average the
numerator has an $\ell^2$ part and an $\ell$-free part $N^\nu$. The
$\ell^2$ part is proportional to $\gamma^\nu$ with a *constant*
coefficient: by §6.6 and §6.7,
$\gamma^\rho\slashed{\ell}\gamma^\nu\slashed{\ell}\gamma_\rho \to
\frac{\ell^2}{4}\gamma^\rho\gamma_\alpha\gamma^\nu\gamma^\alpha\gamma_\rho
= \frac{\ell^2}{4}\,\gamma^\rho(-2\gamma^\nu)\gamma_\rho = \ell^2\gamma^\nu$.
Applying the master integrals ($n=3$: the $\ell^0$ part is finite, the
$\ell^2$ part log divergent):

$$\boxed{\;\Lambda^\nu(k) = \int_0^1\!\!\mathrm{d}u\int_0^{1-u}\!\!\mathrm{d}v
 \left[\frac{P^\nu(k)}{\Delta_{\rm in}(k)}
 + \frac{e^2}{8\pi^2}\left(L_{UV} - \log\Delta_{\rm in}(k)\right)
 \gamma^\nu\right]\;}$$

with $P^\nu(k) = -\frac{e^2}{16\pi^2}N^\nu(k)$ — exactly the
$\Lambda^\nu$ structure of the NLO section
($c = e^2/8\pi^2$), whose explicit $N^\nu$ the pipeline
(`pixi run iic-sympy`) constructs as matrices.

**KK "reduced diagram" subtraction.** Setting $k = 0$ and $q = 0$
(both electron legs on shell) in the *integrand* turns it into the
$(u,v)$-integrand of the one-loop $F_1(0)$, i.e. of
$\delta F_1(0;\lambda)$ (§7.1) — the vertex subgraph collapses to
$\gamma^\nu\times$scalar by §2.2 at $q=0$. Define the pointwise
subtraction

$$\Lambda_R^\nu(k) \equiv \Lambda^\nu(k)
 - \left[\frac{L_a(u,v)}{D_0}
 + \frac{e^2}{8\pi^2}\left(L_{UV}-\log D_0\right)\right]\gamma^\nu
 \quad\text{under the }(u,v)\text{ integral},\qquad
 \int\!\!\int \left[\cdots\right]\mathrm{d}u\,\mathrm{d}v = \delta F_1(0;\lambda)\,\gamma^\nu.$$

$L_{UV}$ cancels *pointwise* (same constant $c$), so $\Lambda_R^\nu$
is UV finite before any outer integration — this is the subtraction
Karplus–Kroll perform "within the diagram", and the reason the NLO
section's IIa and IIc integrands are finite. The log-ratio
$\log(\Delta_{\rm in}/D_0)$ is rationalized with the same $\xi$-identity
as in §7.2.

**IIa version** (both electron legs internal, $p-k$ and $p'-k$; the
subgraph dresses the *external* vertex, so the open index is $\mu$ and
the external momentum $q$ enters the loop):

$$\Lambda^\mu_{\rm IIa}(k) = -\,ie^2\int\frac{\mathrm{d}^4k_2}{(2\pi)^4}\,
 \frac{\gamma^\rho\left(\slashed{p}'-\slashed{k}-\slashed{k}_2+m\right)
 \gamma^\mu\left(\slashed{p}-\slashed{k}-\slashed{k}_2+m\right)\gamma_\rho}
 {\left[k_2^2-\lambda^2\right]\left[(p'-k-k_2)^2-m^2\right]
 \left[(p-k-k_2)^2-m^2\right]},$$

with parameters $u,v$ on the two fermion lines,
$s = u(p'-k)+v(p-k)$ and (verified symbolically)

$$\Delta_{\rm in}^{\rm IIa} = s^2 - u\,(p'-k)^2 - v\,(p-k)^2
 + (u+v)m^2 + (1-u-v)\lambda^2,
 \qquad \hat b_{\rm IIa} = (u+v)(1-u-v)$$

for the $k^2$ coefficient. The same pointwise subtraction
$L(u,v)\gamma^\mu$ applies (same on-shell limit, same
$\delta F_1(0;\lambda)$).

## 8. The five two-loop integrals

These are the order-$e^5$ 1PI cores of §5, written with the rules of
§4.4. Overall factors: four internal vertices $(-ie)^4 = e^4$; the
external vertex's $-ie$ is factored out into the definition of
$\Gamma^\mu$ as always. Momentum routings match the pipeline scripts
exactly. All expressions are for the vertex matrix
$\delta\Gamma^\mu$; the number that the NLO section calls
$\mu_{\rm X}$ is its projection (§2.3) in units $(\alpha/\pi)^2$,
including the stated subtraction and mirror factor.

### 8.1 Diagram I (crossed ladder)

![Diagram I](figures/g2-nlo-I.svg)

Photons $k_a$ (vertices 1,4 counting from the incoming electron) and
$k_b$ (vertices 2,5), external photon at vertex 3. Four fermion
propagators ($i^4 = 1$), two photon propagators ($(-i)^2 = -1$):

$$\delta\Gamma^\mu_{\rm I} = -\,e^4
 \int\frac{\mathrm{d}^4k_a}{(2\pi)^4}\frac{\mathrm{d}^4k_b}{(2\pi)^4}\;
 \frac{\gamma^\rho\left(\slashed{p}'-\slashed{k}_b+m\right)
 \gamma^\nu\left(\slashed{p}'-\slashed{k}_a-\slashed{k}_b+m\right)
 \gamma^\mu\left(\slashed{p}-\slashed{k}_a-\slashed{k}_b+m\right)
 \gamma_\rho\left(\slashed{p}-\slashed{k}_a+m\right)\gamma_\nu}
 {\left[k_a^2-\lambda^2\right]\left[k_b^2-\lambda^2\right]
  \left[(p-k_a)^2-m^2\right]\left[(p-k_a-k_b)^2-m^2\right]
  \left[(p'-k_a-k_b)^2-m^2\right]\left[(p'-k_b)^2-m^2\right]}.$$

(Read the numerator right to left along the electron line: emit $k_a$
($\gamma_\nu$), emit $k_b$ ($\gamma_\rho$), absorb $q$ ($\gamma^\mu$),
absorb $k_a$ ($\gamma^\nu$), absorb $k_b$ ($\gamma^\rho$) — the
crossing is the index pattern
$\nu\ldots\rho\ldots\mu\ldots\nu\ldots\rho$.) No subdiagram is a
vertex or self-energy part, so there is nothing to subtract; the
$k_a$- and $k_b$-subintegrals are separately convergent by power
counting (each one-loop subintegral sees one photon and three fermion
propagators: superficial degree $4-2-3 = -1$), the overall degree is
$8-4-4 = 0$, and the log lives in
the $\gamma^\mu$ ($F_1$) structure only — the $F_2$ projection is UV
finite as it stands, and IR finite (no $\log\lambda$) in the KK scheme.
$$\mu_{\rm I} = \text{projection of }\delta\Gamma^\mu_{\rm I}.$$

### 8.2 Diagram IIa (ladder)

![Diagram IIa](figures/g2-nlo-IIa.svg)

Outer photon $k$ (vertices 1,5), inner photon $k_2$ (vertices 2,4)
nested around the external vertex 3 — i.e. the inner vertex part
$\Lambda^\mu_{\rm IIa}(k)$ of §7.4 inserted at the external vertex of
an LO-type outer loop:

$$\delta\Gamma^\mu_{\rm IIa} = -\,ie^2
 \int\frac{\mathrm{d}^4k}{(2\pi)^4}\;
 \frac{\gamma^\nu\left(\slashed{p}'-\slashed{k}+m\right)\,
 \Lambda^\mu_{\rm IIa}(k)\,
 \left(\slashed{p}-\slashed{k}+m\right)\gamma_\nu}
 {\left[k^2-\lambda^2\right]\left[(p'-k)^2-m^2\right]
 \left[(p-k)^2-m^2\right]} ,$$

which expanded is the double integral

$$\delta\Gamma^\mu_{\rm IIa} = -\,e^4\int\!\!\int
 \frac{\mathrm{d}^4k\,\mathrm{d}^4k_2}{(2\pi)^8}\,
 \frac{\gamma^\nu\left(\slashed{p}'-\slashed{k}+m\right)
 \gamma^\rho\left(\slashed{p}'-\slashed{k}-\slashed{k}_2+m\right)
 \gamma^\mu\left(\slashed{p}-\slashed{k}-\slashed{k}_2+m\right)
 \gamma_\rho\left(\slashed{p}-\slashed{k}+m\right)\gamma_\nu}
 {\left[k^2-\lambda^2\right]\left[k_2^2-\lambda^2\right]
 \left[(p'-k)^2-m^2\right]\left[(p'-k-k_2)^2-m^2\right]
 \left[(p-k-k_2)^2-m^2\right]\left[(p-k)^2-m^2\right]} .$$

**KK subtraction:** replace $\Lambda^\mu_{\rm IIa} \to
\Lambda^\mu_{\rm IIa} - L(u,v)\gamma^\mu$ pointwise (§7.4); the
subtracted piece, $-L(u,v)\times$(LO diagram integrand), factorizes and
integrates to $-\delta F_1(0;\lambda)\cdot F_2^{(2)}(0)$. The ladder is
its own mirror image: **no doubling**.
$$\mu_{\rm IIa} = \text{projection with the pointwise subtraction}.$$
This — with §7.4's $\Delta_{\rm in}^{\rm IIa}$ — is the complete
defining integral for the NLO section's TODO diagram IIa
(`code/g2_iia.py` implements exactly this construction: the three
pieces (a) $P^\mu/\Delta_{\rm in}$, (b) $-L/D_0$, (c) the
$\log(\Delta_{\rm in}/D_0)$ term via the $\xi$-identity).

### 8.3 Diagram IIc (corner)

![Diagram IIc](figures/g2-nlo-IIc.svg)

Outer photon $k$ (vertices 1,4), inner photon $k_2$ (vertices 3,5)
around the outer photon's $p'$-side attachment, external photon at
vertex 2 — i.e. $\Lambda^\nu(k)$ of §7.4 inserted at an *internal*
vertex:

$$\delta\Gamma^\mu_{\rm IIc} = -\,ie^2
 \int\frac{\mathrm{d}^4k}{(2\pi)^4}\,
 \frac{\Lambda^\nu(k)\,
 \left(\slashed{p}'-\slashed{k}+m\right)\gamma^\mu
 \left(\slashed{p}-\slashed{k}+m\right)\gamma_\nu}
 {\left[k^2-\lambda^2\right]\left[(p'-k)^2-m^2\right]
 \left[(p-k)^2-m^2\right]}$$

$$= -\,e^4\int\!\!\int\frac{\mathrm{d}^4k\,\mathrm{d}^4k_2}{(2\pi)^8}\,
 \frac{\gamma^\rho\left(\slashed{p}'-\slashed{k}_2+m\right)
 \gamma^\nu\left(\slashed{p}'-\slashed{k}-\slashed{k}_2+m\right)
 \gamma_\rho\left(\slashed{p}'-\slashed{k}+m\right)
 \gamma^\mu\left(\slashed{p}-\slashed{k}+m\right)\gamma_\nu}
 {\left[k^2-\lambda^2\right]\left[k_2^2-\lambda^2\right]
 \left[(p'-k_2)^2-m^2\right]\left[(p'-k-k_2)^2-m^2\right]
 \left[(p'-k)^2-m^2\right]\left[(p-k)^2-m^2\right]}$$

(note the outer photon propagator carries $k^2$ itself in this
routing, as in `code/g2_iic.py`). **KK subtraction:** the same
pointwise $-L(u,v)\gamma^\nu$; the subtracted cross term again
integrates to $-\delta F_1(0;\lambda)\,F_2^{(2)}(0)$ — this is the NLO
section's piece (b), which factorizes and carries the entire
$-\log\lambda$ of IIc. **Mirror:** the reflected diagram (subgraph at
the $p$-side attachment) contributes equally: factor 2.
$$\mu_{\rm IIc} = 2\times\text{projection of the subtracted diagram},$$

i.e. the sum of the NLO section's pieces (a), (b), (c), whose
integrands $f_a, f_b, f_c$ already include this factor 2.

### 8.4 Diagram IId (self-energy insertion) + mass counterterm

![Diagram IId](figures/g2-nlo-IId.svg)
![mass counterterm insertion](figures/g2-nlo-deltam.svg)

The raw insertion (routing of `code/g2_iid.py`: fermion momenta $k$,
$k' = k+q$, photon $k-p$; the SE sits on the $k$ line, giving the
doubled propagator; five fermion propagators $i^5 = i$, two photons
$(-i)^2 = -1$, vertices $e^4$):

$$\delta\Gamma^\mu_{\rm IId,raw} = -\,ie^2
 \int\frac{\mathrm{d}^4k}{(2\pi)^4}\;
 \frac{\gamma^\nu\left(\slashed{k}'+m\right)\gamma^\mu
 \left(\slashed{k}+m\right)\,
 \Sigma_{\rm loop}(k)\,
 \left(\slashed{k}+m\right)\gamma_\nu}
 {\left[(k-p)^2-\lambda^2\right]\left[k'^2-m^2\right]
 \left[k^2-m^2\right]^2},$$

with $\Sigma_{\rm loop}$ the integral of §7.2 (substituting it makes
this the explicit double integral; the collected prefactor is
$(-ie)^2\,i^3\,(-i)\times(-i\Sigma$-blob$) = -ie^2\,\Sigma$ from two
outer vertices, three fermion and one photon propagator, and the
blob). The $\delta m$ counterterm diagram
is the same string with $-i\Sigma_{\rm loop} \to +i\delta m$, i.e. it
replaces $\Sigma_{\rm loop} \to \Sigma_{\rm loop} - \delta m$. The KK
(= full on-shell) scheme also removes the $\delta Z_2$ part:

$$\delta\Gamma^\mu_{\rm IId} = \text{same string with }
 \Sigma_R(k) = \Sigma_{\rm loop}(k) - \delta m
 - \left(\slashed{k}-m\right)\delta Z_2 .$$

The removed $\delta Z_2$ piece **collapses** exactly (this is the
mechanism behind the assembly theorem): since
$\tilde S(k)\left(\slashed{k}-m\right) = i$,

$$\tilde S(k)\left[-i\left(\slashed{k}-m\right)\delta Z_2\right]\tilde S(k)
 = \delta Z_2\;\tilde S(k),$$

an insertion of the $\delta Z_2$ term on an internal line just
multiplies the LO diagram by $\delta Z_2$. **Mirror** (SE on the $k'$
line): factor 2.
$$\mu_{\rm IId} = 2\times\text{projection with }\Sigma_R
 = \int f_{\rm rat} + \int f_{\rm log}\ \text{of the NLO section.}$$

### 8.5 Diagram IIe (vacuum polarization insertion)

![Diagram IIe](figures/g2-nlo-IIe.svg)

The VP loop of §7.3 between two photon propagators of the LO string
(loop $(-1)$ and trace already inside $i\Pi^{\mu\nu}$):

$$\delta\Gamma^\mu_{\rm IIe} = e^2
 \int\frac{\mathrm{d}^4k}{(2\pi)^4}\;
 \frac{\gamma^\alpha\left(\slashed{k}'+m\right)\gamma^\mu
 \left(\slashed{k}+m\right)\gamma^\beta}
 {\left[k'^2-m^2\right]\left[k^2-m^2\right]}\;
 \frac{-ig_{\alpha\rho}}{k_\gamma^2-\lambda^2}\;
 i\Pi^{\rho\sigma}(k_\gamma)\;
 \frac{-ig_{\sigma\beta}}{k_\gamma^2-\lambda^2},
 \qquad k_\gamma = k-p$$

(prefactor $(-ie)^2 i^2 = e^2$ from the two outer vertices and two
fermion propagators; the photon factors are displayed). Renormalization
replaces $\Pi \to \hat\Pi$ (the $\Pi(0)$ part is the charge
renormalization of the LO diagram, §3.3/§9). Because
$\hat\Pi(k_\gamma^2) \sim k_\gamma^2$ as $k_\gamma^2\to0$, the diagram
is IR finite and $\lambda$ can be set to zero here. Using
transversality and the spectral representation (§7.3),

$$\frac{-ig_{\alpha\rho}}{k_\gamma^2}\;
 i\left(k_\gamma^2 g^{\rho\sigma} - k_\gamma^\rho k_\gamma^\sigma\right)
 \hat\Pi(k_\gamma^2)\;
 \frac{-ig_{\sigma\beta}}{k_\gamma^2}
 = -ig_{\alpha\beta}\,\frac{\hat\Pi(k_\gamma^2)}{k_\gamma^2}
 + \left(k_{\gamma\alpha}k_{\gamma\beta}\ \text{terms}
   \Rightarrow F_1\text{-structures only, Sec. 2.4}\right)$$

$$= -ig_{\alpha\beta}\int_{4m^2}^\infty\mathrm{d}t\,
 \frac{\rho(t)}{k_\gamma^2-t} + \dots,$$

i.e. **diagram IIe = the LO diagram with a photon of mass$^2$ $t$,
folded with $\rho(t)$**. Since §7.1 gives that diagram's $F_2(0)$ as
$(\alpha/\pi)K(t)$ *linearly in the photon propagator*,

$$F_2^{\rm IIe}(0) = \int_{4m^2}^\infty\mathrm{d}t\;\rho(t)\,
 \frac{\alpha}{\pi}\,K(t)
 \quad\Longrightarrow\quad
 \mu_{\rm IIe} = \int_4^\infty\mathrm{d}t\,\frac{1}{3t}
 \left(1+\frac{2}{t}\right)\sqrt{1-\frac{4}{t}}\;K(t)
 \qquad(m=1),$$

exactly the NLO section's starting integral (evaluated there to
$\frac{119}{36}-\frac{\pi^2}{3}$).

## 9. Renormalization and the assembly theorem

Collecting §3.3, §5, §7 and §8, with everything expressed in the
physical $m$ (via $\delta m$) and physical $\alpha$ (via $\hat\Pi$), the
ab-initio order-$\alpha^2$ anomalous moment is

$$a_e\big|_{\alpha^2} = \Big\{F_2\left[\mathrm{I}\right]
 + F_2\left[\mathrm{IIa}\right] + 2F_2\left[\mathrm{IIc}\right]
 + 2F_2\left[\mathrm{IId};\,\Sigma-\delta m\right]
 + F_2\left[\mathrm{IIe};\,\hat\Pi\right]\Big\}
 \;+\; \delta Z_2\,F_2^{(2)}(0),$$

every term a fully specified integral. The NLO section instead
computes the KK-subtracted per-diagram contributions (write $a_{\rm X}$
for the absolute $F_2(0)$ contribution; the NLO section's numbers are
$\mu_{\rm X} = a_{\rm X}/(\alpha/\pi)^2$):

$$a_{\rm IIa} = F_2\left[\mathrm{IIa}\right] - \delta F_1(0)\,F_2^{(2)},\qquad
 a_{\rm IIc} = 2\left(F_2\left[\mathrm{IIc}\right]
 - \delta F_1(0)\,F_2^{(2)}\right),$$

$$a_{\rm IId} = 2\left(F_2\left[\mathrm{IId};\,\Sigma-\delta m\right]
 - \delta Z_2\,F_2^{(2)}\right),\qquad
 a_{\rm I} = F_2\left[\mathrm{I}\right],\qquad
 a_{\rm IIe} = F_2\left[\mathrm{IIe};\,\hat\Pi\right]$$

(§8.2–8.4: the pointwise vertex subtractions integrate to
$\delta F_1(0)F_2^{(2)}$ per insertion, and the on-shell $\delta Z_2$
subtraction collapses to $\delta Z_2F_2^{(2)}$ per line). Summing and
comparing with the master formula:

$$\sum_{\rm X}a_{\rm X}
 = a_e\big|_{\alpha^2}
 - \delta Z_2F_2^{(2)}
 - 3\,\delta F_1(0)\,F_2^{(2)} - 2\,\delta Z_2\,F_2^{(2)}
 = a_e\big|_{\alpha^2}
 - 3\left(\delta F_1(0) + \delta Z_2\right)F_2^{(2)} .$$

By the Ward identity (§4.3) $\delta F_1(0) = -\delta Z_2$, so the last
term vanishes **exactly**:

$$\boxed{\;\mu_{\rm I} + \mu_{\rm IIa} + \mu_{\rm IIc}
 + \mu_{\rm IId} + \mu_{\rm IIe}
 = \frac{a_e\big|_{\alpha^2}}{(\alpha/\pi)^2} = A_2\;}$$

This is the precise sense in which the five scheme-dependent
per-diagram integrals of the NLO section — three vertex subtractions
but only *one* net LSZ cross term, reconciled by $Z_1 = Z_2$ — sum to
the scheme-independent $A_2$, and it derives the NLO section's
renormalization bookkeeping (its "vertex renormalization cross term"
and "caveat on per-diagram comparisons") rather than assuming it. The
infrared side is organized the same way: the $\log\lambda$ of
$\delta F_1(0;\lambda)$ enters IIa and IIc through their subtractions
and $\delta Z_2(\lambda)$ enters IId through its subtraction; the NLO
section verifies explicitly that only the advertised
$\mp\frac12\log(\lambda^2/m^2)$ of IIc/IId survives and cancels in the
sum.

## 10. Summary: what remains to be computed

The Green function
$\langle\Omega|T\psi A^\mu\bar\psi|\Omega\rangle$, expanded by
Gell-Mann–Low + Dyson + Wick, contains at order $e^5$ exactly 10800
complete contractions in 16 topology classes with total weight
$10800/5! = 90$ (`code/wick_enum.py`). Of these:

* disconnected, forward, and tadpole classes are zero or cancel
  against $\langle0|S|0\rangle$;
* class V vanishes by Furry's theorem (proved via charge conjugation);
* external-leg classes are absorbed by LSZ into $Z_2$;
* external-potential classes III/IV are absorbed by on-shell charge
  renormalization ($\hat\Pi(0)=0$);
* the remaining **five 1PI cores** — with the on-shell/KK subtractions
  derived in §7–§9 and the projector of §2.3 — are the integrals
  $\mu_{\rm I}, \mu_{\rm IIa}, \mu_{\rm IIc}, \mu_{\rm IId},
  \mu_{\rm IIe}$, and their sum is exactly $A_2$.

Their evaluation — Feynman-parameter reduction, Fortran quadrature,
and sequential analytic integration — is the subject of the NLO
section: the massive-photon kernel $K(t)$ (§7.1 here = Stage 0 there),
$\mu_{\rm IIe}$ (§8.5 = its dispersion integral), $\mu_{\rm IId}$
(§7.2/§8.4 = its $f_{\rm rat}/f_{\rm log}$), $\mu_{\rm IIc}$
(§7.4/§8.3 = its pieces (a),(b),(c)), $\mu_{\rm IIa}$ (§7.4/§8.2, same
machinery, no mirror), and $\mu_{\rm I}$ (§8.1, the genuine two-loop
integral). Nothing in that section's list of integrals falls outside
this derivation, and nothing derived here requires an integral not
listed there.
