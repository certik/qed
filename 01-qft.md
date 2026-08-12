# Standard Model

# Introduction

The aim of these (work in progress) notes is to use the Standard Model
of particle physics to derive all equations in quantum mechanics (and
quantum field theory) that we need for our research.

We start by deriving the electroweak Standard Model from the
$SU(2)\times U(1)$ symmetry and couple other (standard) assumptions in
the quantum field theory. After that, we only want to derive things and
make nonrelativistic limits or other approximations in order to derive
everything else in quantum mechanics. In particular we show how to
derive the Dirac and Schrödinger equations (as a low energy limit). We
then show some particular ways to solve those equations, like
perturbation theory, scattering theory, \...

The goal is to have a complete theory on about 30 or 40 pages and then
lots of examples (arbitrarily long), that use the theory (but do not
develop new ideas), so that one can learn how the theory works from the
examples. For instance, one can ask \"why is there the term $(\mathbf{p}-e\mathbf{A})^2$ in the Schrödinger equation for electromagnetic field,
why this and not something else?\" or \"why is there the
$\boldsymbol{\sigma}\cdot\mathbf{B}$ term in the Pauli equation?\", to find the
answer, one just finds the Pauli equation in the theory and then looks
at the derivation, so in this case one quickly finds that it follows
from the minimal coupling in QED, e.g. it's the easiest way how
electron-foton interaction can be coupled, e.g. the $U(1)$ symmetry.
Nice thing about QFT is that one can find really nice geometrical
reasons why things are that way and not some other way (just open any
advance book on QFT), but the problem is that basically nowhere is some
easy (but correct) translation of those results to regular QM, so that
everything fits into just couple dozens pages, so that it can serve as a
reference.

The advantage of this top-down approach is that it is easy to see where
things come from and also to understand exactly what approximations one
is using when dealing with any equation in QM. However, as is well-known
in physics, to be a good physicist one has to understand all the
approaches, e.g. both top-down and bottom-up and all other approaches to
QM and QFT, because there are no two approaches that would be 100%
equivalent, so one has to use the right approach for the particular
problem. So these notes do not aspire to be the right way to teach QM,
but rather to serve as a reference to get quickly oriented and to find
the equations to start from.

::: index
Standard Model
:::

# Standard Model

::: index
Electroweak Standard Model
:::

## Electroweak Standard Model

Lagrangian with a global $SU(2)\times U(1)$ symmetry:

$$\mathcal{L}=i\bar L^{(l)}\gamma_\mu\partial^\mu L^{(l)}+i\bar l_R \gamma_\mu\partial^\mu l_R +\frac{1}{2} \partial_\mu\Phi^*\partial^\mu\Phi-m^2\Phi^*\Phi-{1\over4} \lambda(\Phi^*\Phi)^2 -h_e\bar L^{(l)} \Phi e_R - \text{h.c.}$$

where $l=e,\mu,\tau$ and $a=1,2$, $l_{L,R} = \frac{1}{2}(1\mp\gamma_5)l$
and

$$\begin{aligned}
L^{(l)} = \left( \begin{array}{c} \nu_{(l)L} \\ l_L \end{array} \right)
\end{aligned}$$

Local $SU(2)\times U(1)$ symmetry: This consists of two things. First
changing the partial derivatives to covariant ones:

$$\partial^\mu \to D^\mu =\partial^\mu-{i\over2}g\tau_k A_k^\mu - {i\over2}g'YB^\mu$$

and second adding the kinetic terms

$$-{1\over4}F^a_{\mu\nu}F^{a\mu\nu}-{1\over4}B_{\mu\nu}B^{\mu\nu}$$

of the vector gauge particles to the lagrangian.

$$F^a_{\mu\nu} = \partial_\mu A^a_\nu-\partial_\nu A^a_\mu+ g\epsilon^{abc}A^b_\mu A^c_\nu$$

$$B_{\mu\nu} = \partial_\mu B_\nu-\partial_\nu B_\mu$$

$$\begin{aligned}
\Phi = e^{{i\over v}\pi^a(x)\tau^a} \left( \begin{array}{c} 0 \\ {1\over\sqrt{2}}(v+H(x)) \end{array} \right)
\end{aligned}$$

This breaks the gauge invariance. The $\partial^\mu\pi^a$ are going to
be added to $A^a_\mu$ so we can set $\pi_a = 0$ now.

::: index
Higgs boson
:::

### Higgs Terms

$$\mathcal{L}_{Higgs}= \frac{1}{2} \partial_\mu\Phi^*\partial^\mu\Phi-m^2\Phi^*\Phi-{1\over4} \lambda(\Phi^*\Phi)^2$$

Plugging in the covariant derivatives and $\Phi$ in U-gauge (symmetry
breaking):

$$\mathcal{L}_{Higgs} = {1\over2}\Phi^+(\overleftarrow{\partial}_\mu+igA^a_\mu {\tau^a\over 2} + ig'YB_\mu) (\overrightarrow{\partial}^\mu+igA^{a\mu} {\tau^a\over 2} + ig'YB^\mu)\Phi -\lambda(\Phi^+\Phi-{v^2\over2})^2=$$

$$= \Phi^+_U(\overleftarrow{\partial}_\mu+igA^a_\mu {\tau^a\over 2} + ig'YB_\mu) (\overrightarrow{\partial}^\mu+igA^{a\mu} {\tau^a\over 2} + ig'YB\mu)\Phi_U -\lambda(\Phi^+_U\Phi_U-{v^2\over2})^2 =$$

$$= {1\over2}\partial_\mu H\partial^\mu H - \lambda v^2 H^2 - \lambda v H^3 - {1\over 4}\lambda H^4 +$$

$$+{1\over 8}(v+H)^2 \left(2g^2{A^1_\mu+iA^2_\mu\over\sqrt{2}}{A^{1\mu}-iA^{2\mu}\over\sqrt{2}} + (g^2+4Y^2g'^2){gA^3_\mu-2Yg'B_\mu\over\sqrt{g^2+4Y^2g'^2}} {gA^{3\mu}-2Yg'B^\mu\over\sqrt{g^2+4Y^2g'^2}}\right) =$$

$$= {1\over2}\partial_\mu H\partial^\mu H - \lambda v^2 H^2 - \lambda v H^3 - {1\over 4}\lambda H^4 + {1\over 8}(v+H)^2 \left(2g^2W^-_\mu W^{+\mu} + {g^2\over\cos^2\theta_W}Z_\mu Z^\mu\right) =$$

$$= {1\over2}\partial_\mu H\partial^\mu H - \lambda v^2 H^2 +{1\over4}g^2v^2W^-_\mu W^{+\mu}+{g^2v^2\over8\cos^2\theta_W}Z_\mu Z^\mu - \lambda v H^3 - {1\over 4}\lambda H^4 +$$

$$+{1\over2}vg^2W_\mu^-W^{+\mu}H +{g^2\over4\cos\theta_W}vZ_\mu Z^\mu H +{1\over4}g^2W_\mu^-W^{+\mu}H^2 +{g^2\over8\cos\theta_W}Z_\mu Z^\mu H^2$$

Where we put

$$W^{\pm}_\mu = {1\over\sqrt{2}}(A^1_\mu \mp iA^2_\mu)$$

$$Z_\mu = {g\over\sqrt{g^2+4Y^2g'^2}}A^3_\mu- {2Yg'\over\sqrt{g^2+4Y^2g'^2}}B_\mu$$

we defined $\theta_W$ by the relation

$$\cos\theta_W = {g\over\sqrt{g^2+4Y^2g'^2}}$$

so that the expressions simplify a bit, e.g. we now get:

$$\sin\theta_W = {2Yg'\over\sqrt{g^2+4Y^2g'^2}}$$

$$Z_\mu= \cos\theta_W A^3_\mu - \sin\theta_W B_\mu$$

$$g^2+4Y^2g'^2 = {g^2\over\cos^2\theta_W}$$

::: index
Yukawa terms
:::

### Yukawa terms

$$\mathcal{L}_{Yukawa} = -h_e \bar L \Phi e_R - \text{h.c.}= -h_e \bar L \Phi_U e_R - \text{h.c.}=$$

$$=-{1\over\sqrt{2}}h_e(v+H)(\bar e_L e_R + \bar e_R e_L)= -{1\over\sqrt{2}}h_e(v+H)\bar ee=$$

$$=-{1\over\sqrt{2}}h_ev\bar ee-{1\over\sqrt{2}}h_e\bar eeH$$

The term $\bar L \Phi e_R$ is $U(1)$ (hypercharge) invariant, so

$$-Y_L+Y+Y_R=0$$

::: index
leptons
:::

### Leptonic Terms

$$\mathcal{L}=i\bar L\gamma^\mu\partial_\mu L+i\bar e_R \gamma^\mu\partial_\mu e_R \to$$

$$\to i\bar L\gamma^\mu(\partial_\mu-igA^a_\mu{\tau^a\over2}-ig'Y_LB_\mu) L +i\bar e_R \gamma^\mu(\partial_\mu-ig'Y_RB_\mu) e_R =$$

$$= i\bar L\gamma^\mu\partial_\mu L+i\bar e_R \gamma^\mu\partial_\mu e_R +g\bar L\gamma^\mu{\tau^a\over2}LA^a_\mu +g'Y_L\bar L\gamma^\mu LB_\mu +g'Y_R\bar e_R \gamma^\mu e_R B_\mu =$$

$$= i\bar L\gamma^\mu\partial_\mu L+i\bar e_R \gamma^\mu\partial_\mu e_R +{g\over\sqrt{2}}(\bar \nu_L\gamma^\mu e_L W^+_\mu + \text{h.c.}) +{1\over2}g\bar L\gamma^\mu\tau^3L A^3_\mu +g'Y_L\bar L\gamma^\mu LB_\mu +g'Y_R\bar e_R \gamma^\mu e_R B_\mu =$$

$$= i\bar \nu_L\gamma^\mu\partial_\mu \nu_L+i\bar e \gamma^\mu\partial_\mu e +{g\over\sqrt{2}}(\bar \nu_L\gamma^\mu e_L W^+_\mu + \text{h.c.}) +{1\over2}g\bar\nu_L\gamma^\mu\nu_LA^3_\mu -{1\over2}g\bar e_L\gamma^\mu e_LA^3_\mu$$

$$+g'Y_L\bar\nu_L\gamma^\mu\nu_LB_\mu +g'Y_L\bar e_L\gamma^\mu e_LB_\mu +g'Y_R\bar e_R \gamma^\mu e_R B_\mu =$$

$$= i\bar \nu_L\gamma^\mu\partial_\mu \nu_L+i\bar e \gamma^\mu\partial_\mu e +{g\over\sqrt{2}}(\bar \nu_L\gamma^\mu e_L W^+_\mu + \text{h.c.})$$

$$+\left[ (\frac{1}{2} g\sin\theta_W+Y_Lg'\cos\theta_W)\bar\nu_L\gamma^\mu\nu_L +(-\frac{1}{2} g\sin\theta_W +Y_Lg'\cos\theta_W)\bar e_L\gamma^\mu e_L +Y_Rg'\cos\theta_W\bar e_R\gamma^\mu e_R \right]A_\mu$$

$$+\left[ (\frac{1}{2} g\cos\theta_W-Y_Lg'\sin\theta_W)\bar\nu_L\gamma^\mu\nu_L +(-\frac{1}{2} g\cos\theta_W -Y_Lg'\sin\theta_W)\bar e_L\gamma^\mu e_L -2Y_Lg'\sin\theta_W\bar e_R\gamma^\mu e_R \right]Z_\mu$$

Where we substituted new fields $Z_\mu$ and $A_\mu$ for the old ones
$A^3_\mu$ and $B_\mu$ using the relation:

$$Z_\mu= \cos\theta_W A^3_\mu - \sin\theta_W B_\mu$$

$$A_\mu= \sin\theta_W A^3_\mu + \cos\theta_W B_\mu$$

The angle $\theta_W$ must be the same as in the Higgs sector, so that
the field $Z_\mu$ is the same. We now need to make the following
requirement in order to proceed further:

$$Y = -Y_L$$

This follows for example by requiring that neutrinos have zero charge,
i.e. setting $\frac{1}{2} g\sin\theta_W + Y_L g'\cos\theta_W = 0$ and substituting
for $\theta_W$ from the definition (see the Higgs terms), from which
one gets $Y=-Y_L$. From $-Y_L+Y+Y_R=0$ we now get

$$Y_R = 2Y_L$$

it now follows:

$$\frac{1}{2} g\sin\theta_W+Y_Lg'\cos\theta_W = 0$$

$$-\frac{1}{2} g\sin\theta_W +Y_Lg'\cos\theta_W = -g\sin\theta_W$$

$$Y_Rg'\cos\theta_W = -g\sin\theta_W$$

$$\tan\theta_W = -2Y_L {g'\over g}$$

and the Lagrangian can be further simplified:

$$\mathcal{L}= i\bar\nu_L\gamma^\mu\partial_\mu\nu_L+i\bar e\gamma^\mu\partial_\mu e +{g\over\sqrt{2}}(\bar \nu_L\gamma^\mu e_L W^+_\mu + \text{h.c.})$$

$$-g\sin\theta_W(\bar e_L\gamma^\mu e_L+\bar e_R\gamma^\mu e_R) A_\mu$$

$$+{g\over\cos\theta_W}\left[ \frac{1}{2} \bar\nu_L\gamma^\mu\nu_L +(-\frac{1}{2} + \sin^2\theta_W)\bar e_L\gamma^\mu e_L +\sin^2\theta_W\bar e_R\gamma^\mu e_R \right]Z_\mu=$$

$$= i\bar\nu_L\gamma^\mu\partial_\mu\nu_L+i\bar e \gamma^\mu\partial_\mu e +{g\over2\sqrt{2}}(\bar \nu\gamma^\mu (1-\gamma_5) e W^+_\mu + \text{h.c.}) -g\sin\theta_W\bar e\gamma^\mu e A_\mu$$

$$+{g\over2\cos\theta_W}\left[ \bar\nu\gamma^\mu(1-\gamma_5)\nu +\bar e\gamma^\mu (-\frac{1}{2}+2\sin^2\theta_W+\frac{1}{2}\gamma_5) e \right]Z_\mu$$

Where we used the relations $\bar\nu_L\gamma^\mu e_L=\frac{1}{2}\bar\nu\gamma^\mu(1-\gamma_5)e$ and $\bar\nu_R\gamma^\mu e_R=\frac{1}{2}\bar\nu\gamma^\mu(1+\gamma_5)e$ .

::: index
gauge
:::

### Gauge terms

$$\mathcal{L}_{Gauge} = -{1\over4}F^a_{\mu\nu}F^{a\mu\nu} -{1\over4}B_{\mu\nu}B^{\mu\nu}=$$

$$= -{1\over4}(\partial_\mu A^a_\nu-\partial_\nu A^a_\mu+g\epsilon^{abc} A^b_\mu A^c_\nu)(\partial^\mu A^{a\nu}-\partial^\nu A^{a\mu}+g\epsilon^{ajk} A^{j\mu} A^{k\nu}) -{1\over4}B_{\mu\nu}B^{\mu\nu}=$$

$$= -{1\over4}\partial_\mu A^a_\nu\partial^\mu A^{a\nu} -{1\over4}B_{\mu\nu}B^{\mu\nu} -{1\over2}(\partial_\mu A^a_\nu-\partial_\nu A^a_\mu)g\epsilon^{abc} A^{b\mu} A^{c\nu} -{1\over4}g^2\epsilon^{abc}\epsilon^{ajk}A^b_\mu A^c_\nu A^{k\mu} A^{l\nu} =$$

$$= -{1\over2}W^-_{\mu\nu}W^{+\mu\nu} -{1\over4}A_{\mu\nu}A^{\mu\nu} -{1\over4}Z_{\mu\nu}Z^{\mu\nu} -g[(\partial_\mu A^1_\nu-\partial_\nu A^1_\mu)A^{2\mu}A^{3\nu}+ \text{cycl. perm. (123)}]$$

$$-{1\over4}g^2[(A^a_\mu A^{a\mu})(A^b_\nu A^{b\nu})- (A^a_\mu A^a_\nu)(A^{b\mu} A^{b\nu})]=$$

$$= -{1\over2}W^-_{\mu\nu}W^{+\mu\nu} -{1\over4}A_{\mu\nu}A^{\mu\nu} -{1\over4}Z_{\mu\nu}Z^{\mu\nu} -g[A^1_\mu A^2_\nu \overleftrightarrow{\partial}^\mu A^{3\nu}+ \text{cycl. perm. (123)}]$$

$$-{1\over4}g^2[(A^a_\mu A^{a\mu})(A^b_\nu A^{b\nu})- (A^a_\mu A^a_\nu)(A^{b\mu} A^{b\nu})] =$$

$$= -{1\over2}W^-_{\mu\nu}W^{+\mu\nu} -{1\over4}A_{\mu\nu}A^{\mu\nu} -{1\over4}Z_{\mu\nu}Z^{\mu\nu} -ig(W^0_\mu W^-_\nu\overleftrightarrow{\partial}^\mu W^{+\nu}+ \text{cycl. perm. (0-+)})$$

$$-g^2[ \frac{1}{2}(W^+_\mu W^{-\mu})^2 -\frac{1}{2}(W^+_\mu W^{+\mu})(W^-_\nu W^{-\nu}) +(W^0_\mu W^{0\mu})(W^+_\nu W^{-\nu}) -(W^-_\mu W^+_\nu)(W^{0\mu} W^{0\nu})=$$

$$= -{1\over2}W^-_{\mu\nu}W^{+\mu\nu} -{1\over4}A_{\mu\nu}A^{\mu\nu} -{1\over4}Z_{\mu\nu}Z^{\mu\nu} +\mathcal{L}_{WW\gamma}+L_{WWZ}+L_{WW\gamma\gamma}+L_{WWWW}+L_{WWZZ}+L_{WWZ\gamma}$$

Where $W^0_\mu = A^3_\mu=\cos\theta_W Z_\mu + \sin\theta_W A_\mu$ and:

$$\mathcal{L}_{WW\gamma}=-ig\sin\theta_W(A_\mu W^-_\nu\overleftrightarrow{\partial}^\mu W^{+\nu} + \text{cycl. perm. ($A$ $W^-$ $W^+$)})$$

$$\mathcal{L}_{WWZ}=-ig\cos\theta_W(Z_\mu W^-_\nu\overleftrightarrow{\partial}^\mu W^{+\nu}+\text{cycl. perm. ($Z$ $W^-$ $W^+$)})$$

$$\mathcal{L}_{WW\gamma\gamma}=-g^2\sin^2\theta_W(W^-_\mu W^{+\mu}A_\nu A^\nu- W^-_\mu A^\mu W^+_\nu A^\nu)$$

$$\mathcal{L}_{WWWW}=\frac{1}{2} g^2(W^-_\mu W^{-\mu}W^+_\nu W^{+\nu} -W^-_\mu W^{+\mu}W^-_\nu W^{+\nu} )$$

$$\mathcal{L}_{WWZZ}=-g^2\cos^2\theta_W(W^-_\mu W^{+\mu}Z_\nu Z^\nu -W^-_\mu Z^\mu W^+_\nu Z^\nu )$$

$$\mathcal{L}_{WWZ\gamma}=g^2\sin\theta_W\cos\theta_W(-2W^-_\mu W^{+\mu}A_\nu Z^\nu+W^-_\mu Z^\mu W^+_\nu A^\nu+W^-_\mu A^\mu W^+_\nu Z^\nu)$$

::: index
pair: GWS; Lagrangian
:::

### GWS Lagrangian

Plugging everything together we get the GWS Lagrangian:

$$\mathcal{L} = {1\over2}\partial_\mu H\partial^\mu H - \lambda v^2 H^2 +{1\over4}g^2v^2W^-_\mu W^{+\mu}+{g^2v^2\over8\cos^2\theta_W}Z_\mu Z^\mu - \lambda v H^3 - {1\over 4}\lambda H^4 +$$

$$+{1\over2}vg^2W_\mu^-W^{+\mu}H +{g^2\over4\cos\theta_W}vZ_\mu Z^\mu H +{1\over4}g^2W_\mu^-W^{+\mu}H^2 +{g^2\over8\cos\theta_W}Z_\mu Z^\mu H^2$$

$$-{1\over\sqrt{2}}h_ev\bar ee-{1\over\sqrt{2}}h_e\bar eeH$$

$$-{1\over2}W^-_{\mu\nu}W^{+\mu\nu} -{1\over4}A_{\mu\nu}A^{\mu\nu} -{1\over4}Z_{\mu\nu}Z^{\mu\nu} +\mathcal{L}_{WW\gamma}+L_{WWZ}+L_{WW\gamma\gamma}+L_{WWWW}+L_{WWZZ}+L_{WWZ\gamma}$$

$$+i\bar\nu_L\gamma^\mu\partial_\mu\nu_L+i\bar e \gamma^\mu\partial_\mu e +{g\over2\sqrt{2}}(\bar \nu\gamma^\mu (1-\gamma_5) e W^+_\mu + \text{h.c.}) -g\sin\theta_W\bar e\gamma^\mu e A_\mu$$

$$+{g\over2\cos\theta_W}\left[ \bar\nu\gamma^\mu(1-\gamma_5)\nu +\bar e\gamma^\mu (-\frac{1}{2}+2\sin^2\theta_W+\frac{1}{2}\gamma_5) e \right]Z_\mu$$

$$+ (e, \nu_e, h_e \leftrightarrow \mu, \nu_\mu, h_\mu) + (e, \nu_e, h_e \leftrightarrow \tau, \nu_\tau, h_\tau)$$

The free parameters are $g$, $\theta_W$, $v$, $\lambda$, $h_e$,
$h_\mu$, $h_\tau$.

::: index
particle mass
:::

### Particle Masses

The particle masses are deduced from the terms

$$\mathcal{L} = -{1\over2}m_H^2 H^2 +m_W^2 W^-_\mu W^{+\mu} +{1\over2}m_Z^2 Z_\mu Z^\mu -m_e\bar ee +\cdots$$

comparing to the above:

$$\mathcal{L} = -\lambda v^2 H^2 +{1\over4}g^2v^2W^-_\mu W^{+\mu} +{g^2v^2\over8\cos^2\theta_W}Z_\mu Z^\mu -{1\over\sqrt{2}}h_ev\bar ee +\cdots$$

we get

$$m_W = \frac{1}{2} g v$$

$$m_Z = {gv\over2\cos\theta_W}={m_W\over\cos\theta_W}$$

$$m_H = v\sqrt{2\lambda}$$

$$m_e = {1\over\sqrt{2}}h_ev$$

$$m_\mu = {1\over\sqrt{2}}h_\mu v$$

$$m_\tau = {1\over\sqrt{2}}h_\tau v$$

Note that those are the bare masses (e.g. in order to obtain the real
mesaured masses of the particles, one has to renormalize them by
calculating the higher order corrections given by the loop diagrams).

::: index
quarks
:::

### Parameters of the Standard Model

The free parameters are $g$, $\theta_W$, $v$, $\lambda$, then
three masses of the charged leptons $h_e$, $h_\mu$, $h_\tau$, six
quark masses and four parameters of the CKM mixing matrix, which gives
4 + 3 + 6 + 4 = 17 free parameters (if one allows for three neutrino
masses and the corresponding four mixings parameters, one gets 17 + 3 +
4 = 24 free parameters).

They can be traded for other physical parameters (see below), but their
numerical values are not predicted by the theory, so they have to be
measured and their experimental values are approximately:

$$g = 0.631$$

$$\theta_W = 28.67^\circ$$

$$v = 246.218 {\rm\,GeV}$$

$$0.2 < \lambda < 4.0$$

$$h_e = 2.929\cdot 10^{-6} {\rm\,eV}$$

$$h_\mu = 6.065\cdot 10^{-4} {\rm\,eV}$$

$$h_\tau = 1.021\cdot 10^{-2} {\rm\,eV}$$

All the parameters have been measured quite exactly, except $\lambda$.

Other physical constants can then be calculated using the formulas:

$$m_W = \frac{1}{2} g v = 77.7 {\rm\, GeV}$$

$$m_Z = {m_W\over\cos\theta_W} = 88.6 {\rm\, GeV}$$

$$m_H = v\sqrt{2\lambda} = \mbox{from }150 {\rm\,GeV}\mbox{ to }700 {\rm\,GeV}$$

$$m_e = {1\over\sqrt{2}}h_ev = 511{\rm\,KeV}$$

$$m_\mu= {1\over\sqrt{2}}h_\mu v = 105.6{\rm\,MeV}$$

$$m_\tau= {1\over\sqrt{2}}h_\tau v = 1.777{\rm\,GeV}$$

$$G_F = {1\over\sqrt{2} v^2} = (1.16639 \pm 0.00001) \times 10^{-5}
    {\rm\, GeV^{-2}}$$

$$e = g \sin\theta_W = 0.3$$

$$\alpha = {1\over 4\pi} g^2 \sin^2 \theta_W \doteq {1\over 137}$$

Code:

    >>> from math import pi, sin, cos, sqrt
    >>> eV = 1
    >>> KeV = 1e3
    >>> MeV = 1e6
    >>> GeV = 1e9
    >>> g = 0.631
    >>> theta_W = 28.67 * pi / 180
    >>> v = 246.218 * GeV
    >>> h_e = 2.935 * 1e-6 * eV
    >>> h_mu = 6.065 * 1e-4 * eV
    >>> h_tau = 1.021 * 1e-2 * eV
    >>> g*v/2 / GeV
    77.681779
    >>> g*v/2/cos(theta_W) / GeV
    88.5365869768
    >>> h_e * v / sqrt(2) / KeV
    510.99059521630568
    >>> h_mu * v / sqrt(2) / MeV
    105.59311618353983
    >>> h_tau * v / sqrt(2) / GeV
    1.7775856821664329
    >>> 1./sqrt(2)/v**2 / (1e-5 * GeV**-2)
    1.1663943402665491
    >>> g*sin(theta_W)
    0.30273118431564783
    >>> 1. / (g**2*sin(theta_W)**2/(4*pi))
    137.11833915409719

### Quarks

$$\mathcal{L}_{fermion}+\!\!= \sum_{q=d,s,b}i\bar L_0^{(q)}\gamma^\mu\partial_\mu L_0^{(q)} +\sum_{q=d,u,s,c,b,t}i\bar q_{0R}\gamma^\mu\partial_\mu q_{0R}$$

$$\mathcal{L}_{Yukawa}+\!\!= -\sum_{q=d,s,b\atop q'=d,s,b}h_{qq'}i\bar L_0^{(q)}\Phi q_{0R}'+\text{h.c.} -\sum_{q=d,s,b\atop q'=u,c,t}\tilde h_{qq'}i\bar L_0^{(q)}\tilde\Phi q_{0R}'+\text{h.c.}$$
