# Magnetic moment of an electron

In this section we derive the order-$\alpha$ correction to the magnetic
moment of an electron.

We start by computing the electron vertex function for the process
$\gamma(q)\to e^+(p) + e^-(p')$:

$$i M = i e^2 \left(\bar u(p')\Gamma^\mu(p', p)u(p)\right) {1\over q^2}
    \left(\bar u(k') \gamma_\mu u(k)\right)$$

where $k$ corresponds to some heavy target. If $A_\mu^{\mathrm{cl}}$ is a
fixed classical potential, we get:

$$i M 2\pi \delta(p^{0'} - p^0) =
    -i e \bar u(p')\Gamma^\mu(p', p)u(p) A_\mu^\mathrm{cl}$$

Using general arguments (Lorentz invariance, parity-conservation, Ward
identity) we can always write $\Gamma^\mu$ as:

$$\Gamma^\mu(p', p) = \gamma^\mu F_1(q^2) + {i\sigma^{\mu\nu} q_\nu \over
    2m} F_2(q^2)$$

where $F_1$ and $F_2$ ar unknown functions of $q^2 = (p'-p)^2 = -2p'\cdot p + 2m^2$ called form factors. As we will see below, in the
lowest order we get $F_1 = 1$ and $F_2 = 0$.

We can calculate the amplitude for elastic Coulomb scattering of a
nonrelativistic electron from a region of nonzero electrostatic
potential by setting $A_\mu^{\mathrm{cl}}(x)=(\phi(\mathbf{x}), 0)$, then:

$$A_\mu^\mathrm{cl}(q)=(2\pi\delta(q^0)\tilde\phi({\mathbf{q}}), 0)$$

$$i M 2\pi \delta(p^{0'} - p^0) =
    -i e \bar u(p')\Gamma^0(p', p)u(p) 2\pi\delta(q^0)\tilde\phi({\mathbf{q}})$$

$$i M = -i e \bar u(p')\Gamma^0(p', p)u(p) \tilde\phi({\mathbf{q}})$$

If the electrostatic field is very slowly varying over a large (even
macroscopic) region, $\tilde\phi(\mathbf{q})$ will be concentrated about
$\mathbf{q} = 0$, then we can take the limit $\mathbf{q}\to 0$:

$$i M = -i e \bar u(p')\Gamma^0(p', p)u(p) \tilde\phi({\mathbf{q}})$$

$$i M = -i e \bar u(p')
    \left(\gamma^0 F_1(q^2) + {i\sigma^{0\nu} q_\nu \over 2m} F_2(q^2)
        \right)u(p) \tilde\phi({\mathbf{q}})$$

$$i M = -i e \bar u(p') \gamma^0 u(p) F_1(0) \tilde\phi({\mathbf{q}})$$

$$i M = -i e 2m\xi^{'\dag}\xi F_1(0) \tilde\phi({\mathbf{q}})$$

$$i M = -i \left( e F_1(0) \tilde\phi({\mathbf{q}})\right) 2m\xi^{'\dag}\xi$$

This corresponds to the Born approximation for scattering from a
potential

$$V({\mathbf{x}}) = e F_1(0) \phi({\mathbf{x}})$$

Thus $F_1(0)$ is the electric charge of the electron, in units of
$e$. Since $F_1(0) = 1$ already in the first order of perturbation
theory, radiative corrections to $F_1(q^2)$ must vanish at
$q^2=0$.

Now we calculate the scattering from a static vector potential by
setting $A_\mu^{\mathrm{cl}}(x)=(0, \mathbf{A}^{\mathrm{cl}}(\mathbf{x}))$, then:

$$A^\mu_\mathrm{cl}(q)=(0, 2\pi\delta(q^i)\tilde A^i_\mathrm{cl}({\mathbf{q}}))$$

$$i M 2\pi \delta(p^{'i} - p^i) =
    i e \bar u(p')\Gamma^i(p', p)u(p) 2\pi\delta(q^i)\tilde A^i_\mathrm{cl}({\mathbf{q}})$$

$$i M = i e \bar u(p')\Gamma^i(p', p)u(p) \tilde A^i_\mathrm{cl}({\mathbf{q}})$$

$$i M = i e \bar u(p')\left(\gamma^i F_1(q^2) + {i\sigma^{i\nu} q_\nu \over 2m} F_2(q^2) \right)
u(p) \tilde A^i_\mathrm{cl}({\mathbf{q}})$$

In the limit $q\to 0$ this becomes:

$$i M = i e
    2m\xi^{'\dag}\left(-i\epsilon^{ijk}{q^j\sigma^k\over 2m}(F_1(0) + F_2(0)) \right)\xi
\tilde A^i_\mathrm{cl}({\mathbf{q}})$$

$$i M = -i e
    2m\xi^{'\dag}\left(-{\sigma^k\over 2m}(F_1(0) + F_2(0)) \right)\xi
\left(-i\epsilon^{ijk}q^j\tilde A^i_\mathrm{cl}({\mathbf{q}})\right)$$

$$i M = -i e
    2m\xi^{'\dag}\left(-{\sigma^k\over 2m}(F_1(0) + F_2(0)) \right)\xi
    \tilde B^k({\mathbf{q}})$$

$$i M = -i \left(-{e\over m} (F_1(0) + F_2(0))
    2m\xi^{'\dag}{\sigma^k\over 2}\xi
    \tilde B^k({\mathbf{q}})\right)$$

where

$$\tilde B^k({\mathbf{q}}) =
\left(-i\epsilon^{ijk}q^j\tilde A^i_\mathrm{cl}({\mathbf{q}})\right)$$

is the Fourier transform of the magnetic field produced by $\mathbf{A}^{\mathrm{cl}}(\mathbf{x})$.

This corresponds to the Born approximation for scattering from a
potential

$$V({\mathbf{x}}) = -{e\over m} (F_1(0) + F_2(0))
    \xi^{'\dag}{\sigma^k\over 2}\xi
    B^k({\mathbf{x}})$$

$$V({\mathbf{x}}) = -{e\over m} (F_1(0) + F_2(0))
    \xi^{'\dag}{\boldsymbol{\sigma}\over 2}\xi\cdot {\mathbf{B}}({\mathbf{x}})$$

$$V({\mathbf{x}}) = -<{\boldsymbol{\mu}}>\cdot {\mathbf{B}}({\mathbf{x}})$$

where

$$<{\boldsymbol{\mu}}> = {e\over m} (F_1(0) + F_2(0)) \xi^{'\dag}{\boldsymbol{\sigma}\over 2}\xi$$

$$<{\boldsymbol{\mu}}> = g {e\over 2m} {\mathbf{S}}$$

where

$$g = 2(F_1(0) + F_2(0))$$

$${\mathbf{S}} = \xi^{'\dag}{\boldsymbol{\sigma}\over 2}\xi$$

The coefficient $g$ is called the Landé g-factor, and since the
leading order of perturbation theory gives $F_2(0)=0$ (and we know
that $F_1(0)=1$ to all orders), we get:

$$g = 2(F_1(0) + F_2(0)) = 2 + 2F_2(0) = 2 + O(\alpha)$$

This is the standard prediction of the Dirac equation. The anomalous
magnetic moment is then:

$$a_e = {g - 2\over 2} = F_2(0)$$

To calculate that, we need to evaluate the one-loop correction to the
vertex function, so we start by deriving the appropriate Green function
for the process $\gamma(q) + e^+(p) \to e^+(p')$:

$$|i\rangle = a^{r\dag}_{\mathbf{q}} b^{t\dag}_{\mathbf{p}} |\Omega\rangle$$

$$|f\rangle = b^{s\dag}_{\mathbf{p}'} |\Omega\rangle$$

$$\langle f|i \rangle =\langle\Omega| b^s_{\mathbf{p}'} a^{r\dag}_{\mathbf{q}}
     b^{t\dag}_{\mathbf{p}} |\Omega\rangle =$$

$$=\langle\Omega|T b^s_{\mathbf{p}'} a^{r\dag}_{\mathbf{q}}
         b^{t\dag}_{\mathbf{p}} |\Omega\rangle =$$

$$=\langle\Omega|T
         \bar u^s({\mathbf{p}'}){1\over\tilde S(p')}\tilde \psi(p')
         \epsilon_\mu^{r*}({\mathbf{q}}){q^2\over i} \tilde A^\mu(-q)
         \tilde{\bar\psi}(-p){1\over\tilde S(-p)}u^t({\mathbf{p}})
         |\Omega\rangle =$$

$$=\bar u^s({\mathbf{p}'}){1\over\tilde S(p')}
        \epsilon_\mu^{r*}({\mathbf{q}}){q^2\over i}
         \langle\Omega|T
         \tilde \psi(p')
         \tilde A^\mu(-q)
         \tilde{\bar\psi}(-p)
         |\Omega\rangle{1\over\tilde S(-p)}u^t({\mathbf{p}}) =$$

$$=\bar u^s({\mathbf{p}'}){1\over\tilde S(p')}
        \epsilon_\mu^{r*}({\mathbf{q}}){q^2\over i}
         \tilde G(p, p', q)
         {1\over\tilde S(-p)}u^t({\mathbf{p}}) =$$

where:

$$\tilde G(p, p', q) = \langle\Omega|T \tilde \psi(p') \tilde A^\mu(-q)
         \tilde{\bar\psi}(-p)
         |\Omega\rangle$$

is the interacting Green function for the Lagrangian $-\lambda \bar e \gamma^\mu e A_\mu$. In the first order:

$$\tilde G(p, p', q) = \langle\Omega|T \tilde\psi(p') \tilde A^\mu(-q)
         \tilde{\bar\psi}(-p)
         |\Omega\rangle =$$

$$= \int \mathrm{d}^4 x \langle0|T \tilde\psi(p') \tilde A^\mu(-q)
         \tilde{\bar\psi}(-p)
         (-\lambda)\bar e(x) \gamma^\rho e(x) A_\rho(x)
         |0\rangle =$$

$$= (-\lambda)\int \mathrm{d}^4 x \mathrm{d}\hat p'\mathrm{d}\hat q\mathrm{d}\hat p
        e^{i\hat p'p' - \hat q q
        -\hat pp}
        \langle0|T \psi(\hat p') A^\mu(\hat q)
         {\bar\psi}(\hat p)
         \bar e(x) \gamma^\rho e(x) A_\rho(x)
         |0\rangle =$$

$$= (-\lambda)\int \mathrm{d}^4 x \mathrm{d}\hat p'\mathrm{d}\hat q\mathrm{d}\hat p
        e^{i\hat p'p' - \hat q q
        -\hat pp}
        D^\mu_\rho(\hat q-x) S(\hat p' - x)\gamma^\rho S(\hat p-x)
        =$$

$$= (-\lambda)(2\pi)^4\delta(p'-q-p)
        \tilde D^\mu_\rho(q) \tilde S(p')\gamma^\rho \tilde S(p)$$

so the amplitude is:

$$\langle f|i \rangle=\bar u^s({\mathbf{p}'}){1\over\tilde S(p')}
    \epsilon_\mu^{r*}({\mathbf{q}}){q^2\over i}
 (-\lambda)(2\pi)^4\delta(p'-q-p)
    \tilde D^\mu_\rho(q) \tilde S(p')\gamma^\rho \tilde S(p)
     {1\over\tilde S(-p)}u^t({\mathbf{p}}) =$$

$$=(-\lambda)(2\pi)^4\delta(p'-q-p)\epsilon_\mu^{r*}({\mathbf{q}})
        u^s({\mathbf{p}'})\gamma^\mu u^t({\mathbf{p}})$$

and we got $\Gamma^\mu = \gamma^\mu$, so $F_1=1$ and $F_2=0$ in the
lowest order. In the next order we get:

$$\tilde G(p, p', q)
    = (-\lambda)(2\pi)^4\delta(p'-q-p)
        \tilde D^\mu_\rho(q) \tilde S(p')\delta\Gamma^\rho \tilde S(p)$$

$$\delta\Gamma^\mu =
    \int {\mathrm{d}^4 k\over (2\pi)^4} \tilde D_{\nu\rho}(k-p)
        (-ie\gamma^\nu)
        \tilde S(k')
        \gamma^\mu
        \tilde S(k)
        (-ie\gamma^\rho)$$

Now we can write:

$$\bar u(p')\Gamma^\mu(p', p) u(p) =
    \bar u(p')(\gamma^\mu + \delta\Gamma^\mu) u(p)$$

$$\bar u(p')\delta\Gamma^\mu(p', p) u(p) =
    \int {\mathrm{d}^4 k\over (2\pi)^4} \tilde D_{\nu\rho}(k-p)
        \bar u(p')
        (-ie\gamma^\nu)
        \tilde S(k')
        \gamma^\mu
        \tilde S(k)
        (-ie\gamma^\rho)
        u(p) =$$

$$=
    \int {\mathrm{d}^4 k\over (2\pi)^4} {-ig_{\nu\rho}\over (k-p)^2 +i\epsilon}
        \bar u(p')
        (-ie\gamma^\nu)
        {i(\not{k}' + m)\over k'^2-m^2 +i\epsilon}
        \gamma^\mu
        {i(\not{k} + m)\over k^2-m^2 +i\epsilon}
        (-ie\gamma^\rho)
        u(p) =$$

$$= 2ie^2\int {\mathrm{d}^4 k\over (2\pi)^4}
    {\bar u(p') \left(
        \not{k} \gamma^\mu \not{k}' + m^2\gamma^\mu - 2m(k+k')^\mu
        \right) u(p) \over
    ((k-p)^2 + i\epsilon)(k'^2 - m^2 + i\epsilon)(k^2-m^2+i\epsilon)
        }=$$

$$= \cdots =$$

$$= 2i e^2 \int {\mathrm{d}^4 l\over (2\pi)^4} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    {2\over D^3} \bar u(p') \left(
    \gamma^\mu (-\frac{1}{2} l^2+ (1-x)(1-y)q^2 + (1-4z+z^2)m^2)
        + {i\sigma^{\mu\nu}q_\nu\over 2m} (2m^2 z(1-z))
    \right)u(p) =$$

$$= {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    \bar u(p') \left(
    \gamma^\mu \left[\log {z \Lambda^2\over\Delta} + {1\over\Delta}
        \left((1-x)(1-y)q^2 + (1-4z+z^2)m^2\right)\right]
        + {i\sigma^{\mu\nu}q_\nu\over 2m}\left[{1\over\Delta}2m^2 z(1-z)
            \right] \right)u(p)$$

where

$$k' = k + q$$

$$D = l^2 - \Delta + i\epsilon$$

$$\Delta = -xyq^2 + (1-z)^2 m^2 > 0$$

So the expressions for the form factors are:

$$F_1(q^2) = 1 + {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    \left[\log {z \Lambda^2\over\Delta} + {1\over\Delta}
        \left((1-x)(1-y)q^2 + (1-4z+z^2)m^2\right)\right]
        +O(\alpha^2)$$

$$F_2(q^2) = {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
     \left[{1\over\Delta}2m^2 z(1-z) \right]
        +O(\alpha^2) =$$

$$= {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
     \left[2m^2 z(1-z)\over m^2(1-z)^2 - q^2 xy \right]
        +O(\alpha^2)$$

$F_1$ contains both ultraviolet and infrared divergencies. To cure the
infrared divergence, we add a term $\mu^2 z$ to $\Delta$. To cure the
ultraviolet divergence, we make the substitution:

$$F_1(q^2) \to F_1(q^2) - \delta F_1(0)$$

where $\delta F_1$ is the first order (in $\alpha$) correction to
$F_1$ (i.e. $F_1 = 1 + \delta F_1 + O(\alpha^2)$):

$$\delta F_1(0) = {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    \left[\log {z \Lambda^2\over\Delta (q^2=0)} + {1\over\Delta (q^2=0)}
        (1-4z+z^2)m^2\right]$$

so the corrected $F_1$ is:

$$F_1(q^2) = 1 + {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    \left[\log {z \Lambda^2\over\Delta} + {1\over\Delta}
        \left((1-x)(1-y)q^2 + (1-4z+z^2)m^2\right)+\right.$$

$$\left.-\log {z \Lambda^2\over\Delta (q^2=0)} - {1\over\Delta (q^2=0)}
    (1-4z+z^2)m^2\right]
    +O(\alpha^2) =$$

$$= 1 + {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
    \left[\log {m^2 (1-z)^2\over m^2(1-z)^2 - q^2 x y} +
        \left((1-x)(1-y)q^2 + (1-4z+z^2)m^2\over
        m^2(1-z)^2 - q^2 x y +\mu^2z
        \right)+\right.$$

$$\left.-{(1-4z+z^2)m^2\over m^2 (1-z)^2 + \mu^2 z}\right]
    +O(\alpha^2)$$

Neither the ultraviolet nor the infrared divergence affects
$F_2(q^2)$, so we just set $q=0$:

$$F_2(0) = {\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
     \left[2m^2 z(1-z)\over m^2(1-z)^2 \right] +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} x \,\mathrm{d} y \,\mathrm{d} z\,
    \delta(x+y+z-1)
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} y \int_0^1 \,\mathrm{d} z\,
    \theta(1-(1-y-z))\theta((1-y-z)-0)
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} y \int_0^1 \,\mathrm{d} z\,
    \theta(y+z)\theta(1-y-z)
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} y \int_0^1 \,\mathrm{d} z\,
    \theta(1-y-z)
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} z \int_0^{1-z} \,\mathrm{d} y
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} z (1-z)
     {2 z\over 1-z} +O(\alpha^2) =$$

$$={\alpha\over 2\pi} \int_0^1 \mathrm{d} z 2z + O(\alpha^2) =$$

$$= {\alpha\over 2\pi} + O(\alpha^2)$$

Thus we get the correction to the $g$-factor of the electron:

$$a_e = {g - 2\over 2} = F_2(0) = {\alpha\over 2\pi} \approx 0.0011614$$

Code:

    >>> from math import pi
    >>> alpha = 1/137.035999049
    >>> a_e = alpha / (2*pi)
    >>> a_e
    0.0011614097331824923

Experiments give $a_e = 0.00115965218073\pm 0.00000000000028$
([arXiv:1412.8284](http://arxiv.org/abs/1412.8284), eq. (1)).

Higher order corrections from QED can also be calculated:

$$a_e = A_1 \left({\alpha\over \pi}\right) +
      A_2 \left({\alpha\over \pi}\right)^2 +
      A_3 \left({\alpha\over \pi}\right)^3 +
      A_4 \left({\alpha\over \pi}\right)^4 + \cdots$$

we already know that $A_1 = \frac{1}{2}$. See for example
[hep-ph/9410248](http://arxiv.org/abs/hep-ph/9410248) for the expression
for $A_2$:

$$A_2 = \frac{197}{144} + \frac{3}{4} \zeta\left(3\right) - \frac{1}{2}
    \pi^{2} \operatorname{log}\left(2\right) + \frac{1}{12} \pi^{2} =$$

$$= -0.328478965579\dots$$

Code:

    >>> from sympy import zeta, S, log
    >>> A_2 = S(197)/144 + zeta(2)/2 + 3*zeta(3)/4 - 3*zeta(2) * log(2)
    >>> A_2.n()
    -0.328478965579194

See [hep-ph/9602417](http://arxiv.org/abs/hep-ph/9602417) for the
$A_3$ term:

$$A_3 = \frac{28259}{5184} - \frac{215}{24} \zeta\left(5\right)
+ \frac{100}{3} \left(\sum_{n=1}^{\infty} \frac{1}{2^{n} n^{4}} -
  \frac{1}{24} \pi^{2} \operatorname{log}^{2}\left(2\right) + \frac{1}{24}
  \operatorname{log}^{4}\left(2\right)\right) +$$

$$+\frac{139}{18}
\zeta\left(3\right) - \frac{298}{9} \pi^{2}
\operatorname{log}\left(2\right) + \frac{83}{72} \pi^{2}
\zeta\left(3\right) + \frac{17101}{810} \pi^{2} -
\frac{239}{2160} \pi^{4} =$$

$$= 1.181241456\dots$$

Code:

    >>> from sympy import pi, zeta, S, log, sum, var, oo
    >>> var("n")
    n
    >>> a4 = sum(1/(2**n * n**4), (n, 1, oo))
    >>> A_3 = 83*pi**2*zeta(3)/72 - 215*zeta(5)/24 + 100*(a4 + log(2)**4/24 - \
    ...         pi**2*log(2)**2/24)/3 - \
    ...         239*pi**4/2160 + 139*zeta(3)/18 - 298 * pi**2 * log(2)/9 + \
    ...         17101 * pi**2 / 810 + S(28259)/5184
    >>> A_3.n()
    1.18124145658720

Higher terms are only known numerically. The $A_4$ and $A_5$ terms
can be found in [arXiv:1412.8284](http://arxiv.org/abs/1412.8284):

$$A_4 = -1.912 98 (84)$$

$$A_5 = 7.795 (336)$$

We can now sum $a_e$ up to a given order by the following script:

    from sympy import pi, zeta, S, log, summation, var, oo
    var("n")
    a4 = summation(1/(2**n * n**4), (n, 1, oo))
    A1 = S(1)/2
    A2 = S(197)/144 + zeta(2)/2 + 3*zeta(3)/4 - 3*zeta(2) * log(2)
    A3 = 83*pi**2*zeta(3)/72 - 215*zeta(5)/24 + 100*(a4 + log(2)**4/24 - \
            pi**2*log(2)**2/24)/3 - \
            239*pi**4/2160 + 139*zeta(3)/18 - 298 * pi**2 * log(2)/9 + \
            17101 * pi**2 / 810 + S(28259)/5184
    A4 = -1.91298
    A5 = 7.795
    alpha = 1/137.035999049
    a_e_exp = 0.00115965218073
    a_e_exp_err = 0.00000000000028
    a_e_other = 0.00000000000448
    A = [A1, A2, A3, A4, A5]
    a_e= []
    for i in range(len(A)):
        a_e.append((A[i]*(alpha/pi)**(i+1)).n())
    print "========== ================"
    print "Order      $a_e$"
    print "========== ================"
    for i in range(len(A)):
        print "%d          %16.14f" % (i+1, sum(a_e[:i+1]))
    print "Other      %16.14f" % a_e_other
    print "Total      %16.14f" % (sum(a_e) + a_e_other)
    print "Experiment %16.14f" % a_e_exp
    print "Difference %16.14f" % abs(sum(a_e) + a_e_other - a_e_exp)
    print "Exp. err   %16.14f" % a_e_exp_err
    print "========== ================"

and obtain the following table:

| Order      | $a_e$            |
| ---------- | ---------------- |
| 1          | 0.00116140973318 |
| 2          | 0.00115963742812 |
| 3          | 0.00115965223232 |
| 4          | 0.00115965217663 |
| 5          | 0.00115965217716 |
| Other      | 0.00000000000448 |
| Total      | 0.00115965218164 |
| Experiment | 0.00115965218073 |
| Difference | 0.00000000000091 |
| Exp. err   | 0.00000000000028 |

The \"Other\" line are contributions from the dependence on the muon and
tau particle masses, the hadronic vacuum-polarization, the hadronic
light-by-light-scattering and the electroweak contribution (see
[arXiv:1412.8284](http://arxiv.org/abs/1412.8284)). The \"Difference\"
line is the difference from the theory (the \"Total\" line) and
experiment. The \"Exp. err\" line is the experimental error.

At this level of accuracy, the uncertainty of the exact value of
$\alpha$ is the primary cause of the difference from experiment, and
one can use this result to predict a more accurate value for $\alpha$,
assuming that QED and the standard model are valid.
