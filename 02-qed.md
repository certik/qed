# Quantum Electrodynamics (QED)

## Local Gauge Invariance

We use a metric with signature +2 in this section.

The Dirac equation for an electron is:

$$\mathcal{L}=\bar\psi(i\hbar c\gamma^\mu \partial_\mu-mc^2)\psi$$

Physical quantities like a charge density ($\bar\psi\psi$) or a current
($\bar\psi\gamma^\mu\psi$), are all invariant if we add a local phase
$\Lambda(x)$ to the field (this is called a local U(1) gauge
transformation):

$$\psi(x) \to e^{iq\Lambda(x)/\hbar} \psi(x)$$

$$\bar\psi(x) \to \bar \psi(x) e^{-iq\Lambda(x) / \hbar}$$

Where $q$ is a parameter that measures the strength of the phase
transformation (this will be later interpreted as a charge, for example
for electrons $q=-\|e\|$) and $\hbar$ is the Planck constant. And so
we require that the Lagrangian is also invariant under the local gauge
transformation, because there is no experiment that would change if this
local gauge transformation is applied on the wave functions. By putting
this gauge transformation into the Lagrangian density, we otain:

$$\mathcal{L} \to \bar\psi e^{-iq\Lambda(x) / \hbar}
    (i\hbar c\gamma^\mu \partial_\mu-mc^2)
e^{iq\Lambda(x) / \hbar} \psi =$$

$$= \bar\psi(i\hbar c\gamma^\mu (\partial_\mu + iq\partial_\mu\Lambda(x) /
    \hbar) -mc^2)\psi$$

The reason the Lagrangian is not invariant is due to the derivative,
which does not transform covariantly under a local gauge transformation:

$$\bar\psi \partial_\mu \psi
\to
\bar\psi e^{-iq\Lambda(x) / \hbar} \partial_\mu
    e^{iq\Lambda(x) / \hbar} \psi =$$

$$= \bar\psi (\partial_\mu
    + iq\partial_\mu\Lambda(x) / \hbar) \psi \neq \bar\psi \partial_\mu \psi$$

In order to make the derivative transform covariantly (and thus the
Lagrangian gauge invariant), we have to introduce a gauge field, in this
case a vector field $A_\mu(x)$, as follows:

[$$D_\mu = \partial_\mu-{i\over \hbar}qA_\mu$$]{label="covariant-derivative"}

and the field $A_\mu$ must transform as $A_\mu \to A_\mu + \partial_\mu\Lambda(x)$. At this level, we are free to choose either plus
or minus sign in `covariant-derivative`{.interpreted-text role="eq"},
since the sign change can be absorbed in the definition of the $A_\mu$
field without loss of generality (if we change the sign, the field
transformation then changes to $A_\mu \to A_\mu - \partial_\mu\Lambda(x)$). In the +2 metric signature we chose a minus
sign, so that $A_\mu$ coincides with the usual definition of the
electromagnetic 4-potential:

$$D_\mu = \partial_\mu-{i\over \hbar}qA_\mu$$

$$-i\hbar D_\mu = -i\hbar \partial_\mu - qA_\mu$$

$$m\hat v_\mu = \hat p_\mu - qA_\mu$$

$$m\hat{\mathbf{v}} = \hat {\mathbf{p}} - q{\mathbf{A}}$$

With signature -2, we must choose a plus sign and the identification
goes as follows:

$$D_\mu = \partial_\mu+{i\over \hbar}qA_\mu$$

$$i\hbar D_\mu = i\hbar \partial_\mu - qA_\mu$$

$$m\hat v_\mu = \hat p_\mu - qA_\mu$$

$$m\hat{\mathbf{v}} = \hat {\mathbf{p}} - q{\mathbf{A}}$$

And we obtain the same final equation. So the kinematic momentum is
equal to canonical momentum minus charge times the gauge field. The last
expression is independent of a metric signature, and that is what is
e.g. in the kinetic term of a Schrödinger or Pauli equation (with the
minus sign in $\hat{\mathbf{p}} - q\mathbf{A}$). We derive the non-relativistic
limit rigorously later, but it gives the same result. At this level we
just have to make sure we choose the correct sign in
`covariant-derivative`{.interpreted-text role="eq"}, depending on the
metric signature, otherwise we would get the electromagnetic 4-potential
with the opposite sign (the sign of $A_\mu$ is ultimately just a
convention, but later we want to get the same equations as everybody
else).

Another unrelated convention is in choosing the sign of the parameter
$q$. We have choosen it to coincide with an electric charge (negative
for electrons). Some authors choose $q$ to be positive for electrons,
then one must flip the sign in `covariant-derivative`{.interpreted-text
role="eq"}.

We will continue using the +2 signature in the rest of the section.

The operator $D_\mu = \partial_\mu-{i\over \hbar}q A_\mu$ is called a
covariant derivative, because it does not change a form (is invariant)
under a local gauge transformation:

$$\bar\psi D_\mu \psi = \bar\psi (\partial_\mu-{i\over \hbar}qA_\mu) \psi$$

$$\to
\bar\psi e^{-iq\Lambda(x) / \hbar} (\partial_\mu-{i\over \hbar}q(A_\mu
    + \partial_\mu \Lambda(x))) e^{iq\Lambda(x) / \hbar} \psi =$$

$$= \bar\psi (\partial_\mu
    - {i\over \hbar}q A_\mu - {i\over \hbar}q \partial_\mu \Lambda(x)
    + iq\partial_\mu\Lambda(x) / \hbar) \psi =$$

$$= \bar\psi (\partial_\mu - {i\over \hbar}qA_\mu) \psi
= \bar\psi D_\mu \psi$$

Then the Lagrangian

[$$\mathcal{L}=\bar\psi(i\hbar c\gamma^\mu (\partial_\mu-iqA_\mu / \hbar)-mc^2)\psi$$]{label="lag_inv"}

is also gauge invariant:

$$\mathcal{L} \to \bar\psi e^{-iq\Lambda(x) / \hbar}
    (i\hbar c\gamma^\mu (\partial_\mu-iqA_\mu / \hbar
        -iq\partial_\mu\Lambda(x) / \hbar)
    -mc^2) e^{iq\Lambda(x) / \hbar} \psi =$$

$$= \bar\psi(i\hbar c\gamma^\mu (\partial_\mu -iqA_\mu / \hbar
    -iq\partial_\mu\Lambda(x) / \hbar + iq\partial_\mu\Lambda(x) / \hbar)
    -mc^2)\psi =$$

$$= \bar\psi(i\hbar c\gamma^\mu (\partial_\mu -iqA_\mu / \hbar) -mc^2)\psi$$

The Lagrangian `lag_inv`{.interpreted-text role="eq"} can also be
written as:

$$\mathcal{L}=\bar\psi(i\hbar c\gamma^\mu (\partial_\mu-{i \over \hbar} qA_\mu)
    -mc^2)\psi =$$

$$= \bar\psi(i\hbar c\gamma^\mu \partial_\mu-mc^2)\psi
    + qc\bar\psi\gamma^\mu \psi  A_\mu$$

We can see that the condition of a local gauge invariance requires an
interaction with a vector field $A_\mu$. Now we need to add the kinetic
term for the field $A_\mu$:

$$-{1\over4}F_{\mu\nu}F^{\mu\nu}$$

The mass term $\frac{1}{2} m^2 A_\mu A^\mu$ is not gauge invariant, and so
we have to set $m=0$. Here is the full Lagrangian:

$$\mathcal{L}= \bar\psi(i\hbar c\gamma^\mu \partial_\mu-mc^2)\psi
    + q c\bar\psi\gamma^\mu \psi  A_\mu -{1\over4}F_{\mu\nu}F^{\mu\nu}$$

This is a Lagrangian for an electron and a massless vector boson
(photon) of spin 1. We can introduce a current $j^\mu = c\bar\psi\gamma^\mu\psi$, then the Lagrangian density becomes:

$$\mathcal{L}= \bar\psi(i\hbar c\gamma^\mu \partial_\mu-mc^2)\psi
    + q j^\mu  A_\mu -{1\over4}F_{\mu\nu}F^{\mu\nu}$$

For an electron, we can set $q=-e$, where $e$ is the elementary
charge ($e$ is positive).

## QED Lagrangian

We use a metric with signature -2 in this section.

The QED Lagrangian density is

$$\mathcal{L}=\bar\psi(i\hbar c\gamma^\mu D_\mu-mc^2)\psi-{1\over4}F_{\mu\nu}F^{\mu\nu}$$

where

$$\begin{aligned}
\psi=\left( \begin{array}{c} \psi_1 \\ \psi_2 \\ \psi_3 \\ \psi_4 \\ \end{array}\right)
\end{aligned}$$

and we must choose a plus sign in
`covariant-derivative`{.interpreted-text role="eq"} since we use the -2
signature:

$$D_\mu=\partial_\mu+{i\over \hbar}eA_\mu$$

$e$ is the charge (negative for electrons $e=-\|e\|$).

$$F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu$$

is the electromagnetic field tensor. It's astonishing, that this simple
Lagrangian can account for all phenomena from macroscopic scales down to
something like $10^{-13}\,\mathrm{cm}$. So it's not a surprise that Feynman,
Schwinger and Tomonaga received the 1965 Nobel Prize in Physics for such
a fantastic achievement.

Plugging this Lagrangian into the Euler-Lagrange equation of motion for
a field, we get:

$$(i\hbar c\gamma^\mu D_\mu-mc^2)\psi=0$$

$$\partial_\nu F^{\nu\mu}=-ec\bar\psi\gamma^\mu\psi$$

The first equation is the Dirac equation in the electromagnetic field
and the second equation is a set of Maxwell equations ($\partial_\nu F^{\nu\mu}=-e j^\mu$) with a source $j^\mu=c\bar\psi\gamma^\mu\psi$, which
is a 4-current comming from the Dirac equation.
