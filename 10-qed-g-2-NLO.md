# Magnetic moment of an electron (NLO)

In this section we derive the next-to-leading-order (NLO) order-$\alpha^2$
correction to the magnetic moment of an electron. This has been first computed
in 1957 by Petermann and Sommerfield.

See for example
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

