### [Differentiable Monte Carlo Ray Tracing through Edge Sampling](https://dl.acm.org/doi/pdf/10.1145/3272127.3275109)
*Tzu-Mao Li | 2018 | $\approx$ 500 cites*

* **HHMC (2015)** does not does not take geometric discontinuities into account (this paper does).
* "We assume no point light sources and no perfectly specular surfaces" $-$ sounds like bruh.

**Their problem**: differentiate $\int_M f(x, y)dx$ with respect to $y$.  
We cannot do $\frac{\partial}{\partial y}\int_M f(x, y)dx = \int_M \frac{\partial}{\partial y} f(x, y)dx$,
because f(x, y) is not necessarily differentiable by $y$.

Example: $\int_M f(x, y)dx = \int_0^1 sgn(x - y)dx = \int_0^ydx = y$ is differentiable by y (but $sgn(\cdot)$ is not).

**Our (Langevin) problem**: given $\int_M f(x) dx$, differentiate $f(x)$.  
Discontinuity of $f(x)$ is a problem too, but not directly related.

Is this even applicable to Langevin?