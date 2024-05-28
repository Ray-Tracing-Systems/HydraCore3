### Paper: [Differentiable Monte Carlo Ray Tracing through Edge Sampling](https://dl.acm.org/doi/pdf/10.1145/3272127.3275109)
*Tzu-Mao Li | 2018 | ≈ 500 cites*

* **HHMC (2015)** does not does not take geometric discontinuities into account (this paper does).
* "We assume no point light sources and no perfectly specular surfaces" $-$ sounds like bruh.

**Their problem**: differentiate $\int_M f(x, y)dx$ with respect to $y$.  
We cannot do $\frac{\partial}{\partial y}\int_M f(x, y)dx = \int_M \frac{\partial}{\partial y} f(x, y)dx$,
because f(x, y) is not necessarily differentiable by $y$.

Example: $\int_M f(x, y)dx = \int_0^1 sgn(x - y)dx = \int_0^ydx = y$ is differentiable by y (but $sgn(\cdot)$ is not).

**Our (Langevin) problem**: given $\int_M f(x) dx$, differentiate $f(x)$.  
Discontinuity of $f(x)$ is a problem too, but not directly related.

Is this even applicable to Langevin?


### Testing how [EnzymeAD](https://enzyme.mit.edu/) deals with $\log(x \cdot y)$
Idea: it is faster to differentiate $\log(x) + \log(y)$ than $\log(x \cdot y)$.
Experiments show that Enzyme **does not** run any relevant optimizations.

Run `log_test.cpp` in `langevin_specific` folder to see the results.
Example command (use your own `ClangEnzyme` path):
```
clang++-17 log_test.cpp -fplugin=../../Enzyme_103/enzyme/build/Enzyme/ClangEnzyme-17.so -O3 && ./a.out
```

Output:
```
Differentiate `log(x_1 * ... * x_100)`
Elapsed time: 880ms
Gradient: ( ... )

Differentiate `log(x_1) + ... + log(x_100)`
Elapsed time: 70ms
Gradient: ( ... )

```

### Paper: [3D Gaussian Splatting as Markov Chain Monte Carlo](https://arxiv.org/abs/2404.09591)
*Shakiba Kheradmand на стажировке в гугле | 2024*

Contributions:
* **"we rethink Gaussians in Gaussian Splatting as samples from a distribution that represents the 3D
scene, which allows simplification"** — OK
* **"we improve robustness to initialization"** — initialization algorithm does not change, robustness comes from next point
* **"we propose to simply add samples and resample by drawing samples from existing Gaussian
locations under the SGLD paradigm"** — didn't crack the SGLD paradigm; however, instead of splitting / deleting gaussians
they teleport gaussians with respect to the distribution of existing gaussians 
(how? need to check [this](https://epubs.siam.org/doi/10.1137/21M1425062))
* **"we introduce L1 regularization on opacity and scale to encourage using fewer Gaussians and
increase efficiency and performance"** — good, because
    - Heuristics (split / delete) -> theory (regularization, bad gaussians become too small and are teleported)
    - Smaller gaussians means better performance
* **"we provide higher rendering quality when the heuristics are a limiting factor"** — k

Still not quite clear: 
* **Adding Gaussians**: why add gaussians if we can teleport bad ones? If it is important, why limit to 5%?
* SGLD paradigm
