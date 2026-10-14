# Post-processing tools

- `extrapolate.size_extrapolation(L, y, n=1)`: fits $y(L)$ with a
  polynomial of order `n` in $1/L$ and returns the value at
  $L\to\infty$. Use it for energy densities, gaps and order parameters
  from finite chains.

![Size extrapolation of the energy per site](img/post-processing_extrapolation.png#only-light)
![Size extrapolation of the energy per site](img/post-processing_extrapolation-dark.png#only-dark)

The figure shows the ground-state energy per site of open Heisenberg
chains, from $L=16$ to $L=128$, as a function of $1/L$. On an open chain
the edges add an energy of order one, so $E_0/L=e_\infty+a/L+b/L^2$. The
quadratic fit (`n=2`) gives $e_\infty=-0.443143$, which differs from the
Bethe-ansatz value $1/4-\ln 2=-0.443147$ by $4\times10^{-6}$. The linear
fit (`n=1`) misses the curvature and is more than 15 times worse. The
inset zooms on $1/L=0$.

- `functionfit.fit(f, x0)`: minimizes a function `f` with several random
  starts of the Powell method. Use it for fits such as the
  [Calabrese-Cardy formula](entanglement.md#central-charge).
- `analyticcontinuation.imag2real(zs, us)`: Pade continuation of a
  function known at complex frequencies `zs` to the real axis, for
  example the output of `submode="CVMimag"` ([Choosing a correlator method](correlator-methods.md)). It needs the
  module `dmrgpy.padetk`, which is not shipped, and raises an error
  without it.
- `reconstruct.reconstruct_distribution`: maximum-entropy reconstruction of a positive spectral
  function from its moments, used by `submode="maxent"` and
  `get_distribution_maxent`. It needs the module `dmrgpy.maxenttk`, which
  is not shipped, and raises an error without it.
