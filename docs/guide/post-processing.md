# Post-processing tools

- `extrapolate.size_extrapolation(L, y, n=1)`: fits $y(L)$ with a
  polynomial of order `n` in $1/L$ and returns the value at
  $L\to\infty$. Use it for energy densities, gaps and order parameters
  from finite chains.
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
