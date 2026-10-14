# Density matrices and distributions

## Reduced density matrix

The reduced density matrix of one site,

$$\rho_i=\mathrm{Tr}_{j\neq i}|\mathrm{GS}\rangle\langle\mathrm{GS}|,$$

contains every local observable of that site, and its eigenvalues give the
site entanglement entropy ([Entanglement](entanglement.md)).

```python
rho = sc.get_rdm(i=0)     # d x d matrix, d the local dimension
```

It is computed from the MPS, so it is not available with `mode="ED"`,
where it raises an error.

## Operator distributions

For a Hermitian operator $X$, the expectation
value $\langle X\rangle$ is only the mean. The full distribution of the
measured values of $X$ in the ground state is

$$P(x)=\langle\mathrm{GS}|\,\delta(X-x)\,|\mathrm{GS}\rangle,\qquad
\int P(x)\,dx=1.$$

For a conserved charge this is the full counting statistics. A sharp
$P(x)$ means that the ground state has a well-defined value of $X$, and a
broad one means large fluctuations. It is computed with the kernel
polynomial method, as the [dynamical correlators](dynamical-correlators.md):

```python
X = sum(sc.Sz[i] for i in range(n))      # total magnetization
x, p = sc.get_distribution(X=X, delta=0.1)
mus, shift, scale = sc.get_distribution_moments(X=X, delta=0.1)
```

`delta` sets the resolution through the number of Chebyshev moments.
Unlike in [Dynamical correlators](dynamical-correlators.md), it is not calibrated as a broadening, so compare
distributions computed with the same `delta`. `get_distribution` also
works with `mode="ED"`. `get_distribution_moments` returns the raw
Chebyshev moments of the rescaled operator and has no ED version. A
maximum-entropy reconstruction from the moments $\langle X^k\rangle$,
`distribution.get_distribution_maxent`, needs the module
`dmrgpy.maxenttk`, which is not shipped, and raises an error without it.
