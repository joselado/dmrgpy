# Correlators and distributions

On this page we measure the ground state beyond a two-point correlator:
four-point correlators, the two-particle density matrix, and the full
distribution of an operator. Two-point correlators and the equal-time
structure factor are on the [Ground state](ground-state.md#static-correlators)
page, and the single-particle correlation matrix, whose eigenvalues are
the natural orbitals, is under [Entanglement](entanglement.md#single-particle-correlation-matrix).

## Four-point correlators

For a fermionic chain `fc`, the four-point tensor
$\langle c_i^\dagger c_jc_k^\dagger c_l\rangle$ for all $i,j,k,l$ is

```python
t = fc.get_gs().get_four_correlation_tensor()   # shape (n, n, n, n)
```

It holds every density-density correlator $\langle n_in_k\rangle$ (at
$j=i$, $l=k$) and every pair correlator, so it is the input of a
two-particle reduced density matrix or of a cumulant. The keyword
`ctmode` chooses the method (`"batched"`, `"sweep"`, `"fold"`, `"full"`,
`"explicit"`). The default chooses the fastest one for the backend, and
`"explicit"` works on every backend. A single four-point correlator is
also `fc.vev(fc.Cdag[i]*fc.C[j]*fc.Cdag[k]*fc.C[l])`, which is cheaper
when only a few are needed.

The two-particle density matrix
$\langle c_i^\dagger c_j^\dagger c_lc_k\rangle$, with the pairs $i<j$ and
$k<l$ as row and column, is

```python
from dmrgpy.entropytk.correlationentropy import get_highorder_correlation_matrix
d2 = get_highorder_correlation_matrix(fc)
```

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
`distribution.get_distribution_maxent(sc, X=..., n=6)`, gives the smoothest
distribution with the first `n` moments. With few moments it shows the
envelope of the distribution, not its separate peaks.

![Magnetization distribution of the Ising chain](img/correlators_distribution.png#only-light)
![Magnetization distribution of the Ising chain](img/correlators_distribution-dark.png#only-dark)

The figure shows $P(m)$ for the magnetization per site
$m=\frac1n\sum_iS^z_i$ in the transverse-field Ising chain
$H=-\sum_iS^z_iS^z_{i+1}-h\sum_iS^x_i$, which is critical at $h=1/2$.
A small field $10^{-4}\sum_iS^z_i$ selects one of the two ordered states.
In the ordered phase ($h=0.2$), $P(m)$ is one sharp peak close to the
exact bulk value $m=0.489$. At the critical field it is broad and flat,
and in the paramagnet ($h=1$) it is one peak at $m=0$. On the left, with
$n=12$, the DMRG curves lie on the ED curves (dashed lines). There the
peaks are separated, because $m$ only takes the values $k/n$. On the
right, with $n=24$, these values are closer than the resolution set by
`delta`, so the distribution looks smooth.

??? example "Code of this figure (docs/figures/correlators_distribution/)"

    ```python title="compute.py"
    --8<-- "correlators_distribution/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "correlators_distribution/plot.py"
    ```

The many-body density of states, $\rho(E)=\sum_n\delta(E-E_n)$, is not a
distribution of this kind: `get_distribution` computes a ground-state
expectation value, not a trace over all the states, and there is no
method for $\rho(E)$. For a small chain the full ED spectrum from
`get_excited(mode="ED", n=...)` ([Excited states](excited-states.md))
gives it directly.
