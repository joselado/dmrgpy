# Fidelity susceptibility

Let the Hamiltonian depend on a parameter, $H(\lambda)=H_0+\lambda H_1$.
The fidelity susceptibility measures how fast the ground state changes
with $\lambda$. In perturbation theory it is

$$\chi(\lambda)=\sum_{n\neq0}\frac{|\langle 0|H_1|n\rangle|^2}{(E_n-E_0)^2+\delta^2}.$$

At a quantum phase transition the gap closes and $\chi$ has a peak that
grows with system size. So a scan of $\chi(\lambda)$ locates a critical
point without knowing the order parameter. Here `sc` is a spin chain:

```python
from dmrgpy import fidelity
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
chi = fidelity.get_fidelity(sc, h0, h1, 1.0)               # default
chi_pt = fidelity.get_fidelity(sc, h0, h1, 1.0, fmode="PT")  # formula above
```

The default, `fmode="derivative"`, does not use the formula above. It
computes the second derivative of the ground-state overlap by finite
differences, so it needs only ground states. If the ground state is
degenerate within `delta`, it uses the overlap matrix of the degenerate
states and returns the sum of $\chi$ over these states. `fmode="PT"` sums
over the `n` lowest states. Increase `n`
until $\chi$ stops changing. The keywords are:

- `dl`: the finite-difference step in $\lambda$ (default 0.1). Use a
  smaller value when $\chi$ changes fast.
- `n`, `delta`: for the default estimator, the number of low states
  searched for a degenerate ground state, and the energy tolerance for
  calling two states degenerate. For `fmode="PT"`, `n` is the number of
  states in the sum (default 2) and `delta` is the $\delta$ of the
  formula above.
- `mode="ED"`: exact diagonalization, for small chains.

![Fidelity susceptibility of the Ising chain](img/fidelity_ising.png#only-light)
![Fidelity susceptibility of the Ising chain](img/fidelity_ising-dark.png#only-dark)

The figure shows open Ising chains, $H=\sum_i\sigma^z_i\sigma^z_{i+1}+\lambda\sum_i\sigma^x_i$,
computed with `dl=0.02`. On the left is $\chi/n$ for four lengths. The
peak grows with $n$ and moves towards the critical point $\lambda=1$: it
is at $\lambda=0.83$, 0.86, 0.91 and 0.94 for $n=12$, 16, 24 and 32. For
$n=12$ the ED curve (dashed line) covers the DMRG one. On the right, the
difference of the two longest chains removes the contribution of the
edges and gives $\chi/n$ of the bulk. From about $\lambda=1.1$ it lies on the
exact result of the infinite chain, $1/[16\lambda^2(\lambda^2-1)]$. On the
ordered side the exact result is $1/[16(1-\lambda^2)]$, and the two chains
are still too short to reach it: the difference crosses it at
$\lambda=0.85$ and lies below it at $\lambda=0.8$. The exact free-fermion
solution of the open chain shows that at $\lambda=0.8$ the difference
reaches the bulk value only for $n$ of about 50. Close to $\lambda=1$ the correlation length is longer than
the chains, so they do not reach the infinite-chain result there.

We start the scan at $\lambda=0.8$. In the ordered phase the two lowest
states of an open chain have an energy splitting that decays with $n$.
For $n=32$ it is $6\times10^{-4}$ at $\lambda=0.8$ and drops below the
default `delta` near $\lambda=0.75$. Below that point the default
estimator returns $\chi$ summed over the two states. A smaller `delta`
does not help: DMRG does not separate two states this close in energy,
and the overlaps are then not meaningful.
