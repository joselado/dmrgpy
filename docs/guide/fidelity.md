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
differences, so it needs only ground states. It also handles a degenerate
ground state. `fmode="PT"` sums over the `n` lowest states. Increase `n`
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
computed with `dl=0.02`. On the left is $\chi/n$ for four lengths. For the
longer chains the peak grows with $n$ and moves towards the critical
point $\lambda=1$: it is at $\lambda=0.86$, 0.91 and 0.94 for $n=16$, 24
and 32. For $n=12$ the ED curve (dashed line)
covers the DMRG one. On the right, the difference of the two longest
chains removes the contribution of the edges and gives $\chi/n$ of the
bulk. Away from $\lambda=1$ it lies on the exact result of the infinite
chain, $\chi/n=1/[16(1-\lambda^2)]$ for $\lambda<1$ and
$1/[16\lambda^2(\lambda^2-1)]$ for $\lambda>1$. Close to $\lambda=1$ the
correlation length is longer than the chains, so they do not reach the
infinite-chain result there. We start the scan at $\lambda=0.85$: deeper
in the ordered phase, the two ground states of a long open chain become
degenerate within the tolerance `delta`. The default estimator then
treats them as one degenerate manifold, and the number it returns is no
longer $\chi$ of one ground state.
