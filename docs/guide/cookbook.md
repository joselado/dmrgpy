# Cookbook

Each recipe starts from a chain `sc` of `n` sites, built as in [Ground state](ground-state.md)
with the Hamiltonian given.

## Central charge of the critical Ising chain

For
$H=\sum_iS_i^zS_{i+1}^z+\tfrac12\sum_iS_i^x$, at the critical field, the
entanglement entropy of the open chain follows the Calabrese-Cardy formula
with $c=1/2$ ([Entanglement](entanglement.md)):

```python
sc.maxm = 200          # critical chains need a large bond dimension
wf = sc.get_gs()
c = wf.get_CFT_central_charge()
```

## Haldane gap of the spin-1 chain

The open chain has four low states
from the two spin-1/2 edge states. The bulk gap is the fifth state
([Excited states](excited-states.md)):

```python
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
# ... Heisenberg h, set_hamiltonian(h) ...
es = sc.get_excited(n=6)
gap = es[4] - es[0]
```

## Gaps of a Hubbard chain

Here `fc` is a `Spinful_Fermionic_Chain`
with a Hubbard Hamiltonian. `get_gap()` is the gap to the first excited
state, and `get_charge_gap(d=...)` is $E_0(N+d)-E_0(N)$ ([Excited states](excited-states.md)):

```python
fc.get_gap()
fc.get_charge_gap(d=1)
fc.get_charge_gap(d=2)
```

## Fidelity susceptibility across the Ising transition

See [Fidelity susceptibility](fidelity.md).

```python
import numpy as np
from dmrgpy import fidelity
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
lambdas = np.linspace(0.5, 1.5, 11)
chis = [fidelity.get_fidelity(sc, h0, h1, lam) for lam in lambdas]
```

## Dynamical structure factor $S(q,\omega)$

See [Dynamical correlators](dynamical-correlators.md). We compute it from the
site-resolved correlators and a Fourier transform,
$S(q,\omega)=\frac1n\sum_{ij}e^{iq(i-j)}S_{ij}(\omega)$:

```python
es = np.linspace(-0.5, 3, 200)
q = np.pi
Sqw = 0
for i in range(n):
    for j in range(n):
        x, y = sc.get_dynamical_correlator(submode="KPM", es=es, delta=0.1,
                                           name=(sc.Sz[i], sc.Sz[j]))
        Sqw = Sqw + np.exp(1j*q*(i-j))*y/n
```

## Extrapolation to the thermodynamic limit

See [Post-processing tools](post-processing.md). Here we extrapolate the
energy per site computed on several lengths `Ls`:

```python
from dmrgpy import extrapolate
e_inf = extrapolate.size_extrapolation(Ls, e_per_site, n=1)
```
