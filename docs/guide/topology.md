# Topological invariants

## Many-body Berry phase

Let the Hamiltonian depend on a parameter
$\phi$ that goes around a closed loop, for example a flux $\phi$ through a
ring, entering as a twist $e^{i\phi}$ on one bond. The discrete Berry
phase of the ground state is

$$\gamma=\arg\prod_{k}\langle\psi(\phi_k)|\psi(\phi_{k+1})\rangle,$$

with the last state connected to the first one. With inversion symmetry
$\gamma$ is quantized to $0$ or $\pi$, and the two values are two
different phases, for example the two dimerizations of a spin chain. We
compute it by solving the ground state at each $\phi_k$ and multiplying
the overlaps. Here `h_of(phi)` returns the twisted Hamiltonian of `sc`:

```python
import numpy as np
phis = np.linspace(0, 2*np.pi, 20, endpoint=False)
wfs = []
for phi in phis:
    sc.set_hamiltonian(h_of(phi))
    wfs.append(sc.get_gs())
prod = np.prod([wfs[i].dot(wfs[(i+1) % len(wfs)]) for i in range(len(wfs))])
gamma = np.angle(prod)        # 0 or pi
```

See `examples/topological/zak_phase` for the dimerized Heisenberg chain.

## Single-particle Berry phase

For a Bloch Hamiltonian $H(k)$, the Zak
phase of the occupied bands is the phase of the Wilson loop,

$$\gamma=\arg\det\prod_k U(k,k+\delta k),\qquad
U_{mn}(k,k+\delta k)=\langle u_m(k)|u_n(k+\delta k)\rangle,$$

with $m,n$ running over the occupied bands. It gives the polarization of
the band. Here `hk(k)` returns the Bloch matrix, with $k$ from 0 to 1:

```python
from dmrgpy import topology
gamma = topology.berry_phase_matrix(hk, nk=40)   # the lower half of the bands is occupied
```

## Fermion parity

For a fermionic chain the total parity is

$$P=\Big\langle\prod_i(1-2n_i)\Big\rangle.$$

It is $+1$ or $-1$ for a state with an even or odd number of fermions. In
a chain with superconducting pairing the particle number is not
conserved, but the parity is. A change of the ground-state parity as a
function of a parameter is a level crossing between the two sectors, as
in a Kitaev chain or a magnetic impurity in a superconductor. Here `fc` is
a `Fermionic_Chain`:

```python
wf = fc.get_gs()
p = wf.get_fermionic_parity()                     # +1 or -1
p = wf.get_fermionic_parity(fpmode="iterative")   # same number, other route
```

See `examples/topological/parity` for $P$ as a function of the chemical
potential.
