# DMRGPY

DMRGPY is a Python library for quasi-one-dimensional quantum many-body
systems: spin chains, spinless and spinful fermions, bosons and
parafermions. Ground states, excitations, dynamical correlators, real-time
evolution and finite temperature are computed with matrix product states.
Small systems can also be solved with exact diagonalization (ED), so every
DMRG result can be checked on a short chain.

This guide explains how to compute physical quantities. For each
capability it gives the quantity, its formula, the code and the keywords
that matter. The code itself is on
[GitHub](https://github.com/joselado/dmrgpy), and the whole guide is also
available as one [PDF](https://joselado.github.io/dmrgpy/user_guide.pdf).

## Installation

```bash
pip install dmrgpy
```

This installs ED and a DMRG code written in pure Python. It needs no
compiler and works on Linux, Mac and Windows. The faster C++ backend
(ITensor) is compiled from a clone of the repository with
`python install.py`; the
[README](https://github.com/joselado/dmrgpy#how-to-install) explains how.

## A first calculation

The ground-state energy of a spin-1/2 Heisenberg chain of 12 sites, with
DMRG and with ED:

```python
from dmrgpy import spinchain
n = 12
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
print(sc.gs_energy())            # DMRG
print(sc.gs_energy(mode="ED"))   # exact diagonalization
```

Both print $E_0\approx-5.1420906$. Every chain class has the same methods,
and most of them accept `mode="ED"` in the same way.

![Dynamical spin structure factor of the Heisenberg chain](img/dynamical-correlators_sqw.png#only-light)
![Dynamical spin structure factor of the Heisenberg chain](img/dynamical-correlators_sqw-dark.png#only-dark)

A Heisenberg chain with $n=32$ gives the dynamical spin structure factor
$S(q,\omega)$, the quantity measured by inelastic neutron scattering. The
weight fills the two-spinon continuum between its exact edges (dashed). See
[Dynamical correlators](dynamical-correlators.md).

??? example "Code of this figure (docs/figures/dynamical-correlators_sqw/)"

    ```python title="compute.py"
    --8<-- "dynamical-correlators_sqw/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "dynamical-correlators_sqw/plot.py"
    ```

## Backends

The keyword `itensor_version=` of a chain chooses the DMRG code. The
results are the same; only the speed and a few features differ.

- `3`: ITensor v3 in C++, the default when it is compiled. The fastest.
- `"python"`: the pure-Python DMRG, the default after `pip install`. It
  can also run on a GPU ([GPU](gpu.md)).
- `2`: ITensor v2 in C++. Some features are missing.
- `"julia_live"`: ITensors.jl in Julia, after
  `pip install "dmrgpy[julia]"`.

`mode="ED"` uses exact diagonalization instead of DMRG, for small chains.

## Contents of the guide

- [Getting started](models.md): the chain classes, their local operators,
  and how to write a Hamiltonian.
- [Ground state](ground-state.md): energies, expectation values,
  correlators, conserved sectors, excited states and entanglement.
- [Phases](topology.md): Berry phases, fermion parity, fidelity
  susceptibility and mean-field decoupling.
- [Dynamics](dynamical-correlators.md): dynamical correlators, spectral
  functions, real-time evolution and STM spectra.
- [Finite temperature](finite-temperature.md): purification, METTS and
  ED.
- [Infinite chains](infinite-chains.md): iDMRG and VUMPS directly in the
  thermodynamic limit.
- [More](cookbook.md): short recipes, post-processing tools and
  performance.
