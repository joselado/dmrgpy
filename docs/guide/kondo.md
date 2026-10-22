# STM and Kondo spectra

On this page we compute the differential conductance $dI/dV$ measured
by an STM tip above a magnetic atom or a spin chain on a surface. The
electrons tunnel from the tip to the sample, and on the way they exchange
spin and energy with the spins. We follow the perturbation theory of
Ternes, New J. Phys. 17, 063016 (2015), arXiv:1505.04430, to third order
in the coupling. The tip couples to one site of the chain. Here the chain
is a single spin-1/2 in a magnetic field, with $H$ in eV:

```python
import numpy as np
from dmrgpy import spinchain
muB, g, B = 5.7883818066e-5, 2.0, 10.0    # eV/T, g-factor, field in T
sc = spinchain.Spin_Chain(["S=1/2"])
sc.set_hamiltonian(g*muB*B*sc.Sz[0])
eV = np.linspace(-4e-3, 4e-3, 161)        # bias in eV
eV, dIdV = sc.get_kondo_spectrum(eV, site=0, Jrho_s=-0.05, U=0.25, T=1.0)
```

## Second order

Second order (`order=2`) is the Fermi golden rule for spin-flip and
elastic tunneling,

$$\frac{dI}{dV}(eV)\propto\sum_{i,f}p_i\Big[\tfrac12|\langle f|S_-|i\rangle|^2+\tfrac12|\langle f|S_+|i\rangle|^2+|\langle f|S_z|i\rangle|^2\Big]\,\Theta(eV-\epsilon_{if})+4U^2,$$

summed over both tunneling directions. Here $p_i$ is the Boltzmann weight
of eigenstate $i$, $\epsilon_{if}=E_f-E_i$, and $\Theta$ a step broadened
by temperature. Each excitation of the spins gives a step in $dI/dV$ at
$eV=\pm\epsilon_{if}$. This is inelastic tunneling spectroscopy (IETS),
and the step positions measure the spin excitation energies, for example
the Zeeman energy.

## Third order

Third order (`order=3`, the default) adds two terms. Both sum over
every eigenstate $m$ as a virtual intermediate state. The Kondo term
involves three spin matrix elements and a logarithmic function
$F(\epsilon,T)$. It gives a peak at zero bias, which splits into two
peaks at the Zeeman steps in a field. The potential-interference term,
present when `U` is not zero, makes the spectrum asymmetric in bias. The
second-order and Kondo terms are even in $eV$, and the potential term is
odd in $eV$.

At the parameters of the paper ($S=1/2$, $B=0$, $T=1$ K,
$J\rho_s=-0.05$, $\omega_0=20$ meV), the zero-bias peak divided by $2\pi$
is 1.13 for $U=0$ and 1.39 for $U=0.25$, as in the paper's figures. In a
field the code departs from the paper in one place: the exchange diagram
of the Kondo term has its logarithm at $F(eV-(\epsilon_f-\epsilon_m))$,
which is what an exact second-order T-matrix gives, and not at the
argument of the paper's eq. 25. The two agree at zero field for $S=1/2$
and at zero bias. In the cases tested with a field or anisotropy, they
differ by 5 to 20 per cent of the third-order term. The keywords are:

- `site`: the site under the tip.
- `Jrho_s`: the dimensionless Kondo coupling $J\rho_s$, in the sign
  convention of the paper (which uses $-0.05$).
- `U`: the ratio of potential scattering to spin scattering.
- `T`: the temperature in kelvin. `kB` (default in eV/K) converts it to
  the energy unit of $H$.
- `omega0`, `Gamma0`: the band cutoff and the lifetime broadening of the
  Kondo function $F$.
- `T0`: the tunneling strength. It only sets the overall scale.
- `order`: 2 or 3.

![STM spectra of magnetic atoms and chains](img/kondo_dIdV.png#only-light)
![STM spectra of magnetic atoms and chains](img/kondo_dIdV-dark.png#only-dark)

Left: the second-order spectrum of an $S=2$ atom with
$H=D(S^z)^2+E[(S^x)^2-(S^y)^2]$, $D=-1.55$ meV and $E=0.31$ meV (Fe on
Cu$_2$N), at $T=0.5$ K. The steps sit at the excitation energies of the
spin multiplet, 0.18, 3.90 and 5.76 meV (dotted lines), which come from
the exact $5\times5$ Hamiltonian. Middle: the third-order spectrum of an
$S=1/2$ with $J\rho_s=-0.05$, $U=0$ and $T=1$ K, divided by $2\pi$
to give units of $e^2T_0^2/h$. At $B=0$ there is a Kondo peak at zero bias. Its height,
1.137, and its value at $\pm4$ mV, 0.885, match the paper's 1.13 and 0.886
(dotted lines). In a field the peak splits. The spectrum drops in the gap
$|eV|<g\mu_BB$, and the Kondo term gives two peaks just outside the
Zeeman steps. Right: the $T=0$ spectrum of a 3-site chain
from `mode="DMRG"` agrees with `mode="ED"` to 0.5 per cent (see
[ED and DMRG](#ed-and-dmrg)).

??? example "Code of this figure (docs/figures/kondo_dIdV/)"

    ```python title="compute.py"
    --8<-- "kondo_dIdV/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "kondo_dIdV/plot.py"
    ```

## ED and DMRG

`mode="ED"` (the default) diagonalizes the full
Hamiltonian, since every eigenstate enters as a virtual state. It works
at any temperature, on small chains. `mode="DMRG"`
works at $T=0$ only and never needs excited states. The second-order and
potential terms come from the [dynamical correlator](dynamical-correlators.md), and the
Kondo term from a real-time evolution of a three-point function. It runs
on `itensor_version=3` and `"python"`, on chains of at least 3 sites.
Here `sc` is a 3-site chain with couplings of order $10^{-3}$ eV:

```python
es = np.linspace(-3e-3, 3e-3, 800)
eV, dIdV = sc.get_kondo_spectrum(eV, Jrho_s=-0.05, U=0.0, T=0.0, order=3,
               omega0=2e-3, Gamma0=2e-3, mode="DMRG", submode="KPM",
               delta=2e-5, es=es,
               dt2=150.0, n_t2_half=20, dtau=150.0, n_tau_half=20)
```

These keywords have no default and must be chosen for the problem:

- `es`: the frequency grid of the dynamical correlator. For `order=2` it
  must cover the transitions below $\max|eV|$. With `U` not zero it must
  cover the whole spin excitation spectrum of the site, plus a few
  `delta`. The grid may be non-uniform.
- `delta`, `submode`: the broadening and the method of the dynamical
  correlator.
- `dt2`, `n_t2_half`, `dtau`, `n_tau_half`: the two time grids of the
  Kondo term. They must resolve the time scales $1/\omega_0$ and
  $1/\Gamma_0$, so the default $\omega_0$ and $\Gamma_0$ need very long
  grids. Too coarse a grid gives a wrong result without any error, so
  compare with `mode="ED"` on a small chain first.
- `n_gs`: the number of degenerate ground states. With a degenerate
  ground state, `mode="ED"` averages over the multiplet, while DMRG
  returns one member. `n_gs=g` makes DMRG average over `g` states too.

With `submode="TASK"`, the default on `itensor_version=3` and
`"python"`, there is no frequency or time grid at all. On
`itensor_version=2` the default stays `"KPM"` with the two-time route
above. The tangent-space Krylov method of
[Choosing a correlator method](correlator-methods.md) gives the ground
state and the excited states that the spin of the site reaches, as
discrete energies, and the exact sums of `mode="ED"` run over them:

```python
eV, dIdV = sc.get_kondo_spectrum(eV, Jrho_s=-0.05, U=0.0, T=0.0, order=3,
               omega0=2e-3, Gamma0=2e-3, mode="DMRG", submode="TASK",
               nkry=100)
```

The result is exact when `maxm` holds the whole Hilbert space (with
`sc.cutoff = 0`) and `nkry` is large enough; with a smaller `maxm` the
excited states are one-site changes of the ground state, which moves the
steps slightly. It takes only `nkry` (default 100, increase until the
spectrum stops changing) and `n_gs`.

See `examples/kondo/kondo_third_order_dmrg_VS_ED`, which is the right
panel of the figure in [Third order](#third-order), and
`examples/kondo/kondo_third_order_task_VS_ED`. The tip couples to one
site only; several sites with independent couplings are not implemented.

## Orbital-resolved IETS of a magnetic atom

`dmrgpy.atom` describes a single atom with five $d$ orbitals instead of a
chain. `generate_atom` builds a `Spinful_Fermionic_Chain` whose sites are
the orbitals, with a crystal field `tij`, spin-orbit coupling `soc`, the
interactions `U` and `J`, a magnetic field `B` acting as
$2\mathbf S+\mathbf L$, and `Ne` electrons. Two IETS signals are
available:

```python
from dmrgpy import atom
orbs = ["dxy", "dyz", "dz2", "dxz", "dx2y2"]
tij = np.diag([0.0, 0.05, 0.6, 0.05, 0.6])           # crystal field, eV
fc = atom.generate_atom(orbs=orbs, tij=tij, Ne=2, soc=0.03, U=4., J=0.5)
es = np.linspace(0., 0.1, 301)
x, y = atom.get_spinflip(fc, es=es, iorb=2, delta=1e-3)
x, y = atom.get_orbital_cotunneling(fc, es=es, iorb=2, delta=1e-3)
```

- `get_spinflip`: the spin-flip signal, $\sum_aS_{aa}(\omega)$ built from
  $S^a-\langle S^a\rangle$ on orbital `iorb`.
- `get_orbital_cotunneling`: the orbital cotunneling signal, from
  $c^\dagger_{\texttt{iorb}\,\sigma}c_{j\sigma'}$ summed over the other
  orbitals $j$ and both spins.

`iorb` is the orbital the tip couples to, for example $d_{z^2}$ for a tip
on top of the atom. With a crystal field, different orbitals give
different weights and different peak positions. Both functions use ED.
The low multiplet of an atom is often nearly degenerate, so the remarks
of [Degenerate ground states](correlator-methods.md#degenerate-ground-states) on `dex` and `T` apply here. See
`examples/kondo/atom_iets_orbital_resolved`.
