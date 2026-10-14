# Dynamical correlators

On this page we compute zero-temperature dynamical correlators: the
dynamical spin structure factor, the single-particle spectral function, and
the correlator of any pair of operators $(A,B)$.

## The quantity

For a pair $(A,B)$, `get_dynamical_correlator` returns the
Lehmann spectral density

$$S_{AB}(\omega)=\sum_n M_n\,\frac{\delta/\pi}{(\omega-\Delta_n)^2+\delta^2}
=\frac{i}{2\pi}\Big[G^R_{AB}(\omega)-G^A_{AB}(\omega)\Big],$$

with $M_n=\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle$,
$\Delta_n=E_n-E_0$, and the resolvents
$G^{R/A}_{AB}(\omega)=\langle\mathrm{GS}|A\,(\omega+E_0\pm i\delta-H)^{-1}B|\mathrm{GS}\rangle$.
The frequency $\omega$ is measured from the ground-state energy $E_0$, and
$\delta$ is the broadening. When every $M_n$ is real, $S_{AB}$ is real and
equal to $-\frac1\pi\mathrm{Im}\,G^R_{AB}(\omega)$. This is always the case
for a Hermitian pair $A=B^\dagger$, where $S_{AB}\ge0$, and also for real
operators on a real Hamiltonian. In general $M_n$ is complex, so the
returned array is complex. A converged calculation satisfies the sum rule

$$\int d\omega\,S_{AB}(\omega)=\langle\mathrm{GS}|AB|\mathrm{GS}\rangle,$$

which is a simple check of a calculation. For EX and SECTOR ([Choosing a correlator method](correlator-methods.md)) the
missing weight is the part of the spectrum outside the states kept. Here $|\mathrm{GS}\rangle$ is the
ground state of the chain, or the state set with `set_gs()`.

All the examples use the Heisenberg chain `sc` of [Ground state](ground-state.md).

```python
import numpy as np
es = np.linspace(-1, 4, 300)                      # frequencies
(x, y) = sc.get_dynamical_correlator(name=(sc.Sz[0], sc.Sz[0]),
                                     es=es, delta=0.1)
print(np.sum(y.real)*(x[1]-x[0]))                 # sum rule: <Sz_0 Sz_0> = 0.25
```

The keywords that every method shares are:

- `name=(A, B)`: the operator pair, built from the chain operators. A
  string with two sites also works, `name="ZZ", i=0, j=3`, as do `"XX"`,
  `"YY"`, `"+-"`, `"cdc"` and `"densitydensity"`.
- `es`: the frequency grid.
- `delta`: the broadening. A pole becomes a line of full width $2\delta$.
- `submode`: the method, see [Choosing a correlator method](correlator-methods.md). The default is
  `"KPM"`.
- `mode="ED"`: exact diagonalization, for small chains and for checking.
  It accepts the submodes `"KPM"`, `"ED"`, `"INV"`, `"EX"`, `"CVM"`,
  `"ROOTN"` and `"TD"`.

The default `es` and `delta` change from one method to another, so it is
best to always give them. A KPM line has the width $2\delta$ at the
center of the many-body spectrum, and down to about 0.7 of it for low-energy
excitations of a long chain. It is not a Lorentzian: its peak is about 1.5
times higher than a Lorentzian of the same width and weight. So compare
peak heights only between calculations with the same method.

## Dynamical spin structure factor

The dynamical spin structure factor is

$$S^{zz}_{ij}(\omega)=\sum_n\langle\mathrm{GS}|S^z_i|n\rangle\langle n|S^z_j|\mathrm{GS}\rangle\,\delta(\omega-\Delta_n).$$

The local one, $i=j$, is what inelastic spin-flip spectroscopy with an STM
measures on a single magnetic atom. Its Fourier transform

$$S(q,\omega)=\frac1n\sum_{i,j}e^{iq(i-j)}S^{zz}_{ij}(\omega)$$

is what inelastic neutron scattering measures. We compute it from all the
site pairs:

```python
es = np.linspace(0, 4, 100)
qs = np.linspace(0, 2*np.pi, 41)
sqw = np.zeros((len(qs), len(es)))
for i in range(n):
    for j in range(n):
        (x, y) = sc.get_dynamical_correlator(name=(sc.Sz[i], sc.Sz[j]),
                                             es=es, delta=0.1)
        sqw += np.real(np.exp(1j*qs*(i-j))[:, None]*y[None, :])/n
```

For a Heisenberg chain the weight is largest at $q=\pi$, the
antiferromagnetic wave vector, and its lower edge follows the two-spinon
dispersion $\frac{\pi}{2}|\sin q|$ for long chains.

![Two-spinon continuum of the Heisenberg chain](img/dynamical-correlators_sqw.png#only-light)
![Two-spinon continuum of the Heisenberg chain](img/dynamical-correlators_sqw-dark.png#only-dark)

The figure shows $S^{zz}(q,\omega)$ of an open Heisenberg chain with
$n=32$, with KPM and $\delta=0.15$. On an open chain the momentum states
are standing waves, so we use one operator per momentum,
$S^z_q=\sqrt{2/(n+1)}\sum_i\sin(q(i+1))S^z_i$ with $q=\pi k/(n+1)$, and one
KPM run per $q$. The weight fills the two-spinon continuum between the
exact edges $\frac{\pi}{2}J|\sin q|$ and $\pi J|\sin(q/2)|$ (dashed). It is
largest just above the lower edge and near $q=\pi$. The right panel shows
two cuts, with the lower edge as a dotted line. Their integrals match
$\langle S^z_qS^z_q\rangle$ from `vev`. At every $q$ this sum rule holds
to a few parts in $10^5$.

Two methods give the spin structure factor directly:

- `sc.get_spin_spectral_function(i, j=None, es=..., delta=..., nex=20)`:
  returns $S^{zz}_{ij}+(S^{+-}_{ij}+S^{-+}_{ij})/2$, which equals
  $S^{xx}+S^{yy}+S^{zz}$. With `return_poles=True` it also returns the
  three channels separately, so we can see which spin channel each
  excitation belongs to. `connected=True` removes the elastic peak
  $\langle S^z_i\rangle\langle S^z_j\rangle$ at $\omega=0$. It uses the
  [SECTOR method](correlator-methods.md), so it needs `itensor_version=3` or `"python"`.
  The sector states are shared by all the pairs $(i,j)$, so a full
  $S(q,\omega)$ map costs about as much as one site. This is the better
  choice on a long chain. See `examples/dynamical_correlator/sector_spin_structure_factor`.
- `sc.get_full_SS_correlator(mode="ED", i=0, j=None, es=..., delta=...)`:
  returns $\frac13\sum_\alpha S^{\alpha\alpha}_{ij}(\omega)$, averaged over a
  degenerate ground-state manifold. This is what an unpolarized measurement
  on an isotropic magnet sees. ED only.

## Single-particle spectral function

For a fermionic chain the single-particle spectral function is

$$A_{ij}(\omega)=\sum_n\langle\mathrm{GS}|c_i|n^{N+1}\rangle\langle n^{N+1}|c_j^\dagger|\mathrm{GS}\rangle\,\delta(\omega-\omega_n^+)
+\sum_m\langle\mathrm{GS}|c_j^\dagger|m^{N-1}\rangle\langle m^{N-1}|c_i|\mathrm{GS}\rangle\,\delta(\omega-\omega_m^-),$$

with $\omega_n^+=E_n^{N+1}-E_0^N-\mu$ and
$\omega_m^-=E_0^N-E_m^{N-1}-\mu$. The first term is electron addition, at
$\omega>0$, and the second is electron removal, at $\omega<0$. The local
one, $A_{ii}(\omega)$, is the local density of states measured by STM
($dI/dV$); its Fourier transform $A(k,\omega)$ is measured by ARPES
(removal part) and inverse photoemission (addition part).

Here `fc` is a `Fermionic_Chain` and `hc` a
`Spinful_Fermionic_Chain_Native`:

```python
fc.set_conserved_sector(Nf=n//2)                  # the N-particle ground state
(x, A, info) = fc.get_spectral_function(i=0, es=np.linspace(-4, 4, 400),
                                        delta=0.1, nex=20, return_poles=True)
print(info["mu"], info["gap"])                    # mu and the charge gap

hc.set_conserved_sector(Nf=n, Sz=0)               # half-filled Hubbard chain
(x, A) = hc.get_spectral_function(i=0, spin="up", es=np.linspace(-6, 6, 400),
                                  delta=0.2, nex=20)
```

The chemical potential $\mu=(E_0^{N+1}-E_0^{N-1})/2$ is subtracted by
default, so the gap is centered at $\omega=0$. Pass `shift=None` to keep
the raw frequencies. `info` also gives the charge gap
$E_0^{N+1}+E_0^{N-1}-2E_0^N$ and the two halves separately. This method
needs `itensor_version=3` or `"python"`. See
`examples/dynamical_correlator/sector_spectral_function_hubbard`, where the
Mott gap opens with $U$.

On other backends, or for a broad continuum, we build the two halves with
`get_dynamical_correlator`. The second operator acts first on the ground
state, so $(c_i,c_i^\dagger)$ is the addition part and
$(c_i^\dagger,c_i)$ the removal part:

```python
es = np.linspace(0, 4, 200)
(x, yp) = fc.get_dynamical_correlator(name=(fc.C[0], fc.Cdag[0]), es=es, delta=0.1)
(x, yh) = fc.get_dynamical_correlator(name=(fc.Cdag[0], fc.C[0]), es=es, delta=0.1)
w = np.concatenate([-es[::-1], es])               # removal at w<0, addition at w>0
A00 = np.concatenate([yh.real[::-1], yp.real])
```

Here the frequency is measured from $E_0^N$ in both halves, without the
$\mu$ shift. If the chain is not in a sector, $|\mathrm{GS}\rangle$ is the
global ground state, at whatever particle number it has.

![Mott gap of the Hubbard chain](img/dynamical-correlators_spectral_function.png#only-light)
![Mott gap of the Hubbard chain](img/dynamical-correlators_spectral_function-dark.png#only-dark)

The left panel shows $A_{ii}(\omega)$ at a bulk site of a half-filled
Hubbard chain with $n=12$, from the two KPM correlators above with
$\delta=0.25$. We write the interaction as
$U\sum_i(n_{i\uparrow}-\frac12)(n_{i\downarrow}-\frac12)$, so the global
ground state is at half filling and $\mu=0$. At $U=0$ the result lies on
the exact free-fermion spectrum (dashed), drawn with the same KPM kernel. As $U$ grows, the single band
splits into a lower and an upper Hubbard band. The right panel shows the
charge gap from the ground states of the sectors $N-1$, $N$ and $N+1$. It
approaches the Lieb-Wu gap of the infinite chain (dashed) as $n$ grows,
and the finite chains stay above it because of the level spacing. For
these sectors we use $U\sum_in_{i\uparrow}n_{i\downarrow}$, and the
chemical potential comes out exactly at $U/2$, as particle-hole symmetry
requires.

## Any pair of operators

`name=(A, B)` accepts any operators built from the chain operators, for
example a bond operator:

```python
A = hc.Cdagup[0]*hc.Cup[1]
(x, y) = hc.get_dynamical_correlator(name=(A, A.get_dagger()), es=es, delta=0.1)
```

The pair may also be operators already converted with `toMPO()` (section
2). This avoids building them again in a scan over parameters. It works
with `mode="ED"` and with the DMRG submodes `"EX"`, `"CVM"` and `"ROOTN"`;
the other submodes raise an error and need the operator itself.
