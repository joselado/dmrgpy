# Ground state

On this page we compute the ground state of a chain and measure it. All
the examples use a spin-1/2 Heisenberg chain:

```python
from dmrgpy import spinchain
n = 20
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
```

## Ground-state energy

The ground state is the state with the lowest
energy,

$$E_0=\min_{|\psi\rangle}\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}.$$

```python
e0 = sc.gs_energy()   # E0
wf = sc.get_gs()      # the ground state |GS> as an MPS
```

DMRG finds this minimum within matrix product states of bond dimension
`sc.maxm`. The result is exact when `maxm` is large enough, and an upper
bound on $E_0$ otherwise. The keywords that control accuracy are:

- `sc.maxm`: the bond dimension. Increase it until $E_0$ stops changing.
  Gapped chains need a small `maxm`; critical chains and wide ladders need
  a large one.
- `sc.nsweeps`: the number of DMRG sweeps.
- `gs_energy(maxde=1e-4)`: doubles `maxm` until the energy fluctuation per
  site ([Energy fluctuation](#energy-fluctuation)) is smaller than `maxde`.
- `mode="ED"`: exact diagonalization, for small chains and for checking
  DMRG.

The chain stores the ground state. A second call to `gs_energy()` returns
it without solving again. The sweep schedule and warm starts are on
[Controlling the solver](solver.md).

![Energy convergence with the bond dimension](img/ground-state_convergence.png#only-light)
![Energy convergence with the bond dimension](img/ground-state_convergence-dark.png#only-dark)

The figure shows the energy error and the energy fluctuation per site of
a Heisenberg chain with $n=16$, as a function of `maxm`. The reference is
ED. Both fall fast with `maxm`. The energy error follows $\delta E^2$
(dashed line). From `maxm=64` on, both stay at the numerical precision.

??? example "Code of this figure (docs/figures/ground-state_convergence/)"

    ```python title="compute.py"
    --8<-- "ground-state_convergence/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "ground-state_convergence/plot.py"
    ```

## Expectation values

For any operator $O$ built from the chain
operators, as $H$ is,

$$\langle O\rangle=\langle\mathrm{GS}|O|\mathrm{GS}\rangle,\qquad
\langle O^n\rangle=\langle\mathrm{GS}|O^n|\mathrm{GS}\rangle.$$

```python
mz = [sc.vev(sc.Sz[i]).real for i in range(n)]   # magnetization profile
e2 = sc.vev(h, npow=2)                           # <H^2>
```

Common profiles have their own methods, each returning one value per site:

- `get_density()`: the occupation $\langle n_i\rangle$, summed over spin.
  Fermionic chains.
- `get_density_fluctuation()`: $\langle n_i^2\rangle-\langle n_i\rangle^2$.
  Fermionic chains.
- `get_onsite_pairing()`: $\langle c_{i\uparrow}c_{i\downarrow}\rangle$.
  Spinful fermionic chains.
- `get_magnetization()`: $\langle S^x_i\rangle$, $\langle S^y_i\rangle$,
  $\langle S^z_i\rangle$. A $3\times n$ array on spin chains, and an
  $n\times 3$ array on spinful fermionic chains.

## Static correlators

A two-point correlator is the expectation value of
a product of two operators:

```python
czz = [sc.vev(sc.Sz[0]*sc.Sz[i]).real for i in range(n)]       # <Sz_0 Sz_i>
cdc = fc.vev(fc.Cdag[0]*fc.C[5])   # <c^dag_0 c_5> on a fermionic chain fc
```

The equal-time structure factor is its Fourier transform,
$S(q)=\sum_{i,j}e^{iq(i-j)}\langle O_iO_j\rangle$. This is the
quantity measured by neutron or X-ray diffraction, and it is the
frequency integral of the [dynamical correlators](dynamical-correlators.md).

![Spin correlations in spin-1/2 and spin-1 chains](img/ground-state_correlators.png#only-light)
![Spin correlations in spin-1/2 and spin-1 chains](img/ground-state_correlators-dark.png#only-dark)

The figure shows $|\langle S^z_iS^z_{i+r}\rangle|$ in two open Heisenberg
chains with $n=64$, measured from site $i=16$. The spin-1/2 chain is
gapless, and its correlations decay as a power law, $\sqrt{\ln r}/r$
(dashed line). The small zigzag between even and odd $r$ is an effect of
the open ends. The spin-1 chain has the Haldane gap, and its correlations decay
exponentially, as $e^{-r/\xi}/\sqrt r$ with $\xi=6.03$ (dotted line). A
free fit of the DMRG data gives $\xi=6.09$. In the spin-1 chain, the two
end sites are spin 1/2. They bind the edge states of the open chain, so
that the ground state is unique.

??? example "Code of this figure (docs/figures/ground-state_correlators/)"

    ```python title="compute.py"
    --8<-- "ground-state_correlators/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "ground-state_correlators/plot.py"
    ```

## Energy fluctuation

The energy variance tells us how close the state is
to an eigenstate of $H$:

$$\delta E=\sqrt{\langle H^2\rangle-\langle H\rangle^2}.$$

```python
de = sc.gs_energy_fluctuation()
```

For an exact eigenstate $\delta E=0$. In DMRG, $\delta E$ comes from the
truncation of the bond dimension, so it is a direct check of `maxm`. It is
a stricter check than the energy itself: roughly, the energy error is of
order $\delta E^2$, while the error of other observables is of order
$\delta E$. The figure in
[Ground-state energy](#ground-state-energy) shows this on a Heisenberg
chain.

## Fixed particle number or fixed magnetization

Without any constraint,
DMRG returns the global ground state, at whatever particle number or total
$S^z$ it has. When $H$ conserves a charge $Q$, we can instead ask for the
ground state of one sector,

$$E_0(q)=\min_{\hat Q|\psi\rangle=q|\psi\rangle}
\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}.$$

Here `fc` is a `Fermionic_Chain` and `hc` a `Spinful_Fermionic_Chain_Native`.

```python
fc.set_conserved_sector(Nf=6)          # exactly 6 fermions
sc.set_conserved_sector(Sz=0)          # total Sz = 0
hc.set_conserved_sector(Nf=8, Sz=0)    # Hubbard chain: both
fc.set_conserved_sector()              # back to the full Hilbert space
```

`Sz` is given in units of $1/2$, so `Sz=2` means $S^z_\mathrm{tot}=1$.
The available charges are `Nf` for fermions, `Sz` for spins, both for the
native Hubbard chain, and `Nb` for bosons. Parafermion chains have none.

With sectors we compute the addition spectrum $E_0(N)$, the charge gap
$E_0(N+1)+E_0(N-1)-2E_0(N)$, and the magnetization curve, each point being
a true ground state of its own sector. The same calls work with
`mode="ED"`. On large chains a sector is also faster, by about 4 times at
$n=60$.

![Magnetization curves from sector energies](img/ground-state_sectors.png#only-light)
![Magnetization curves from sector energies](img/ground-state_sectors-dark.png#only-dark)

The figure shows the magnetization curve of two open chains with $n=32$.
In a field $h$, the ground state is in the sector $M=S^z_\mathrm{tot}$ that
minimizes $E_0(M)-hM$. So the magnetization jumps from $M$ to $M+1$ at
$h_M=E_0(M+1)-E_0(M)$, and one solve per sector gives the whole curve. For
the XX chain ($S^xS^x+S^yS^y$) the steps follow the free-fermion curve
$\arcsin(h/J)/\pi$ (dashed line), and every $E_0(M)$ agrees with free
fermions on the same chain within $3\times10^{-10}$. The Heisenberg chain saturates
close to $h_s=2J$; on this open chain the last step is at
$J(1+\cos(\pi/n))=1.995\,J$.

??? example "Code of this figure (docs/figures/ground-state_sectors/)"

    ```python title="compute.py"
    --8<-- "ground-state_sectors/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "ground-state_sectors/plot.py"
    ```

Inside a sector, every operator must conserve the charge. An operator that
changes it, for example $c_i$ at fixed $N_f$ or $S^x_i$ at fixed $S^z$,
raises an error. To apply such an operator to the sector ground state,
first convert the state to an ordinary MPS:

```python
fc.set_conserved_sector(Nf=6)
fc.gs_energy()
fc.promote_to_dense()                    # same state, no sector
wf = fc.applyoperator(fc.C[3], fc.wf0)   # c_3 |GS_N>, now with 5 fermions
```

Two ground states with $N$ and $N-1$ fermions, both converted this way, can
be compared directly. This gives for example the photoemission weight
$Z_i=|\langle\mathrm{GS}_{N-1}|c_i|\mathrm{GS}_N\rangle|^2$
(see `examples/groundstate/sector_promotion_to_dense`).

!!! note "Backends"
    Sectors work on `itensor_version=3`, `itensor_version="python"` and
    `mode="ED"`. They are not available on `itensor_version=2`, on
    `"julia_live"`, with METTS and for infinite chains; these raise an error.
