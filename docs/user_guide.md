# DMRGPY Physics User Guide

This guide explains how to compute physical quantities with DMRGPY. For each
capability it gives the quantity, its formula, the code and the keywords that
matter. It does not explain how the code works inside; for that, see
`docs/documentation.md`.

Every chain class (`Spin_Chain`, `Fermionic_Chain`, `Bosonic_Chain`, ...) has
the same methods. Most calls accept `mode="DMRG"` or `mode="ED"`, so a DMRG
result on a small chain can be checked against exact diagonalization.

## Contents

1. [Physical models and Hilbert spaces](#1-physical-models-and-hilbert-spaces)
2. [Building a Hamiltonian and observables](#2-building-a-hamiltonian-and-observables)
3. [Ground-state properties](#3-ground-state-properties)
4. [Excited states and energy gaps](#4-excited-states-and-energy-gaps)
5. [Entanglement and quantum information](#5-entanglement-and-quantum-information)
6. [Dynamical (frequency-dependent) correlators](#6-dynamical-frequency-dependent-correlators)
7. [Real-time dynamics: quenches](#7-real-time-dynamics-quenches)
8. [Density of states](#8-density-of-states)
9. [Finite temperature](#9-finite-temperature)
10. [Topological invariants](#10-topological-invariants)
11. [Mean-field decoupling](#11-mean-field-decoupling)
12. [Fidelity susceptibility and quantum phase transitions](#12-fidelity-susceptibility-and-quantum-phase-transitions)
13. [Ground-state degeneracy](#13-ground-state-degeneracy)
14. [Reduced density matrices and operator distributions](#14-reduced-density-matrices-and-operator-distributions)
15. [Post-processing tools](#15-post-processing-tools)
16. [Worked-example cookbook](#16-worked-example-cookbook)
17. [STM/Kondo tunneling spectra (third-order perturbation theory)](#17-stmkondo-tunneling-spectra-third-order-perturbation-theory)
18. [Infinite chains (iDMRG)](#18-infinite-chains-idmrg)
19. [Running the pure-Python backend on a GPU](#19-running-the-pure-python-backend-on-a-gpu)
20. [Performance: BLAS threads](#20-performance-blas-threads)

## 1. Physical models and Hilbert spaces

Every model is a chain of $n$ local Hilbert spaces,
$\mathcal H=\bigotimes_{i=1}^n\mathcal H_i$. The chain class fixes the local
space and gives lists of local operators, indexed by site. The Hamiltonian
and the observables are then built from these operators (section 2).

| Chain class | Local Hilbert space | Operators |
|---|---|---|
| `spinchain.Spin_Chain(sites)` | spin $S=\tfrac12$ to $3$, one label per site (`"S=1/2"`, `"S=1"`, ...) | `Sx`, `Sy`, `Sz` |
| `fermionchain.Fermionic_Chain(n)` | spinless fermion, empty or occupied | `C`, `Cdag`, `N`, Jordan-Wigner string `F` |
| `fermionchain.Majorana_Chain(n)` | $n$ Majorana modes, built from $n/2$ fermions | `G`, with $\gamma_{2j}=(c_j+c_j^\dagger)/\sqrt2$, $\gamma_{2j+1}=i(c_j-c_j^\dagger)/\sqrt2$ |
| `fermionchain.Spinful_Fermionic_Chain(n)` | spin-1/2 fermion: $0,\uparrow,\downarrow,\uparrow\downarrow$ | `Cup`, `Cdagup`, `Cdn`, `Cdagdn`, `Nup`, `Ndn`, `Ntot`, `Sx`, `Sy`, `Sz`, `Delta` |
| `fermionchain.Spinful_Fermionic_Chain_Native(n)` | the same, one 4-state site per orbital | the same as above |
| `bosonchain.Bosonic_Chain(n, maxnb=...)` | boson with $0$ to $\text{maxnb}_i-1$ particles (default `maxnb=4`) | `A`, `Adag`, `N`, projectors `D[i][k]` $=|k\rangle\langle k|$ |
| `bosonchain.SpinBoson_Chain(sites)` | boson and spin sites mixed, labels `"B"`, `"B6"`, `"S=1/2"`, ... | `A`, `Adag`, `N`, `D` on boson sites; `Sx`, `Sy`, `Sz` on spin sites |
| `parafermionchain.Parafermionic_Chain(n, Z=3)` | $\mathbb Z_N$ clock model, $N=2,3,4$ | clock and shift `Sig`, `Tau`, parafermions `Chi`, `Psi` |
| `mixedchain.Mixed_Spin_Fermion_Chain(sites)` | spin sites and spinful-fermion sites mixed, labels `"1/2"`, `"1"`, ..., `"F"` | `Sx`, `Sy`, `Sz` everywhere; `Cup`, `Cdagup`, `Cdn`, `Cdagdn`, `Nup`, `Ndn`, `Ntot`, `Delta` on fermion sites |

The spin operators of a spinful fermion are
$S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$,
$S^+_i=c^\dagger_{i\uparrow}c_{i\downarrow}$, and the onsite pairing is
$\Delta_i=\tfrac12c_{i\uparrow}c_{i\downarrow}$. In the mixed chains, an
operator that does not exist on a site (for example `Cup` on a spin site)
is zero.

`Spinful_Fermionic_Chain` stores each orbital as two spinless sites, and
`Spinful_Fermionic_Chain_Native` as one 4-state site. The physics and the
sign conventions are the same. In practice the native class is not faster,
so `Spinful_Fermionic_Chain` is the default choice.
`Mixed_Spin_Fermion_Chain` is for models with a true local moment next to
conduction electrons, such as the Kondo lattice
(`examples/fermion_models/mixed_spin_fermion_chain`).

Some classes do not run on every backend:

- `Spinful_Fermionic_Chain_Native`, `Mixed_Spin_Fermion_Chain` and
  `SpinBoson_Chain`: `itensor_version=3` and `"python"` only.
- `SpinBoson_Chain`: no `mode="ED"`. Use `Bosonic_Chain` for a purely
  bosonic model.
- `Bosonic_Chain` with `maxnb` different from 4: not on
  `itensor_version=2`, it raises an error.

Every constructor takes `itensor_version=` and the chain settings as
keywords. A keyword has the same effect as an assignment after
construction, so

```python
from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["S=1/2"]*20, maxm=60, nsweeps=20)
```

is the same as setting `sc.maxm = 60` and `sc.nsweeps = 20` afterwards. A
keyword that is not a setting (for example `hubbard=` or `hamiltonian=`)
raises an error.

## 2. Building a Hamiltonian and observables

A Hamiltonian or an observable is written as its operator expression, with
the operator lists of section 1. For the spin-1/2 Heisenberg chain,

$$H=J\sum_{i=1}^{n-1}\mathbf S_i\cdot\mathbf S_{i+1}
=J\sum_{i=1}^{n-1}\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+S_i^zS_{i+1}^z\right),$$

```python
from dmrgpy import spinchain
n = 30
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
```

Any sum of products of local operators is allowed, with any range and
complex coefficients. This covers single-ion anisotropy $D(S^z_i)^2$,
Dzyaloshinskii-Moriya terms $\mathbf D\cdot(\mathbf S_i\times\mathbf S_{i+1})$,
biquadratic exchange $(\mathbf S_i\cdot\mathbf S_{i+1})^2$, hopping
$t\sum_i(c_i^\dagger c_{i+1}+\mathrm{h.c.})$, the Hubbard term
$U\sum_in_{i\uparrow}n_{i\downarrow}$, and so on. Section 16 gives complete
models. `h.get_dagger()` returns $h^\dagger$, so `h + h.get_dagger()` adds
the Hermitian conjugate. It is best to write $H$ in units where its largest
coefficient is of order 1: a few routes, such as real-time evolution, use
absolute thresholds.

**Shortcuts on fermionic chains.** These methods take a function `f(i,j)`
of two sites and add the whole family of terms:

- `fc.set_hoppings(f)`: $\sum_{ij}f(i,j)\,c_i^\dagger c_j$. On a
  `Spinful_Fermionic_Chain` the indices run over the spinless sites
  (up and down alternate).
- `fc.set_hoppings_spinful(f)`: the same hopping for both spins, with `f`
  indexed by orbital. Spinful chains.
- `fc.set_hubbard(f)`: $\sum_{ij}f(i,j)\,n_in_j$, with
  $n_i=n_{i\uparrow}+n_{i\downarrow}$ on spinful chains. Note that the
  $i=j$ term is $f(i,i)\,n_i^2=f(i,i)(n_i+2n_{i\uparrow}n_{i\downarrow})$.

```python
from dmrgpy import fermionchain
fc = fermionchain.Fermionic_Chain(4)
fc.set_hoppings(lambda i,j: 1.0 if abs(i-j)==1 else 0.0)
fc.set_hubbard(lambda i,j: 2.0 if abs(i-j)==1 else 0.0)
```

**Operators.** On an operator `A`, `A.simplify()` collects equal terms,
and `A.is_hermitian()`, `A.is_antihermitian()` and `A.is_zero()` answer by
reading the terms. A `True` answer is exact. A `False` answer can be wrong
when the property depends on an identity of the local space, such as
$S^xS^x=1/4$ for spin 1/2. `sc.is_hermitian(A)` does not have this
limitation, since it also checks numerically on the chain.

`sc.toMPO(A)` builds the matrix product operator of `A` once. Two such
operators can be added, subtracted and multiplied by numbers
(`A + 2*B - A`). Not available on `"julia_live"`. To combine Hamiltonians,
it is simpler to sum the operator expressions before calling
`set_hamiltonian`.

**States and operators.** The chain also gives the basic operations on
states. Here `wf` is a state, for example `wf = sc.get_gs()`, and `A` an
operator:

| Call | Returns |
|---|---|
| `sc.overlap(a, b)` | $\langle a|b\rangle$ |
| `sc.aMb(a, A, b)` | $\langle a|A|b\rangle$ |
| `sc.applyoperator(A, wf)`, or `A*wf` | $A|\psi\rangle$ |
| `sc.applyinverse(A, wf)` | $A^{-1}|\psi\rangle$, solved iteratively |
| `sc.exponential(A, wf)` | $e^{A}|\psi\rangle$ |
| `sc.summps(a, b)` | $|a\rangle+|b\rangle$ |
| `sc.scale_mps(x, wf)`, or `wf*x` | $x|\psi\rangle$ |
| `sc.trace(A)`, `sc.inverse_trace(A)` | $\mathrm{Tr}\,A$, $\mathrm{Tr}\,A^{-1}$ |
| `sc.operator_norm(A)` | an estimate of $\lVert A\rVert$ |
| `sc.is_zero_operator(A)` | whether $A=0$ |
| `wf.dot(wf2)`, `wf.norm()` | $\langle\psi|\psi_2\rangle$, $\lVert\psi\rVert$ |
| `wf.normalize()` | $|\psi\rangle/\lVert\psi\rVert$ |
| `wf.get_conjugate()` | the complex conjugate state |

A state remembers the backend that made it (an MPS for DMRG, a vector for
`mode="ED"`), and these operations run on that backend. Give `mode="ED"`
to `overlap` and `aMb` when the states come from ED. On DMRG,
`exponential` needs a Hermitian or an anti-Hermitian `A`, so
`sc.exponential(-1j*t*h, wf)` is a real-time step. For a general `A` use
`mode="ED"` on a small chain.

## 3. Ground-state properties

In this section we compute the ground state of a chain and measure it. All
the examples use a spin-1/2 Heisenberg chain:

```python
from dmrgpy import spinchain
n = 20
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
        for i in range(n-1))
sc.set_hamiltonian(h)
```

**Ground-state energy.** The ground state is the state with the lowest
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
  site (see below) is smaller than `maxde`.
- `mode="ED"`: exact diagonalization, for small chains and for checking
  DMRG.

The chain stores the ground state. A second call to `gs_energy()` returns
it without solving again.

**Expectation values.** For any operator $O$ built from the chain
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
  $\langle S^z_i\rangle$ as a $3\times n$ array. Spin and spinful fermionic
  chains.

**Static correlators.** A two-point correlator is the expectation value of
a product of two operators:

```python
czz = [sc.vev(sc.Sz[0]*sc.Sz[i]).real for i in range(n)]       # <Sz_0 Sz_i>
cdc = fc.vev(fc.Cdag[0]*fc.C[5])   # <c^dag_0 c_5> on a fermionic chain fc
```

The equal-time structure factor is its Fourier transform,
$S(q)=\sum_{i,j}e^{iq(i-j)}\langle O_iO_j\rangle$. This is the
quantity measured by neutron or X-ray diffraction, and it is the
frequency integral of the dynamical correlators of section 6.

**Energy fluctuation.** The energy variance tells us how close the state is
to an eigenstate of $H$:

$$\delta E=\sqrt{\langle H^2\rangle-\langle H\rangle^2}.$$

```python
de = sc.gs_energy_fluctuation()
```

For an exact eigenstate $\delta E=0$. In DMRG, $\delta E$ comes from the
truncation of the bond dimension, so it is a direct check of `maxm`. It is
a stricter check than the energy itself: roughly, the energy error is of
order $\delta E^2$, while the error of other observables is of order
$\delta E$.

**Fixed particle number or fixed magnetization.** Without any constraint,
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

Sectors work on `itensor_version=3`, `itensor_version="python"` and
`mode="ED"`. They are not available on `itensor_version=2`, on
`"julia_live"`, with METTS and for infinite chains; these raise an error.

## 4. Excited states and energy gaps

In this section we compute the low-energy spectrum: excited states, gaps,
non-Hermitian and generalized eigenproblems, and effective Hamiltonians of
a low-energy manifold. Here `sc` is the Heisenberg chain of section 3.

**Low-lying spectrum.** The lowest eigenvalues $E_0\le E_1\le\dots$ of $H$
and the gap $\Delta=E_1-E_0$ are

```python
es = sc.get_excited(n=4)                 # [E0, E1, E2, E3]
gap = sc.get_gap()                       # E1 - E0
es, ws = sc.get_excited_states(n=4)      # energies and states
```

The list has one entry per eigenstate, so a degenerate level appears once
for each state of the multiplet. If the ground state is degenerate,
`get_gap()` is zero, and we must ask for enough states to reach the first
excited level. When DMRG splits a multiplet slightly, identify the states
by a quantum number (for example $\langle\mathbf S^2\rangle$) and not by
their position in the list.

DMRG finds each excited state as the ground state of $H$ plus a penalty
on the states already found. This search can converge to a wrong state.
When the energy variance of a state is much larger than that of the ground
state, a warning is printed, and the state should be checked with more
sweeps, a larger `maxm` or `mode="ED"` on a smaller chain. The keywords are:

- `n`: the number of states.
- `scale`: the weight of the penalty, in units of the bandwidth of $H$
  (default 10).
- `purify=True` (default, `get_excited_states`): diagonalizes $H$ again in
  the space of the states found, which separates nearly degenerate states.

For example, an open spin-1 Heisenberg chain has four low states from the
two spin-1/2 edge modes (a singlet and a triplet), and the Haldane gap is
the gap to the next state:

```python
es = s1.get_excited(n=6)    # s1: an open spin-1 Heisenberg chain
haldane_gap = es[4] - es[0]
```

**Gap to another particle number.** For a conserved charge $A$, for
example the particle number $N=\sum_in_i$, `get_charge_gap(d=...)` gives
the energy to change $\langle A\rangle$ by $d$ from its ground-state value
$A_0$. It solves

$$H_\lambda=H+\lambda\big(A-A_0-d\big)^2,$$

increasing $\lambda$ until the constraint holds, and returns
$\Delta(d)=E_0(A_0+d)-E_0(A_0)$.

```python
gap2 = fc.get_charge_gap(d=2)    # fc: a fermionic chain; add two particles
```

`d=2` is the cost of adding a pair, and `d=-2` of removing one.
`Fermionic_Chain`, `Spinful_Fermionic_Chain` and the native spinful chain
have this method. The conserved sectors of section 3 give the same energies
directly, $E_0(N)$ for each $N$.

**Non-Hermitian Hamiltonians.** When $H\neq H^\dagger$ (gain and loss,
complex or non-reciprocal hopping), the ground state is the eigenvalue with
the smallest real part. Eigenstates come in left and right pairs,

$$H|\psi_R\rangle=E|\psi_R\rangle,\qquad
H^\dagger|\psi_L\rangle=E^*|\psi_L\rangle,\qquad
\langle\psi_L|\psi_R\rangle=1.$$

`gs_energy()` detects a non-Hermitian $H$ and uses non-Hermitian DMRG
(NH-DMRG). Here `fc` is a fermionic chain with a non-Hermitian Hamiltonian,
for example a staggered imaginary potential $i\gamma(-1)^in_i$:

```python
e, psil, psir = fc.nhdmrg()     # E, left and right eigenvectors
e = fc.gs_energy()              # the same E; stores psir as the ground state
```

The energy is not a variational bound, so the result is checked with the
residuals $\lVert H|\psi_R\rangle-E|\psi_R\rangle\rVert$ and
$\lVert H^\dagger|\psi_L\rangle-E^*|\psi_L\rangle\rVert$. The keywords are:

- `tol`: the largest relative residual accepted (default `1e-4`).
- `ntries`: how many times to restart from a random state when the residual
  is too large (default 5).

Excited states of a non-Hermitian $H$ come from `get_excited_states(n)` and
`get_excited(n)`, which use an Arnoldi method on matrix product states. It
is less accurate than NH-DMRG, so check it against `mode="ED"` on a small
chain (`examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi`). NH-DMRG works on
`itensor_version=2`, `3`, `"python"` and `"julia_live"`.

**Generalized eigenvalue problem.** For a Hermitian positive definite
operator $A$,

$$H|\psi\rangle=\lambda A|\psi\rangle .$$

```python
A = fc.Id + 0.5*sum(fc.N)           # for example
lam = fc.gs_energy_generalized(A)   # smallest lambda
```

DMRG solves this directly, with no inverse of $A$. A non-Hermitian $H$ is
also accepted, and then $\lambda$ is the complex eigenvalue with the
smallest real part. Not available on `itensor_version=2`, it raises an
error.

After this call the chain stores the generalized eigenstate as its ground
state, and `fc.e0` is its energy
$\langle\psi|H|\psi\rangle/\langle\psi|\psi\rangle$ (with left and right
states for a non-Hermitian $H$), while $\lambda$ is in
`fc.lam_generalized`. Every later calculation (excited states,
correlators) uses this state. So call it last, or call `fc.restart()` and
`fc.gs_energy()` to get back the true ground state.

**Effective low-energy Hamiltonians.** Take the $n$ lowest eigenstates and
the projector $P=\sum_{k<n}|\psi_k\rangle\langle\psi_k|$ on them. We write
the projected Hamiltonian as a combination of projected operators $O_a$,

$$PHP\simeq\sum_aJ_a\,PO_aP,$$

and fit the couplings $J_a$ by least squares. This gives, for example, the
spin exchange of a Hubbard model at strong coupling,
$J\approx4t^2/U$:

```python
# fc: a Spinful_Fermionic_Chain(2) Hubbard dimer
J = fc.get_heff(operators=[fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1],
                           fc.Sz[0]*fc.Sz[1]], n=4)    # one J per operator

from dmrgpy import effectivehamiltonian
coef = effectivehamiltonian.get_effective_hamiltonian_couplings(fc, n=4)
latex = effectivehamiltonian.get_effective_hamiltonian(fc, n=4)
```

`get_heff` fits the operators given. `get_effective_hamiltonian_couplings`
fits all products $S^\alpha_iS^\beta_j$ of the chain spin operators and
returns a dictionary, and `get_effective_hamiltonian` returns the same fit
as a LaTeX string. The keywords are:

- `n`: the number of states in the manifold. It must not cut a degenerate
  multiplet, otherwise the fit is not meaningful and it raises an error.
- `operators`: the operators of the fit. They must be linearly independent
  in the manifold, otherwise the couplings are not unique and it raises an
  error.
- `tol`: couplings smaller than this are set to zero.
- `method`: `"single"` (default) fits $PO_aP\,PO_bP$, and `"full"` fits
  $P\,O_aO_b\,P$. They differ when the manifold is not closed under $O_a$.
- `mode="ED"`: exact eigenstates, for small systems.

See `examples/groundstate/effective_hamiltonian`, which fits the twisted
exchange of a Hubbard dimer with a Peierls phase.

## 5. Entanglement and quantum information

In this section we measure the entanglement of a state. Here `wf` is a
ground state, for example `wf = sc.get_gs()` for the Heisenberg chain of
section 3.

**Entanglement entropy of a cut.** We cut the chain between sites $i$ and
$i+1$. The Schmidt decomposition of the state is
$|\psi\rangle=\sum_\alpha\lambda_\alpha|\alpha\rangle_L|\alpha\rangle_R$,
and the von Neumann entropy of the left part is

$$S_i=-\sum_\alpha\lambda_\alpha^2\log\lambda_\alpha^2
=-\mathrm{Tr}\,\rho_L\log\rho_L,\qquad
\rho_L=\mathrm{Tr}_R|\psi\rangle\langle\psi|.$$

```python
s = wf.get_bond_entropy(i, i+1)     # cut between sites i and i+1
s = wf.get_entropy(b)               # the same, cut between sites b-1 and b
```

In a gapped chain $S_i$ saturates with the size of the block, and in a
critical chain it grows as $\log\ell$. It also tells us the bond dimension
we need, roughly $\chi\gtrsim e^{S}$. The bond entropy needs an MPS, so it
is not available with `mode="ED"`.

**Entropy of one site or two sites.** `wf.get_site_entropy(i)` and
`wf.get_pair_entropy(i, j)` give the entropy of the reduced density matrix
of site $i$, or of the two sites $i$ and $j$ (not necessarily neighbors),
with the rest of the chain. They are available for spin-1/2 chains and
spinless fermionic chains. With `mode="ED"`, call them on the chain:
`sc.get_site_entropy(wf, i)`, `sc.get_pair_entropy(wf, i, j)`.

**Mutual information.** For two sites $i$ and $j$,

$$I(i,j)=S_i+S_j-S_{ij},$$

```python
mi = wf.get_mutual_information(i, j)
```

$I(i,j)$ measures all the correlations, classical and quantum, between the
two sites, so it shows which sites are coupled in the state even when no
term of $H$ couples them. Same availability as the site entropy.

**Central charge.** For a critical open chain of length $L$, the entropy
of a cut at position $\ell$ follows the Calabrese-Cardy formula,

$$S(\ell)=\frac c6\log\!\left[\frac{2L}{\pi}\sin\frac{\pi\ell}{L}\right]+\text{const}.$$

The central charge $c$ identifies the conformal field theory, for example
$c=1/2$ for the critical Ising chain and $c=1$ for the XX and Heisenberg
chains. `get_CFT_central_charge()` fits this formula to the entropy of
every bond:

```python
c = wf.get_CFT_central_charge()
```

On an open chain the entropy also has a term that alternates between even
and odd bonds, and the fit does not remove it. So the result converges
slowly with $L$: for the spin-1/2 XX chain it gives about 0.54 at $L=20$
and 0.74 at $L=40$. Use long chains, or fit the bond entropies yourself
with only even (or only odd) $\ell$.

**Single-particle correlation matrix.** For fermions, the correlation
matrix is

$$C_{ij}=\langle c_i^\dagger c_j\rangle .$$

Its eigenvalues $n_\alpha\in[0,1]$ are the occupations of the natural
orbitals. For a single Slater determinant every $n_\alpha$ is 0 or 1.
Interactions move them away from 0 and 1, and the correlation entropy
measures how far the state is from a single Slater determinant:

$$S=-\sum_\alpha\Big[n_\alpha\log n_\alpha+(1-n_\alpha)\log(1-n_\alpha)\Big].$$

Here `fc` is a fermionic chain with its Hamiltonian set.

```python
cm = fc.get_correlation_matrix()          # C_ij
na = fc.get_correlation_eigenvalues()     # n_alpha
s = fc.get_correlation_entropy()          # S
orb = fc.get_correlated_orbitals()        # |orbital|^2, most correlated first
rho = fc.get_correlated_density(n=2)      # weight of the 2 most correlated orbitals
sd = fc.get_gs().get_correlation_entropy_density()   # S resolved by site
```

The most correlated orbitals are those with $n_\alpha$ closest to $1/2$.
The keywords are:

- `T`: temperature of `get_correlation_matrix` (default 0). Finite $T$ is
  computed with exact eigenstates, for small chains.
- `dmmode`: how $C_{ij}$ is computed (`"fast"`, `"full"`, `"explicit"`,
  `"simple"`). The default chooses the fastest one for the state.
  `"explicit"` works on every backend.

**Two-particle correlators.** The four-point tensor
$\langle c_i^\dagger c_jc_k^\dagger c_l\rangle$ for all $i,j,k,l$ is

```python
t = fc.get_gs().get_four_correlation_tensor()   # shape (n, n, n, n)
```

The keyword `ctmode` chooses the method (`"batched"`, `"sweep"`, `"fold"`,
`"full"`, `"explicit"`). The default chooses the fastest one for the
backend, and `"explicit"` works on every backend.
The two-particle density matrix
$\langle c_i^\dagger c_j^\dagger c_lc_k\rangle$, with the pairs $i<j$ and
$k<l$ as row and column, is

```python
from dmrgpy.entropytk.correlationentropy import get_highorder_correlation_matrix
d2 = get_highorder_correlation_matrix(fc)
```

## 6. Dynamical (frequency-dependent) correlators

In this section we compute zero-temperature dynamical correlators: the
dynamical spin structure factor, the single-particle spectral function, and
the correlator of any pair of operators $(A,B)$.

**The quantity.** For a pair $(A,B)$, `get_dynamical_correlator` returns the
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

which is a simple check of a calculation. For EX and SECTOR (below) the
missing weight is the part of the spectrum outside the states kept. Here $|\mathrm{GS}\rangle$ is the
ground state of the chain, or the state set with `set_gs()`.

All the examples use the Heisenberg chain `sc` of section 3.

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
- `submode`: the method, see "Choosing a method" below. The default is
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

### Dynamical spin structure factor

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

Two methods give the spin structure factor directly:

- `sc.get_spin_spectral_function(i, j=None, es=..., delta=..., nex=20)`:
  returns $S^{zz}_{ij}+(S^{+-}_{ij}+S^{-+}_{ij})/2$, which equals
  $S^{xx}+S^{yy}+S^{zz}$. With `return_poles=True` it also returns the
  three channels separately, so we can see which spin channel each
  excitation belongs to. `connected=True` removes the elastic peak
  $\langle S^z_i\rangle\langle S^z_j\rangle$ at $\omega=0$. It uses the
  SECTOR method below, so it needs `itensor_version=3` or `"python"`.
  The sector states are shared by all the pairs $(i,j)$, so a full
  $S(q,\omega)$ map costs about as much as one site. This is the better
  choice on a long chain. See `examples/dynamical_correlator/sector_spin_structure_factor`.
- `sc.get_full_SS_correlator(mode="ED", i=0, j=None, es=..., delta=...)`:
  returns $\frac13\sum_\alpha S^{\alpha\alpha}_{ij}(\omega)$, averaged over a
  degenerate ground-state manifold. This is what an unpolarized measurement
  on an isotropic magnet sees. ED only.

### Single-particle spectral function

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

### Any pair of operators

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

### Choosing a method

The `submode` keyword selects how $S_{AB}(\omega)$ is computed. All the
methods return the same quantity, defined above.

- `"KPM"` (default): the Kernel Polynomial Method. It expands the spectrum
  in Chebyshev polynomials of $H$ and gives the whole frequency range at
  once. The cost grows as $1/\delta$. Use it for a first look and for broad
  spectra.
- `"CVM"`: the correction-vector method. It solves one linear system per
  frequency, so it is accurate at each $\omega$ but the cost grows with
  the number of points in `es`. Use it to resolve a narrow window or a
  sharp peak. `"CVM_explicit"` computes the same resolvent in a different
  way, and is useful as a cross-check.
- `"ROOTN"`: the root-$N$ Krylov correction vector of Nocera and Alvarez
  ([arXiv:2204.03165](https://arxiv.org/abs/2204.03165)). It applies the
  resolvent in `N` smaller steps, each in a Krylov space of size `nkry`.
  Use it for high frequencies, near the top of the spectrum, where a
  single correction vector needs a large bond dimension. It is the most
  expensive method.
- `"TD"`: real-time evolution of $C(t)=\langle\mathrm{GS}|A(t)B(0)|\mathrm{GS}\rangle$
  followed by a Fourier transform with the damping $e^{-\delta t}$. Use it
  for fine resolution and when time evolution is cheap (short-range
  Hamiltonians, small entanglement growth). The time propagator is the one
  of section 7.
- `"TDZ"`: the same, but the evolution follows a complex-time contour (Cao,
  Lu, Stoudenmire and Parcollet,
  [arXiv:2311.10909](https://arxiv.org/abs/2311.10909)). The bond dimension
  grows much more slowly than in `"TD"`, so it reaches longer times.
- `"EX"`: an explicit sum over the lowest `nex` excited states found by
  DMRG. Use it when a few low-energy states carry the physics, for example
  in a small gapped system.
- `"SECTOR"`: an explicit sum over the lowest `nex` eigenstates of the
  quantum-number sector that $B|\mathrm{GS}\rangle$ belongs to ($N+1$ for
  $c^\dagger$, $S^z+1$ for $S^+$). Use it for sharp low-energy
  excitations: the poles come out at converged energies with exact weights.
  It needs `itensor_version=3` or `"python"`. `get_spectral_function` and
  `get_spin_spectral_function` are built on it.
- `"ED"` and `"INV"` (`mode="ED"` only): the exact Lehmann sum and the
  exact resolvent. These are the references to check the other methods.
- `"maxent"` and `"CVMimag"`: need the separate modules `dmrgpy.maxenttk`
  and `dmrgpy.padetk`, which are not distributed with dmrgpy. Without them
  they raise an error.

As a rule: KPM for the full spectrum, CVM or TD for a narrow window at
high resolution, TDZ when TD needs too large a bond dimension, EX or
SECTOR when a few states are enough.

```python
pair = (sc.Sz[0], sc.Sz[0])
(x, y) = sc.get_dynamical_correlator(submode="CVM", name=pair, es=es, delta=0.1)
(x, y) = sc.get_dynamical_correlator(submode="TD", name=pair, es=es, delta=0.1, dt=0.05)
(x, y) = sc.get_dynamical_correlator(submode="EX", name=pair, es=es, delta=0.1, nex=20)
(x, y) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=pair, es=es, delta=0.1)
```

**Keywords of each method.** The accuracy of every method also depends on
the ground state itself, so `maxm` and `nsweeps` (section 3) come first.

- `sc.kpmmaxm` (KPM): the bond dimension of the Chebyshev vectors. Increase
  it until the spectrum stops changing.
- `sc.kpm_scale` (KPM, default 0.7): how much of $[-1,1]$ the spectrum of
  $H$ fills after rescaling. Do not go below 0.5; the calculation raises an
  error when the expansion becomes unstable.
- `sc.kpm_n_scale` (KPM, default 1): a positive integer that multiplies the
  number of Chebyshev moments, for a line sharper than $2\delta$.
- `kernel=` (KPM): how the moments become a spectrum, `"jackson"`
  (default), `"lorentz"`, `"plain"` or `"hodc"`. The last is the high-order
  kernel of [arXiv:2512.03149](https://arxiv.org/abs/2512.03149), with
  `hodc_order` (default 6) and `hodc_eta` (default `delta`). It is more
  accurate where the spectrum is a smooth continuum, but it can be slightly
  negative and is more sensitive to a small `kpmmaxm`. See
  `examples/dynamical_correlator/hodc_VS_jackson_kernel`.
- `sc.kpm_energy_truncate=True` (KPM, `itensor_version=3` and `"python"`):
  energy truncation of Holzner et al.,
  [PRB 83, 195115 (2011)](https://doi.org/10.1103/PhysRevB.83.195115). It
  allows `kpm_scale` well below 0.5, so that more moments are spent on the
  low-energy window. It is slower per moment. The settings are
  `kpm_truncate_dK` (increase until converged, 30 is a good value),
  `kpm_truncate_nsweeps` and `kpm_truncate_threshold`. See
  `examples/dynamical_correlator/dynamical_correlator_kpm_energy_truncation`.
- `sc.cvm_maxm` and `sc.cvm_tol` (CVM): the bond dimension of the
  correction vector and the relative tolerance of its solver (default
  $10^{-5}$). If a frequency point warns that it did not converge,
  increase `cvm_maxm`.
- `N` and `nkry` (ROOTN, defaults 8 and 20): the number of steps and the
  Krylov dimension. Larger `N` helps at high frequency when `nkry` is
  small.
- `dt` and `damping_periods` (TD and TDZ): the time step, and the total
  time $T$ in units of $1/\delta$ (default 6). The resolution is set by
  $T$.
- `damping` and `predict` (TD): the time window, `"exp"` (default, a
  Lorentzian line), `"gaussian"` (faster decaying tails, line 1.2 times
  wider) or `"parzen"`, and linear prediction of $C(t)$ to longer times
  (default `True`), which makes the line sharper and more accurate when
  the simulated time is short.
- `alpha0` and `n_max` (TDZ, defaults 0.1 and 4): the angle of the time
  contour and the order of the expansion back to real time. A larger
  `alpha0` needs a smaller bond dimension but a larger `n_max`, at most 4.
- `nex` (EX and SECTOR): the number of states in the sum. With
  `return_poles=True`, SECTOR also returns the poles and their weights,
  and `info["captured"]` is the fraction of the weight
  $\langle\mathrm{GS}|B^\dagger B|\mathrm{GS}\rangle$ they account for. A
  warning appears below 0.9; increase `nex` then.

**Reusing the KPM moments.** The Chebyshev moments are the expensive part
of KPM. We can compute them once and build the spectrum with several
kernels:

```python
from dmrgpy import kpmdmrg
mus, emin, emax, scale, npol, delta = sc.get_dynamical_correlator_moments(
        name=(sc.Sz[0], sc.Sz[0]), delta=0.1)
(x, y) = kpmdmrg.dynamical_correlator_from_moments(mus, emin, emax, scale,
        npol, es, kernel="hodc", delta=delta)
```

### Degenerate ground states

With `mode="ED", submode="ED"`, a degenerate ground state is handled by
averaging the correlator over all the states with excitation energy below
`dex` (default $10^{-5}$), with equal weights. This is a hard cutoff, so in
a scan where a level splitting crosses `dex` (a Zeeman field, for example)
the result jumps. A warning appears when a level is close to `dex`. For a
scan, use a small temperature instead, `T=0.01`, which weights the states
with Boltzmann factors and changes smoothly (section 9). See
`examples/kondo/atom_iets_orbital_resolved`.

### Non-Hermitian Hamiltonians

For a non-Hermitian $H$ (section 4) the ground state is a pair of left and
right eigenstates, and the correlator is
$\langle\psi_L|A(z)B|\psi_R\rangle$. `submode="KPM"` then uses the
non-Hermitian KPM of [PRL 130, 100401](https://doi.org/10.1103/PhysRevLett.130.100401),
and needs `E_max`, an upper bound on the spectral radius of $H$. It
recomputes the moments at every frequency, so it is much slower. `"CVM"`,
`"CVM_explicit"` and `"EX"` also work; the other submodes raise an error.
Non-Hermitian KPM works with `mode="ED"`, `itensor_version=3` and
`"python"`. See `examples/non_hermitian/nhkpm_v3_VS_ED`.

## 7. Real-time dynamics: quenches

In this section we evolve a state in real time,
$|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle$, and measure an observable along
the way. The typical case is a quantum quench: we prepare the ground state
of a Hamiltonian $H_0$, change the Hamiltonian to $H_1$ at $t=0$, and follow
how observables oscillate, relax or thermalize. Here `h0` and `h1` are two
Hamiltonians of the chain `sc`, for example a staggered field and the
Heisenberg chain of section 3.

```python
from dmrgpy import timedependent
sc.set_hamiltonian(h0)
wf0 = sc.get_gs()                          # ground state of H0
sc.set_hamiltonian(h1)                     # quench to H1
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                            nt=200, dt=0.02, wf=wf0)
```

This gives $\langle\psi(t)|S^z_0|\psi(t)\rangle$ at the times
$t=0,\,dt,\,2dt,\dots$. The keywords are:

- `nt`, `dt`: the number of time steps (default 1000) and the time step
  (default 0.01). A smaller `dt` is more accurate.
- `wf`: the initial state. Without it, the ground state of the chain is
  used.
- `h`: the Hamiltonian of the evolution, if not the one of the chain.
- `mode="ED"`: exact evolution, the reference for small chains.
- `return_wf=True` (DMRG): also return the final state.

The accuracy of the DMRG evolution is set by `sc.maxm` and `dt`. The
entanglement grows with time after a quench, so the bond dimension needed
grows too. When the result stops agreeing with a larger `maxm`, the
simulation has reached its maximum time.

**Local quench.** `evolution_ABA` first applies an operator $A$ to the
state, for example a spin flip or an added electron, and then measures $B$:

$$\langle\psi|A^\dagger e^{iHt}Be^{-iHt}A|\psi\rangle.$$

Here we flip up the spin at the edge of the Heisenberg chain and follow
how the magnetization arrives at site 3:

```python
Sp0 = sc.Sx[0] + 1j*sc.Sy[0]               # S^+ on site 0
(ts, sz3) = timedependent.evolution_ABA(sc, A=Sp0, B=sc.Sz[3], nt=200, dt=0.02)
```

It takes the same keywords as `evolve_and_measure`, and starts from the
ground state of the chain unless `wf` is given. Note that $A$ is applied as
it is, so the state $A|\psi\rangle$ is not normalized.

**Evolved states.** To get the states themselves, `sc.exponential(-1j*t*h,
wf)` returns $e^{-iHt}|\psi\rangle$ for one time (section 2), and
`timeevolution.evolve_WF(h, wf, ts=ts)` returns a list of states, one for
each time in `ts`. Note that `evolve_WF` applies $e^{+iht}$, so pass `-h`
for the physical evolution. With the states, we compute any overlap, for
example the return amplitude $\langle\psi(0)|\psi(t)\rangle$.

**Choosing the propagator.** `sc.tevol_method` selects the algorithm:

- `"TDVP"` (default): two-site TDVP. It works for any Hamiltonian and
  grows the bond dimension as needed. On `itensor_version=2` it is
  replaced by `"MPO"`.
- `"TDVP_GSE"`: one-site TDVP with global subspace expansion (Yang and
  White, [PRB 102, 094315 (2020)](https://arxiv.org/abs/2005.06104)). Use
  it when the initial state has a small bond dimension, for example a
  product state. The bond dimension grows only during the first
  `sc.tdvp_gse_sweeps` steps (default 3); the other settings are
  `sc.tdvp_gse_krylov_order` (default 3) and `sc.tdvp_gse_cutoff`
  (default $10^{-8}$).
- `"TEBD"`: second-order Trotter evolution with two-site gates. Faster
  than TDVP, but only for nearest-neighbor Hamiltonians; otherwise it
  raises an error. A `Spinful_Fermionic_Chain_Native` Hubbard chain is
  nearest-neighbor, the interleaved `Spinful_Fermionic_Chain` is not.
- `"AUTO"`: `"TEBD"` when the Hamiltonian is nearest-neighbor, `"TDVP"`
  otherwise.
- `"MPO"`: a second-order Taylor expansion of $e^{-iH\,dt}$ applied as an
  MPO. The only option on `itensor_version=2`, and less accurate than TDVP
  at the same bond dimension.

```python
sc.tevol_method = "TEBD"
```

The same propagator is used by the `"TD"` and `"TDZ"` dynamical
correlators of section 6. See `examples/time_evolution/time_evolution_quench`.

## 8. Density of states

"Density of states" can mean two different quantities, and only the second
one has a method.

**Many-body density of states.** This is the density of all the
eigenstates of $H$,

$$\rho(E)=\sum_n\delta(E-E_n).$$

There is no method for it. `get_distribution` (section 14) is not a
replacement: it computes a ground-state expectation value
$\langle\mathrm{GS}|B\,\delta(X)\,A|\mathrm{GS}\rangle$, not a trace over
all the states. For a small chain, the full ED spectrum from
`get_excited(mode="ED", n=...)` (section 4) gives $\rho(E)$ directly.

**Local density of states.** For fermions, the density of states measured
by STM ($dI/dV$) is the local single-particle spectral function

$$A_i(\omega)=\sum_n|\langle n^{N+1}|c_i^\dagger|\mathrm{GS}\rangle|^2\delta(\omega-\omega_n^+)
+\sum_m|\langle m^{N-1}|c_i|\mathrm{GS}\rangle|^2\delta(\omega-\omega_m^-),$$

with electron addition at $\omega>0$ and electron removal at $\omega<0$.
This is $A_{ii}(\omega)$ of section 6, and we compute it with
`get_spectral_function(i=i)`, or from the two correlators
$(c_i,c_i^\dagger)$ and $(c_i^\dagger,c_i)$ as shown there. On a spinful
chain, pass `spin="up"` or `spin="dn"` and add the two.

## 9. Finite temperature

In this section we compute thermal averages,

$$\langle O\rangle_T=\frac{\mathrm{Tr}\,[e^{-\beta H}O]}{\mathrm{Tr}\,e^{-\beta H}},
\qquad\beta=1/T.$$

There are three methods: purification, METTS and exact diagonalization.

**Purification.** Each physical site gets an ancilla site. The product of
singlets between each site and its ancilla is a pure state whose reduced
density matrix on the physical chain is $\rho\propto\mathbb 1$, the
$T=\infty$ state. Evolving it in imaginary time with the physical $H$,

$$|\Psi(\beta)\rangle=e^{-\beta H/2}|\Psi(0)\rangle,$$

and tracing out the ancillas gives $\rho(T)\propto e^{-\beta H}$.

```python
from dmrgpy import thermal
n = 6
tc = thermal.Thermal_Spin_Chain(["S=1/2"]*n, T=0.5)
h = sum(tc.Sx[i]*tc.Sx[i+1] + tc.Sy[i]*tc.Sy[i+1] + tc.Sz[i]*tc.Sz[i+1]
        for i in range(n-1))
tc.set_hamiltonian(h)
tc.get_gs()                                   # builds |Psi(beta)>
E = tc.MBChain.vev(h).real                    # thermal energy
czz = tc.MBChain.vev(tc.Sz[0]*tc.Sz[1]).real  # thermal <Sz_0 Sz_1>
```

`tc.Sx`, `tc.Sy`, `tc.Sz` are the operators of the physical sites, and
`tc.MBChain` is the doubled chain that holds the purified state. Static
expectation values on `tc.MBChain` are thermal averages. A dynamical
correlator on `tc.MBChain` is not the thermal one, so use METTS or ED for
those. `T=0` gives the ordinary ground state. The imaginary-time evolution
is accurate to about $10^{-6}$ relative at the defaults. The keywords are:

- `tc.anneal_step`: the size of each imaginary-time step, in units of the
  spectral width (default 2). Smaller is more accurate and slower.
- `tc.anneal_order`: the order of the expansion of each step (default 8).
- `tc.MBChain.maxm`: the bond dimension of the purified state. It needs
  to be larger at intermediate temperatures.

See `examples/finite_temperature/thermal_purification_VS_exact`.

**METTS.** Minimally entangled typical thermal states (Stoudenmire and
White, arXiv:1002.1305) sample product states $|i\rangle$ with a Markov
chain. Each one is evolved to
$|\phi_i\rangle\propto e^{-\beta H/2}|i\rangle$ and then collapsed to a new
product state. The plain average of $\langle\phi_i|O|\phi_i\rangle$ over
the chain is the thermal average. The states stay weakly entangled, so
no ancillas are needed. Here `sc` is the Heisenberg chain of section 3:

```python
T = 0.5
mean, err = sc.metts_vev(sc.Sz[0]*sc.Sz[1], T, nsamples=300, nwarmup=30)
res = sc.metts_vev([sc.Sz[0], sc.Sz[1], h], T)   # list of (mean, err)
```

Being a Monte Carlo method, it returns a mean and a statistical error.
The error assumes independent samples, so it is too small when samples
are correlated. Measuring several operators in one call is much cheaper
than one call per operator, since they share the samples. The keywords
are:

- `nsamples`, `nwarmup`: the number of samples, and the number of first
  steps discarded for equilibration.
- `dbeta_half_step`: the imaginary-time step (default 0.05).
- `basis_ops`: the two bases used in turn for the collapse (default
  `("Sz","Sx")`). Alternating them reduces autocorrelation.
- `njobs`: runs independent Markov chains in parallel processes. Only
  for `itensor_version="python"`.

METTS works with `itensor_version=3` and `"python"`. It is not available
on `itensor_version=2`, where it raises an error.

The same sampling gives real-time thermal correlators
$C_{AB}(t)=\langle A(t)B\rangle_T$ (arXiv:2405.18484). Each sample
$|\psi_i\rangle$ gives $\langle\psi_i|e^{iHt}Ae^{-iHt}B|\psi_i\rangle$ by
real-time evolution, and the average is $C_{AB}(t)$:

```python
ts, c, err = sc.metts_dynamical_correlator((sc.Sz[0], sc.Sz[0]), T,
                 nt=100, dt=0.1, nsamples=200, nwarmup=30)
```

The result is in the time domain, on $t=0,\Delta t,\dots,(n_t-1)\Delta t$.
For a spectrum, apply a window and a Fourier transform yourself.

**Exact diagonalization.** For small chains we sum over the eigenstates,
$\langle O\rangle_T=\sum_n p_n\langle n|O|n\rangle$ with
$p_n=e^{-E_n/T}/Z$. Here `sc` is a short chain, `n` of about 10 sites:

```python
sc.vev(sc.Sz[0]*sc.Sz[1], mode="ED", T=0.5)
sc.get_correlation_matrix(T=0.5, operators=[sc.Sz[i] for i in range(n)])
sc.get_dynamical_correlator(mode="ED", submode="ED", T=0.5,
                            name=(sc.Sz[0], sc.Sz[0]))
```

The sum uses the full spectrum when the Hilbert space has at most about
2000 states. Above that, only the `n=` lowest states are used, and an
error is raised if the last one still has a large Boltzmann weight. The
ED dynamical correlator is the exact reference for
`metts_dynamical_correlator` (see
`examples/finite_temperature/dynamical_metts_VS_ED`).

## 10. Topological invariants

**Many-body Berry phase.** Let the Hamiltonian depend on a parameter
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

**Single-particle Berry phase.** For a Bloch Hamiltonian $H(k)$, the Zak
phase of the occupied bands is the phase of the Wilson loop,

$$\gamma=\arg\det\prod_k U(k,k+\delta k),\qquad
U_{mn}(k,k+\delta k)=\langle u_m(k)|u_n(k+\delta k)\rangle,$$

with $m,n$ running over the occupied bands. It gives the polarization of
the band. Here `hk(k)` returns the Bloch matrix, with $k$ from 0 to 1:

```python
from dmrgpy import topology
gamma = topology.berry_phase_matrix(hk, nk=40)   # the lower half of the bands is occupied
```

**Fermion parity.** For a fermionic chain the total parity is

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

## 11. Mean-field decoupling

For a spin Hamiltonian with exchange couplings $J_{ij}$, the
self-consistent mean-field (Weiss-field) decoupling replaces the exchange
by a field on each site,

$$H_{\rm MF}=p\,H_{\rm exch}+H_{\rm onsite}+\sum_i\mathbf h_i\cdot\mathbf S_i,
\qquad\mathbf h_i=(1-p)\sum_jJ_{ij}\langle\mathbf S_j\rangle.$$

We solve $H_{\rm MF}$, compute the new $\langle\mathbf S_i\rangle$, rebuild
$\mathbf h_i$, mix the old and new values, and repeat until the
magnetization stops changing. With $p=0$ this is pure mean-field theory.
With $0<p<1$ a fraction $p$ of the exchange is kept and treated exactly,
together with the field. With $p=1$ it is the original model. The
couplings are read from the Hamiltonian of the chain, so any one- and
two-site spin Hamiltonian works, including anisotropic and long-range
exchange. Here `sc` is the Heisenberg chain of section 3:

```python
from dmrgpy import meanfield
n = sc.ns
m0 = [[0, 0, 0.5*(-1)**i] for i in range(n)]       # staggered start
mf = meanfield.spinchain_meanfield(sc, p=0.0, m0=m0, mix=0.5)
m = mf.get_magnetization()       # self-consistent <S_i>, a 3 x n array
e = mf.gs_energy()               # mean-field energy
```

The function returns a new chain whose Hamiltonian is the converged
$H_{\rm MF}$. This Hamiltonian includes the constant
$-(1-p)\sum_{ij}J_{ij}\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$,
so the exchange is not counted twice. At $p=0$ its `gs_energy()` is the
expectation value of the original Hamiltonian in the mean-field state.
The keywords are:

- `p`: the fraction of the exchange kept exactly.
- `m0`: the starting magnetization, one 3-vector per site. Random by
  default. For an antiferromagnet start from a staggered pattern, since
  simple mixing can oscillate otherwise.
- `mix`: the weight of the new magnetization in each step (default 0.9).
  Lower it if the loop oscillates.
- `maxerror`, `maxite`: the convergence tolerance and the maximal number
  of iterations. After `maxite` the loop stops with a warning.
- `mixmode`: `"default"` (linear mixing) or `"broyden"`.
- `mode="ED"`: solve each step with exact diagonalization.

`meanfield.decompose_spin_hamiltonian(sc)` returns the onsite fields
`b[i,a]`, the exchange `J[i,j,a,b]` and the constant read from the
Hamiltonian. A Hamiltonian that is not of this form, for example with a
three-site term, raises an error.

## 12. Fidelity susceptibility and quantum phase transitions

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

## 13. Ground-state degeneracy

A degenerate ground state signals a broken symmetry or a topological
phase. Examples are the two ordered states of an Ising chain and the four
edge states of an open Haldane chain. We count the levels close to $E_0$,

$$g\approx\sum_i\exp\!\left[-\left(\frac{(E_i-E_0)^2}{\delta}\right)^2\right],$$

adding excited states $E_i$ until the count stops growing:

```python
g = sc.get_gs_degeneracy()             # close to an integer
g = sc.get_gs_degeneracy(delta=1e-4)   # a stricter window
```

A level counts when $|E_i-E_0|$ is smaller than about $\sqrt{\delta}$,
which is 0.1 at the default `delta=1e-2`. A result that is not close to
an integer means that a level sits at the edge of the window. On a finite
chain a degeneracy is often split by an amount that decays with the
length, so choose `delta` relative to that splitting. `mode="ED"` uses
exact diagonalization.

## 14. Reduced density matrices and operator distributions

**Reduced density matrix.** The reduced density matrix of one site,

$$\rho_i=\mathrm{Tr}_{j\neq i}|\mathrm{GS}\rangle\langle\mathrm{GS}|,$$

contains every local observable of that site, and its eigenvalues give the
site entanglement entropy (section 5).

```python
rho = sc.get_rdm(i=0)     # d x d matrix, d the local dimension
```

It is computed from the MPS, so it is not available with `mode="ED"`,
where it raises an error.

**Operator distributions.** For a Hermitian operator $X$, the expectation
value $\langle X\rangle$ is only the mean. The full distribution of the
measured values of $X$ in the ground state is

$$P(x)=\langle\mathrm{GS}|\,\delta(X-x)\,|\mathrm{GS}\rangle,\qquad
\int P(x)\,dx=1.$$

For a conserved charge this is the full counting statistics. A sharp
$P(x)$ means that the ground state has a well-defined value of $X$, and a
broad one means large fluctuations. It is computed with the kernel
polynomial method, as the dynamical correlators of section 6:

```python
X = sum(sc.Sz[i] for i in range(n))      # total magnetization
x, p = sc.get_distribution(X=X, delta=0.1)
mus, shift, scale = sc.get_distribution_moments(X=X, delta=0.1)
```

`delta` sets the resolution through the number of Chebyshev moments.
Unlike in section 6, it is not calibrated as a broadening, so compare
distributions computed with the same `delta`. `get_distribution` also
works with `mode="ED"`. `get_distribution_moments` returns the raw
Chebyshev moments of the rescaled operator and has no ED version. A
maximum-entropy reconstruction from the moments $\langle X^k\rangle$,
`distribution.get_distribution_maxent`, needs the module
`dmrgpy.maxenttk`, which is not shipped, and raises an error without it.

## 15. Post-processing tools

- `extrapolate.size_extrapolation(L, y, n=1)`: fits $y(L)$ with a
  polynomial of order `n` in $1/L$ and returns the value at
  $L\to\infty$. Use it for energy densities, gaps and order parameters
  from finite chains.
- `functionfit.fit(f, x0)`: minimizes a function `f` with several random
  starts of the Powell method. Use it for fits such as the
  Calabrese-Cardy formula of section 5.
- `analyticcontinuation.imag2real(zs, us)`: Pade continuation of a
  function known at complex frequencies `zs` to the real axis, for
  example the output of `submode="CVMimag"` (section 6). It needs the
  module `dmrgpy.padetk`, which is not shipped, and raises an error
  without it.
- `reconstruct.reconstruct_distribution`: maximum-entropy reconstruction of a positive spectral
  function from its moments, used by `submode="maxent"` and
  `get_distribution_maxent`. It needs the module `dmrgpy.maxenttk`, which
  is not shipped, and raises an error without it.

## 16. Worked-example cookbook

Each recipe starts from a chain `sc` of `n` sites, built as in section 3
with the Hamiltonian given.

**Central charge of the critical Ising chain.** For
$H=\sum_iS_i^zS_{i+1}^z+\tfrac12\sum_iS_i^x$, at the critical field, the
entanglement entropy of the open chain follows the Calabrese-Cardy formula
with $c=1/2$ (section 5):

```python
sc.maxm = 200          # critical chains need a large bond dimension
wf = sc.get_gs()
c = wf.get_CFT_central_charge()
```

**Haldane gap of the spin-1 chain.** The open chain has four low states
from the two spin-1/2 edge states. The bulk gap is the fifth state
(section 4):

```python
sc = spinchain.Spin_Chain(["S=1" for i in range(n)])
# ... Heisenberg h, set_hamiltonian(h) ...
es = sc.get_excited(n=6)
gap = es[4] - es[0]
```

**Gaps of a Hubbard chain.** Here `fc` is a `Spinful_Fermionic_Chain`
with a Hubbard Hamiltonian. `get_gap()` is the gap to the first excited
state, and `get_charge_gap(d=...)` is $E_0(N+d)-E_0(N)$ (section 4):

```python
fc.get_gap()
fc.get_charge_gap(d=1)
fc.get_charge_gap(d=2)
```

**Fidelity susceptibility across the Ising transition** (section 12):

```python
import numpy as np
from dmrgpy import fidelity
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
lambdas = np.linspace(0.5, 1.5, 11)
chis = [fidelity.get_fidelity(sc, h0, h1, lam) for lam in lambdas]
```

**Dynamical structure factor $S(q,\omega)$** (section 6), from the
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

**Extrapolation to the thermodynamic limit** (section 15), here of the
energy per site computed on several lengths `Ls`:

```python
from dmrgpy import extrapolate
e_inf = extrapolate.size_extrapolation(Ls, e_per_site, n=1)
```

## 17. STM/Kondo tunneling spectra (third-order perturbation theory)

In this section we compute the differential conductance $dI/dV$ measured
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

**Second order** (`order=2`) is the Fermi golden rule for spin-flip and
elastic tunneling,

$$\frac{dI}{dV}(eV)\propto\sum_{i,f}p_i\Big[\tfrac12|\langle f|S_-|i\rangle|^2+\tfrac12|\langle f|S_+|i\rangle|^2+|\langle f|S_z|i\rangle|^2\Big]\,\Theta(eV-\epsilon_{if})+4U^2,$$

summed over both tunneling directions. Here $p_i$ is the Boltzmann weight
of eigenstate $i$, $\epsilon_{if}=E_f-E_i$, and $\Theta$ a step broadened
by temperature. Each excitation of the spins gives a step in $dI/dV$ at
$eV=\pm\epsilon_{if}$. This is inelastic tunneling spectroscopy (IETS),
and the step positions measure the spin excitation energies, for example
the Zeeman energy.

**Third order** (`order=3`, the default) adds two terms. Both sum over
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

**ED and DMRG.** `mode="ED"` (the default) diagonalizes the full
Hamiltonian, since every eigenstate enters as a virtual state. It works
at any temperature, on small chains. `mode="DMRG"`
works at $T=0$ only and never needs excited states. The second-order and
potential terms come from the dynamical correlator of section 6, and the
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

See `examples/kondo/kondo_third_order_dmrg_VS_ED`. The tip couples to one
site only; several sites with independent couplings are not implemented.

### Orbital-resolved IETS of a magnetic atom

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
of section 6 on `dex` and `T` apply here. See
`examples/kondo/atom_iets_orbital_resolved`.

## 18. Infinite chains (iDMRG)

An infinite chain repeats a unit cell of `n_uc` sites forever in both
directions. It describes the thermodynamic limit directly, with no edges
and no finite-size extrapolation. The state is a uniform matrix product
state, and the results are densities per site and correlators at any
distance.

**Building the Hamiltonian.** Each operator exists in three copies per
site `i` of the unit cell: `SzC[i]` in the central cell, `SzR[i]` in the
next cell and `SzL[i]` in the previous cell. A term is written once, and
it is repeated in every cell:

```python
from dmrgpy import infinitechain
ic = infinitechain.Infinite_Spin_Chain(["1/2"])          # n_uc = 1
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h)
e = ic.gs_energy()      # energy per site
```

A two-site cell describes for example a dimerized chain:

```python
ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
h = (1.0*(ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1])
     + 0.4*(ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0]))
ic.set_hamiltonian(h)
```

For longer couplings, `ic.get_operator(name, i, group=c)` gives site `i`
of the cell at distance `c`, so `c=-1,0,1` are `L`, `C`, `R`. For example
a J1-J2 chain on a one-site cell:

```python
h = 0
for op in ("Sx", "Sy", "Sz"):
    o0 = ic.get_operator(op, 0)
    h = h + o0*ic.get_operator(op, 0, group=1)        # J1
    h = h + 0.4*o0*ic.get_operator(op, 0, group=2)    # J2
```

At $J_2=0.4$ this chain is dimerized, so its ground state has a period of
two sites. A one-site cell cannot hold it and the solver does not
converge. Use the same couplings on a two-site cell for this phase.

Fermionic chains use `Infinite_Many_Body_Chain` with a list of site
codes: `0` for a spinless fermion (`C`, `Cdag`, `N`) and `1` for a
spinful site (`Cup`, `Cdagup`, `Cdn`, `Cdagdn`, `Nup`, `Ndn`). The
Jordan-Wigner string is added automatically, also for hoppings that skip
sites. A term with an odd number of fermion operators raises an error.

**Two solvers.** `ic.gs_method` selects the algorithm:

- `"vumps"` (the default): variational uniform MPS (arXiv:1701.07035).
  It finds the best state of bond dimension `maxm` directly in the
  thermodynamic limit. It works for any unit cell and any range of
  coupling.
- `"idmrg"`: infinite DMRG, which grows a finite chain from the center
  until the energy per site converges. Only `n_uc` of 1 or 2.

Both run on `itensor_version="python"` and `3`. `itensor_version=2` and
`"julia_live"` raise an error. The keywords are:

- `ic.maxm`, `ic.cutoff`: the bond dimension (default 30) and the
  truncation cutoff. Increase `maxm` until the results stop changing.
- `ic.maxiter`, `ic.etol`: the maximal number of iterations and the
  energy tolerance.
- `ic.vumps_nrestarts`: the number of random starts of VUMPS (default 4).
- `ic.noise`, `ic.noise_iters`: for iDMRG. A particle-conserving model has
  the empty and the full state as exact eigenstates, and iDMRG can get
  stuck in them. The noise (default $10^{-4}$) breaks this, and is applied
  only while the state is a product state.

After the solve, `ic.converged` says whether the tolerance was reached.
If it was not, `gs_energy()` also gives a `RuntimeWarning` with the
reason. VUMPS starts from a random state, so on a hard model
call `gs_energy()` a few times and keep the lowest energy. A simple
product state, such as the Neel state, gives an upper bound to compare
with. For a gapless
model the convergence of iDMRG is slow, so `converged` can be `False`
while the energy is already accurate to several digits. Entanglement
entropies are not available for infinite chains.

**Static correlators.**

```python
ic.vev("Sz", 0)                    # <Sz> on site 0 of the cell
ic.correlator("Sz", 0, "Sz", r)    # <Sz_0 Sz_r>, r in sites, r >= 0
ic.correlator("Cdag", 0, "C", r)   # fermionic chain, with the string
```

They work with both solvers and both backends. A pair with an odd number
of fermion operators raises an error. A good check of `maxm` needs no
reference value: the sum of the bond correlators of one unit cell must be
equal to `n_uc` times the energy per site.

**Excitations.** In an infinite chain the excitations form bands $E(k)$.
We use the quasiparticle ansatz: one tensor $B$ inserted on top of the
ground state, in a plane wave of momentum $k$. This gives the dispersion
of the lowest modes, as measured by inelastic neutron scattering:

```python
ic.excitation_energies(k, n=1)   # n lowest energies at momentum k
ic.excitation_gap()              # minimum over k
```

$k$ is in radians per unit cell. It needs `gs_method="vumps"`, `n_uc` of
1 or 2, and couplings to the neighboring cell only.

**Dynamical structure factor.** Each quasiparticle branch $a$ contributes
a delta peak to $S(k,\omega)$ with weight
$w_a=|\langle k,a|O_k|\mathrm{GS}\rangle|^2$:

$$S(k,\omega)=\sum_aw_a\,\delta(\omega-E_a(k)).$$

```python
e, w = ic.spectral_weights("Sz", k, n=2)
e, w, total = ic.spectral_weights("Sz", k, n=2, return_total=True)
ks, es, S = ic.dynamical_structure_factor("Sz", delta=0.05)  # map in (k, w)
```

The peaks are exact in momentum and energy. They do not contain the
continuum of two or more particles. `total` is the static structure
factor $\sum_re^{ikr}(\langle O_0O_r\rangle-\langle O\rangle^2)$, so
`w.sum()/total` is the fraction of the weight in the computed branches.
When it is close to 1, a single-mode picture is valid at that $k$. In
`dynamical_structure_factor`, `delta` is a Lorentzian broadening for
plotting. These two methods need `itensor_version="python"` and
`gs_method="vumps"`, and a parity-even operator. See
`examples/idmrg/dynamical_structure_factor_tfim`.

**Local gap estimate.** `ic.local_excitation_gap()` takes the converged
two-site block of iDMRG and returns the gap between its two lowest
states. It is cheap, but less accurate than `excitation_gap()` (10 per
cent too high on the XX chain). On fermionic models it can be a charge
excitation and shifts with a uniform chemical potential. Use it only as a
rough check. It needs `gs_method="idmrg"`. `window=w`, with
`itensor_version="python"` and `n_uc=1`, adds `w` sites on each side of
the block and improves the estimate.

**Dynamical correlators.** Two methods compute
$\langle O_i(t)O_j\rangle$ and its spectrum in an infinite chain:

```python
es, ys = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=16,
             window_chain_kwargs=dict(maxm=30, nsweeps=15),
             delta=0.3, es=np.linspace(-1, 5, 200))
ic.gs_method = "idmrg"
ks, es, Skw = ic.td_dynamical_correlator("Sz", 0, "Sz", n_window=14,
             dt=0.05, nt=200, maxdim=60, x_values=range(-6, 7),
             ks=np.linspace(-np.pi, np.pi, 41), delta=0.1)
```

- `kpm_finite`: the KPM correlator of section 6 on a finite open window of
  `n_window` unit cells. It has finite-size and edge effects, so converge
  it in `n_window`. A small `delta` needs many moments, and then the
  edges matter even for a long window. The window always runs on
  `itensor_version="python"`.
- `td_dynamical_correlator`: real-time evolution of a window whose two
  ends are attached to the converged infinite environment
  (arXiv:1804.09163), so it has no edges. One evolution gives every
  distance `x` and every momentum `k`. `n_window` and `x_values` must be
  converged separately. It needs `gs_method="idmrg"`. `connected=True`
  (default) subtracts $\langle O_i\rangle\langle O_j\rangle$. It returns
  the one-sided Fourier transform, so take the real part of `Skw` for the
  spectral function.

See `examples/idmrg/td_dynamical_correlator` for both.

**Lower-level tools.** These act on the converged state of
`itensor_version="python"` (`ic._result` for iDMRG, `ic._vumps_result`
for VUMPS) and have no method on the chain:

- `pyitensor.idmrg.apply_mpo` and `pyitensor.vumps.apply_mpo`: apply a
  bounded periodic operator, for example a product of one-site operators
  or a gate, and truncate. Not for a Hamiltonian, which is extensive.
- `pyitensor.idmrg.imps_overlap`: the overlap per site of two infinite
  states, 1 in magnitude for the same state.
- `pyitensor.idmrg.imps_sum` and `pyitensor.vumps.imps_sum`: the direct
  sum of two states. Two normalized ground states give a cat state, which
  this cannot represent, and raise an error.

See `examples/idmrg/vumps_apply_mpo`, `imps_overlap` and `vumps_imps_sum`.

## 19. Running the pure-Python backend on a GPU

The `itensor_version="python"` backend can keep its tensors on a GPU. It
needs `jax` with CUDA support. The switch is global and is set before
building a chain:

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")       # "numpy", the default, goes back to the CPU
print(backend.device_info())     # for example "jax: cuda:0"
```

The rest of the script does not change, and the results agree with the
CPU ones to about $10^{-11}$. A GPU helps only when the bond dimension is
large: below about 120 it is slower than the CPU, and at 240 it is about
5 times faster. So it pays for real-time evolution, where the bond
dimension grows, and rarely for a 1D ground state. Infinite chains and
non-Hermitian DMRG run but do not get faster. For a script that does one
calculation and exits, `backend.set_pad_bonds(maxm)` fixes all tensor
shapes so that each kernel is compiled only once. With padding on, the
default `backend.set_jit("auto")` also fuses the kernels. For a KPM
correlator,
set `kpmmaxm` equal to `maxm`. See `docs/gpu_cpu_performance.md` for the
measurements.

## 20. Performance: BLAS threads

DMRG makes many small matrix operations, and on small matrices several
BLAS threads are slower than one. It matters most when the machine is
shared, for example on a cluster node or with several runs in parallel,
where a calculation can become 3 to 17 times slower. Pin one thread
before Python imports numpy:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 myscript.py
```

Inside a script, `blasthreads.limit(1)` does the same for one block. It
needs the package `threadpoolctl`, and without it does nothing and gives a
warning:

```python
from dmrgpy import blasthreads
with blasthreads.limit(1):
    e0 = sc.gs_energy()
```

dmrgpy never changes the number of threads by itself.
