# Controlling the solver

On this page we control the DMRG solve itself: the sweep schedule, the
state it starts from, and the state the chain keeps. The physics of the
ground state is on the [Ground state](ground-state.md) page; here `sc` is
its Heisenberg chain.

## Sweep schedule

Each sweep moves through the chain and updates two sites at a time,
keeping at most `sc.maxm` Schmidt states on every bond. The settings are:

- `sc.maxm` (default 30): the bond dimension. The one that sets the
  accuracy; increase it until the result stops changing.
- `sc.nsweeps` (default 15): the number of sweeps. A gapped chain
  converges in a few sweeps, a critical or frustrated one needs more.
- `sc.cutoff` (default $10^{-12}$): Schmidt weights below this are
  discarded, so a weakly entangled bond uses fewer than `maxm` states.
- `sc.noise` (default $10^{-7}$): a small density-matrix perturbation
  during the first half of the sweeps, which helps the state leave a
  local minimum, for example a product state that is an exact
  eigenstate. It is zero in the final sweeps, so it does not bias the
  result. Raise it to $10^{-5}$ when a solve is stuck, set it to 0 for a
  model that converges without it.
- `sc.bond_ramp` (default `True`): the first half of the sweeps grow the
  bond dimension geometrically from `sc.bond_ramp_start` (default 10) to
  `maxm`, and the second half run at `maxm`. The early sweeps mostly
  locate the right state and gain little from a large bond dimension, so
  this makes the solve cheaper without changing the result. Set it to
  `False` to run every sweep at `maxm`; `sc.bond_ramp_fraction` (default
  0.5) is the fraction of the sweeps spent growing.
- `sc.verbose = True`: prints the progress of each sweep, to see whether
  `nsweeps` is enough.

These settings are also constructor keywords,
`Spin_Chain(sites, maxm=60, nsweeps=20)`. The chain remembers the solved
state: `gs_energy()` and `get_gs()` return it without a new solve as long
as the Hamiltonian and the settings above do not change, and
`sc.restart()` forgets it.

## One-site updates with bond expansion

The default sweep updates two sites at a time. This lets a bond grow, but
the local problem has dimension $d^2D^2$, so it becomes expensive for a
large local dimension $d$, for example $d=4$ in a Hubbard chain. A one-site
update costs only $dD^2$, but alone it cannot grow a bond: from a product
state it stays a product state. Controlled bond expansion (Gleis, Li and
von Delft, [PRL 130, 246402 (2023)](https://arxiv.org/abs/2207.14712))
solves this. Before each one-site update, the bond gets the few new states
on which $H|\psi\rangle$ has weight outside the current two-site space,
the part that a two-site update would add. After the update the bond is
cut back to `sc.maxm`. So it converges like two-site DMRG with one-site
local problems.

```python
sc.dmrg_method = "CBE"
e0 = sc.gs_energy()
```

- `sc.dmrg_method` (default `"TWO_SITE"`): `"CBE"` for one-site updates
  with controlled bond expansion. Use it for chains with a large local
  dimension; for spin-1/2 chains the two-site update is usually faster.
- `sc.cbe_delta` (default 0.1): each bond is expanded to $(1+\delta)$`maxm`
  states before the update, as in the paper.

On `itensor_version="python"` it also works with a conserved sector
([Fixed particle number or fixed magnetization](ground-state.md#fixed-particle-number-or-fixed-magnetization)).

!!! note "Backends"
    `"CBE"` exists on `itensor_version="python"`, on
    `itensor_version=3` and on `"julia_live"`; on the last two it does not
    work in a conserved sector. The
    other backends raise an error, and so do the excited-state search and
    the non-Hermitian ground state.

See `examples/groundstate/cbe_ground_state`.

## Starting from a given state

By default a solve starts from a state the backend chooses, a random one
on `itensor_version=3`, `"python"` and `"julia_live"`. `wf0=` starts the
sweeps from a
state we give, for example the ground state of the previous
point of a parameter scan:

```python
import numpy as np
bs = np.linspace(0.0, 0.5, 11)
hst = sum((-1)**i*sc.Sz[i] for i in range(n))     # staggered field
wf = None
es = []
for b in bs:
    sc.set_hamiltonian(h + b*hst)
    es.append(sc.gs_energy(wf0=wf, reconverge=True))
    wf = sc.get_gs()                                # start of the next point
```

`sc.set_initial_wf_guess(wf)` does the same for the next solve. A
warm start follows one branch of states along the scan, which is what we
want when a symmetry-broken state has to be tracked, or when the solve
from a random start lands in a wrong state. It does not make the solve
faster in the cases we tried: with the bond ramp, a solve from a random
start converges in as few sweeps as one from the previous point. And it
can be wrong across a level crossing of a conserved quantity: in a scan
over a uniform field $B\sum_iS^z_i$ the ground state jumps from one
$S^z_\mathrm{tot}$ sector to the next, and a state in the old sector
leaves it only through the noise term, slowly. With a few sweeps the
warm start stays in the old sector with an energy that is too high (on
16 sites, by 0.1 to 0.2 with 2 or 3 sweeps, and still by 0.03 after 15
sweeps on `"python"`). Start from a random state there.

A state can also be taken as it is, without sweeping:

```python
sc.set_hamiltonian(sum((-1)**i*sc.Sz[i] for i in range(n)))
neel = sc.get_gs()                 # Neel product state, as the ground state of a staggered field
sc.set_hamiltonian(h)
sc.set_initial_wf(neel)            # or sc.set_gs(neel)
e_neel = sc.gs_energy()            # <Neel|H|Neel> = -(n-1)/4, no sweep
mz = [sc.vev(sc.Sz[i]).real for i in range(n)]   # measured in the Neel state
```

After `set_initial_wf` or `set_gs`, every method of the chain measures
the state we set, with `gs_energy()` returning its energy
$\langle\psi|H|\psi\rangle$. This gives the variational energy of a trial
state, and the correlators or the time evolution of a prepared state.
The state is normalized when it is set. `sc.set_initial_wf(None)`, or a
new `set_hamiltonian`, returns to the solved ground state.

Two more states are useful here. `sc.random_mps()` is a random
normalized MPS, a start that tests whether a solve depends on where it
begins. `sc.get_gs(best=True, n=3)` solves `n` times from independent
random starts and keeps the lowest energy, for a frustrated or disordered
chain where a single solve can end in an excited state.

## Free-fermion states

For a fermion chain, a natural start is the ground state without the
interaction. This is a Slater determinant: the $N_f$ lowest orbitals of
the single-particle Hamiltonian $H_0=\sum_{ij}h_{ij}c^\dagger_ic_j$ are
filled. A Slater determinant is fixed by its correlation matrix
$\Lambda_{ij}=\langle c^\dagger_ic_j\rangle$, and `get_gaussian_state`
writes it as an MPS. Here we start the Hubbard chain at half filling
from its Fermi sea:

```python
import numpy as np
from dmrgpy import fermionchain
n = 8
fc = fermionchain.Spinful_Fermionic_Chain(n)
fc.maxm = 60                                # bond dimension
h = np.zeros((n,n))
for i in range(n-1): h[i,i+1] = h[i+1,i] = -1.0          # hopping, the same for both spins
H = sum(h[i,j]*(fc.Cdagup[i]*fc.Cup[j] + fc.Cdagdn[i]*fc.Cdn[j])
        for i in range(n) for j in range(n) if h[i,j] != 0)
H = H + 4.0*sum((fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) for i in range(n))   # U=4
fc.set_hamiltonian(H)
wf = fc.get_gaussian_state(h=h, nf=n)       # the U=0 Fermi sea, n electrons
fc.set_initial_wf_guess(wf)
e0 = fc.gs_energy()
```

The start does not fix the particle number: when another filling has a
lower energy, DMRG leaves the start and goes there. Here the
particle-hole symmetric form of $U$ keeps the ground state at half
filling; otherwise we fix the filling with a chemical potential or a
penalty $\lambda(N-N_f)^2$ in $H$.

The correlation matrix of the MPS, `fc.vev(fc.Cdag[i]*fc.C[j], wf=wf)`,
is $\Lambda_{ij}$ of the filled orbitals, up to an error of roughly the
square root of `eigval_cutoff`: on a 40-site chain at half filling it is
2.6e-5 at the default 1e-8 and 1.4e-6 at 1e-10. The truncation at `maxm`
adds to this. The check works at sizes far beyond exact
diagonalization, so a free-fermion state is also a test of the MPS
machinery on a large chain. With `gutzwiller=True` every doubly
occupied orbital is projected out, $\prod_i(1-n_{i\uparrow}n_{i\downarrow})$,
which gives the Gutzwiller wave function of strong coupling:
`fc.vev(fc.Nup[i]*fc.Ndn[i], wf=wf)` is then zero on every site.

A superconductor without interactions has a Gaussian ground state too,
with pairing. With a pairing matrix $\Delta_{ij}=-\Delta_{ji}$ the state
is the ground state of the Bogoliubov-de Gennes (BdG) Hamiltonian
$H_0=\sum_{ij}h_{ij}c^\dagger_ic_j+\frac12\sum_{ij}(\Delta_{ij}c^\dagger_ic^\dagger_j+\mathrm{h.c.})$,
where all the quasiparticle states of negative energy are filled. The
particle number is not fixed, so the chemical potential goes in the
diagonal of `h` and there is no `nf`:

```python
fc = fermionchain.Fermionic_Chain(n)
mu = 2.5
h = -mu*np.eye(n)
D = np.zeros((n,n))
for i in range(n-1):
    h[i,i+1] = h[i+1,i] = -1.0
    D[i,i+1], D[i+1,i] = 1.0, -1.0       # p-wave pairing of the Kitaev chain
wf = fc.get_gaussian_state(h=h, delta=D)
```

The anomalous correlations `fc.vev(fc.C[i]*fc.C[j], wf=wf)` are then
those of the BdG matrix. On a spinful chain an $n\times n$ `delta` is the
singlet pairing $\sum_{ij}\Delta_{ij}c^\dagger_{i\uparrow}c^\dagger_{j\downarrow}$,
so `delta=D0*np.eye(n)` is an s-wave superconductor, and
`gutzwiller=True` works with pairing too. In the topological phase of
the Kitaev chain, with open ends, the two Majorana end modes give two
ground states of almost the same energy; when they are degenerate a
warning says so, and the state is one of the two, of definite parity. The keywords are:

- `h=`, `nf=`: the hopping matrix and the number of fermions. On a
  spinful chain an $n\times n$ matrix is used for both spins; a
  $2n\times 2n$ matrix has mode $2i$ for spin up and $2i+1$ for spin down
  of site $i$.
- `delta=`: the pairing matrix, with `h` and without `nf`.
- `cm=`: the correlation matrix $\Lambda$ instead of `h` and `nf`, or
  the $2n\times 2n$ one $\langle\gamma^\dagger_a\gamma_b\rangle$ with
  $\gamma=(c_0,c^\dagger_0,c_1,c^\dagger_1,\ldots)$ for a state with
  pairing. Its eigenvalues must be 0 or 1, meaning that the state is pure.
- `gutzwiller=True`: the Gutzwiller projection, on a spinful chain.
- `eigval_cutoff` (default 1e-8): the accuracy of the circuit of
  rotations that builds the state. A smaller value gives more gates and
  a more accurate state; 1e-8 is enough for a DMRG start.
- `maxblocksize`: a smaller value gives a smaller bond dimension and a
  less accurate state.
- `maxm=`, `cutoff=`: the truncation (default: the chain's). A warning
  says when it changes the state.
- `mode="ED"`: the same state as an exact vector.

The examples `examples/fermion_models/gaussian_state` and
`examples/fermion_models/gaussian_state_pairing` compare the correlations
of a 40-site Fermi sea and of a 40-site Kitaev chain with the exact ones.

## Several chains

`sc.copy()` returns an independent chain with the same Hamiltonian and
settings, to compare two values of `maxm` on the same model, or to keep a
converged chain while another one is modified. Two chains never share a
state.

Two methods solve an auxiliary problem on a copy of the chain and leave
`sc` untouched:

```python
e_min = sc.lowest_eigenvalue(sum(sc.Sz))   # smallest eigenvalue of any Hermitian operator, here -n/2
W = sc.bandwidth(h)                        # E_max - E_min of H
```

The bandwidth is the energy scale the kernel polynomial method rescales
by ([Choosing a correlator method](correlator-methods.md)).

!!! note "Backends"
    `cutoff`, `noise`, the warm starts and the bond ramp work on every
    backend. `verbose` prints on `itensor_version=2`, `3` and `"python"`.
