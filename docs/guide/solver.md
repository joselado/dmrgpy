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

## Starting from a given state

By default a solve starts from a state the backend chooses, a random one
on `itensor_version=3` and `"python"`. `wf0=` starts the sweeps from a
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
    `cutoff`, `noise` and the warm starts work on every backend. The
    bond ramp runs on `itensor_version=2`, `3` and `"python"`; on
    `"julia_live"` every sweep runs at `maxm`, from a random start of
    bond dimension `bond_ramp_start`. `verbose` prints on
    `itensor_version=2`, `3` and `"python"`.
