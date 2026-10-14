# Infinite chains

An infinite chain repeats a unit cell of `n_uc` sites forever in both
directions. It describes the thermodynamic limit directly, with no edges
and no finite-size extrapolation. The state is a uniform matrix product
state, and the results are densities per site and correlators at any
distance.

## Building the Hamiltonian

Each operator exists in three copies per
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

## Two solvers

`ic.gs_method` selects the algorithm:

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

## Static correlators


```python
ic.vev("Sz", 0)                    # <Sz> on site 0 of the cell
ic.correlator("Sz", 0, "Sz", r)    # <Sz_0 Sz_r>, r in sites, r >= 0
ic.correlator("Cdag", 0, "C", r)   # fermionic chain, with the string
```

They work with both solvers and both backends. A pair with an odd number
of fermion operators raises an error. A good check of `maxm` needs no
reference value: the sum of the bond correlators of one unit cell must be
equal to `n_uc` times the energy per site.

## Excitations

In an infinite chain the excitations form bands $E(k)$.
We use the quasiparticle ansatz: one tensor $B$ inserted on top of the
ground state, in a plane wave of momentum $k$. This gives the dispersion
of the lowest modes, as measured by inelastic neutron scattering:

```python
ic.excitation_energies(k, n=1)   # n lowest energies at momentum k
ic.excitation_gap()              # minimum over k
```

$k$ is in radians per unit cell. It needs `gs_method="vumps"`, `n_uc` of
1 or 2, and couplings to the neighboring cell only.

## Dynamical structure factor

Each quasiparticle branch $a$ contributes
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

## Local gap estimate

`ic.local_excitation_gap()` takes the converged
two-site block of iDMRG and returns the gap between its two lowest
states. It is cheap, but less accurate than `excitation_gap()` (10 per
cent too high on the XX chain). On fermionic models it can be a charge
excitation and shifts with a uniform chemical potential. Use it only as a
rough check. It needs `gs_method="idmrg"`. `window=w`, with
`itensor_version="python"` and `n_uc=1`, adds `w` sites on each side of
the block and improves the estimate.

## Dynamical correlators

Two methods compute
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

- `kpm_finite`: the KPM [dynamical correlator](dynamical-correlators.md) on a finite open window of
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

## Lower-level tools

These act on the converged state of
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
