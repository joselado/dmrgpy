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
model on a one-site cell, iDMRG can end in a dimerized state (see
[Static correlators](#static-correlators)). Then the energy per site
alternates between iterations, so `converged` stays `False` while the
energy is already accurate to several digits. Entanglement
entropies are not available for infinite chains.

![Energy and spin correlator of the infinite Heisenberg chain](img/infinite-chains_energy.png#only-light)
![Energy and spin correlator of the infinite Heisenberg chain](img/infinite-chains_energy-dark.png#only-dark)

The left panel shows the error of the iDMRG energy per site of the
spin-1/2 Heisenberg chain, against the Bethe-ansatz value
$e_\infty=1/4-\ln 2$, as a function of `maxm`. The chain is gapless, so
the error falls as a power of `maxm`, from $3\times10^{-4}$ at `maxm=8`
to $1.6\times10^{-6}$ at `maxm=64`. On the one-site cell it stops there.
A two-site cell gives $5\times10^{-7}$ at `maxm=96` (square). Here
`converged` stays `False` up to `maxm=64`, because the one-site state is
dimerized, while the energy is already accurate. The right panel is described in
[Static correlators](#static-correlators).

??? example "Code of this figure (docs/figures/infinite-chains_energy/)"

    ```python title="compute.py"
    --8<-- "infinite-chains_energy/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "infinite-chains_energy/plot.py"
    ```

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

With iDMRG at finite `maxm`, a gapless chain can end in a dimerized
state: the bonds alternate in strength. This comes from the finite
`maxm`, and the alternation becomes smaller when `maxm` grows. A one-site
cell cannot hold this state. The energy per site then alternates between
iterations, so `converged` stays `False`, and
`ic.correlator("Sz", 0, "Sz", 1)` gives only one of the two bonds. In
this case `gs_energy`, `vev` and `correlator` give a `RuntimeWarning`
with the entanglement entropies of the two bonds. On a two-site cell both
bonds are measured, and their mean is the energy per site.

When VUMPS did not converge, `vev`, `correlator` and the excitations also
give a `RuntimeWarning`. Such a state is not a consistent uniform MPS, and
its correlators can be wrong by 10 per cent while its energy is accurate.

The right panel of the figure in [Two solvers](#two-solvers) shows
$|\langle S^z_0S^z_r\rangle|$ of the Heisenberg chain on a two-site
cell at `maxm=96`, averaged over the two sites. There is no edge, so the
correlator decays as the power law $a\sqrt{\ln r}/r$ of a Luttinger
liquid (dashed, $a$ fitted at $r=8$ to 20). At $r=1$ it gives
$3\langle S^z_0S^z_1\rangle=-0.443147$, which is $e_\infty$. At the
largest $r$ the data fall a little below the power law, because a finite
`maxm` gives a finite correlation length.

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

![Excitation bands and structure factors of infinite chains](img/infinite-chains_excitations.png#only-light)
![Excitation bands and structure factors of infinite chains](img/infinite-chains_excitations-dark.png#only-dark)

The left panel is the transverse-field Ising chain
$H=-\sum\sigma^x\sigma^x-g\sum\sigma^z$ at $g=2.5$, in the paramagnetic
phase, with `maxm=4`. The map is $S^{xx}(k,\omega)$. The VUMPS band
(dashed) lies on the exact free-fermion band
$2\sqrt{1+g^2-2g\cos k}$ (dotted) within $10^{-8}$. The right panel is
the spin-1 Heisenberg chain with `maxm=16`, and the map is
$S^{zz}(k,\omega)$ on a log color scale. The magnon branch has its
minimum at $k=\pi$. There it gives the Haldane gap $0.4098$, against the
reference $0.4105$ (dotted line). Near $k=0$ the weight goes to zero,
because $\sum_jS^z_j$ is conserved. The weak features above the branch
come from higher states of the ansatz. The ansatz does not describe the
two-magnon continuum, which starts at $2\Delta$.

??? example "Code of this figure (docs/figures/infinite-chains_excitations/)"

    ```python title="compute.py"
    --8<-- "infinite-chains_excitations/compute.py"
    ```

    ```python title="plot.py"
    --8<-- "infinite-chains_excitations/plot.py"
    ```

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

## Operators on the state, overlaps and sums

Three methods act on the converged state itself and return the result as
a new chain, so that `vev` and `correlator` measure it. Here `ic` and
`ic2` are two converged chains on the same unit cell, on
`itensor_version="python"`:

```python
flipped = ic.apply_operator(2*ic.SxC[0])   # Pauli X on every site
flipped.vev("Sz", 0)                       # -<Sz> of ic
flipped.eta                                # 1 for a unitary operator
f = ic.overlap(flipped)                    # fidelity per site, 0 to 1
f2 = ic.overlap(ic2)                       # between two solved chains
try:
    cat = ic.summps(ic2)                   # |ic> + |ic2>
except RuntimeError:
    pass                                   # two normalized states, see below
```

- `apply_operator(op)`: `op` is an operator on the central cell, built
  from `ic.SxC`, `ic.SzC[0]*ic.SzC[1]` and so on, applied in every cell at
  once. It must stay inside the cell (an `SxR` factor raises an error),
  since it is repeated cell by cell. The state is normalized afterwards,
  and `eta` is the factor by which the operator changed its norm per
  cell. A gate that entangles neighboring sites grows the bond
  dimension, and the result is refused when truncating it back loses
  weight; pass `maxdim=` larger than `ic.maxm`. The gauge fixing behind
  this holds only for a state solved at a modest `maxm` (4 worked, 8 did
  not, on a free-fermion cell), so apply such a gate to a small state.
- `overlap(other)`: the fidelity per site of two states on the same unit
  cell, $|\langle\phi|\psi\rangle|^{1/N}$ in the limit of $N$ sites. It is 1 for the
  same physical state in any gauge, and below 1 otherwise. Two ground
  states on the two sides of a phase transition have a fidelity that
  drops at the transition, the infinite-chain version of
  [Fidelity susceptibility](fidelity.md).
- `summps(other)`: the direct sum of two states. Tiled to infinity, the
  sum keeps only the branch of larger norm per site, so two normalized
  ground states tie and raise an error: a cat state of two ordered
  states has no uniform MPS representation.

The returned chain is not a ground state, so `gs_energy()` and the
excitations raise on it; `set_hamiltonian()` makes it an ordinary chain
again. See `examples/idmrg/apply_mpo_to_infinite_chain` and
`examples/idmrg/imps_overlap`.
