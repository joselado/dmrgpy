# DMRGPY Physics User Guide

This guide explains what DMRGPY computes from a **physics** point of
view — the quantities each method returns, the formula behind them, and
when to reach for which method. It intentionally says nothing about the
solver backends (DMRG vs ED, C++ vs Python vs Julia) or code layout —
see `docs/documentation.md` for that. Every `Many_Body_Chain` subclass
(`Spin_Chain`, `Fermionic_Chain`, `Bosonic_Chain`, ...) exposes the same
physics-facing API described here, `mode="DMRG"|"ED"` selectable on
almost every call so a result can be cross-checked against an exact
reference on small systems.

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
19. [Performance: BLAS threads](#19-performance-blas-threads)

## 1. Physical models and Hilbert spaces

Every model is a chain of $n$ local Hilbert spaces $\mathcal H=\bigotimes_{i=1}^n \mathcal H_i$, on which a Hamiltonian and observables are built out of local operators.

| Chain class | Local Hilbert space | Key operators |
|---|---|---|
| `spinchain.Spin_Chain` | spin-$S$, $S\in\{\tfrac12,1,\tfrac32,2,\tfrac52,3\}$ per site | $S^x_i,S^y_i,S^z_i$ |
| `fermionchain.Fermionic_Chain` | spinless fermion (occupied/empty) | $c_i,c_i^\dagger,n_i=c_i^\dagger c_i$, Jordan-Wigner string $F_i$ |
| `fermionchain.Majorana_Chain` | Majorana fermion | Majorana operators built from `Fermionic_Chain` |
| `fermionchain.Spinful_Fermionic_Chain` | spin-$\tfrac12$ fermion (4 states: $0,\uparrow,\downarrow,\uparrow\downarrow$), built from two interleaved spinless sites per physical site | $c_{i\sigma},c^\dagger_{i\sigma},n_{i\sigma}$, plus derived $S^x_i,S^y_i,S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$, onsite pairing $\Delta_i=\tfrac12 c_{i\uparrow}c_{i\downarrow}$ |
| `fermionchain.Spinful_Fermionic_Chain_Native` | same physics as `Spinful_Fermionic_Chain`, but on a genuinely 4-dimensional local space (one tensor-network site per physical site, `itensor_version=3` only) | identical operator lists/formulas as `Spinful_Fermionic_Chain` |
| `bosonchain.Bosonic_Chain` | truncated boson Fock space, $n_i\in\{0,\ldots,\text{maxnb}_i-1\}$, per-site dimension `maxnb` (default 4, i.e. up to 3 bosons/site) settable via `Bosonic_Chain(n, maxnb=[...])` | $a_i,a_i^\dagger,n_i$, occupation projectors $\hat n_i^{(k)}=\lvert k\rangle\langle k\rvert$ for $k=0,\ldots,\text{maxnb}_i-1$ (`bc.D[i][k]`, plus `bc.D0`..`bc.D3` when every site has `maxnb`$\,\ge 4$) |
| `parafermionchain.Parafermionic_Chain` | $\mathbb Z_N$ parafermion (clock model), $N\in\{2,3,4\}$ | clock/shift operators $\sigma_i,\tau_i$ and composite parafermion operators $\chi_i,\psi_i$ built as $\tau$-string $\times\sigma_i$ |
| `mixedchain.Mixed_Spin_Fermion_Chain` | mixes genuine spin-$S$ sites and spinful-fermion locations *in the same chain*, one entry per logical location | at a spin location: native $S^x_i,S^y_i,S^z_i$; at a fermion location: $c_{i\sigma},c^\dagger_{i\sigma},n_{i\sigma}$ plus derived $S^x_i,S^y_i,S^z_i,\Delta_i$ as in `Spinful_Fermionic_Chain` |

Spinful fermionic chains are built by *interleaving* two spinless
fermionic sites per physical site (site $2i$ = spin up, site $2i+1$ =
spin down) rather than by a genuinely 4-dimensional local space, so that
the same Jordan-Wigner machinery used for spinless fermions applies
unchanged; `Spinful_Fermionic_Chain` wraps this bookkeeping for you.

`Spinful_Fermionic_Chain_Native` is the alternative built directly on a
genuinely 4-dimensional local space (ITensor v3's own `Electron`/
`Hubbard` site type) instead: one tensor-network site per physical site,
with the same operator lists (`Cup`/`Cdagup`/`Cdn`/`Cdagdn`/`Nup`/`Ndn`/
`Ntot`/`Sx`/`Sy`/`Sz`/`Delta`) and identical physics/sign convention as
`Spinful_Fermionic_Chain` -- the two classes are drop-in equivalent for
any given Hamiltonian, cross-checked to agree exactly under ED and to
DMRG tolerance under `itensor_version=3` (the only DMRG backend this
class wires up). Despite halving the site count, it is *not* generally
faster in practice: see its class docstring
(`fermionchain.py`) for a measured comparison against
`Spinful_Fermionic_Chain` -- two-site DMRG's per-sweep cost is driven by
the local dimension of the two-site block being diagonalized, which
grows faster (dimension $4\times4=16$ per pair of native sites, versus
$2\times2=4$ per pair of interleaved sites) than the site-count halving
saves. The same disadvantage held up under every other regime checked
too: strong on-site coupling, long-range/power-law hopping, two-site
and one-site+subspace-expansion (`"TDVP_GSE"`) real-time evolution, and
the KPM dynamical correlator itself -- see the class docstring for the
full rundown. One case does flip in its favor, though: the 4-point
correlator tensor `<Cdag_i C_j Cdag_k C_l>`
(`mps.MPS.get_four_correlation_tensor()`, §5) is a Python loop of
independent *static* overlaps rather than an iterative two-site search,
so it does not pay the two-site combined-local-dimension penalty above.
Both classes support a `ctmode="full"` C++-accelerated path in addition
to the generic `ctmode="explicit"` one -- `Spinful_Fermionic_Chain_Native`
gets its own (`Chain::four_correlation_tensor_spinful()`, using ITensor's
own automatic Jordan-Wigner insertion on the flavor-resolved operator
names, since ITensor's `ElectronSite` has no bare `"Cdag"`/`"C"` the
plain version needs). Measured (n=3,4,5,6,12 orbitals),
`Spinful_Fermionic_Chain_Native`'s `ctmode="full"` is the fastest of all
four combinations at every size tried, including n=12 (24 flat modes:
~620s vs ~890s for `Spinful_Fermionic_Chain`'s own `ctmode="full"`, a
~30% win). Prefer `ctmode="full"` for this class whenever it's
available (it always is, for `itensor_version=3`). Otherwise prefer
`Spinful_Fermionic_Chain`; no other calculation tried so far makes the
native-site class faster.

`Mixed_Spin_Fermion_Chain` is for models that need a literal local
moment next to a conduction-electron site (e.g. Kondo-lattice-like
Hamiltonians), rather than the large-$U$ two-fermion-site trick
`Spinful_Fermionic_Chain`/`spinfermionchain.py` use to emulate a spin.
Its `sitesin` constructor argument is a list with one entry per logical
location, each either a spin label (`"1/2"`, `"1"`, ... as in
`Spin_Chain`) or a fermion marker (`"F"`); a fermion location expands
internally to a spin-up/spin-down site pair, exactly like
`Spinful_Fermionic_Chain`. All operator lists (`Sx`/`Sy`/`Sz`,
`Cup`/`Cdagup`/`Cdn`/`Cdagdn`/`Nup`/`Ndn`/`Ntot`/`Delta`) are indexed by
*logical* location, not by physical site — the fermion-only operators
read as the literal integer `0` at spin locations, since they have no
meaning there. Currently only `itensor_version=3` (and `"python"`) are
supported; see `mixedchain.py`'s module docstring for why
`itensor_version=2` isn't yet, and `examples/fermion_models/mixed_spin_fermion_chain`
for a worked Kondo-lattice example.

`Bosonic_Chain(n, maxnb=[...])` takes a per-site local dimension list
(defaulting to `[4]*n`, i.e. up to 3 bosons/site); ED always honors it
exactly (`pyboson/boson.py`). On the DMRG side, `itensor_version=3` and
`itensor_version="python"` both thread a non-default `maxnb` through to
the tensor-network site itself (encoded as the site type code
$100+\text{maxnb}_i$, see `mpscpp3/get_sites.h`/`extra/bosonfour.h` and,
on the pure-Python side, `pyitensor/sites/boson.py`'s `get_boson_site()`
factory) — `itensor_version=2` and the Julia backend still only
understand the single fixed 4-level boson site regardless of what
`maxnb` requests, so a non-default `maxnb` should be run under
`itensor_version=3` (the default) or `"python"` for DMRG/ED results to
actually agree; see `examples/boson_models/boson_maxnb_v3_VS_ED`.

## 2. Building a Hamiltonian and observables

Any Hamiltonian or observable is written exactly as its second-quantized
or spin-operator expression, using operator lists indexed by site, e.g.
for the spin-$\tfrac12$ Heisenberg chain

$$H=\sum_{i=1}^{n-1} J\,\mathbf S_i\cdot\mathbf S_{i+1}=\sum_{i=1}^{n-1} J\left(S_i^xS_{i+1}^x+S_i^yS_{i+1}^y+S_i^zS_{i+1}^z\right)$$

```python
from dmrgpy import spinchain
spins = ["S=1/2" for i in range(30)]
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(29):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
```

Any sum of products of these operators is a valid `MultiOperator`, so
essentially any local (or long-range) Hamiltonian on the given lattice
can be built term by term: single-ion anisotropy $D(S_i^z)^2$,
Dzyaloshinskii-Moriya terms $\mathbf D\cdot(\mathbf S_i\times\mathbf
S_{i+1})$, biquadratic exchange $(\mathbf S_i\cdot\mathbf S_{i+1})^2$,
hopping $t\sum_i(c_i^\dagger c_{i+1}+\text{h.c.})$, Hubbard interaction
$U\sum_i n_{i\uparrow}n_{i\downarrow}$, and so on — see §16 for concrete
Hamiltonians.

**Algebra on already-built operators.** `sc.toMPO(h)` compiles a
`MultiOperator` into a `StaticOperator` (an already-built matrix product
operator, `itensor_version` 2, 3, and `"python"` only — not
`"julia_live"` yet). Two `StaticOperator`s can then be combined directly
with `+`, `-`, unary `-`, and scalar `*`/`/`, without going back through
the symbolic `MultiOperator` form:

```python
A = sc.toMPO(sc.Sz[0])
B = sc.toMPO(sc.Sx[0]*sc.Sx[1])
C = A + 2*B - A          # still a StaticOperator
```

This is a compressed direct sum at the tensor-network level (like
ITensorMPS.jl's `+(::MPO, ::MPO)`), useful for combining operators that
only exist as already-built `StaticOperator`s (e.g. two independently
constructed products or exponentials); for the common case of combining
Hamiltonians before ever building an MPO, summing the underlying
`MultiOperator`s directly (as `h = h + ...` above) remains preferred.

## 3. Ground-state properties

**Ground-state energy**

$$E_0=\langle\mathrm{GS}|H|\mathrm{GS}\rangle=\min_{|\psi\rangle}\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}$$

```python
e0 = sc.gs_energy()          # E0
wf = sc.get_gs()             # the |GS> wavefunction itself
```

**The DMRG sweep schedule: bond-dimension ramp.** `sc.maxm` is the
*target* bond dimension, not the one every sweep runs at. By default
(`sc.bond_ramp = True`) the ground-state sweep schedule spends the first
`sc.bond_ramp_fraction` of its `sc.nsweeps` sweeps growing the bond
dimension geometrically from `sc.bond_ramp_start` up to `sc.maxm`, and
holds it at `sc.maxm` for the rest:

$$m_i=\Big\lceil m_\mathrm{start}\Big(\frac{m_\mathrm{max}}{m_\mathrm{start}}\Big)^{i/n_r}\Big\rfloor\ (i<n_r),\qquad m_i=m_\mathrm{max}\ (i\ge n_r),\qquad n_r=\lfloor n_\mathrm{sweeps}\,f_\mathrm{ramp}\rfloor$$

Two-site DMRG costs $\mathcal O(m^3)$ per bond, while the early sweeps
mostly just locate the right variational subspace and gain little from a
large $m$ — so running them small makes them nearly free, and the
expensive full-`maxm` sweeps then start from an already-good state. At
the default `bond_ramp_fraction = 0.5` the second half of the schedule —
and hence the returned energy — always runs at the full `sc.maxm`, so the
ramp is a pure scheduling change: the same ground state, less time.

```python
sc.bond_ramp = True             # on by default
sc.bond_ramp_start = 10         # bond dimension of the first sweep
sc.bond_ramp_fraction = 0.5     # fraction of the sweeps spent ramping
sc.bond_ramp_noise_decay = 0.1  # noise decay per ramping sweep
sc.bond_ramp = False            # ...or restore a flat schedule at sc.maxm
```

The noise term (`sc.noise`, White's density-matrix perturbation, which
keeps DMRG from stalling in a local minimum) is tied to the same
schedule: it starts at `sc.noise`, decays by `sc.bond_ramp_noise_decay`
per ramping sweep, and is switched off entirely once the schedule reaches
`sc.maxm`, so the final, converged sweeps are noise-free. With
`bond_ramp = False` the original schedule is used instead — full
`sc.noise` for the first half of the sweeps, none for the second.

A *warm* start is never truncated by the ramp: when `gs_energy()` is
re-entered with a wavefunction already in hand (`set_initial_wf`, or
simply a previous `gs_energy()` call's own solution), the ramp's starting
bond dimension is floored at whatever that state already carries, so a
re-run can only improve the energy. The ramp applies to the ground-state
solve on all three DMRG backends (`itensor_version` 2, 3 and `"python"`);
`"julia_live"` keeps its own schedule. See
`examples/groundstate/bond_dimension_ramp` for a 30-site inhomogeneous
Heisenberg--Hubbard chain timing ramp against flat.

**Targeting a quantum-number sector.** By default DMRG searches the whole
Hilbert space and returns the global ground state, at whatever particle
number or total magnetization that happens to have.
`set_conserved_sector` instead confines the entire calculation --
starting state, every sweep, the returned wavefunction -- to one sector
of a conserved quantity:

$$E_0(q)=\min_{|\psi\rangle\,:\,\hat Q|\psi\rangle=q|\psi\rangle}\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle},\qquad [\hat Q,H]=0$$

```python
fc.set_conserved_sector(Nf=6)      # exactly 6 particles
sc.set_conserved_sector(Sz=0)      # total Sz = 0
hc.set_conserved_sector(Nf=8, Sz=0)  # native Hubbard chain: both at once
fc.set_conserved_sector()          # no arguments -> back to the full space
```

Which quantities are available follows from the chain's site types:

| chain | conserved quantities |
|---|---|
| `Fermionic_Chain`, `Spinful_Fermionic_Chain` (Jordan--Wigner) | `Nf` |
| `Spinful_Fermionic_Chain_Native` (native Hubbard sites) | `Nf`, `Sz`, or both |
| every spin chain (`Spin_Chain`, any $S$) | `Sz` |
| `Bosonic_Chain` | `Nb` |
| parafermion chains | none |

`Sz` is in ITensor's integer $2S^z$ units, so `Sz=0` is $S^z_\mathrm{tot}=0$
and `Sz=2` is $S^z_\mathrm{tot}=1$. Asking for only part of what a site
type offers means exactly what it says: a native Hubbard chain with `Nf`
alone fixes the particle number while leaving $S^z$ free, so spin-flip
terms stay legal.

This is the direct route to an addition spectrum $E_0(N)$, the charge gap
$E_0(N{+}1)+E_0(N{-}1)-2E_0(N)$, or a magnetization curve — each point a
genuine ground-state energy of its own sector, rather than something
extracted by tuning a chemical potential or a field. Its cost relative to the
unconstrained search depends on how big the solve is: the quantum numbers
restore the block sparsity the default mode gives up, but they also add
per-block bookkeeping, and only the former scales. Measured on Heisenberg
chains (BLAS pinned to one thread), sector mode is **0.6x** — i.e. slower
— at $n=20$, `maxm=60`; **2.0x** faster at $n=40$, `maxm=100`; and
**4.2x** faster at $n=60$, `maxm=200`, always for the same energy. Expect
a small penalty on toy systems and a real speedup on the ones that
actually cost something.

Two consequences to be aware of, both reported rather than left to be
discovered. First, every operator built on the chain while a sector is
set must itself conserve the requested quantities — a non-conserving one
raises a `ValueError` naming it, whether it is in the Hamiltonian, in a
`vev`, or a dynamical-correlator vertex:

```python
sc.set_conserved_sector(Sz=0)
sc.vev(sc.Sz[0]*sc.Sz[1])   # fine
sc.vev(sc.Sx[0])            # ValueError: ... changes ... by Sz=2 ...
```

A Hamiltonian written as $S^xS^x+S^yS^y+S^zS^z$ is fine even though no
single term of it conserves $S^z$: that expansion is recognized and its
$S^+S^+$/$S^-S^-$ strings cancel exactly. Second, the sector invalidates
the Hamiltonian and ground state already held by the backend, so the next
`gs_energy()` re-solves from scratch.

Sector targeting requires `itensor_version=3` and DMRG. A sector-mode
chain deliberately raises rather than falling back to ED (or to any other
backend), since a solver without quantum numbers would silently answer
with the *global* ground state instead.

What works under a sector: ground state and excited states, expectation
values and static correlators, entanglement entropies, real-time
evolution and dynamical correlators via the default
`tevol_method="TDVP"` (and `"TDVP_GSE"`) — all subject to the
conserving-operator rule above. What does not, and says so: METTS
(`metts_vev`, `metts_dynamical_correlator`), `tevol_method="TEBD"`, and
the infinite-chain algorithms (iDMRG, VUMPS). Those assemble tensors by
hand in a way that is dense-only, so they refuse rather than produce a
wrong number. See
`examples/groundstate/sector_targeted_groundstate` for the addition
spectrum and charge gap of a $t$--$V$ chain, checked against ED.

**Leaving a sector without losing the state.** The conserving-operator
rule above rules out exactly the quantities one usually wants *from* a
fixed-$N$ ground state: $c_i|\mathrm{GS}\rangle$, a one-body density
matrix built from it, a photoemission weight, a pairing quench. Clearing
the sector with `set_conserved_sector()` does not help on its own — it
also throws the state away, so the next `gs_energy()` re-solves
unconstrained and answers with the global ground state.
`promote_to_dense()` is the other option: it leaves sector mode while
*keeping* the state, converted exactly to its dense equivalent on the
chain's ordinary site indices.

```python
fc.set_conserved_sector(Nf=6)          # sweeps confined to 6 particles
fc.gs_energy()
fc.applyoperator(fc.C[3], fc.wf0)      # ValueError: c changes Nf by -1
fc.promote_to_dense()                  # same state, ordinary indices
wf = fc.applyoperator(fc.C[3], fc.wf0) # now legal
fc.vev(sum(fc.N), wf=wf)               # -> 5
```

Nothing is re-solved, truncated or approximated: a QN-conserving MPS is
the same wavefunction as its dense counterpart, only stored
block-sparsely, and promotion just scatters the blocks back into full
tensors. What it costs is the block sparsity itself, so promote *after*
the expensive sweeps, not before. The Hamiltonian and the band-edge
caches are rebuilt on the next call that needs them, while the
ground-state energy and wavefunction are kept — a bare `gs_energy()`
afterwards therefore still returns the sector's energy rather than
re-solving; call `restart()` if an unconstrained re-solve is what you
want. A wavefunction Python is already holding is not reached by
`promote_to_dense()` and needs `wf = fc.promote_mps(wf)` of its own.

Promotion always rebases onto the chain's *own* site indices, kept from
construction, so states promoted out of different sectors at different
times remain comparable:

```python
fc.set_conserved_sector(Nf=n//2)   ; fc.gs_energy() ; fc.promote_to_dense()
wf_N = fc.wf0.copy()
fc.set_conserved_sector(Nf=n//2-1) ; fc.gs_energy() ; fc.promote_to_dense()
wf_Nm = fc.wf0.copy()
Z = abs(fc.overlap(wf_Nm, fc.applyoperator(fc.C[i], wf_N)))**2  # photoemission weight
```

`promote_to_dense()` is `itensor_version=3` only, like
`set_conserved_sector` itself, and does nothing if no sector is set. See
`examples/groundstate/sector_promotion_to_dense` for the $t$--$V$ chain's
CDW profile, one-body density matrix and site-resolved photoemission
weight $Z_i=|\langle \mathrm{GS}_{N-1}|c_i|\mathrm{GS}_N\rangle|^2$, all
computed this way and checked against sector-restricted ED.


**Expectation values and moments.** For any operator $O$ built the same
way as $H$,

$$\langle O\rangle=\langle\mathrm{GS}|O|\mathrm{GS}\rangle,\qquad \langle O^n\rangle=\langle\mathrm{GS}|O^n|\mathrm{GS}\rangle$$

```python
mz = [sc.vev(sc.Sz[i]).real for i in range(n)]        # local magnetization profile
e2 = sc.vev(h, npow=2)                                  # <H^2>, e.g. for fluctuations
```

**Energy fluctuation** (a measure of how sharply the DMRG/ED state is an eigenstate, and physically the variance of $H$ in the prepared state):

$$\delta E=\sqrt{\langle H^2\rangle-\langle H\rangle^2}$$

```python
de = sc.gs_energy_fluctuation()
```

**Static two-point correlators.** `sc.vev(sc.Sz[0]*sc.Sz[i])` gives
$\langle S^z_0 S^z_i\rangle$ directly; `correlator.get_correlator`
provides shorthand names for common correlators over a list of site
pairs, e.g. `"SS"` for the full dot product

$$\langle\mathbf S_i\cdot\mathbf S_j\rangle=\langle S_i^xS_j^x\rangle+\langle S_i^yS_j^y\rangle+\langle S_i^zS_j^z\rangle$$

or fermionic correlators like `"cdc"` ($\langle c_i^\dagger c_j\rangle$), `"density"`/`"densitydensity"` ($\langle n_in_j\rangle$), and pairing correlators `"delta"`/`"deltadeltad"`. These static correlators are the equal-time limit of the dynamical correlators in §6, and are what you Fourier-transform to get an equal-time structure factor $S(q)=\sum_{i,j}e^{iq(i-j)}\langle O_iO_j\rangle$.

## 4. Excited states and energy gaps

**Low-lying spectrum.**

```python
es = sc.get_excited(n=4)          # [E0, E1, E2, E3]
gap = sc.get_gap()                # E1 - E0
```

Both `mode="DMRG"` and `mode="ED"` return one entry per eigen*state*, so a
degenerate level appears once per member in both and the two agree index by
index. `get_gap()` is therefore ~0, correctly, whenever the ground level is
degenerate — ask for enough levels that the multiplet and the state you want
both fit (with a three-fold degenerate ground state, the first genuine
excitation is `n=4`). Degenerate members converge to slightly different
energies under DMRG, so when that spread approaches the splitting you are
after, identify states by a quantum number (e.g. $\langle S^2\rangle$)
rather than by position in the list.

`get_excited_states(n, purify=True)` additionally returns the
wavefunctions; `purify=True` re-diagonalizes $H$ in the Gram-Schmidt
orthogonalized subspace spanned by the raw excited-state MPS, correcting
for near-degenerate states that DMRG's iterative solvers can otherwise
mix.

Excited states beyond the ground state are found with an overlap-penalty
method (each state is optimized against $H+w\sum_k|\psi_k\rangle\langle\psi_k|$
for the states already found), which can occasionally converge to a
spurious, non-eigenstate stationary point instead of the true excited
state (`itensor_version` 2, 3, and `"python"` are all susceptible). Every
call to `get_excited_states`/`get_excited` on a DMRG backend checks each
returned state's energy fluctuation $\langle H^2\rangle-\langle
H\rangle^2$ against the ground state's own (both already computed as a
byproduct of the search) and emits a `UserWarning` if a state's
fluctuation is far above that reference — a cheap, always-on sanity check,
though not a fix; a warned-about state's energy and wavefunction should
be treated with caution (e.g. cross-checked against `mode="ED"` on a
smaller system, or recomputed with a larger `scale`/more sweeps).

**Sector (charge) gaps.** For a conserved quantity $A$ with $[H,A]=0$
(e.g. total particle number $\hat N=\sum_i n_i$, or total $S^z$), the gap
to the lowest state with $\langle A\rangle$ shifted by $d$ from the
ground-state sector is obtained by adding a Lagrange-multiplier penalty
and increasing $\lambda$ until the constraint is satisfied:

$$H_\lambda=H+\lambda\big(A-\langle A\rangle_0-d\big)^2,\qquad \Delta_A(d)=E_\lambda-E_0$$

```python
gap_charge = fc.get_charge_gap(d=2)   # gap to add/remove a pair of particles
```

This is exactly how a superconducting/charging gap is extracted in a
finite fermionic chain: $d=2$ probes the energy cost of adding a pair of
particles (Cooper-pair-like excitation) rather than a single particle,
which avoids odd/even-parity finite-size artifacts.

**Example: the Haldane gap.** For the spin-1 Heisenberg chain, the
famous Haldane gap between the (nearly four-fold degenerate, due to
fractionalized $S=\tfrac12$ edge states) ground-state manifold and the
first bulk excitation is

```python
es = sc.get_excited(n=6)
haldane_gap = es[4] - es[0]     # skip the 4 near-degenerate edge states
```

**Non-Hermitian Hamiltonians.** For $H\neq H^\dagger$ (complex hopping
amplitudes, gain/loss terms, PT-symmetric models, ...) the "ground
state" convention throughout dmrgpy is the eigenvalue with the smallest
real part. Eigenvalues come in left/right eigenpairs,

$$H\,|\psi_R\rangle=E\,|\psi_R\rangle,\qquad
H^\dagger|\psi_L\rangle=E^{*}|\psi_L\rangle,\qquad
\langle\psi_L|\psi_R\rangle=1 .$$

On every session backend (ITensor v3, ITensor v2, and the pure-Python
engine — `itensor_version` 2, 3 or `"python"`), `gs_energy()` solves
this with a genuine non-Hermitian DMRG (NH-DMRG) — a port of
[ITensorNHDMRG.jl](https://github.com/tipfom/ITensorNHDMRG.jl) in its
default configuration: independent Arnoldi solves of the two-site
eigenproblem for $H$ and $H^\dagger$ ("onesided" solver, targeting the
smallest real part), with both MPS truncated by the same isometry
obtained from the Hermitian average $\rho=(\rho_L+\rho_R)/2$ of the
left/right reduced density matrices (the "fidelity" algorithm of
Yamamoto et al., [PRB 105, 205125
(2022)](https://doi.org/10.1103/PhysRevB.105.205125)). The full eigenpair
is available directly:

```python
e, psil, psir = sc.nhdmrg()     # E, left and right eigenvector MPS
sc.gs_energy()                  # same E; stores psir as the ground state
```

Because a non-Hermitian "energy" is not a variational bound, `nhdmrg()`
certifies convergence through both eigen-residuals
$\lVert H|\psi_R\rangle-E|\psi_R\rangle\rVert$ and
$\lVert H^\dagger|\psi_L\rangle-E^{*}|\psi_L\rangle\rVert$ and re-runs
from a fresh random state (up to `ntries` times) if either stalls;
`krylovdim`/`restarts` tune the per-bond Arnoldi effort. The MPS Arnoldi
route (`get_excited_states`, now used for non-Hermitian *excited* states,
$n\ge 2$) remains available on every backend, but is typically several
orders of magnitude less accurate than NH-DMRG at comparable cost — see
`examples/non_hermitian/nhdmrg_VS_ED_VS_arnoldi`, which cross-checks
NH-DMRG on all three ported backends against exact diagonalization and
the Arnoldi route on an interacting fermionic chain with a staggered
imaginary potential.

`itensor_version="julia_live"` also implements `nhdmrg()`, but is not a
port: it calls the real ITensorNHDMRG.jl package, with two adaptations
applied in `mpsjulialive/nhdmrg.jl` (see its header for the full
derivation of each). First, that package's "adjoint" sweep is against
$H^{T}$ rather than $H^\dagger$, so its left vector solves the transpose
eigenvalue equation and is complex-conjugated here to match dmrgpy's
convention. Second, it does not tie its left solve to whichever
eigenvalue its right solve picked, so a spectrum with a
complex-conjugate pair tied for the smallest real part (the generic
PT-symmetric/Hatano-Nelson situation) could return a left and a right
vector belonging to *different* eigenvalues; this is broken by re-solving
against $e^{i\theta}H$, which leaves every eigenvector untouched and
rotates the spectrum just enough to separate the pair's real parts. The
rotation can itself re-target a different eigenvalue when
$|\mathrm{Im}\,\lambda|$ is large, so the untied run's eigenvalue is kept
as an anchor and a tie-break result is accepted only if it reproduces it;
otherwise the attempt is failed rather than answered. Neither adaptation is observable on a complex-
*symmetric* $H$, which is what most textbook non-Hermitian models (and
every other non-Hermitian example in this repository) happen to be --
see `examples/non_hermitian/nhdmrg_julia_asymmetric_VS_ED` for a chain
with asymmetric hopping, where both are. The biorthogonal pair $|\psi_R\rangle,|\psi_L\rangle$
is also what feeds the non-Hermitian dynamical correlator,
`get_dynamical_correlator(submode="KPM")` for $H\neq H^\dagger$ — see §6.

Two independent MPS Arnoldi implementations are available, both
matrix-free (they only ever apply $H$ to a wavefunction, never build a
matrix) and both usable in ED mode or DMRG mode on any backend.
`dmrgpy.mpsalgebra.lowest_energy_non_hermitian_iram` (`algebra/
arpacktk.py`) is an Implicitly Restarted Arnoldi Method (IRAM), adapted
from [ARPACK](https://bitbucket.org/chaoyang2013/arpack)'s
`znaupd`/`znaup2`/`znaitr`/`znapps` — the same exact-shift
polynomial-filter restart ARPACK's own `eigs`-style solvers use,
re-derived here for dmrgpy's own MPS/ED wavefunction objects instead of
flat arrays (no BLAS/LAPACK Fortran dependency). It is the **default**
MPS Arnoldi solver (`dmrgpy.mpsalgebra.lowest_energy_non_hermitian`, and
the route `get_excited_states` uses for non-Hermitian $H$ with $n\ge2$,
`excited_states_non_hermitian` in `excited.py`), since it compresses and
reuses its existing Krylov subspace instead of rebuilding it from
scratch on every restart, needing noticeably fewer $H\,|\psi\rangle$
applications to reach the same tolerance on most spectra.
`dmrgpy.mpsalgebra.lowest_energy_non_hermitian_arnoldi` (`algebra/
arnolditk.py`, a restarted Arnoldi method with explicit Rayleigh-Ritz
reseeded restarts) remains available directly for comparison — the two
can trade places on hard (near-degenerate) spectra, see
`examples/non_hermitian/arnoldi_vs_iram_benchmark`, which benchmarks
both head to head (Op-count, wall time, accuracy) in both ED and DMRG
mode — and is still used internally for Krylov operators with a
projector baked in (`fermionchain.py`'s `Spinon_Chain.get_gs`, which
IRAM's interface doesn't (yet) support).

`arpacktk.py` also implements ARPACK's shift-invert mode
(`dmrgpy.mpsalgebra.mpsiram_shift_invert`, ARPACK's mode 3,
$\mathrm{OP}=(A-\sigma I)^{-1}$), used to find the eigenvalues closest to
a target energy $e$ rather than an extremal one — `degeneracy.py`'s
`eigenvalue_degeneracy` (used by `gs_degeneracy(mode="DMRG")` on
non-Hermitian Hamiltonians) is built on it. Since dmrgpy has no exact
MPO inverse — `self.applyinverse` is itself only an iterative
correction-vector solve — this departs from ARPACK's own mode-3
assumption that $\mathrm{OP}$'s eigenvalues are *exactly*
$1/(\lambda-\sigma)$: only the IRAM restart *direction* (which Krylov
directions the implicit shifts filter away) comes from
$\mathrm{OP}$'s own cheap Hessenberg matrix (`which="LM"`, since the
largest $|\mathrm{OP}\text{-eigenvalue}|$ is the $H$-eigenvalue closest
to $e$); convergence and the reported eigenpairs are instead recomputed
every outer iteration from $H$'s own exact (still cheap, $O(\text{ncv}^2)$)
representation on the Krylov basis $\mathrm{OP}$ builds — the same
accommodation arnolditk's own `mode="ShiftInv"` path already made.

`arpacktk.py` also implements ARPACK's mode 2 — the generalized
eigenproblem $A|\psi\rangle=\lambda M|\psi\rangle$ for a Hermitian
positive-definite $M$
(`dmrgpy.mpsalgebra.mpsiram_generalized`/`generalized_excited_states`,
$\mathrm{OP}=M^{-1}A$, $B=M$). Unlike mode 1/3, the Krylov basis must be
orthonormal in the $M$-weighted inner product
$\langle u,v\rangle_M=\langle u|M|v\rangle$ rather than the plain one, so
building it uses a separate routine
(`arnoldi_extend_generalized`) rather than adding a conditional $M$ to
the mode-1/3 one. $\mathrm{OP}$ again goes through
`self.applyinverse` (the same approximate-inverse caveat as mode 3, here
with the trivial shift $\sigma=0$) — but unlike mode 3, $\mathrm{OP}$'s
own Hessenberg-matrix eigenvalues stay directly meaningful (there is no
shift to undo), so convergence/selection follow the same pattern as
plain `mpsiram`, just built with the $M$-inner product. Verified against
`scipy.linalg.eigh`'s generalized Hermitian-definite eigensolver in both
ED and DMRG mode.

**Generalized-eigenvalue DMRG.** For the same problem
$H|\psi\rangle=\lambda A|\psi\rangle$ ($A$ Hermitian positive definite),
`gs_energy_generalized(A)` solves it with a genuine DMRG sweep instead of
a Krylov method, needing no approximate operator inverse at all (unlike
`mpsiram_generalized` above, which must invert $M$ iteratively since
dmrgpy has no exact MPO inverse):

```python
lam = fc.gs_energy_generalized(a_metric_operator)   # smallest lambda
```

The trick is a self-consistent Lagrange multiplier: minimizing
$\langle\psi|H|\psi\rangle$ subject to the metric normalization
$\langle\psi|A|\psi\rangle=1$ has stationarity condition
$(H-\lambda A)|\psi\rangle=\mu|\psi\rangle$ for multiplier $\lambda$ --
the *ordinary* ($\mu$,$\psi$) eigenproblem of the plain-normalized shifted
operator $H-\lambda A$, exactly what a standard two-site DMRG sweep
already finds. At $\mu=0$ this is precisely
$H|\psi\rangle=\lambda A|\psi\rangle$, so each outer iteration (i) builds
the MPO $H-\lambda A$ from the current $\lambda$ estimate, (ii) runs one
ordinary DMRG sweep against it, then (iii) updates $\lambda$ to the
freshly-swept state's generalized Rayleigh quotient
$\langle\psi|H|\psi\rangle/\langle\psi|A|\psi\rangle$ -- one outer
iteration per `Sweeps` schedule entry, so bond dimension ramps exactly as
an ordinary `gs_energy()` run's own schedule does. $A=\mathrm{Id}$
reduces this exactly to plain ground-state DMRG.

Implemented on the pure-Python (pyitensor) backend
(`itensor_version="python"`, i.e. after `chain.setup_python()`) --
`dmrgpy.pyitensor.dmrg.dmrg_generalized` is the underlying routine --,
on compiled ITensor v3 (`itensor_version=3`, i.e. after
`chain.setup_cpp(version=3)`), where `Chain::gs_energy_generalized`
(`mpscpp3/chain_session.h`) runs the identical algorithm against ITensor
v3's own `dmrg()`/`Sweeps`/`sum()` instead of pyitensor's hand-rolled
two-site sweep, and on the live Julia backend
(`itensor_version="julia_live"`, i.e. after `chain.setup_julia()`), where
`get_gs_generalized` (`mpsjulialive/generalized.jl`) runs it once more
against ITensorMPS.jl's own `dmrg()`/`Sweeps`/`add()`. `itensor_version=2`
(mpscpp2) doesn't have this session method yet and raises
`NotImplementedError`. See
`examples/groundstate/dmrg_generalized_benchmark`, which heads all three
solvers up against each other on the same interacting-fermion-chain test
problem: needing no approximate inverse, both DMRG routes are far more
accurate than ARPACK mode 2 (whose own correction-vector solve also gets
more expensive, and less accurate, as the chain grows -- two orders of
magnitude slower and $10^{9}$ times less accurate on an 8-site chain),
and v3 is itself consistently ~2-4x faster than pyitensor at these sizes.
`julia_live` matches both DMRG routes to machine precision and, once
warm, lands between them: at n=8, v3 0.16s, `julia_live` 0.20s,
pyitensor 0.50s. Its *first* call in a session is dominated by Julia's
JIT compilation (~46s), so a single-point timing of this backend measures
compilation, not the solver -- read its second and later rows.

**Non-Hermitian generalized-eigenvalue DMRG.** `gs_energy_generalized`
also accepts a non-Hermitian $H$: `self.hamiltonian` is checked for
Hermiticity, and a non-Hermitian one is transparently dispatched to a
dedicated solver (nhdmrg.py's `nhdmrg_generalized`) instead of raising --
same calling convention, same "smallest real part" convention as
non-Hermitian `gs_energy()`/NH-DMRG (§4) above:

```python
lam = fc.gs_energy_generalized(a_metric_operator)   # complex lambda, smallest Re
```

The metric $A$ must still be Hermitian positive definite (mirroring
`mpsiram_generalized`'s own $M$ precondition -- and, following ARPACK's
own convention there, the *primary* operator needs no Hermiticity
assumption at all, only the metric does). The self-consistent
Lagrange-multiplier trick generalizes directly: since $A$ is Hermitian,
$(\lambda A)^\dagger=\bar\lambda A$, so minimizing (in the NH-DMRG sense
-- smallest real part, not a variational bound) subject to the metric
*biorthogonal* normalization $\langle\psi_L|A|\psi_R\rangle=1$ gives
stationarity condition
$(H-\lambda A)|\psi_R\rangle=\mu|\psi_R\rangle$,
$(H-\lambda A)^\dagger|\psi_L\rangle=\bar\mu|\psi_L\rangle$ -- the
*ordinary* ($\langle\psi_L|\psi_R\rangle=1$-normalized) NH-DMRG
eigenproblem of the shifted pair $(H-\lambda A,\,H^\dagger-\bar\lambda A)$,
exactly what one ordinary NH-DMRG sweep already finds. At $\mu=0$ this is
precisely $H|\psi_R\rangle=\lambda A|\psi_R\rangle$, so each outer
iteration (i) builds $H-\lambda A$ and $H^\dagger-\bar\lambda A$ from the
current (complex) $\lambda$ estimate, (ii) runs one ordinary NH-DMRG
sweep against them, then (iii) updates $\lambda$ to the freshly-swept
pair's generalized *biorthogonal* Rayleigh quotient
$\langle\psi_L|H|\psi_R\rangle/\langle\psi_L|A|\psi_R\rangle$ in place of
the Hermitian case's plain one. $A=\mathrm{Id}$ reduces this exactly to
plain NH-DMRG.

Implemented on the pure-Python (pyitensor) backend
(`dmrgpy.pyitensor.nhdmrg.nhdmrg_generalized`), on compiled ITensor v3
(`Chain::nhdmrg_generalized`, `mpscpp3/chain_session.h`, a line-for-line
port against this file's own hand-rolled two-site sweep) and on the live
Julia backend (`get_gs_generalized_nhdmrg`,
`mpsjulialive/generalized.jl`, the same outer loop wrapped around real
ITensorNHDMRG.jl sweeps) -- `mpscpp2`
(`itensor_version=2`) still raises `NotImplementedError` for a
non-Hermitian $H$, no analogous session method there. Since
`nhdmrg_generalized`'s two-site sweep never calls ITensor v3's own
`dmrg()` (it is hand-rolled directly against a restarted Arnoldi solve
and manual ITensor contractions, unlike the Hermitian path above), it
does *not* need the Hermitian path's short-chain guard: chains shorter
than 3 sites work fine here on `itensor_version=3`. See
`examples/non_hermitian/nhdmrg_generalized_benchmark`, which heads all
three routes up against each other on a non-Hermitian
interacting-fermion-chain test problem (a staggered imaginary on-site
potential): ARPACK mode 2 needs no adaptation at all for a non-Hermitian
primary operator (its $M$-positive-definite precondition was always the
only one), so this is a fair like-for-like comparison, and the same
pattern as the Hermitian benchmark holds -- needing no approximate
inverse, both DMRG routes are far more accurate and (as the chain grows)
up to two orders of magnitude faster than ARPACK mode 2. `julia_live` also matches to machine precision here, sitting
just ahead of pyitensor once warm (n=8: v3 1.16s, `julia_live` 2.26s,
pyitensor 2.42s) with the same first-call JIT caveat as the Hermitian
benchmark. Unlike the
Hermitian case, v3 isn't consistently faster than pyitensor here (roughly
on par at the sizes benchmarked) -- NH-DMRG's per-bond cost already pays
for *two* Arnoldi solves (right block and its adjoint) regardless of
backend, narrowing the compiled-vs-pure-Python gap relative to plain
ground-state DMRG's single local diagonalization per bond.

**Caveat (both Hermitian and non-Hermitian `gs_energy_generalized`).**
Afterward, `fc.wf0`/`fc.e0` hold the eigenstate/eigenvalue of the
*generalized* problem, not a plain eigenstate of `fc.hamiltonian` alone --
every other method that reads `fc.wf0` as an ordinary ground state
(`get_excited_states()`, any dynamical/KPM correlator, ...) has no way to
tell the difference and will silently build on the wrong reference state
if called afterward. Call `gs_energy_generalized()` as the last step of a
calculation, or recompute a genuine ground state (`gs_energy()`/`nhdmrg()`)
first if you need one of those other methods too.

**Effective low-energy Hamiltonians.** Given the $n$ lowest eigenstates
$\{|\psi_k\rangle\}$ and the projector $P=\sum_k|\psi_k\rangle\langle\psi_k|$
onto the manifold they span, the projected Hamiltonian $PHP$ can be fitted
as a real linear combination of the same projections of a chosen set of
operators $O_a$:

$$PHP\simeq\sum_a J_a\,PO_aP,\qquad J_a=\arg\min_{J}\Big\|PHP-\sum_a J_aPO_aP\Big\|$$

and the fitted $J_a$ are the effective couplings.

```python
# fit against an explicit list of operators -> one coefficient each
J = fc.get_heff(operators=[fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1],
                           fc.Sz[0]*fc.Sz[1]], n=4, mode="ED")

# ...or against all pairwise products of the chain's own Sx/Sy/Sz
from dmrgpy import effectivehamiltonian
coef = effectivehamiltonian.get_effective_hamiltonian_couplings(fc, n=4)
latex = sc.get_effective_hamiltonian(n=1)   # the same fit, as a latex string
```

`n` selects the manifold, `tol` is the magnitude below which a fitted
coupling is dropped, and `operators` overrides the default basis
(`get_projection_operators`, every $S^x/S^y/S^z$ of the chain).
`method="single"` builds $PO_aPO_bP$ — the product of the two *projected*
operators, one projection per operator — while `method="full"` builds
$P(O_aO_b)P$, the projection of the product, at one projection per
*pair*. These are the same operator only when the manifold is closed
under $O_a$; they are not two implementations of one thing.

Two ways this calculation can be meaningless, both of which are now
refused rather than returned:

- **`n` must not cut a degenerate multiplet.** The retained states are
  then not a symmetry-invariant subspace and the fitted couplings are
  physically meaningless (on a Hubbard dimer, whose low manifold is a
  singlet plus a three-fold triplet, `n=3` gave couplings ~50x the
  correct ones). A `ValueError` names the degeneracy and the values of
  `n` that do cut cleanly.
- **The operator basis must be linearly independent on that manifold.**
  Otherwise the coefficients can be shifted by any null-space vector
  without changing the fit, so no unique set of couplings exists; the
  raw candidate set is typically very much overcomplete (37 candidates
  spanning a rank-16 space on the dimer). A `ValueError` reports the
  rank and the number of null directions.

The fit itself is a linear least-squares solve, so it is exact and
reproducible; the eigenstates are Löwdin-orthonormalized first when they
need it, which DMRG's overlap-penalty excited states generally do. See
`examples/groundstate/effective_hamiltonian`, which recovers the exact
$J\cos2\phi$/$J\sin2\phi$ twist of a Hubbard dimer's exchange under a
Peierls phase $\phi$.


## 5. Entanglement and quantum information

**Entanglement entropy of a real-space bipartition.** Cutting the chain
at bond $(i,i+1)$, the Schmidt decomposition of the ground state is
$|\mathrm{GS}\rangle=\sum_\alpha\lambda_\alpha|\alpha\rangle_L\otimes|\alpha\rangle_R$
and the von Neumann entanglement entropy of either half is

$$S_{i}=-\sum_\alpha\lambda_\alpha^2\log\lambda_\alpha^2=-\mathrm{Tr}\big[\rho_L\log\rho_L\big],\qquad \rho_L=\mathrm{Tr}_R|\mathrm{GS}\rangle\langle\mathrm{GS}|$$

```python
s = wf.get_bond_entropy(i, i+1)
```

`get_site_entropy(i)`/`get_pair_entropy(i,j)` instead build a reduced
density matrix from local *projectors* (e.g. $S^z=\pm\tfrac12$ for
spin-$\tfrac12$, occupation number for fermions) rather than a full MPS
bond cut, which lets you ask about the entanglement of a single site or
a pair of (possibly non-adjacent) sites with the rest of the system.

**Mutual information** between two sites (or subsystems) $i,j$:

$$I(i,j)=S_i+S_j-S_{ij}$$

```python
mi = wf.get_mutual_information(1, 2)
```

quantifies total (classical + quantum) correlation between $i$ and $j$,
and is often used to map out effective couplings or emergent
degrees of freedom (e.g. between the two edge spins of a Haldane chain).

**CFT central charge.** For a critical (gapless) chain, the
entanglement entropy of a cut at position $\ell$ in an open chain of
length $L$ follows the Calabrese-Cardy formula:

$$S(\ell)=\frac{c}{6}\log\!\left[\frac{2L}{\pi a}\sin\frac{\pi\ell}{L}\right]+\text{const.}$$

DMRGPY computes $S(\ell)$ at every bond and least-squares fits this
formula to extract the central charge $c$ — the universal number that
identifies the underlying conformal field theory (e.g. $c=\tfrac12$ for
the critical Ising chain, $c=1$ for a free boson / XX chain):

```python
c = wf.get_CFT_central_charge()
```

**Single-particle correlation matrix and orbital entanglement.** For
fermionic systems, the correlation matrix

$$C_{ij}=\langle c_i^\dagger c_j\rangle$$

is directly accessible (`sc.get_correlation_matrix()`), and diagonalizing
it, $C=U\,n\,U^\dagger$, gives natural-orbital occupations
$n_\alpha\in[0,1]$ (`get_correlated_orbitals`) and an orbital-resolved
entanglement entropy

$$S=-\sum_\alpha\Big[n_\alpha\log n_\alpha+(1-n_\alpha)\log(1-n_\alpha)\Big]$$

(`get_correlation_entropy`). `get_highorder_correlation_matrix` and
`get_four_correlation_tensor` give the corresponding two-particle
correlators $\langle c_i^\dagger c_j^\dagger c_l c_k\rangle$ and
$\langle c_i^\dagger c_j c_k^\dagger c_l\rangle$.

`get_four_correlation_tensor(ctmode=...)` has three implementations:
`ctmode="explicit"` (backend-agnostic Python loop of `vev()`s, always
available), `ctmode="full"` (native per-element AutoMPO build — C++ for
`itensor_version` `2`/`3`, pure Python for `"python"` — builds and
applies an independent AutoMPO per $(i,j,k,l)$ tuple), and
`ctmode="sweep"` (`itensor_version` `3` or `"python"`, plain
non-native-spinful fermionic sites only) — a single-sweep,
environment-reuse implementation following the algorithmic idea of
[ITensorCorrelators.jl](https://github.com/ITensor/ITensorCorrelators.jl):
rather than rebuilding a fresh MPO for every tuple, it reuses partial
tensor-network contractions across the whole $(N,N,N,N)$ tensor. Agrees
with `ctmode="full"` and ED to machine precision / solver tolerance on
both backends, and is substantially faster: at $n=12$, 9.3s → 1.2s under
`itensor_version=3` and 37.7s → 2.0s under `"python"` (measured
single-threaded — see §19 on BLAS threads), i.e. roughly 8x and 19x,
against a `ctmode="full"` that is slower still.

Those numbers reflect a fix worth knowing about if you read the older
docstrings. The tensor splits into pairwise-distinct-index entries (handled
by the sweep) and repeated-index entries, and the latter used to fall back
to the same per-tuple AutoMPO build `ctmode="full"` uses. There are only
$O(n^3)$ of those against the sweep's $O(n^4)$, which is why they were
described as subdominant — but a smaller count times a far more expensive
per-tuple cost still dominates, and measured, that fallback was 96% of the
whole runtime at $n=12$. Every one of those operators is a product of four
$C^\dagger/C$ factors on at most 3 distinct sites, i.e. *local*, so both
backends now fold them directly instead of compiling an MPO over the whole
chain. `accelerate` (default `True`, accepted on every `ctmode`) still only
affects the repeated-index entries for `ctmode="sweep"` — the
pairwise-distinct sweep has no equivalent conjugate-pair saving to skip
(see either backend's own docstring) — so don't expect the usual ~2x win
`accelerate` gives `ctmode="full"`. See
`examples/staticcorrelators/four_correlation_tensor_sweep_VS_full`,
`Chain::four_correlation_tensor_sweep` (`mpscpp3/chain_session.h`) and
`pyitensor/chain.py`'s port for the algorithm.

`ctmode="fold"` (`itensor_version="python"` or `3`) evaluates every tuple as
a *local* operator fold — no MPO is built for any of them. It is
flavour-agnostic: it reads each mode's `(operator name, site)` off the
chain's own `C`/`Cdag`, so it covers spinless chains and
`Spinful_Fermionic_Chain_Native` alike. That matters because native spinful
sites previously had only `ctmode="explicit"` under `"python"`, which builds
an MPO and sweeps the whole chain per tuple; `fold` is exact against it to
machine precision and 8–12x faster (5 sites / 10 modes: 9.8s → 0.9s), and is
now the default there. The C++ backend has the same port, where it replaces a
`ctmode="full"` that was itself *slower* than the pure-Python fold: 2.8s → 0.7s
at 5 sites, a 4x win and an algorithm-beats-language result worth remembering. It does not reuse environments *across* tuples the way
`ctmode="sweep"` does, so for a plain spinless chain prefer `"sweep"`.

`ctmode=None` (the default) auto-selects the fastest method actually
available for the wavefunction's backend/chain type: `"sweep"` whenever
it applies, else `"full"` whenever it applies (any `itensor_version` in
`2`/`3`/`"python"`, or `Spinful_Fermionic_Chain_Native` under
`itensor_version=3` via its own `four_correlation_tensor_spinful()`),
else the always-correct `"explicit"` fallback (e.g. for
`itensor_version="julia_live"`, which has no `"full"`/`"sweep"`
implementation). Passing a `ctmode` explicitly is still a hard request —
it raises rather than silently falling back if that method isn't
available for the wavefunction at hand. `Spinful_Fermionic_Chain_Native`
always needs `ctmode="full"` (or the default, which resolves to it):
the sweep has no native-spinful-site counterpart yet.

## 6. Dynamical (frequency-dependent) correlators

The central quantity is a retarded correlator (Green's function) between
two operators $A,B$,

$$G_{AB}(\omega)=\langle\mathrm{GS}|A\,\frac{1}{\omega-H+E_0+i\delta}\,B|\mathrm{GS}\rangle,\qquad S_{AB}(\omega)=-\frac1\pi\,\mathrm{Im}\,G_{AB}(\omega)$$

where $\delta$ is a small broadening. Choosing $A=B=S^z_i$ at the same
site gives the local dynamical spin structure factor $S^{zz}_{ii}(\omega)$
(what a local probe like NMR/ESR couples to); choosing $A=S^z_i$,
$B=S^z_j$ at different sites and Fourier-transforming over $i-j$ gives
the momentum-resolved dynamical structure factor $S(q,\omega)$ measured
in inelastic neutron scattering. All five submodes below compute
$G_{AB}(\omega)$ (or equivalently $S_{AB}(\omega)$); they differ only in
*how*, and therefore in what energy range/resolution/cost trade-off they
offer:

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                      name=(sc.Sz[0], sc.Sz[0]))
```

**`submode="KPM"` — Kernel Polynomial Method.** $H$ is linearly rescaled
into $\tilde H=(H-b)/a$ so its full spectrum lies in $[-1,1]$
(`kpm_scale`, default 0.7, sets how much margin is left inside $[-1,1]$),
then $G_{AB}$'s spectral function is expanded in Chebyshev polynomials
$T_m(\tilde H)$:

$$\mu_m=\langle\mathrm{GS}|A\,T_m(\tilde H)\,B|\mathrm{GS}\rangle,\qquad S_{AB}(\tilde\omega)=\frac{1}{\pi\sqrt{1-\tilde\omega^2}}\left[\mu_0+2\sum_{m=1}^{N-1}g_m\,\mu_m\,T_m(\tilde\omega)\right]$$

with Jackson-kernel damping coefficients $g_m$ suppressing Gibbs
ringing from the truncation at $N$ moments (`kpm_n_scale` scales $N$,
i.e.\ the energy resolution, relative to the rescaled bandwidth). This is
the default, general-purpose method: robust, works across the whole
spectrum at once, cost grows only linearly with the number of moments.
The band edges entering the rescaling are obtained variationally: the
lower edge reuses the ground-state energy, and the upper edge runs a
deliberately reduced-effort DMRG on $-H$ (few sweeps at modest bond
dimension) — it is only a spectral *bound*, protected by the
`kpm_scale` margin, not a physical result, and a variational
underestimate only shrinks the number of moments. If the bound is ever
too tight for the chosen `kpm_scale`, the moment recursion detects the
resulting exponential divergence and aborts with an explicit error
rather than returning a silently wrong spectrum.

**Energy truncation (`itensor_version="python"` and `3`) — narrowing the KPM window
below the full bandwidth.** The default `kpm_scale` rescales $H$ so its
*entire* many-body spectrum fits in $[-1,1]$, which is always safe but
wastes resolution: a local operator's correlator usually has real
spectral weight only over a width $W_A$ much smaller than the full
bandwidth $W$, so most of the Chebyshev expansion's $N$ moments are
"spent" on frequency regions the correlator never visits. Choosing a
narrower window (smaller `kpm_scale`) concentrates the same $N$ onto the
physically relevant region instead — but without a safeguard, the
moment recursion's own numerical noise lets a little high-energy weight
leak in, and since Chebyshev polynomials are unbounded outside $[-1,1]$,
that leakage grows exponentially over the recursion (exactly the failure
the `kpm_scale`-too-tight error above is detecting, just deliberately
triggered by choosing a small window on purpose rather than by an
inaccurate band-edge estimate). *Energy truncation* (Holzner,
Weichselbaum, McCulloch & von Delft, "Chebyshev matrix product state
approach for spectral functions", [PRB 83, 195115
(2011)](https://doi.org/10.1103/PhysRevB.83.195115), Sec. III-B) is the
fix: after every Chebyshev vector is formed, it is swept site by site,
and at each site a small ($\le$ `kpm_truncate_dK`) Krylov subspace of
that site's local effective Hamiltonian is diagonalized and any
component with rescaled energy $|\varepsilon|\ge$ `kpm_truncate_threshold`
is projected out, before the vector is used for a moment or fed into the
next recursion step. This is a *precautionary* measure, not an exact
projector (a finite Krylov dimension only approximates the local
spectrum), so a small residual above threshold is expected and does not
vanish with more `kpm_truncate_nsweeps` — matching the original paper's
own characterization of it (its Sec. V.C).

```python
sc.kpm_scale = 0.3            # narrower than the safe ~0.5 floor
sc.kpm_energy_truncate = True # enable energy truncation (default: False)
sc.kpm_truncate_dK = 30       # per-site Krylov subspace dimension
sc.kpm_truncate_nsweeps = 10  # number of truncation sweeps per vector
sc.kpm_truncate_threshold = 1.0  # rescaled-energy cutoff (paper's eps_P)
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=(sc.Sz[0], sc.Sz[0]))
```

`kpm_truncate_dK` is a *convergence* knob, not a tuning one: the
projection it defines is exact once $d_K$ reaches the dimension of a
site's local space, so raise it until the reconstructed spectrum stops
moving and then stop (larger $d_K$ only costs time, since the per-site
Krylov build is where energy truncation's whole overhead lives). The
paper's own Table I recommends $d_K=30$. Both backends
build that subspace with full re-orthogonalization; a single
Gram-Schmidt pass is not enough at $d_K\sim 30$, where the Krylov
vectors converge toward the site's locally dominant eigendirection and
a singly-orthogonalized basis silently ceases to be orthonormal (which
used to corrupt pyitensor's projection and move a 6-site Heisenberg
chain's spectral peak from $0.48\,J$ to $1.06\,J$).

Enabling `kpm_energy_truncate` also switches the rescaling convention
itself: instead of centering the window on the full bandwidth's
*midpoint*, it anchors it at the ground state $E_0$ (the paper's own Eq.
21b), placing $E_0$ at one edge of the window and $E_0+W_s$
($W_s=(E_{\max}-E_{\min})\cdot$ `kpm_scale`) at the other. This matters
because a dynamical correlator's Chebyshev vectors are built by acting
$A$/$B$ on the ground state, so the physically relevant region sits just
above $E_0$, not around the spectrum's geometric middle — the
midpoint-centered window would otherwise clip the ground state itself out
of the window before ever reaching a genuinely useful, narrower regime.
Available for `itensor_version="python"` and `itensor_version=3`; there
is no `itensor_version=2` port (mpscpp2 has no equivalent machinery to
build a per-site local effective Hamiltonian from, unlike mpscpp3's
`LocalMPO`/`diagHermitian`). Requesting it on `itensor_version=2` is
silently a no-op (the existing, always-safe rescaling is used instead).
The v3 port (`mpscpp3/chain_session.h`'s
`kpm_dynamical_correlator_truncated()`) is a wholly independent method
from `kpm_dynamical_correlator()` — a deliberate design choice so the
existing, always-safe v3 KPM path is never touched by this feature — not
a branch inside the same function the way the pyitensor port is; both
implement the identical algorithm and agree on the same physical answer
(see `test_kpm_energy_truncation_v3_accuracy.py`'s cross-backend check).
Available only for the `(A, B)`-operator dynamical correlator, not the
lower-level "arbitrary operator" KPM (`general_kpm`/`kpm_wfa_wfb`), which
has no ground-state reference to anchor to — setting the flag has no
effect on those, on any backend (they expand an already-bounded operator
rescaled into $[a,b]$ by `scale_operator()`, which needs no safeguard). See
`examples/dynamical_correlator/dynamical_correlator_kpm_energy_truncation`
for a worked example (including the divergence this fixes, reproduced on
purpose), `src/dmrgpy/pyitensor/kpm_energy_truncation.py` for the
pyitensor implementation, and `mpscpp3/chain_session.h`'s own
`kpm_dynamical_correlator_truncated`/`kpm_energy_truncate` comments for
the v3 one.

Performance note: energy truncation is a controlled *slowdown*, not a
speedup — it buys resolution/feasibility (a window that would otherwise
diverge), not raw speed. Measured directly (a 4-site chain, `kpm_scale`
narrowed from the safe 0.7 to 0.65): the pyitensor backend costs
~5.5–12.5x more per Chebyshev moment (`kpm_truncate_dK`/`nsweeps` of
10/3 vs. the paper's own recommended 30/10, respectively) than the
untruncated path, and v3's native port ~15–32x more per moment for the
same settings on this small system — the per-site Krylov-subspace cost
scales with bond dimension, so this overhead grows quickly with system
size (measured directly at 8 sites: ~245x per moment at
`kpm_truncate_dK=30`/`nsweeps=10`). Narrowing the window does reduce the
moment count needed for a given frequency resolution, but that reduction
only outweighs the per-moment cost once the correlator's own spectral
width is genuinely much smaller than the full bandwidth (the paper's own
premise for large, extensive systems) — for a small test system the two
effects roughly cancel or net out to a slowdown, not a speedup.

**`submode="CVM"` — correction-vector method.** Instead of a global
polynomial expansion, this solves directly for the correction vector at
one frequency $\omega$ at a time, via the positive-definite linear system

$$\big[(H-\omega-E_0)^2+\eta^2\big]\,x_c=-\eta\,B|\mathrm{GS}\rangle$$

(solved by conjugate gradient in MPS form), from which

$$x=(\omega+E_0+i\eta-H)^{-1}B|\mathrm{GS}\rangle=i\,x_c+\frac{H-\omega-E_0}{\eta}\,x_c,\qquad G_{AB}(\omega)=\langle\mathrm{GS}|A|x\rangle$$

Here $\eta$ (`delta`) is the artificial broadening that regularizes the
resolvent at a real frequency. CVM is more accurate at a single targeted
frequency/energy window (e.g. zooming in on a sharp resonance) than a
global KPM expansion, at the cost of re-solving the linear system for
every $\omega$ on the requested grid. Two things keep that per-$\omega$
cost down. First, everything $\omega$-independent is computed once per
sweep instead of once per point — in particular the right-hand side
$-\eta B|\mathrm{GS}\rangle$. Second, the CG loop terminates early once
the bond-dimension truncation (`cvm_maxm`) puts a floor under the
reachable residual — past that floor the iteration cannot improve the
tracked best solution any further (the truncated recurrence in fact
diverges), so the solver returns the best correction vector seen
instead of burning the full `cvm_nit` iteration budget (`cvm_patience`
sets how many iterations without meaningful improvement conclude the
floor is reached, and `cvm_blowup` how far past the best residual the
running one may diverge before stopping). Neither feature
changes the answer. Each point reports its CG iteration count and best
residual; if that residual stalls far above `cvm_tol`, the correction
vector is not converged at this bond dimension and the fix is a larger
`cvm_maxm`, not more iterations. (Warm-starting each point's CG from
the neighboring point's correction vector was tried and measured to
*hurt* — truncated CG from a nearby-but-wrong start can stagnate at a
much worse residual than from the cold start — so each point is solved
independently.)

**`submode="ROOTN"` — root-$N$ Krylov-space correction vector.**
Implements Nocera & Alvarez, "Root-$N$ Krylov-space correction-vectors
for spectral functions with the density matrix renormalization group"
([arXiv:2204.03165](https://arxiv.org/abs/2204.03165)). Rather than
building the correction vector $x(\omega+i\eta)=(\omega-H+E_0+i\eta)^{-1}B|\mathrm{GS}\rangle$
in one shot, it is built as $N$ sequential fractional-power steps,

$$x^{p/N}(\omega+i\eta)=\Big(\frac{1}{\omega-H+E_0+i\eta}\Big)^{p/N}B|\mathrm{GS}\rangle,\qquad p=1,\dots,N,$$

where each step re-seeds a Krylov (Lanczos) subspace of dimension `nkry`
with the *previous* step's vector, tridiagonalizes $H$ in that subspace,
and applies the $1/N$-power resolvent in the resulting eigenbasis before
handing the result forward as the next step's seed. `N=1` reduces to the
"conventional" Krylov-space correction vector (Nocera, PRE 2016) that the
paper compares against. In the paper's MPS/DMRG setting, `nkry`'s role is
played by the bond dimension $m$, and building the correction vector in
$N$ smaller steps lets the entanglement grow gradually instead of needing
a single very large $m$ at high target frequency; here, `nkry` directly
caps the Lanczos subspace size at each step, and the same qualitative
effect is observed: at a fixed, small `nkry`, increasing $N$ reduces the
error against the exact answer, especially near the top of the many-body
bandwidth. With `nkry` equal to the full Hilbert space dimension, root-$N$
reproduces the exact answer for any $N$, since the Krylov subspace is
then exact — a useful self-consistency check of the recursion itself,
independent of the Krylov truncation.

```python
(x, y) = sc.get_dynamical_correlator(mode="ED", submode="ROOTN",
                                      name=(sc.Sz[0], sc.Sz[0]),
                                      N=8, nkry=20)
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="ROOTN",
                                      name=(sc.Sz[0], sc.Sz[0]),
                                      N=8, nkry=20) # itensor_version in (2,3,"python")
```

Two implementations exist, both cross-checked to agree at machine
precision on small chains (see
`examples/dynamical_correlator/dynamical_correlator_rootn_ED` and
`examples/dynamical_correlator/dynamical_correlator_rootn_v3_VS_ED`):

- `mode="ED"` (`src/dmrgpy/algebra/rootn.py`): the recursion above against
  the exact ED Hamiltonian, with the Krylov subspace built from plain
  numpy vectors.
- `mode="DMRG"` (`src/dmrgpy/rootndmrg.py`, `itensor_version` in
  `(2,3,"python")`): the *same* recursion, but with the Krylov subspace
  built out of *global* MPS vectors — each Lanczos step is a truncated
  MPO application (`self.toMPO(self.hamiltonian)*v`) rather than a
  local, per-bond update. This is deliberately *not* how the paper's own
  Appendix algorithm works (a multi-target state-averaged DMRG sweep,
  jointly representing $|\mathrm{GS}\rangle$, $B|\mathrm{GS}\rangle$,
  $\mathrm{Re}(x)$, $\mathrm{Im}(x)$ in one block MPS compressed together
  at every bond, so the correction vector's local tensor is built
  directly from the ground state's own environment): a first attempt at
  a more literal, per-bond local-sweep implementation (in
  `mpscpp3/chain_session.h`, using `LocalMPO` and a Lanczos-based
  fractional-power update at each bond, modeled on `applyExp`/TDVP) gave
  numerically *wrong* results (sign-flipping, unstable) when cross-checked
  against exact ED — reapplying the local update bond-by-bond on a single
  self-referential MPS does not correctly realize "apply $f(H)$ once"
  globally, unlike TDVP (a well-defined local Trotter step) or
  ground-state DMRG (repeated local energy minimization provably
  converges to the global minimum). The global-Krylov approach here
  avoids that failure mode entirely by only ever using already-tested
  whole-MPS primitives (truncated MPO application, inner products, MPS
  addition) — the same primitives `submode="CVM"`'s conjugate gradient
  already relies on — at the cost of not reproducing the paper's own
  bond-dimension bookkeeping across the four channels, and of every
  Lanczos step costing a full MPO application rather than a cheap local
  tensor contraction.

**`submode="TD"` — time-dependent DMRG.** Real-time evolution gives the
correlator directly in the time domain,

$$C(t)=\langle\mathrm{GS}|A(t)B(0)|\mathrm{GS}\rangle=e^{iE_0t}\langle\mathrm{GS}|A\,e^{-iHt}\,B|\mathrm{GS}\rangle$$

which is then windowed with a damping/taper factor $w_\delta(t)$ (see
`damping` below; the default is an exponential
$w_\delta(t)=e^{-\delta t}$, equivalent to a Lorentzian broadening of
width $\delta$ in frequency) and Fourier transformed,

$$S_{AB}(\omega)=\frac1\pi\,\mathrm{Re}\!\int_0^{T}\!dt\;e^{i\omega t}\,C(t)\,w_\delta(t)$$

The total simulated time $T$ (`damping_periods`/$\delta$) must be long
enough that the damping has suppressed truncation ringing by $t=T$.
Frequency resolution is set by $T$ (via the usual $\Delta\omega\sim
1/T$ time-frequency uncertainty), so this method is best when you want
fine resolution over a *narrow* frequency window, at the cost of a real
dynamical simulation (bond dimension grows with entanglement generated
during the evolution). The propagator used for this evolution is the
same `sc.tevol_method` selector described in §7, including
`tevol_method="TEBD"` — cheaper per step than TDVP whenever the
Hamiltonian is strictly nearest-neighbor, e.g.\ a
`Spinful_Fermionic_Chain_Native` Hubbard chain (the standard interleaved
`Spinful_Fermionic_Chain` is *not* nearest-neighbor after Jordan-Wigner
threading, so TEBD raises `NotImplementedError` there). A script that
wants the `"TEBD"` speedup on whichever of the two representations turns
out nearest-neighbor, without hardcoding that choice per model, can use
`tevol_method="AUTO"` instead (§7) to fall back to `"TDVP"` automatically
on the interleaved representation rather than raising.

*Sharpening `"TD"`'s lineshape (`damping`, `predict`).* An exponential
window gives $S_{AB}(\omega)$ an exact Lorentzian lineshape, whose
$1/(\omega-\omega_0)^2$ tail is much heavier away from a peak than
`"KPM"`'s default Jackson-kernel reconstruction (see `docs/
td_dynamical_correlator_sharpening_plan.md` for the full design,
literature, and the empirical comparison behind the default below). Two
independent, composable knobs address this:

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="TD",
        name=(sc.Sz[0], sc.Sz[0]), delta=0.05,
        damping="exp",                # default; "gaussian"/"parzen" also available
        predict=True, lp_order=None, lp_extend_factor=10)  # predict=True is the default
```

- `damping` selects $w_\delta(t)$: `"exp"` (default), the exponential
  above; `"gaussian"`, $w_\delta(t)=e^{-(\delta t)^2/2}$, whose Fourier
  transform decays as $e^{-\omega^2}$ (far faster than the Lorentzian's
  algebraic tail), at the cost of a slightly wider FWHM at the same
  $\delta$ ($2\sqrt{2\ln2}\,\delta\approx2.35\delta$ vs the Lorentzian's
  $2\delta$); `"parzen"`, a taper with compact support that vanishes
  (with zero derivative) exactly at $t=T$, which instead targets the
  Gibbs ringing from truncating $C(t)$ at a finite $T$, independent of
  the peak-broadening tradeoff above.
- `predict` (default `True`) extrapolates the raw, undamped $C(t)$ via
  linear prediction (`dynamicstk.linearprediction.linear_predict_extend`,
  following White & Affleck and Barthel, White & Schollwöck,
  [arXiv:0901.2342](https://arxiv.org/abs/0901.2342)) before any
  windowing: an autoregressive model of order `lp_order` (default
  `None`, auto-picked as $\min(20,\max(4,\lfloor n_t/10\rfloor))$ so it
  stays safe even for a short simulation) is fit to the tail of $C(t)$
  and used to synthesize `lp_extend_factor` times as many additional
  samples, so the same real TDVP/TEBD simulation yields a sharper
  spectral function than windowing alone — the standard fix in the
  DMRG-dynamics literature for finite-simulated-time resolution loss,
  since it needs no additional entanglement growth. Pass `predict=False`
  to recover the exact pre-existing behavior. `lp_fit_start_fraction`
  (default 0.5) skips the corresponding leading fraction of $C(t)$ before
  fitting, and `lp_max_pole_radius` (default 1.0) reflects any fitted
  pole outside the unit circle back onto it, since the physical
  correlator only has weight for poles on or inside it and a slightly
  unstable fit would otherwise diverge under extrapolation.

`damping="exp"` paired with `predict=True` is the default combination
because it empirically gave the narrowest, best-centered peak among
every combination tried on a test system with a known exact gap:
pairing `predict=True` with `"gaussian"` instead came out *worse* (its
larger intrinsic FWHM at fixed $\delta$ partly cancels prediction's own
narrowing), and pairing it with `"parzen"` came out statistically tied
with plain `"exp"`, not better (`docs/td_dynamical_correlator_sharpening_plan.md`
has the numbers). Note that neither knob changes the *asymptotic*
algebraic decay of the tail far from any resonance — `"KPM"`'s far tail
still decays much faster, since its Jackson-kernel-damped Chebyshev
reconstruction has no Lorentzian/algebraic tail to begin with; these
knobs sharpen `"TD"`'s peaks and suppress its near-tail ringing, not
close that specific gap.

Both `damping` and `predict`/`lp_*` are also accepted by `submode="TDZ"`
and by the internal `S(k,\omega)` reduction (`sxt_to_skomega`), since all
three share the same windowing/FFT tail -- but their own defaults are
unchanged (`damping="exp"`, `predict=False`), since the empirical
comparison above was only run for `"TD"`.

**`submode="TDZ"` — complex-time evolution (Cao, Lu, Stoudenmire &
Parcollet, arXiv:2311.10909).** Real-time evolution (`"TD"` above) grows
entanglement, so the MPS bond dimension needed for a given accuracy
grows with the simulated time $T$. TDZ instead evolves along a complex
time contour

$$z(t,\alpha_0)=\int_0^t e^{-i\alpha_0 f(t')}\,dt',\qquad f(t)=e^{-t\omega_0},\qquad \omega_0=2\pi/t_{\max}$$

Since $\mathrm{Im}\,z(t,\alpha_0)<0$ for $\alpha_0>0$, this progressively
damps high-energy content as it evolves, so the bond dimension needed
for a given accuracy grows far more slowly than under real-time
evolution alone (the original paper reports $\chi\sim20$–$30$ vs
$\chi\sim500$–$700$ for comparable accuracy on the Anderson impurity
model). The true real-time ($\alpha_0=0$) correlator is then recovered
order by order via a perturbative Taylor expansion in $\alpha_0$ around
the simulated contour,

$$C(t,0)\approx\phi^{(0)}(t,\alpha_0)+\sum_{n=1}^{n_{\max}}g^{(n)}(t,\alpha_0)$$

where $\phi^{(n)}(t,\alpha_0)=\langle H^n B|\mathrm{GS}\rangle\cdot
|\psi(t,\alpha_0)\rangle$ (precomputed once per $n$, reused as a fixed
overlap target at every time step) and $g^{(n)}$ are explicit
combinatorial expressions in $\phi^{(1..n)}$ and the pure contour
integrals $J^{(n)}(t,\alpha_0)=-i\,\partial^n_{\alpha_0}z(t,\alpha_0)$
(see the paper's Appendix B; this implementation hardcodes $n\le4$,
which the paper finds always suffices for $\alpha_0\lesssim0.3$). The
reconstructed $C(t,0)$ is then windowed/Fourier-transformed exactly as
in `"TD"`.

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="TDZ",
                                      name=(sc.Sz[0], sc.Sz[0]),
                                      alpha0=0.1, n_max=4, dt=0.05)
```

`alpha0` is the contour angle parameter (larger reduces the bond
dimension needed further, but requires a larger `n_max` to reconstruct
the real axis accurately); `n_max` (≤4) is the reconstruction order;
`dt`/`tmax`/`nt` set the underlying time step/duration exactly as in
`"TD"`. Uses two-site TDVP when available (`itensor_version` 3 or
`"python"`, `tevol_method="TDVP"` or `"AUTO"` — `"AUTO"` reduces to plain
`"TDVP"` here rather than trying `"TEBD"` first, since a complex/imaginary
time step has no TEBD counterpart on any backend to try, see §7's own
note — the paper's own setup), one-site TDVP with global subspace
expansion (`tevol_method="TDVP_GSE"`, see §7 — same `itensor_version`
support as `"TDVP"`), or falls back to the MPO-Taylor propagator otherwise
(`tevol_method="MPO"` or `"TEBD"`, or `itensor_version=2`, which has no
TDVP) — the same TDVP-vs-Taylor choice `"TD"` already makes. Current
scope: only the "greater" branch of the correlator is computed (the same
simplification `"TD"` itself already makes), so this is best used the
same way as `"TD"`: high-resolution work in a narrow frequency window,
now reachable at a lower bond-dimension cost for a given simulated time.

**`submode="EX"` — exact diagonalization in a truncated DMRG subspace.**
Builds $A$, $B$, $H$ explicitly in the subspace spanned by the lowest
`nex` DMRG excited states, then evaluates the exact Lehmann sum in that
subspace,

$$G_{AB}(\omega)=\sum_{n=1}^{n_{ex}}\frac{\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle}{\omega-E_n+i\delta}$$

i.e.\ a small, explicit sum over poles at the computed excited-state
energies $E_n$, each with residue given by the transition matrix
elements. Cheap and exact *within* the truncated subspace; only as good
as how many/which excited states were computed. The `nex` excited-state
MPS are not assumed to be orthonormal (`dcex.py` builds their own overlap
matrix and solves a generalized eigenvalue problem $Hc=eSc$ rather than
assuming $S=1$), which is important in practice since Gram-Schmidt over
bond-truncated MPS is itself only approximate.

**`submode="maxent"` — maximum-entropy reconstruction.** Reconstructs a
positive-definite spectral function from a finite set of moments
$\langle(H-E_0)^k\rangle$ using a maximum-entropy method
(`distribution.get_distribution_maxent`), rather than a Chebyshev
expansion — useful when positivity of the reconstructed $S(\omega)$
matters more than matching KPM's polynomial-expansion artifacts.

**`submode="KPM"` for non-Hermitian Hamiltonians.** When $H\neq H^\dagger$
(§4), `submode="KPM"` automatically routes to a different algorithm, a
port of the non-Hermitian Kernel Polynomial Method (NH-KPM) of
[NHKPM.jl](https://github.com/GUANGZECHEN/NHKPM.jl) ([Phys. Rev. Lett.
130, 100401](https://doi.org/10.1103/PhysRevLett.130.100401)):

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                      name=(sc.Sz[0], sc.Sz[0]), E_max=10)
```

The ground state is now the biorthogonal pair $|\psi_R\rangle,|\psi_L\rangle$
from NH-DMRG (§4) rather than a single self-dual $|\mathrm{GS}\rangle$,
and the correlator computed is
$\langle\psi_L|A(z)\,B|\psi_R\rangle$. Because the spectrum of $H$ is
complex, the ordinary Chebyshev recursion (which needs a *real* rescaled
spectrum in $[-1,1]$) does not apply; NH-KPM instead expands the
frequency-shifted operator $\tilde H(z)=(z\,\mathbb{1}-H)/E_{\max}$ using
a *coupled* forward/adjoint recursion built from both $\tilde H(z)$ and
$\tilde H(z)^\dagger$ (see `src/dmrgpy/algebra/kpm.py`'s
`get_mu_n_nh`/`spec_from_moments_nh`, ported line-for-line from the
reference's `get_vn_NH`/`get_spec_kpm_NH`). The key practical consequence
is that, unlike the Hermitian case, the moments depend on $z$ itself, so
they are recomputed from scratch at every requested frequency rather than
amortized once over the whole spectrum — this mirrors the reference
algorithm's own cost profile, and is why NH-KPM is noticeably more
expensive per frequency point than the Hermitian `"KPM"` path. `E_max`
(an upper bound on the spectral radius of $H$) must be supplied
explicitly: unlike the Hermitian case's variational band-edge estimate
(see above), there is no automatic estimator yet for a non-Hermitian
spectral bound. Implemented so far for the ED backend, `itensor_version=3`,
and `itensor_version="python"`; `itensor_version=2` raises
`NotImplementedError`. See `examples/non_hermitian/nhkpm_v3_VS_ED`
(ED vs `itensor_version=3`, machine-precision agreement on a small
interacting fermionic chain with a staggered imaginary potential) and
`examples/non_hermitian/nhkpm_python_VS_v3_timing` (`itensor_version=3`
vs `"python"` on a non-uniform hopping/non-uniform imaginary-onsite-energy
chain, same machine-precision agreement — the pure-Python backend runs
roughly 2x slower than v3 for this workload, since NH-KPM's
per-frequency moment recursion is far more matvec-heavy than the
Hermitian KPM path).

**The `dex` cutoff — which states count as "the" ground state
(`mode="ED", submode="ED"`).** The $T=0$ ED sum above starts from the
ground state, but a degenerate ground manifold has no single one. The ED
path therefore takes every eigenstate whose excitation energy lies below
a tolerance `dex` (default $10^{-5}$) as the initial manifold and averages
the correlator over it with **equal weights** $1/n_{\rm ex}$. That is the
right thing to do for a genuinely degenerate manifold, and only for that
case: `dex` is a hard step, so a level at $0.99\,\texttt{dex}$ contributes
in full while one at $1.01\,\texttt{dex}$ contributes nothing at all.

Consequently `dex` must be chosen relative to the *splittings being
resolved*, not to the broadening $\delta$ or the frequency grid. The
pathological case is a parameter sweep — e.g.\ a Zeeman field whose
splitting grows through `dex` — where the result jumps discontinuously as
each level crosses the threshold, and the sweep partially averages over
the very splitting it was meant to resolve. No single `dex` works at both
ends of such a sweep: it must be large enough to keep the whole manifold
at zero field and small enough to exclude the split states at finite
field. DMRGPY emits a `RuntimeWarning` whenever an eigenvalue lies within
a factor of three of `dex` on either side, i.e.\ exactly when the answer
depends on where the cutoff was placed.

The clean way out of that regime is to stop using a cutoff at all and
weight the manifold physically instead: pass a small `T` (see §9), which
switches to the Boltzmann-weighted Lehmann sum. It reproduces the
equal-weight degenerate average smoothly as $T\to0$ — at exact degeneracy
the two agree to machine precision — while varying continuously, not in
steps, once the manifold splits. Use `dex` for a one-off spectrum at a
fixed parameter, and `T` for anything swept.
See `examples/kondo/atom_iets_orbital_resolved`.

**Choosing a method:** KPM (default) for a first look at the full
spectrum; CVM or TD when you need high resolution in a specific,
narrow frequency window; TDZ instead of TD when that window also needs
long simulated times/low frequencies, where TD's real-time bond-dimension
growth becomes limiting; EX when a handful of excited states already
capture the physics (e.g. a small gapped system); maxent when you want a
guaranteed-positive reconstruction from limited moment data (e.g.\
combined with finite-temperature ED, see §9). For a non-Hermitian $H$,
`"KPM"` is currently the only submode with a genuine biorthogonal
implementation (see above); the other submodes fall back to a
correction-vector method that assumes $A^\dagger=B$.

## 7. Real-time dynamics: quenches

Beyond frequency-domain correlators, DMRGPY directly simulates real-time
unitary evolution $|\psi(t)\rangle=e^{-iHt}|\psi(0)\rangle$ and measures
an observable along the way:

```python
from dmrgpy import timedependent
wf0 = sc.get_gs()                          # prepare |GS> of H0
sc.set_hamiltonian(h1)                     # quench to a different H
(ts, sz) = timedependent.evolve_and_measure(sc, operator=sc.Sz[0],
                                             nt=200, dt=1e-2, wf=wf0)
```

giving $\langle\psi(t)|S^z_0|\psi(t)\rangle$ as a function of time after
a Hamiltonian quench $H_0\to H_1$ — the standard non-equilibrium quantum
quench setup: prepare the ground state of one Hamiltonian (e.g.\ a
symmetry-broken/Néel-ordered $H_0$), then let it evolve under a
different $H_1$ (e.g.\ the isotropic Heisenberg point) and watch
observables relax, oscillate, or thermalize. `evolution_ABA` similarly
lets you apply an operator $A$ as an instantaneous local quench (e.g.\
flip a spin, or add a particle) before evolving and measuring $B(t)$.
`timeevolution.imaginary_exponential` computes autocorrelation functions
directly, $\langle\psi_0|e^{iHt}|\psi_0\rangle$, without a separate
measurement operator.

**Choosing the propagator: `sc.tevol_method`.** Five options:

- `"TDVP"` (the default) — two-site TDVP, which grows the MPS bond
  dimension via SVD the same way ground-state DMRG does. Used whenever
  `itensor_version` is `3`, `"python"` or `"julia_live"`;
  `itensor_version=2` falls back to `"MPO"` (below) even with this
  default.
- `"TDVP_GSE"` (`itensor_version` 3, `"python"` or `"julia_live"`, same
  support as plain `"TDVP"` above) — one-site TDVP preceded, for the first
  `sc.tdvp_gse_sweeps` steps (default 3), by a *global subspace
  expansion* step: a Krylov subspace $\{\psi,H\psi,H^2\psi,\dots\}$ of
  dimension `sc.tdvp_gse_krylov_order` (default 3) is used to enlarge the
  MPS's bond dimension *without changing the state it represents*, using
  a cutoff `sc.tdvp_gse_cutoff` (default $10^{-8}$) — the scheme of Yang
  & White, [arXiv:2005.06104](https://arxiv.org/abs/2005.06104)/Phys.
  Rev. B 102, 094315 (2020). Most useful when the starting state's bond
  dimension is small (e.g.\ a product-state quench) and one-site TDVP
  alone (which conserves bond dimension exactly) wouldn't be able to grow
  into the entanglement the subsequent evolution generates. On
  `"julia_live"` neither half is a dmrgpy port -- ITensorMPS.jl ships the
  expansion itself as `expand(psi,H; alg="global_krylov")` (citing the
  same paper) and its `tdvp` takes `nsite=1` directly, so
  `mpsjulialive/tdvp.jl` only wires the two together; that route also
  handles a bond-dimension-1 (product-state) start fine, unlike
  `itensor_version=3` (see
  `examples/time_evolution/tdvp_gse_VS_ED_time_evolution`'s own note on
  that edge case).
- `"TEBD"` (`itensor_version` `3`, `"python"`, or `"julia_live"`, and only
  for a strictly nearest-neighbor Hamiltonian — any term touching 3 or
  more distinct sites raises `NotImplementedError` (`"python"`), a
  catchable `RuntimeError` (`3`), or a `juliacall.JuliaError`
  (`"julia_live"`)) — the standard 2nd-order-Trotter, even/odd-bond
  ("brick-wall") algorithm: $e^{-iH\,dt}\approx
  e^{-iH_{\rm odd}\,dt/2}\,e^{-iH_{\rm even}\,dt}\,e^{-iH_{\rm odd}\,dt/2}$,
  where $H_{\rm odd}$/$H_{\rm even}$ sum the bond Hamiltonians on
  odd-/even-indexed bonds (each internally commuting, since every bond in
  one group acts on disjoint sites). Onsite terms are split half-and-half
  onto each site's neighboring bond(s) (full weight at a chain boundary).
  Unlike the TDVP variants, every bond's evolution gate
  $e^{-i\tau h_{\rm bond}}$ is the exact exponential of the *bare* local
  2-site Hamiltonian, built once up front and reused unchanged for every
  time step — no per-step Krylov/Lanczos work at all — so `"TEBD"` is
  typically cheaper per step than TDVP whenever the Hamiltonian qualifies.
  `itensor_version=3` builds the gate via ITensor's own `BondGate`
  primitive (`mpscpp3/tebd.h`); `"python"` exponentiates the bare 2-site
  matrix directly (`pyitensor/tebd.py`'s `TEBDEvolver`); `"julia_live"`
  builds the gate as an ITensor outer product of per-site operators and
  exponentiates via ITensors.jl's own tensor `exp()`
  (`mpsjulialive/tebd.jl`) — all three agree with each other (and with
  exact diagonalization) on a fermionic hopping+onsite benchmark. Unlike
  the other two backends, the Julia port resolves the Jordan-Wigner
  string itself, from scratch, off the *raw* `"C"/"Cdag"` term list
  (the same one the MPO-based methods there already serialize), rather
  than `MultiOperator.to_terms()`'s Jordan-Wigner-*predressed*
  `"A"/"Adag"/"F"` form the other two backends consume — real
  ITensors.jl's builtin `"Fermion"` site type only defines
  `op("C",..)`/`op("Cdag",..)`/`op("F",..)`, not dmrgpy's own `"A"/"Adag"`
  names.
- `"MPO"` — a hand-rolled 2nd-order Taylor expansion of $e^{-iH\,dt}$
  applied as an MPO each step; the only option on `itensor_version=2`
  (which has no TDVP or TEBD at all), and available (if slower/less
  accurate for a given bond dimension) everywhere else too.
- `"AUTO"` (`itensor_version` `3`, `"python"`, or `"julia_live"`, same
  support as `"TEBD"`) — tries `"TEBD"` first and transparently retries
  as plain `"TDVP"` if the Hamiltonian turns out not to be strictly
  nearest-neighbor, instead of raising. This exists because "is my
  Hamiltonian nearest-neighbor" is often a property of *which terms a
  caller adds*, not something pinned down once at chain-construction
  time — a script that starts nearest-neighbor and later grows a
  longer-range term (or is reused across several models) would otherwise
  need its own `tevol_method` bookkeeping to keep getting the cheaper
  `"TEBD"` path whenever it still applies. `"TEBD"` itself deliberately
  stays a hard opt-in that raises rather than silently falling back
  (see the note at the end of this list) — `"AUTO"` is the explicit way
  to ask for the fallback instead, so a caller who really meant `"TEBD"`
  and got a surprising long-range term still finds out. The check costs
  at most one discarded MPO build on the fallback path: both
  `bond_hamiltonians()` implementations (`pyitensor/tebd.py`,
  `mpscpp3/tebd.h`, `mpsjulialive/tebd.jl`) reject a non-nearest-neighbor
  term before touching the wavefunction at all, so retrying with
  `"TDVP"` never discards a partial time evolution. Not the default,
  precisely because it isn't a free win: whether TEBD applies becomes
  something you find out from behavior (a quietly different integrator,
  and thus different truncation-error characteristics) rather than from
  reading `sc.tevol_method` — worth it once a script's Hamiltonian shape
  is genuinely dynamic, not worth it as a blanket default for scripts
  that already know their Hamiltonian is nearest-neighbor and can just
  say `"TEBD"`.

`"julia_live"` implements `"TDVP"`, `"TDVP_GSE"`, `"TEBD"`, and `"AUTO"`;
the legacy `"MPO"` path raises `NotImplementedError` there rather than
silently running plain TDVP instead, so a backend-comparison script can't
quietly end up comparing different integrators.

```python
sc.setup_cpp(version=3)    # or sc.setup_python(), or sc.setup_julia()
sc.tevol_method = "TEBD"
```

```python
sc.tevol_method = "TDVP_GSE"
sc.tdvp_gse_sweeps = 3
sc.tdvp_gse_krylov_order = 3
sc.tdvp_gse_cutoff = 1e-8
```

## 8. Density of states

**Many-body density of states.** The full many-body spectral density

$$\rho(E)=\sum_n\delta(E-E_n)$$

(sum over *all* eigenstates of $H$, not resolved by any operator) has no
ready-made helper -- the closest available tool is
`get_distribution`/`kpmdmrg.general_kpm` (§6), but that computes a
ground-state expectation value $\langle\mathrm{gs}|B\,\delta(X)\,A|\mathrm{gs}\rangle$,
not a full-spectrum trace, so it is not a drop-in replacement for this
quantity.

**Single-particle local density of states.** For fermionic chains, the
physically distinct quantity usually meant by "DOS" (e.g.\ as measured
by scanning tunneling spectroscopy) is the local single-particle
spectral function at site $i$,

$$A_i(\omega)=-\frac1\pi\,\mathrm{Im}\,G^R_{ii}(\omega),\qquad G^R_{ii}(\omega)\ \text{built from}\ \langle c_i(t)c_i^\dagger(0)\rangle\ (\omega>0)\ \text{and}\ \langle c_i^\dagger(t)c_i(0)\rangle\ (\omega<0)$$

i.e.\ the particle-addition and particle-removal (electron/hole)
branches of the local Green's function concatenated across $\omega=0$.
There is no ready-made helper for this (the module this section used to
point to, `fermiondos.py`, relied on a `get_dynamical_correlator` calling
convention that no longer exists and was removed) -- build the two
branches directly with `get_dynamical_correlator` instead, using
`(fc.C[i], fc.Cdag[i])` for the particle-removal branch and
`(fc.Cdag[i], fc.C[i])` for the particle-addition branch, then
concatenate.

## 9. Finite temperature

Thermal (mixed-state, finite-$T$) expectation values are obtained via
**purification**: each physical site is paired with an ancilla site, and
the maximally-entangled state of every physical-ancilla pair (e.g.\ the
singlet-forming Heisenberg coupling $\mathbf S_i^{\rm phys}\cdot\mathbf
S_i^{\rm anc}$ as the "Hamiltonian" preparing that state) is exactly the
purification of the infinite-temperature ($T=\infty$) physical density
matrix $\rho\propto\mathbb 1$. Imaginary-time evolving this purified
state under the *physical* Hamiltonian $H$,

$$|\Psi(\beta)\rangle=e^{-\beta H/2}|\Psi(0)\rangle,\qquad \beta=1/T$$

(stepped via repeated small first-order updates
$|\Psi(\beta+\Delta\beta)\rangle\propto(1-\Delta\beta\,H)|\Psi(\beta)\rangle$,
renormalized at each step) and tracing out the ancillas gives exactly the
thermal (Gibbs) density matrix of the physical chain,
$\rho(T)=\mathrm{Tr}_{\rm anc}|\Psi(\beta)\rangle\langle\Psi(\beta)|\propto e^{-\beta H}$:

```python
from dmrgpy import thermal
tc = thermal.Thermal_Spin_Chain(spins, T=0.1)
```

`T=0` recovers ordinary ground-state DMRG. For small systems (Hilbert
space dimension $\lesssim2000$, `algebra.maxsize`), `get_correlation_matrix(T=...)`
(ED only, `entropytk/correlationentropy.py`'s `get_correlation_matrix_finiteT`)
instead computes the exact thermal average directly by brute-force ED,
$\rho=\sum_nP_n|n\rangle\langle n|$, $P_n=Z^{-1}e^{-(E_n-E_0)/T}$ — a
useful cross-check of the purification approach, following the same
excited-states pattern as `vev(mode="ED", T=...)` below (full spectrum
whenever it's small enough to diagonalize exactly, an explicit `n=` and
a `RuntimeError` on inadequate truncation otherwise). Automatic default
operators (no `operators=` kwarg) are only defined for fermionic chains
today (`self.C`); on a `Spin_Chain`/`Thermal_Spin_Chain` pass an explicit
`operators=[...]` (e.g.\ `[sc.Sz[i] for i in range(n)]`).

**`vev(O, mode="ED", T=...)`** is a second, independent way to get the
same brute-force ED thermal average
$\langle O\rangle=\sum_nP_n\langle n|O|n\rangle$, $P_n=Z^{-1}e^{-(E_n-E_0)/T}$,
for a single operator rather than the whole correlation matrix
(`edtk/edchain.py`'s `EDchain.vev` → `vevtk/thermalvev.py`'s
`thermal_vev_ex`). The sum only ever runs over a finite set of exact
eigenstates obtained from `algebra.lowest_states`; whenever the full
Hilbert space is small enough to diagonalize exactly anyway (the same
`maxsize=2000`-state threshold `algebra.lowest_states` itself uses),
`thermal_vev_ex` defaults to using the *entire* spectrum, so the thermal
average is exact for any `T`. Above that size only a sparse, truncated
set of the lowest states is available (an explicit `n=` kwarg, forwarded
from `vev(...)`, controls how many); `thermal_vev_ex` then checks the
Boltzmann weight of the highest included state and raises `RuntimeError`
rather than silently returning a wrong average if that state still
carries non-negligible weight at the requested `T` — pass a larger `n=`
in that case.

**`metts_vev(O, T, ...)`** is a third, independent finite-temperature
method: METTS (Minimally Entangled Typical Thermal States, E.M.
Stoudenmire and S.R. White, *New J. Phys.* **12**, 055026 (2010),
arXiv:1002.1305), implemented for `itensor_version="python"`
(`pyitensor/metts.py`), `itensor_version=3`
(`mpscpp3/chain_session.h`'s `Chain::metts_vev`, a direct port of the
same algorithm onto real ITensor v3, not an independent
reimplementation), and `itensor_version="julia_live"`
(`mpsjulialive/metts.jl`'s `metts_vev`, a value-level port of the same
algorithm reusing `tdvp.jl`'s `tdvp_step` unchanged for the imaginary-time
evolution) — not for `itensor_version=2` (mpscpp2 has no TDVP module at
all, and METTS needs imaginary-time TDVP). On `itensor_version=3`, a
chain shorter than 3 sites raises `NotImplementedError` rather than
crashing: ITensor v3's two-site TDVP hits the same "LocalOp is default
constructed" abort as its two-site `dmrg()` for such short chains (see
§Architecture's note on that `mpscpp3` bug). Rather than
purification's single, growing-entanglement wavefunction or ED's exact
sum over eigenstates, METTS samples a Markov chain of *unentangled*
classical product states (CPS) $|i\rangle$: each is imaginary-time
evolved by half the inverse temperature,
$|\phi_i\rangle=e^{-\beta H/2}|i\rangle/\lVert e^{-\beta H/2}|i\rangle\rVert$
(via imaginary-time TDVP — `tdvp_step` with a purely real, rather than
imaginary, effective time step), then collapsed back down to a new CPS
$|i'\rangle$ by sampling a definite outcome at each site in turn, with
probability $|\langle i'|\phi_i\rangle|^2$, from a single left-to-right
sweep (no need to ever form the full $2^N$ amplitude vector — a
canonicalized MPS's marginal at one site, conditioned on those already
sampled to its left, is already diagonal). A plain (unweighted) sample
average of $\langle\phi_i|O|\phi_i\rangle$ over the resulting chain then
converges to the true thermal average, because this specific Markov
chain's stationary distribution already carries the correct Boltzmann
weight — no importance reweighting needed. The collapse basis is
alternated between two choices from one sample to the next
(`basis_ops=("Sz","Sx")` by default) since collapsing repeatedly in the
same basis can trap the chain in one symmetry sector for many steps,
inflating autocorrelation (the paper's own finding):

```python
mean, stderr = sc.metts_vev(sc.Sz[0], T, nsamples=300, nwarmup=30,
                             dbeta_half_step=0.05, basis_ops=("Sz","Sx"))
```

`nwarmup` discards that many initial Markov-chain steps before
averaging (equilibration); `dbeta_half_step` sets the imaginary-time
TDVP step size for the $e^{-\beta H/2}$ evolution (split into
$\lceil(\beta/2)/\texttt{dbeta\_half\_step}\rceil$ equal steps). Being a
Monte Carlo method, `metts_vev` returns `(mean, stderr)` rather than an
exact value — `stderr` is a naive i.i.d. estimate and, since consecutive
METTS samples are Markov-correlated, is likely optimistic unless
`dbeta_half_step`/`nwarmup` are generous enough that samples actually
decorrelate. Its main advantage over purification is that the sampled
states stay unentangled classical product states between imaginary-time
evolutions (no ancilla doubling, and bond dimension never has to carry
the entanglement of a single ever-more-thermalized wavefunction) — see
`examples/finite_temperature/metts_VS_exact` for a cross-check against
the exact ED thermal average on a small Heisenberg chain.

`O` can also be a list/tuple of operators, measured together on one
shared sampled Markov chain instead of resampling from scratch for each:

```python
results = sc.metts_vev([sc.Sz[0], sc.Sz[1], H], T, nsamples=300, nwarmup=30)
# results == [(mean_Sz0, stderr_Sz0), (mean_Sz1, stderr_Sz1), (mean_H, stderr_H)]
```

Since the `nwarmup+nsamples` imaginary-time evolutions dominate the cost
of `metts_vev`, not the handful of extra `<phi|O_k|phi>` measurements per
sample, batching several observables this way is far cheaper than calling
`metts_vev` once per operator.

For `itensor_version="python"`, an `njobs` keyword (default 1) runs
`njobs` independent METTS Markov chains in parallel worker processes and
pools their statistics, instead of one longer sequential chain —
`nsamples` is split as evenly as possible across them, each chain gets
its own `nwarmup` equilibration and an independently-seeded RNG, and
their per-chain `(mean, stderr, count)` triples are combined into the
same `(mean, stderr)` a single pooled run over all the raw samples would
give (no raw samples need to cross the process boundary, so this is
exact, not an approximation). Measured directly on a 6-site chain
(`nsamples=200`, `nwarmup=20`): wall time went from 34s at `njobs=1`
down to 14s at `njobs=8` — real but sub-linear speedup, since every
extra chain repeats the full `nwarmup` rather than sharing it. This is
the effective optimization for this backend: its per-sample cost is
dominated by `pyitensor`'s own generic-tensor-engine Python overhead
(confirmed by profiling), not shared BLAS work, so splitting across OS
processes helps where the `kernels.py` JAX/numba contraction-kernel
route does not (see that module's docstring, and
`pyitensor/metts.py`'s own comment, on why numba is measurably *slower*
for METTS specifically — each sample restarts from a fresh
bond-dimension-1 product state, so a run exercises many distinct
contraction shapes rather than reusing one, and the fixed per-shape
compile tax never amortizes). `njobs` is not available for
`itensor_version=3` or `"julia_live"`: each is a single live in-process
session (a C++ object, or a live Julia process) with no per-worker copy a
process pool could hand out; requesting `njobs>1` on either raises rather
than silently falling back to `njobs=1`.

**`metts_dynamical_correlator(name, T, ...)`** extends METTS from static
expectation values to real-time finite-temperature *dynamical*
correlators
$\mathcal{C}_{AB}(t)=\langle A(t)B\rangle_T=\langle e^{iHt}Ae^{-iHt}B\rangle_T$
(Z. Wang, P. McClarty, D. Dankova, A. Honecker and A. Wietek,
"Spectroscopy and complex-time correlations using minimally entangled
typical thermal states", arXiv:2405.18484, Sec. II, "Dynamical METTS
algorithm"), implemented for `itensor_version="python"`, `3`, and
`"julia_live"` (`mpsjulialive/metts.jl`'s `metts_dynamical_correlator`, a
value-level port reusing the same `tdvp_step` `metts_vev` already uses,
now with a purely real time step for the real-time evolution of
$|v_i(t)\rangle$/$|w_i(t)\rangle$ instead of the purely imaginary one used
for sampling). For every METTS sample
$|\psi_i\rangle$ produced by the exact same Markov chain `metts_vev`
already samples (imaginary-time evolution + sequential-sampling
collapse), define $|v_i(0)\rangle=B|\psi_i\rangle$,
$|w_i(0)\rangle=|\psi_i\rangle$, real-time evolve both independently
under $H$ (two-site TDVP), and measure
$\mathcal{C}^i(t)=\langle w_i(t)|A|v_i(t)\rangle$ at each requested time
step. A plain (unweighted) sample average of $\mathcal{C}^i(t)$ over
retained samples converges to $\mathcal{C}_{AB}(t)$, for the same reason
`metts_vev`'s own plain average converges to the thermal average — no
importance reweighting needed:

```python
ts, means, stderrs = sc.metts_dynamical_correlator(
    (sc.Sz[0], sc.Sz[0]), T, nt=100, dt=0.1, nsamples=200, nwarmup=30,
    dbeta_half_step=0.05, basis_ops=("Sz","Sx"))
```

`name` follows the same `(A,B)` convention as
`get_dynamical_correlator`'s own `name=` (a string like `"ZZ"`, or an
explicit `(MultiOperator,MultiOperator)` tuple/list) — `A`,`B` are used
exactly as given, with no dagger applied to either, matching
`get_dynamical_correlator(mode="ED", submode="ED", T=...)`'s own
convention (see below) so the two are directly comparable. `nt`,`dt` set
the (uniformly spaced) real-time measurement grid
$t=0,\Delta t,\dots,(n_t-1)\Delta t$; `nsamples`, `nwarmup`,
`dbeta_half_step`, `basis_ops`, `seed`, `niter`, `njobs` all mean exactly
what they mean for `metts_vev` (same shared Markov chain, same caveats on
`stderrs` being a Markov-correlated, likely-optimistic naive estimate).
`tdvp_niter` separately bounds the Krylov iterations used for the
*real-time* evolution of $|v_i(t)\rangle$/$|w_i(t)\rangle$ specifically
(default 50), independent of `niter`'s bound on the imaginary-time
sampling step (default 30) — the two generally warrant different
settings since $|v_i(t)\rangle$/$|w_i(t)\rangle$ typically become more
entangled over the course of real-time evolution than the METTS samples
$|\psi_i\rangle$ themselves ever do (for `itensor_version="julia_live"`,
`niter`/`tdvp_niter` are accepted for signature parity but silently
ignored, same as `metts_vev`'s own `niter` — ITensorMPS.jl's `tdvp()`
manages its own internal Krylov dimension with no exposed per-step
iteration-count knob). No Fourier transform/windowing is
performed internally: `metts_dynamical_correlator` returns the raw
time-domain samples/statistics, matching `evolution_DC`'s own `(ts,cs)`
convention — apply a window (e.g. a Hann window, as the reference paper
recommends) and FFT separately if a frequency-domain spectral function is
wanted.

For a direct, exact reference to validate against, `get_dynamical_correlator(
mode="ED", submode="ED", T=..., name=...)` computes the finite-temperature
dynamical correlator's spectral function via a full Boltzmann-weighted
Lehmann sum over *every* ED eigenstate
$\mathcal{C}_{AB}(\omega)=\frac{1}{\mathcal{Z}}\sum_{n,m}e^{-\beta E_n}\langle n|A|m\rangle\langle m|B|n\rangle\,[\text{kernel at }\omega=E_m-E_n]$
— the finite-$T$ generalization of the existing T=0 `submode="ED"`
near-degenerate-ground-state sum to every eigenstate weighted by its
exact Boltzmann factor (exact, since it starts from a full dense
diagonalization, unlike `thermal_vev_ex`'s own partial-diagonalization
truncation-safety check, which this doesn't need). See
`examples/finite_temperature/dynamical_metts_VS_ED` (`itensor_version` in
`("python", 3)`) and `examples/finite_temperature/dynamical_metts_julia_VS_ED`
(`itensor_version="julia_live"`) for a cross-check of
`metts_dynamical_correlator` against this exact ED reference, evaluated
directly in the time domain, on a small Heisenberg chain.

## 10. Topological invariants

**Many-body Berry phase.** For a ground state that depends on an
adiabatic parameter $k$ threaded around a closed loop (e.g.\ inserted
flux), the discretized Berry phase

$$\gamma=\arg\prod_{k}\langle\psi(k)|\psi(k+\delta k)\rangle$$

quantized (typically in units of $\pi$) signals a topologically
nontrivial ground state / obstruction to adiabatic continuity — the
many-body generalization of a Zak/Berry phase, computed by running DMRG
at a discrete set of parameter points and chaining ground-state
overlaps.

**Single-particle Berry phase / Wilson loop.** For a translationally
invariant single-particle Hamiltonian $H(k)$ with a set of occupied
bands, the non-Abelian Wilson loop over one Brillouin-zone circuit is

$$W=\det\Big[\textstyle\prod_k U(k,k+\delta k)\Big],\qquad U_{mn}(k,k+\delta k)=\langle u_m(k)|u_n(k+\delta k)\rangle$$

with $m,n$ running over occupied bands; $\gamma=\arg W$ is the
polarization/Zak phase of that band manifold (and, combined with a scan
over a second momentum direction, the ingredient for a Chern number).

## 11. Mean-field decoupling

For a spin Hamiltonian with exchange couplings $J_{ij}$, a
self-consistent mean-field (Weiss-field) decoupling replaces the
two-body exchange term with a one-body field,

$$H_{\rm MF}=\sum_i\mathbf h_i\cdot\mathbf S_i,\qquad \mathbf h_i=(1-p)\sum_jJ_{ij}\langle\mathbf S_j\rangle$$

iterated to self-consistency: solve for $\langle\mathbf S_i\rangle$ under
$H_{\rm MF}$ (by DMRG), rebuild $\mathbf h_i$ from the new expectation
values, mix old/new fields, and repeat until $\max_i|\Delta\langle\mathbf
S_i\rangle|$ falls below a tolerance. The parameter $p\in[0,1]$
interpolates between pure mean-field theory ($p=0$) and keeping a
fraction $p$ of the original many-body exchange treated exactly
alongside the self-consistent field (useful for hybrid
mean-field-plus-fluctuations treatments of, e.g., magnetically ordered
phases where pure MF overestimates order).

```python
from dmrgpy import meanfield
meanfield.spinchain_meanfield(sc, p=0.0)
```

## 12. Fidelity susceptibility and quantum phase transitions

For a Hamiltonian $H(\lambda)=H_0+\lambda H_1$ depending on a tuning
parameter $\lambda$, the fidelity susceptibility measures how sharply
the ground state changes as $\lambda$ varies — and diverges at a quantum
phase transition, where the gap to the first excited state closes:

$$\chi(\lambda)=\sum_{n\neq0}\frac{|\langle 0|H_1|n\rangle|^2}{(E_n-E_0)^2+\delta^2}\qquad(\text{perturbative form, }\texttt{fmode="PT"})$$

`get_fidelity`'s **default** is not this perturbative form but a
non-perturbative estimator (`fmode="derivative"`) that computes $\chi$
directly from finite differences of the ground-state overlap matrix
between nearby $\lambda$ values, which also handles a (near-)degenerate
ground-state manifold via a smooth gauge choice:

```python
from dmrgpy import fidelity
chi = fidelity.get_fidelity(sc, h0, h1, lam, n=3) # fmode="derivative" (default)
chi_pt = fidelity.get_fidelity(sc, h0, h1, lam, n=3, fmode="PT") # perturbative form above
```

Scanning $\lambda$ and plotting $\chi(\lambda)$, a peak (sharpening and
diverging with system size) locates a quantum critical point without
needing to know its universality class in advance.

## 13. Ground-state degeneracy

Exact and near (e.g.\ symmetry-protected, or finite-size-split
topological) ground-state degeneracies are estimated by a narrow,
super-Gaussian-broadened level count around the ground-state energy
(`degeneracy.py`'s `gs_degeneracy_simple`/`eigenvalue_degeneracy`):

$$g(E_0)\approx\sum_i\exp\!\left[-\left(\frac{(E_i-E_0)^2}{\delta}\right)^2\right]$$

Note the quartic falloff in $(E_i-E_0)$ (not a plain Gaussian) — this
makes the window considerably narrower than $\delta$ would suggest for a
standard Gaussian, so `delta` needs to be picked accordingly (typically
larger than the target energy resolution) to count near-degenerate
levels rather than only exactly-degenerate ones.

summed over a growing number of low-lying computed eigenstates $E_i$
until the count converges — a value near an integer $g$ signals a
$g$-fold degenerate (or near-degenerate, at the working precision
$\delta$) ground-state manifold, as expected e.g.\ for the four
edge-state-split ground states of an open Haldane chain, or a
symmetry-broken ordered phase.

```python
from dmrgpy import degeneracy
g = sc.get_gs_degeneracy()
```

## 14. Reduced density matrices and operator distributions

**Reduced density matrix.** The single-site reduced density matrix
$\rho_i=\mathrm{Tr}_{j\neq i}|\mathrm{GS}\rangle\langle\mathrm{GS}|$ is
available directly (`sc.get_rdm(i=0)`), the basic object entanglement
entropies, local observables, and further post-processing (e.g.\ local
susceptibilities) are built from.

**Operator distributions.** More generally, the full probability
distribution of an arbitrary Hermitian operator $X$ (not just its
expectation value) in the ground state,

$$P(x)=\langle\mathrm{GS}|\,\delta(X-x)\,|\mathrm{GS}\rangle$$

is computed via the same KPM machinery as the dynamical correlators of
§6 (`sc.get_distribution`), or reconstructed from a finite set of raw
moments $\langle X^k\rangle$ via maximum entropy
(`get_distribution_maxent`). Useful for e.g.\ full counting statistics of
a conserved charge, or distinguishing a sharply peaked (well-defined
quantum number) ground state from a broadly spread one.

## 15. Post-processing tools

- **Analytic continuation** (`analyticcontinuation.py`): Padé
  continuation of a correlator known on the imaginary/complex-frequency
  axis (e.g.\ from a Matsubara-like or complex-shifted CVM calculation,
  §6's `submode="CVMimag"`) to the real frequency axis, where the
  physical spectral function lives.
- **Function fitting** (`functionfit.py`): a generic multi-start Powell
  minimizer used e.g.\ to fit the Calabrese-Cardy entropy formula in §5.
- **Finite-size extrapolation** (`extrapolate.py`): polynomial
  extrapolation in $1/L$ of a size-dependent quantity $y(L)$ toward the
  thermodynamic limit $L\to\infty$ — standard practice for extracting
  bulk quantities (energy density, order parameters, gaps) from finite
  DMRG chains.
- **Maximum-entropy reconstruction** (`reconstruct.py`): reconstructs a
  positive spectral function from a truncated moment expansion,
  underlying both the `"maxent"` dynamical-correlator submode and
  `get_distribution_maxent`.

## 16. Worked-example cookbook

**Central charge of a critical transverse-field Ising chain**
($H=\sum_iS_i^zS_{i+1}^z+\tfrac12\sum_iS_i^x$ at the critical field,
expected $c=\tfrac12$):

```python
sc.maxm = 200          # larger bond dimension: needed at criticality
wf = sc.get_gs()
print(wf.get_CFT_central_charge())
```

**Haldane gap of the spin-1 chain** (see §4):

```python
spins = ["S=1" for i in range(n)]
sc = spinchain.Spin_Chain(spins)
# ... build Heisenberg h ...
es = sc.get_excited(n=6)
print("Haldane gap:", es[4]-es[0])
```

**Charge gap of a Hubbard chain** (see §4):

```python
print("Single-particle gap:", fc.get_gap())
print("Charge (pair) gap:", fc.get_charge_gap(d=2))
```

**Fidelity susceptibility across the Ising transition** (see §12):

```python
h0 = 4*sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))   # Ising coupling
h1 = 2*sum(sc.Sx[i] for i in range(n))                # transverse field
for lam in lambdas:
    chis.append(fidelity.get_fidelity(sc, h0, h1, lam, n=3))
```

**Momentum- and frequency-resolved dynamical structure factor**
$S(q,\omega)$ (see §6), by combining the site-resolved KPM correlator
with a lattice Fourier transform:

```python
Sqw = {}
for i in range(n):
    for j in range(n):
        x, y = sc.get_dynamical_correlator(submode="KPM", name=(sc.Sz[i], sc.Sz[j]))
        Sqw[(i, j)] = (x, y)          # combine with sum_ij e^{iq(i-j)} S_ij(w) offline
```

## 17. STM/Kondo tunneling spectra (third-order perturbation theory)

`Spin_Chain.get_kondo_spectrum` computes the differential tunneling
conductance $dI/dV(eV)$ of an STM tip coupled to one site of a spin
chain, following the weak-coupling (Kondo-scattering) perturbation
theory of Ternes, *New J. Phys.* **17**, 063016 (2015),
[arXiv:1505.04430](https://arxiv.org/abs/1505.04430). This is a
different observable from the dynamical correlators of §6: instead of a
retarded Green's function of the spin system alone, it is the full
Fermi's-golden-rule tunneling current through tip+spin+sample, expanded
to third order in the tip-sample tunneling amplitude.

Two backends are available via `mode=`:

- `mode="ED"` (default): full exact diagonalization of the chain's
  Hamiltonian (every eigenstate is needed as a possible virtual
  intermediate state, not just the low-energy ones), independent of the
  chain's own `itensor_version`/mode setting. Works at any `T>=0`.
- `mode="DMRG"`: `itensor_version=3` throughout, never diagonalizing
  beyond the ground state — only `T=0` is supported. See "T=0 and the
  DMRG backend" below.

```python
sc = spinchain.Spin_Chain(["1/2"])
sc.set_hamiltonian(g*muB*B*sc.Sz[0])   # Zeeman-split S=1/2 impurity
eV, dIdV = sc.get_kondo_spectrum(eV_grid, site=0, Jrho_s=-0.05, U=0.25,
                                  T=1.0, order=3)
```

**Second order** (`order=2`) is the plain spin-flip/potential-scattering
Fermi golden rule result,

$$\frac{\partial I}{\partial V}(eV)\propto\sum_{i,f}p_i\Big[\tfrac12|\langle f|S_-|i\rangle|^2+\tfrac12|\langle f|S_+|i\rangle|^2+|\langle f|S_z|i\rangle|^2\Big]\,\Theta(eV-\epsilon_{if})+4U^2$$

summed over both tunneling directions, with $p_i$ the Boltzmann
occupation of eigenstate $i$ at temperature $T$ and $\Theta$ a
temperature-broadened step function. This reproduces the textbook
inelastic-tunneling spin-flip steps at $eV=\pm(\epsilon_f-\epsilon_i)$
(e.g. the Zeeman step of a single $S=1/2$ impurity).

**Third order** (`order=3`, the default) adds two corrections that
require summing over *all* eigenstates as virtual intermediate states
$m$ (energy conservation is not required for $m$): a Kondo term (a
Levi-Civita triple product of spin matrix elements $\langle i|S|f\rangle$,
$\langle f|S|m\rangle$, $\langle m|S|i\rangle$, weighted by a
temperature-broadened logarithmic function $F(eV-\epsilon_m,T)$ that
produces the characteristic zero-bias Kondo-like resonance, splitting
into two peaks under a Zeeman field), and — when `U!=0` — a
potential-scattering interference term responsible for a bias-asymmetric
lineshape.

**Both tunneling directions** are summed at every order, since the
measured $dI/dV$ is for the net current $I=I^{t\to s}-I^{s\to t}$. For
unpolarized tip and sample the matrix elements are direction independent,
so writing $g(eV)$ for the $t\to s$ expression the measured terms are

$$\frac{\partial I}{\partial V}\Big|_{\rm 2nd,\,Kondo}=g(eV)+g(-eV),\qquad\frac{\partial I}{\partial V}\Big|_{U\text{-}M}=g(eV)-g(-eV).$$

The relative sign is fixed by the paper's own worked $S=1/2$ example: the
purely Kondo-like processes contribute with the *same* sign in both
directions, while the potential-scattering ones change sign when the
tunneling direction is inverted. So the second-order and third-order
Kondo terms are **even** in bias (the zero-field Kondo resonance sits
exactly at $eV=0$), and the potential-interference term is **odd** (zero
at $eV=0$) — that term is the sole source of bias asymmetry in the
spectrum.

**Scope and known limitations**, worth reading before trusting specific
numbers:

- Only a single chain site couples to the tip (`site=`); the paper's own
  model allows several sites with independent tip couplings, not
  implemented here.
- The paper's own closed-form equations for two numerical building
  blocks — the temperature-broadened step $\Theta(x)$ and the
  temperature-broadened Kondo log function $F(\epsilon,T)$ — do not
  reproduce the physics the paper itself describes for them (checked
  directly: the printed $\Theta(x)$ diverges rather than saturating, and
  the printed $F$ closed form drops its own temperature broadening).
  `kondospectrumtk/stepfunctions.py` uses corrected forms instead,
  re-derived from the paper's own unambiguous defining integrals and
  verified against digitized values from the paper's own figures (see
  that module's docstring, and `examples/kondo/kondo_spectrum_VS_paper/`).
- The potential-interference term's general-spin closed form (`U!=0` in
  `order=3`) is an extrapolation from the paper's own worked $S=1/2$
  example (only that special case is spelled out in closed form in the
  paper). Its overall normalization is nonetheless pinned down, together
  with every other prefactor in the spectrum, by the paper's absolutely
  scaled Figs. 3b/3d — see below.
- **Normalization is checked against absolute figure values.** Figs. 3b
  and 3d of the paper are plotted in absolute units ($e^2T_0^2/h$), so
  their zero-bias peak heights (1.13 at $U=0$, 1.39 at $U=0.25$, for a
  single $S=1/2$ at $B=0$, $T=1\,$K, $J\rho_s=-0.05$, $\omega_0=20\,$meV)
  fix every prefactor at once. `examples/kondo/kondo_spectrum_VS_paper/`
  and `tests/test_kondo_spectrum.py` assert on them directly. In
  particular they fix the spin-average normalization ("SA factor"): the
  paper's spin-averaged transition matrix element (its own
  eq. for $|M_{if}|^2$) is twice the plain electron-spin trace, so the
  third-order Levi-Civita coefficient is $\mathrm{Im}[X]/2$ and the
  elastic potential channel is $4|U|^2$ rather than the bare $|U|^2$ of
  the printed $|\mathcal{M}^{(1)}|^2$ equation — see
  `kondospectrumtk/conductance.py`'s module docstring for the full
  bookkeeping.

**T=0 and the DMRG backend** (`mode="DMRG"`). At $T=0$ only the ground
state is thermally populated, which simplifies both terms enough to
avoid diagonalizing beyond the ground state entirely — the actual
motivation for supporting `mode="DMRG"` in the first place, since DMRG
cannot enumerate excited states the way `mode="ED"` does:

- **Second order** reduces to a $\Theta_0$-weighted (the exact,
  closed-form Heaviside limit of $\Theta$) cumulative integral of the
  ordinary $T=0$ dynamical structure factor
  $S_{\alpha\alpha}(\omega)=\sum_f|\langle f|S_\alpha|\mathrm{GS}\rangle|^2\delta(\omega-\epsilon_{f0})$,
  which `get_dynamical_correlator` already computes (`submode="KPM"` or
  `"CVM"`) without any excited-state enumeration
  (`kondospectrumtk/secondorder_dc.py`).
- **Third order** cannot be reduced to a two-operator correlator (it is
  a three-vertex object, with two *different* intermediate states each
  weighted by a different function, $\Theta_0$ and $F_0$) — instead it
  is built from a Heisenberg three-point function
  $G(t_2,\tau)=\langle\mathrm{GS}|S_l(t_2+\tau)S_k(t_2)S_j(0)|\mathrm{GS}\rangle$,
  obtained via real-time TDVP evolution with a "checkpoint-and-branch"
  construction (evolve $S_j|\mathrm{GS}\rangle$ forward/backward in
  $t_2$, apply $S_k$ at each checkpoint, evolve each branch further in
  $\tau$, overlap with a fixed $S_l|\mathrm{GS}\rangle$ reference at
  every step), then extracted via two closed-form time-domain kernels
  (derived by inverse-Fourier-transforming $\Theta_0$ and $F_0$) rather
  than by evaluating those functions pointwise on a discrete frequency
  grid, which does not converge robustly for this construction — see
  `kondospectrumtk/twotime.py`'s module docstring for the full
  derivation and the numerical pitfalls it was built to avoid ($\Theta_0$'s
  kernel is a Cauchy principal value, computed via an FFT-based Hilbert
  transform for machine-precision accuracy). This is the expensive part:
  cost scales with the number of $t_2$ checkpoints, each its own short
  TDVP trajectory (`kondospectrumtk/dmrgtwotime.py`).
- **Potential-interference term** (`U!=0`, part of `order=3`) is also
  supported: its own $T=0$ limit collapses the excited-state sum to a
  convolution of the *same* $T=0$ dynamical structure factor against the
  $F_0$ kernel instead of $\Theta_0$'s cumulative-sum weighting, so it
  reuses `get_dynamical_correlator` exactly like the second-order term
  above, needing no excited-state enumeration either
  (`kondospectrumtk/potentialdc.py`). Carries the same general-spin
  extrapolation caveat as `conductance.third_order_potential_dIdV` (see
  that function's docstring). All three DMRG-side routines sum over both
  tunneling directions exactly as their `mode="ED"` counterparts do —
  only the closed-form kernels depend on $eV$, so the $s\to t$ direction
  costs nothing beyond evaluating them a second time at $-eV$ (the
  expensive dynamical correlator / $G(t_2,\tau)$ is still built once).

Further `mode="DMRG"` limitations beyond the general ones above: the
second-order term's `es`
frequency-grid parameter has no safe default and must be supplied
explicitly (it needs to cover every relevant transition energy, which is
a property of the chain's spectrum, not of the `eV` sweep range — see
`second_order_dIdV_dc`'s docstring); the third-order term's `dt2`,
`n_t2_half`, `dtau`, `n_tau_half` time-grid parameters likewise have no
safe default and must be supplied explicitly (a grid fine/wide enough for
the default $\omega_0$/$\Gamma_0$ needs $\sim10^5$–$10^6$ $t_2$
checkpoints, each its own real TDVP trajectory — infeasible as a silent
default — while a small, fast default is wildly under-resolved and
returns a finite but silently wrong result instead of erroring, confirmed
directly; see `two_time_kondo_term_dmrg`'s docstring); chains need at
least 3 sites (a 1-site chain hits an internal ITensor v3 error unrelated
to this feature, building the Hamiltonian MPO).

`kondospectrumtk/dmrgtwotime.py` was written against this codebase's
existing, verified DMRG API and validated once a compiled ITensor v3
backend became available: $G(t_2,\tau)$ matches the ED reference to
$\sim10^{-9}$–$10^{-10}$ pointwise, and the swept third-order Kondo term
matches a grid-consistent ED reference to $\sim10^{-10}$
(`tests/test_kondo_spectrum_dmrgtwotime.py`, skipped automatically when
no compiled `itensor_version=3` backend is available). Getting there
surfaced three real bugs, none of which showed up in the ED-only testing
this module was originally written against — see that module's own
docstring for the details (in short: `tdvp_step` silently renormalizes
every step to unit norm, discarding $S_j|\mathrm{GS}\rangle$'s true
amplitude unless corrected for explicitly; a forward/backward
time-stepping bug meant "backward" checkpoints never actually reached
negative times; and a naive per-chunk trapezoidal integral is exactly 0
for the single-$t_2$-point chunks real-time evolution necessarily
produces, silently zeroing the entire term). The second-order term
(`submode="KPM"`) was spot-checked too, agreeing to within a few tens of
percent at thresholds, consistent with the expected
$\delta$-broadening/moment-truncation error on top of what the ED path
already has.

### Orbital-resolved IETS of a magnetic atom

`dmrgpy.atom` describes a *single* multi-orbital magnetic atom rather
than a chain: `generate_atom(orbs, tij, U, J, soc, B, Js, Ne)` builds a
`Spinful_Fermionic_Chain` whose "sites" are the atom's five $d$ orbitals,
with a crystal field `tij`, spin-orbit coupling `soc`, Hund's coupling
`J`, an $-U S^2-J L^2$ interaction, an external field `B` acting as
$2\mathbf{S}+\mathbf{L}$, and the electron count `Ne` fixed by a Lagrange
multiplier. Two inelastic-tunneling (IETS) observables are then available:

- `get_spinflip(fc, es=..., iorb=...)` — the spin-flip contribution,
  $\sum_{a=x,y,z} S_{aa}(\omega)$ built from the fluctuation operators
  $S^a_{i}-\langle S^a_{i}\rangle$ on orbital `iorb`, returned over both
  bias polarities.
- `get_orbital_cotunneling(fc, es=..., iorb=...)` — the orbital
  cotunneling contribution, summing correlators of
  $c^\dagger_{\texttt{iorb}\sigma}c_{j\sigma'}$ over every *other*
  orbital $j$ and both spin channels, again for both bias polarities.

`iorb` selects the orbital the STM tip couples to (e.g.\ $d_{z^2}$ for a
tip directly above the atom), and it matters: with a crystal field on,
different orbitals give different spectral weights *and* different peak
positions. An out-of-range `iorb` raises `ValueError`. Both routines run
in ED by default (`mode="ED"`, `submode="ED"`), so the `dex`/`T` caveat of
§6 applies directly to them — the low-lying multiplet of a magnetic atom
is precisely the near-degenerate manifold that discussion is about. See
`examples/kondo/atom_iets_orbital_resolved`.

> **Note (behaviour change).** Before 2026-08-21 both entry points
> accepted `iorb` and then silently overwrote it with `0`, so every call
> returned orbital 0's spectrum regardless of what was requested. Results
> obtained with `iorb!=0` on an earlier version need regenerating.

## 18. Infinite chains (iDMRG)

Every method above works on a *finite* chain of `n` sites. `Infinite_Many_Body_Chain`/`Infinite_Spin_Chain` (`infinitechain.py`) instead describe a translationally-invariant chain that repeats a fixed-size **unit cell** forever in both directions, and solve it with infinite DMRG (iDMRG, White's growing algorithm generalized to a multi-site unit cell) rather than sweeping a fixed-length system. Currently limited to a 1- or 2-site unit cell (`n_uc<=2`) — a >2-site unit cell raises `NotImplementedError` at construction, see `pyitensor/idmrg.py`'s module docstring for why.

`itensor_version="python"` (default) or `3` (`ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)`) are both supported — `2` and `"julia_live"` have no iDMRG port and also raise `NotImplementedError`. `ic.gs_method` (`"vumps"` by default, on EITHER backend since 2026-08-08 — see below) picks the ground-state solver. Both solvers, on both backends, support `vev`/`correlator` (see "Static correlators" below): the v3 C++ backend's `gs_method="idmrg"` (`mpscpp3/chain_session.h`'s `Chain::idmrg_ground_state`) is a line-by-line port of `pyitensor/idmrg.py`'s own growing algorithm, and the three things the Python side had gained since that port — McCulloch's wavefunction prediction, the gauge-consistent unit-cell extraction described below, and the per-site energy baseline — are ported too, along with the static observables built on them (`Chain::idmrg_onsite_expectation`/`idmrg_two_point_correlator`/`idmrg_local_excitation_gap`). Relative speed depends on the model: v3's own `idmrg_ground_state` is ~2.6-2.7x slower than `"python"`'s for a gapless (critical) chain (the local 2-site solve needs close to its full Krylov dimension every macro-iteration), but can be *faster* for a gapped model, where that same solve converges quickly — see `docs/documentation.md`'s iDMRG section for benchmark numbers.

**Building the Hamiltonian: `L`/`C`/`R`-suffixed operators.** Instead of one operator list per absolute site index, an infinite chain exposes each operator three times per unit-cell site `i` (`i=0..n_uc-1`): `SxC[i]` (site `i` of the *central* cell — ordinary intra-cell use), `SxR[i]` (site `i` of the *next* cell), and `SxL[i]` (site `i` of the *previous* cell, provided purely so a coupling can be phrased in whichever direction reads most naturally). A term may touch at most two *adjacent* unit cells (any subset of `{L,C}` or of `{C,R}`, never both `L` and `R` at once) — `set_hamiltonian` raises `ValueError` otherwise:

```python
from dmrgpy import infinitechain
ic = infinitechain.Infinite_Spin_Chain(["1/2"])       # n_uc=1, uniform chain
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]  # NN Heisenberg
ic.set_hamiltonian(h)
```

A bare, single-site operator (e.g. `ic.SzC[0]`, no product) is also a valid term — a Zeeman-field-style onsite Hamiltonian, either on its own or added alongside bond terms.

A two-site unit cell can express a dimerized (alternating-bond) chain, and a coupling can skip over an intermediate site by reaching straight from `C` to `R`:

```python
ic = infinitechain.Infinite_Spin_Chain(["1/2", "1/2"])
h = (1.0*(ic.SxC[0]*ic.SxC[1] + ic.SyC[0]*ic.SyC[1] + ic.SzC[0]*ic.SzC[1])   # strong intra-cell bond
     + 0.4*(ic.SxC[1]*ic.SxR[0] + ic.SyC[1]*ic.SyR[0] + ic.SzC[1]*ic.SzR[0]))  # weak inter-cell bond
ic.set_hamiltonian(h)
```

**Ground-state energy density.** Since the chain is infinite, the physically meaningful quantity is the energy *per site*, not a total:

```python
density = ic.gs_energy()   # converged energy density (or the best value reached after ic.maxiter)
ic.converged                # True iff the density stabilized below ic.etol
```

`ic.maxm`/`ic.cutoff` cap the bond dimension/SVD truncation exactly like a finite chain's `maxm`/`cutoff`; `ic.maxiter`/`ic.etol` are iDMRG-specific: `maxiter` macro-iterations (each adds one full unit cell to each side) are run, stopping early once the energy-density finite difference between consecutive macro-iterations drops below `etol`. For a *gapless* model (e.g. the uniform Heisenberg chain above) this finite-difference convergence is only power-law in the number of iterations at any fixed `maxm`, not exponential, so it may not trip `etol` within a practical `maxiter` even though the energy density itself is already accurate to several digits — check the returned value directly rather than relying solely on `ic.converged`. A *gapped* model (e.g. the dimerized chain above, or any chain with `j_weak` bonds well below the uniform point) has a finite correlation length and converges both faster and more reliably.

**Fermionic infinite chains.** `Infinite_Many_Body_Chain`'s first argument is a list of site-type codes, so an infinite chain can be fermionic as well as a spin chain: `0` is a spinless fermion site (`C`/`Cdag`/`N`/`F`), `1` a native spinful (Electron/Hubbard) site carrying both flavours at once (`Cup`/`Cdagup`/`Cdn`/`Cdagdn`/`Nup`/`Ndn`, local dimension 4). The native spinful code is what makes a two-orbital model fit the `n_uc<=2` limit — a `Spinful_Fermionic_Chain`-style representation would need two tensor-network sites per orbital, i.e. four for a two-orbital cell:

```python
ic = infinitechain.Infinite_Many_Body_Chain([1, 1])   # two native spinful sites per cell
Cup  = [ic.get_operator("Cup", i, "C") for i in range(2)]
Cdup = [ic.get_operator("Cdagup", i, "C") for i in range(2)]
CupR = [ic.get_operator("Cup", i, "R") for i in range(2)]
ic.set_hamiltonian(Cdup[0]*CupR[0] + ...)
```

Both backends thread the Jordan-Wigner string themselves, *locally* between each term's own two endpoints, rather than from some absolute origin — an infinite chain has no site 1 for a finite-chain-style string to start at. Two consequences for what a Hamiltonian may contain:

- The "at most 2 distinct sites per term" rule is counted in *sites*, not operator factors: a three-operator product like `(N_f - 1/2) * (Cdag_f * C_c)` touches two sites and is fine. A term whose two endpoints have sites strictly between them (e.g. `Cdag` at cell site 0 and `C` at site 0 of the *next* cell, for `n_uc=2`) is also fine — the string across the intervening sites is inserted automatically.
- A term with **odd total fermion parity** (e.g. a bare `Cdag`) is rejected with `ValueError` from `set_hamiltonian`: its string would have to run to infinity in both directions. Such a term is not parity-conserving, so it cannot appear in a physical Hamiltonian anyway.

**Fermionic correlators** work the same way: `ic.correlator("Cdag", 0, "C", r)` returns the physical fermionic correlator, i.e. `<Cdag_0 (prod_{0<k<r} F_k) C_r>` with the Jordan-Wigner string threaded across every site strictly between the two operators — not the stringless product of two bare matrices, which is a different quantity whose error *grows with separation* (so it corrupts exactly the decay rate such a correlator is usually measured for). Supported on all three paths: `itensor_version="python"` with `gs_method="idmrg"` or `"vumps"`, and `itensor_version=3` with `gs_method="vumps"`. The endpoint matrices and the decision of whether a string is open at all come from the same helper that builds the Hamiltonian's own 2-site terms, so a correlator and a Hamiltonian term written with the same operator names cannot disagree about the convention. Parity-even operators (`N`, `Nup`, `Sz`, ...) get no string, exactly as before.

A pair with **odd total fermion parity** (e.g. `("Cdag", "N")`) raises `ValueError`: its string can never close on an infinite chain, and the quantity vanishes identically in any parity-conserving state anyway.

See `examples/idmrg/fermionic_infinite_chain/main.py`, which checks both the energy density (against the exact free-fermion band integral) and `<Cdag_0 C_r>` (against the exact one-body density matrix) on every backend/solver combination, for a chain including a hopping that skips a site so a real string is involved.

**Static correlators.** After `gs_energy()` (called automatically by `vev`/`correlator` if not already run), one- and two-point expectation values of the converged infinite chain are available via the standard infinite-MPS transfer-matrix formalism. All four combinations work: `itensor_version` `"python"` or `3`, each under `gs_method="vumps"` (the default, see below) or `gs_method="idmrg"` (described in this paragraph). The two backends are cross-checked directly against each other (`tests/test_idmrg_correlator_v3.py`, `examples/idmrg/idmrg_correlator_python_VS_v3`), agreeing to ~1e-8 or tighter on a gapped chain:

```python
ic.vev("Sz", 0)                    # <Sz> at site 0 of the unit cell
ic.correlator("Sz", 0, "Sz", r)    # <Sz(0) Sz(0+r)>, r measured in physical sites, r>=0
```

These are reconstructed *after* convergence, from the gauge-consistent unit cell the growing algorithm extracts from a single micro-step's own two-site wavefunction (`pyitensor/idmrg.py`'s `_theta_cell`, `Chain::idmrg_theta_cell` on the C++ side). Fermionic operators are handled: the Jordan-Wigner string is threaded across every site strictly between the two endpoints, so `ic.correlator("Cdag", 0, "C", r)` is the physical fermionic correlator, not a stringless product of two bare matrices (a pair of odd total fermion parity, whose string could never close, raises `ValueError`). Two cheap diagnostics are worth knowing about. `<H_uc>` — the sum of every bond's correlator inside one unit cell — must equal `n_uc * density` exactly for any translationally invariant state, so it is a model-agnostic check on whether a given `maxm`/`maxiter` is enough, needing no reference value. `ic._result.state_overlap` (on the raw `pyitensor.idmrg.IDMRGResult`) is the second: the overlap between each macro-iteration's converged local state and the prediction it started from, which approaches 1 once the state has genuinely stopped changing.

Both were, until recently, routinely far from their ideal values, and correlators built on `n_uc=1` in particular could come out with the wrong *sign*. That is fixed: the growing algorithm now carries the state across iterations with McCulloch's wavefunction prediction, and extracts the unit cell in a single, self-consistent gauge. Measured against exactly solvable references, `<H_uc> - n_uc*density` now lands at `1e-15..1e-9` (it previously missed by up to `0.12`), `state_overlap` reaches `1-1e-13` (it previously plateaued around `0.5-0.65` for `n_uc=1`), and the XX chain's `<Sz>` comes out at `1e-13` against an exact `0`. Correlators still converge more slowly in `maxm` than the energy density does — that is ordinary finite-bond-dimension physics, most visible for a gapless model — so the `<H_uc>` check remains the right thing to run before trusting a number.

**The default ground-state solver: VUMPS** (`ic.gs_method = "vumps"`, the default since 2026-08-08 — Variational Uniform Matrix Product States, Zauner-Stauber et al., arXiv:1701.07035; see `pyitensor/vumps.py`'s own module docstring for the algorithm) — instead of growing a finite window and truncating it down to `maxm` at every step (the `gs_method="idmrg"` growing algorithm above), VUMPS solves directly, in the thermodynamic limit, for the actual `maxm`-dimensional variational optimum (`ic.maxm` sets VUMPS's own target bond dimension `D` here too). Both `itensor_version="python"` (`pyitensor/vumps.py`) and `itensor_version=3` (`mpscpp3/chain_session.h`'s `Chain::vumps_ground_state`, a C++ port of the same algorithm — built from plain dense arrays closed over LAPACK rather than ITensor tensor-network objects, since the bond/physical dimensions this feature targets are always small; see that method's own doc comment) support `gs_method="vumps"`; the two are cross-checked directly against each other to ~1e-10 or tighter on TFIM/Heisenberg at `D=1,2,3` (`tests/test_vumps_v3.py`). Explicitly set `ic.gs_method = "idmrg"` instead for `local_excitation_gap`/`td_dynamical_correlator` (no VUMPS equivalent, see their own sections below) or for VUMPS's documented D>1 convergence-robustness gap (see below) if the growing algorithm's own more battle-tested behavior is preferred:

```python
ic = infinitechain.Infinite_Spin_Chain(["1/2"])
ic.gs_method = "vumps"
ic.maxm = 4        # VUMPS's own target bond dimension D
ic.maxiter = 800   # VUMPS outer-iteration cap per bond dimension in its own D-ramp
ic.vumps_nrestarts = 6   # independent random-restart attempts per bond dimension
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]
ic.set_hamiltonian(h)
density = ic.gs_energy()
```

`vev`/`correlator` also work under `gs_method="vumps"`, on BOTH backends: `pyitensor.vumps.onsite_expectation`/`two_point_correlator` for `itensor_version="python"`, and `Chain::vumps_onsite_expectation`/`vumps_two_point_correlator` (a line-for-line C++ port of the same formula) for `itensor_version=3` — cross-checked directly against each other to ~1e-14 or tighter on TFIM at `D=2,3` (`tests/test_vumps_correlator_v3.py`). Both are computed directly from the converged mixed-gauge `{AC, AR}` rather than `pyitensor.idmrg`'s dominant-right-fixed-point eigenproblem: `AC` is already the exactly-normalized single-(super)site reduced state by construction of the mixed canonical gauge (Vanderstraeten, Haegeman, Verstraete, "Tangent-space methods for uniform matrix product states", arXiv:1810.07006, Eq.(34)), and `AR`'s exact right-orthonormality lets a two-point correlator spanning multiple unit cells close by a direct trace with no eigenproblem either (the mixed-gauge analogue of that same review's Eq.(37)-(39)). Unlike the growing-algorithm's own reconstructed-from-the-last-macro-iteration correlators above, these carry no `maxiter`/`<H_uc>`-self-consistency caveat — VUMPS solves directly at the target bond dimension in the thermodynamic limit, so once `.converged` is `True` the correlator is exact for that converged `{AL,AR,C}`, only limited by the bond dimension `D` itself (same caveat as the energy density, see below). `local_excitation_gap` is the one method that still requires `gs_method="idmrg"` specifically (it re-diagonalizes the growing algorithm's own final 2-site effective Hamiltonian, which has no VUMPS equivalent), on either backend at `window=0` (`Chain::idmrg_local_excitation_gap` for `itensor_version=3`); its `window>0` variant is `itensor_version="python"`-only, being an explicit prototype rather than stable API; conversely `excitation_energies`/`excitation_gap` (below) *require* `gs_method="vumps"` (the default) — they need `VUMPSResult`'s own mixed-gauge `{AL,AR,C,GL,GR}`, which the growing algorithm's `IDMRGResult` (`gs_method="idmrg"`) has no equivalent of.

```python
ic.gs_method = "vumps"
ic.vev("Sz", 0)                    # <Sz> at site 0, from the converged VUMPSResult
ic.correlator("Sz", 0, "Sz", r)    # <Sz(0) Sz(0+r)>, same signature as gs_method="idmrg"
```

**A real, scoped reliability caveat.** `D=1` converges reliably and exactly for every exactly-solvable (product-state-like) model tried (a pure field, a fully-decoupled Heisenberg dimer). Already at `D>1`, though — confirmed directly on the transverse-field Ising model, a simple gapped model with no special critical/symmetry complications — independent calls at the same `D` can land on noticeably different converged energies (most within a percent or two of the exact answer, but occasionally further off), even with `vumps_ground_state`'s own built-in D-ramp, multi-restart, and variational-principle safety-net machinery (see `pyitensor/vumps.py`'s own "Convergence robustness" docstring section for the full, numerically-confirmed account). A caller needing a reliable `D>1` result should call `gs_energy()` a few times independently (e.g. on freshly-constructed chains) and keep the lowest reported density, the same "rerun and take the best" discipline `idmrg_ground_state`'s own unseeded-random-MPS initialization already recommends elsewhere in this section. See `examples/idmrg/vumps_TFIM/main.py` for a worked example sweeping `D` and comparing VUMPS against both `gs_method="idmrg"` and the transverse-field Ising model's own exact (free-fermion) energy density.

Entanglement/entropy are not implemented for infinite chains yet.

**Excited states: the tangent-space/quasiparticle excitation ansatz.** Unlike a finite chain, an infinite chain's excitations form a momentum-resolved band `E(k)`, not a single discrete state, so "the first excited state" here means the standard single-mode/quasiparticle ansatz (Haegeman et al.), built on top of a **VUMPS** ground state's own mixed-gauge `{AL,AR,C}` representation (`ic.gs_method = "vumps"` — required, see above): a tangent-space vector with one excitation tensor `B` (a function of a free matrix `X`) inserted at every unit-cell position, weighted by momentum `k`. `ic.excitation_energies(k, n=1)` returns the lowest `n` excitation energies (above the ground state) at momentum `k` (radians, per unit cell); `ic.excitation_gap(ks=None)` scans `k` (default `numpy.linspace(-pi, pi, 41)`) and returns the minimum — the scalar "gap", mirroring the finite-chain `get_gap()` naming (§4):

```python
ic.gs_method = "vumps"             # required for excitation_energies/excitation_gap
ic.excitation_energies(0.0, n=1)   # lowest excitation energy at k=0
ic.excitation_gap()                # min_k E(k), the scalar gap
```

**Scope.** Both `itensor_version="python"` (`pyitensor/idmrg_excitations.py`) and `itensor_version=3` (`mpscpp3/chain_session.h`'s `Chain::vumps_excitation_energies`, a C++ port of the same algorithm, same dense-array/LAPACK approach as `vumps_ground_state` above) are supported — both require `gs_method="vumps"` (`NotImplementedError` otherwise, or if a different `itensor_version` is used). Any converged bond dimension `D>=1` is supported — including a genuinely entangled ground state (`D>1`, e.g. the transverse-field Ising or a dimerized Heisenberg chain), which used to be an explicit, rejected scope limit here (see `pyitensor/idmrg_excitations.py`'s own module docstring, "History" section, for the eight-pass investigation that limit came from and how it was eventually resolved by rewriting the ansatz from scratch on top of VUMPS's own mixed-gauge state, mirroring MPSKit.jl's own architecture). For `D=1` the computed dispersion matches the exact free-fermion single-magnon dispersion of a field-polarized XX chain to ~14 digits across the whole Brillouin zone (`examples/idmrg/excitation_gap_xx/main.py`); for `D=2` it matches an independently-converged MPSKit.jl transverse-field Ising result to 6 significant figures (`examples/idmrg/excitation_gap_tfim/main.py`), and `H_eff(k)` is Hermitian to machine precision at every `D` tried. `itensor_version=3` is cross-checked directly against `itensor_version="python"` across a full momentum scan, matching to ~1e-10 or tighter on the gapped TFIM case and to a looser ~1e-4..1e-7 on the gapless/critical Heisenberg case (both backends' own non-convex VUMPS restart search can land on slightly different local optima there, not a discrepancy in the ported algorithm itself) — see `tests/test_vumps_excitations_v3.py` and `examples/idmrg/vumps_excitation_v3_VS_python/main.py`. A deliberately constructed longer-range term (spanning more than 2 adjacent unit cells after `n_uc`-grouping) is rejected with `NotImplementedError`, raised from `gs_energy()` itself (`vumps_ground_state`'s own reach-1 check) rather than from the excitation machinery.

**Unit cells of more than 2 sites.** `Infinite_Many_Body_Chain([...])` accepts a cell of any length, and `gs_method="vumps"` (the default) then runs the **sequential multi-site** algorithm — the state is a list of per-site tensors `AL[n]`/`AR[n]`/`C[n]`, and one iteration sweeps the cell solving a one-site eigenproblem `H_AC[n]` at each site and a zero-site `H_C[n]` at each bond (`pyitensor/vumps_ms.py`). Nothing is grouped, so the cost is **linear** in the cell size.

That distinction is the whole point. The older way to support a multi-site cell is to fold it into one supersite of dimension `prod(d_p)` and run the single-site algorithm — exact, but exponential: a 4-site spinful cell is `d_g = 256`. Nietner, Vanhecke, Verstraete, Eisert and Vanderstraeten ([arXiv:2003.01142](https://arxiv.org/abs/2003.01142)) state it directly — *"the cost of a naive application of the VUMPS algorithm would scale exponentially with the size of the unit cell"* — and give as their key property *"a computational effort that scales linearly rather than exponentially in the size of the unit cell"*. The sequential algorithm is the one production codes implement (TeNPy's `SingleSiteVUMPSEngine`/`TwoSiteVUMPSEngine` sweep over the MPS unit cell; MPSKit likewise), and it originates with the multi-site VUMPS of Zauner-Stauber, Vanderstraeten, Fishman, Verstraete and Haegeman ([PRB 97, 045145 (2018)](https://arxiv.org/abs/1701.07035)).

`vev` and `correlator` work at any `n_uc`. Two things do not, and say so rather than failing obscurely: `gs_method="idmrg"` still requires `n_uc <= 2` (its growth loop pairs sublattice `m` with `n_uc-1-m`, which are only genuinely adjacent for `n_uc <= 2`), and so do `excitation_energies`/`excitation_gap`, whose tangent-space ansatz is still written against the grouped single-supersite gauge. `n_uc <= 2` keeps using the grouped VUMPS path, so its values are unchanged; the two agree to machine precision where both apply.

**Product-state traps, and the noise that breaks them.** Every particle-number-conserving Hamiltonian has the vacuum (and the filled state) as an *exact* eigenstate, and an exact eigenstate is an absorbing fixed point of the growing algorithm: the local solve warm-started there returns it immediately, its Schmidt rank is 1, so every subsequent truncation keeps a single singular value and nothing grows the bond dimension back. The run then reports a product state with `converged=True` in a fraction of a second — a silently wrong ground state, not a visible failure. Both backends do this; it is a property of the algorithm, not of either port.

`ic.noise` (default `1e-4`) applies White's density-matrix perturbation to break it, together with a matching random admixture on the local solve's start vector — both halves are needed, since enlarging the basis alone leaves the solve pinned to an exact eigenvector it cannot leave. The schedule is **demand-driven**: noise arms only while the state genuinely is a product state (keyed on the purity of the noise-free reduced state) and runs a few iterations past that, so an already-entangled model never sees it and is bit-for-bit unaffected. `ic.noise_iters` (default 40) caps the total, so a model whose ground state genuinely *is* a product state — a field-polarized chain — stops re-arming it and still ends on a clean tail; such models come back exact. Convergence is never declared while noise is active. Set `ic.noise = 0` to disable the mechanism entirely and recover the previous numerics exactly. See `docs/known_issue_idmrg_product_state_collapse.md` for the measured before/after against an exact band integral.

**A cheaper, cruder alternative: the local superblock gap.** `ic.local_excitation_gap(niter=200)` re-diagonalizes the growing algorithm's own final, converged 2-site effective Hamiltonian — already solved once for the ground state — for its *second*-lowest eigenvalue instead, and returns the difference. This is the direct infinite-chain analogue of a well-known finite-DMRG trick: at the last sweep, ask the same local effective Hamiltonian for its two lowest Ritz pairs rather than just the ground state. It is also the natural place a Lagrange-multiplier/orthogonality-penalty idea (as used by finite DMRG's own dedicated excited-state method, `get_excited`/`get_gap`, §4) would show up here — except no penalty weight is needed: since there is no separate re-sweep (the "state to stay orthogonal to" is just the local ground vector already found, in the very same local Hilbert space), the constraint is enforced exactly via deflation (projecting out the ground vector), which is what a penalty method converges to anyway as its weight → ∞. Unlike `excitation_gap`, this requires `gs_method="idmrg"` (not `"vumps"`, the default) and carries no momentum label — it is a single number, not a dispersion:

```python
ic.local_excitation_gap()   # a single scalar, no k dependence
```

**Accuracy caveat.** This is a genuinely cruder notion of "gap": it reuses the ground state's own `HL`/`HR` environments unmodified (never letting them relax for whatever the second local eigenstate actually represents), so it is *not* guaranteed to match the true minimum-momentum gap the tangent-space ansatz targets. Measured directly on the two cases `excitation_gap` and this method were both cross-checked against: for the exactly-solvable field-polarized XX chain (`D=1`), it comes out ~10% too high (5.5 vs. the exact 5.0); for a genuinely entangled (`D>1`) dimerized Heisenberg chain, it lands within ~0.5% of a large finite open chain's own extrapolated ED gap (see `examples/idmrg/local_excitation_gap/main.py`). Prefer `excitation_gap` for a physically principled answer (it now handles any `D`, see above); `local_excitation_gap` remains useful as a cheap, order-of-magnitude cross-check, or when `gs_method="idmrg"` is what's already been run for other reasons (`vev`/`correlator`). **Both of those calibrations were measured on spin models**, where the charge channel described below does not exist at all; on a fermionic model with `ConserveQNs=false` there is no reason those error bars carry over, and the estimate should be cross-checked against a finite chain, the correlator decay, or an `n(mu)` plateau before it is relied on.

**Why the excited solve shifts rather than only projects.** The second eigenvalue is the lowest one of the *same* stored operator restricted to its ground state's orthogonal complement. Writing that as the bare projector `P = I - |psi0><psi0|` and diagonalizing `P H P` is wrong, and used to make this method return large **negative** "gaps": `P H P` agrees with `H` on the complement, but it also carries `psi0` itself as an eigenvector with eigenvalue *exactly zero*. That is harmless only while the rest of the spectrum lies below zero — the moment the stored operator's own ground eigenvalue is **positive**, zero becomes `P H P`'s smallest eigenvalue and is exactly what a smallest-eigenvalue solver is asked to find (deflation keeps `psi0` out of the Krylov space in exact arithmetic only; rounding regrows the component and a restarted solver locks onto it). The reported gap is then `0 - e0`, i.e. precisely minus the stored superblock energy — reported from a spinful (c,f) Kondo chain at `maxm=16` as −358.003 meV against a stored energy of +0.358002587681, digit for digit. The sign of that stored energy is nobody's choice: the per-site energy baseline leaves it as a small residual boundary term of either sign, which is why the same model was right at `maxm=8` (energy −1.05, negative) and wrong at `maxm=16`. Both backends now use `P H P + sigma |psi0><psi0|` with `sigma` far above the bottom of the spectrum, so `psi0`'s own eigenvalue can never be the answer while the complement is untouched. Both eigenvalues are also re-solved from the stored superblock at the same solver strength rather than the ground one being read back from the growing algorithm, and if the re-solve lands *below* the growing algorithm's own value a `RuntimeWarning` is raised — the returned number is then still a genuine spectral gap of the stored operator, but no longer a gap measured above the state whose observables `vev`/`correlator` report.

**Charge excitations, and why a uniform on-site shift moves this number.** Both backends build every site with `ConserveQNs=false` (`mpscpp3/get_sites.h`; `pyitensor` likewise), so the local superblock spectrum contains particle-number-*changing* excitations against the frozen environment, with nothing to confine the deflated solve to the ground state's own charge sector. On a fermionic model the gap can therefore be set by a charge excitation — and adding a constant `mu` to every on-site energy, which in a gapped plateau leaves the converged state and all its correlators untouched, genuinely changes the answer. Measured exactly (by dense diagonalization of the stored superblock) on a gapped SSH spinless-fermion chain at `maxm=16`: the gap is 0.631321 at `mu=0` and 0.431325 / 0.431324 at `mu=±0.2`, i.e. exactly `|mu|` lower, because the two ±1-electron excitations that are degenerate at `mu=0` split; the stored ground vector is the exact ground state at every point (overlap 1.0000), so no solver is at fault. This is a property of the estimator, not a defect to fix: "the lowest state orthogonal to the ground state" simply is not a fixed-charge quantity here. Read it as an order-of-magnitude cross-check, and cross-check a fermionic gap against a finite-chain calculation, the correlator decay, or an `n(mu)` plateau before relying on it.

**Tightening it further: `window=`.** `ic.local_excitation_gap(window=w, niter=200)` (default `w=0`, exactly the behavior above) grows the local diagonalization block by `w` extra *free* physical sites on each side of the original 2, re-solving both the ground state and the deflated first excited state fresh within this larger block (rather than reusing the growing algorithm's own ground vector) — i.e. it lets the excitation spread across more of the chain instead of only ever living on the same frozen 2 sites, at the cost of an exponentially larger local Hilbert space (`d**(2*w)` more, `d` = the physical dimension). This is *not* the same thing as throwing more Krylov vectors at the original 2-site problem — `local_excitation_gap` solves that one to Krylov convergence already, so more iterations there do not change the answer (with the caveat above that the *ground* eigenvalue it deflates against must be the right one, which is why both are re-solved together). Measured directly: on the field-polarized XX chain the error drops from 10% at `w=0` to 3.8%/2.0%/0.81%/0.44% at `w=1/2/4/6`; on a gapped, genuinely entangled (`D>1`) transverse-field Ising chain, `w=3` (an 8-site local block) matches an 18-site open finite chain's own ED gap to <1%, converging at least as fast as growing the finite chain itself does. The improvement is not free, and its *rate* depends on the model's correlation length: on the S=1 Heisenberg chain (the Haldane gap, correlation length ~6 sites), `w=0/1/2` only move the estimate from 29%/24%/21% too high — still improving, but far more slowly, since a handful of extra sites barely dents a correlation length that long, and the physical dimension `d=3` makes each extra site pair 9x more expensive rather than XX/TFIM's 4x. Only `n_uc=1` is supported for `w>0` — widening needs to know which sublattice position each extra site takes, which is not tracked for `n_uc=2` yet (raises `NotImplementedError`). `w>0` is also `itensor_version="python"`-only: it is a prototype rather than stable API, so it was deliberately left out of the v3 C++ port (which does cover `w=0`).

**Dynamical correlators (finite-window KPM).** `ic.kpm_finite(opname_i, p_i, opname_j, r, n_window, window_chain_kwargs=None, **kwargs)` computes `<opname_i(site p_i) opname_j(site p_i+r)>(omega)` of the infinite chain, reusing the existing finite-chain KPM machinery (`kpmdmrg.get_dynamical_correlator`) unmodified rather than a new Chebyshev-recursion implementation. Named `kpm_finite` (not `get_dynamical_correlator`, the finite-chain method it wraps) to flag up front that it is the finite-window *approximation* described below, not an exact infinite-size calculation. There is no infinite-chain analogue of `vev`/`correlator`'s transfer-matrix formalism for a *dynamical* quantity — the Hamiltonian is extensive/unbounded in the thermodynamic limit, so a literal Chebyshev expansion of the full `H` has no meaning (unlike `apply_mpo`/`imps_sum`'s bounded-operator scope above). Instead, this method builds an ordinary finite, open-boundary chain of `n_window` repeats of the unit cell (Hamiltonian = `h_intra` tiled onto every cell plus `h_inter` tiled onto every adjacent pair of cells — one bond fewer than a periodic ring at the two open ends), places the two operators at the window's *central* unit cell (as far as possible from both ends), and delegates directly to the same KPM code path an ordinary finite `Spin_Chain`/`Fermionic_Chain` uses:

```python
es, ys = ic.kpm_finite("Sz", 0, "Sz", 0, n_window=16,
        window_chain_kwargs=dict(maxm=30, nsweeps=15, kpmmaxm=60),
        delta=0.3, es=np.linspace(-1, 5, 200))
```

`n_window` has no default — see the convergence caveat below. `window_chain_kwargs` is an optional dict of attribute overrides applied to the temporary finite `Many_Body_Chain` (`maxm`, `nsweeps`, `kpmmaxm`, `kpm_scale`, ...), independent of `ic`'s own `maxm`/etc.; remaining `**kwargs` (`delta`, `kernel`, `es`, `deconvolve`, ...) are forwarded to `kpmdmrg.get_dynamical_correlator` unchanged.

**Scope restriction — a finite-window approximation, read before use.** This is *not* an exact infinite-size method: results carry finite-size/open-boundary corrections that must be checked by convergence in `n_window`, exactly as a static `vev`/`correlator` caller would check `maxm`/`etol` convergence of the original iDMRG ground state. One Chebyshev moment corresponds to one application of the (nearest-neighbor) window Hamiltonian, so it can only move information by ~1 site per moment (a Lieb-Robinson-style bound) — but KPM's own moment count scales with the *window's own extensive bandwidth* divided by the requested `delta` (an ordinary finite chain's KPM already has this property, nothing new here), so a genuinely fine `delta` can require a moment count comparable to (or larger than) `n_window` itself, at which point open-boundary reflections contaminate the result regardless of how large `n_window` is. Prefer a coarser `delta`, or check that the correlator has visibly converged with growing `n_window`, for quantitative work (especially near a gapless point, where a fine `delta` is most tempting). Unlike `vev`/`correlator`, this does not need `ic._result` (no dependency on a previously converged `IDMRGResult`, or even on `ic.itensor_version`), so it works regardless of which backend `gs_energy()` itself used. See `examples/idmrg/dynamical_correlator_finite_window/main.py` for a worked example sweeping `n_window`.

**A genuinely infinite-chain dynamical correlator, via real-time TDVP (`td_dynamical_correlator`).** `kpm_finite`'s own open-boundary window has a real error source no amount of `n_window` alone can fix: an open chain's own ground state carries boundary artifacts (e.g. Friedel-oscillation-like features) that contaminate even the *central* region, not just the two edges. `ic.td_dynamical_correlator(opname_i, p_i, opname_j, n_window, dt=0.1, nt=200, x_values=None, maxdim=60, cutoff=1e-10, niter=50, connected=True, **kwargs)` fixes this by capping the window's two ends with the *converged* iDMRG growth environment (`idmrg_ground_state`'s own `HL`/`HR`, already computed during growth and exposed on `IDMRGResult`) instead of plain open boundaries — infinite boundary conditions (IBC), following Milsted/Vanderstraeten et al., "Infinite boundary conditions for response functions and limit cycles in iDMRG" (arXiv:1804.09163) — and evolves the perturbed window in real time via two-site TDVP rather than expanding in Chebyshev moments. Supports both `itensor_version="python"` and `itensor_version=3` (calls `gs_energy()` automatically if needed, like `vev`/`correlator`, though `vev`/`correlator` themselves still require `itensor_version="python"` regardless):

The native ITensor v3 backend (`Chain::td_dynamical_correlator_window`, `mpscpp3/chain_session.h`) reuses the vendored ITensorTDVP library's own boundary-tensor `tdvp(psi,H,t,LH,RH,sweeps,args)` overload directly against a tiled window MPS/MPO — unlike the `"python"` backend, which has to hand-roll its own window-aware TDVP sweep (pyitensor's generic TDVP infers a site's Link via a same-Index chain-neighbor lookup that cannot see a window's extra boundary legs). Two scope differences versus `"python"`: (1) `x_values` may not extend beyond the window's own explicit range (`center+x` must stay within the window, i.e. increase `n_window` instead of relying on padding) — the `"python"` backend pads beyond the window with extra unevolved unit-cell copies, not ported here; (2) as of this writing, `itensor_version=3`'s own `idmrg_ground_state` has a known, pre-existing convergence bug for Hamiltonians with an onsite ("field") term (energy diverging every macro-iteration) — unrelated to `td_dynamical_correlator` itself, but it means a v3 `td_dynamical_correlator` call inherits that limitation for such models; a purely bond-coupled Hamiltonian (e.g. plain Heisenberg) is unaffected.

```python
ks, es, Skw = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=14, dt=0.05, nt=200,
        maxdim=60, cutoff=1e-10, x_values=range(-6, 7),
        ks=np.linspace(-np.pi, np.pi, 41), delta=0.1, window=[-1, 6])
```

`opname_j` is applied at sublattice position `p_i` and evolved forward in time under the window's own ground-state-energy-shifted Hamiltonian (`e^{-i(H-E_GS)t}`, matching this codebase's own established real-time correlator convention, e.g. `mpscpp3::quench_tdvp`'s `Hshift=H-EGS*Id`); `opname_i` is inserted at the shifted position (bra side, not itself evolved) — this is the paper's own headline efficiency result: every `x` (and, via the spatial Fourier transform below, every `k`) comes from this *one* window evolution, not one run per distance the way a naive real-time approach (or `kpm_finite`'s own one-run-per-`r` KPM calls) would need. Cross-checked against an exact non-interacting (free-fermion) reference for the XX chain, on systems far larger than many-body ED could reach — see `docs/documentation.md`'s own architecture-level notes on this method for the full derivation and what that check found. `connected=True` (default) subtracts the disconnected background `<opname_i><opname_j>` before the spatial Fourier transform — turning it off produces a spurious, dominant `k=0` contribution with no discernible dispersion (the raw correlator approaches `<opname_i><opname_j>`, not 0, at large separation). `Skw` (shape `(len(ks), len(es))`) is obtained via a spatial DFT (`S(k,t)=sum_x e^{-ikx}S(x,t)`) followed by the *same* damping/FFT convention (`delta` -> Lorentzian broadening) every other dynamical-correlator submode in this codebase already uses.

**Scope**: this simplifies the paper's own Eq. 7 to `t1=0` (the ground state is perturbed by `opname_j` and evolved forward only; `opname_i` is never itself time-evolved) rather than the full two-branch trick (evolving a *second*, independent window backward in time too, which doubles the accessible total time for the same TDVP cost) — a documented, straightforward follow-up. `n_window` and `x_values` are two *separate* convergence axes to check (not one): `n_window` controls how much environment margin surrounds the perturbation, `x_values` how far the spatial sum/Fourier transform reaches — growing `x_values` together with `n_window` can keep changing results (a slowly decaying connected-correlator tail keeps contributing as the range widens), so converge each independently. See `examples/idmrg/td_dynamical_correlator/main.py` for a worked `n_window`-convergence sweep and a cross-check against `kpm_finite`.

**Applying an operator/gate to the converged chain (advanced).** `pyitensor.idmrg.apply_mpo(result, W_bulk, cutoff=..., maxdim=...)` is the infinite-chain analogue of the finite backends' `applyMPO`: it contracts a periodic MPO onto every site of the converged unit cell and re-canonicalizes/truncates the grown bond dimension back down via the standard two-sided fixed-point infinite-MPS canonicalization procedure, returning a new `pyitensor.idmrg.PeriodicMPS` that `onsite_expectation`/`two_point_correlator` accept exactly like an `IDMRGResult`. There is no `Infinite_Many_Body_Chain`-level wrapper yet, so it is reached by working with `pyitensor.idmrg` directly against `ic._result` (only meaningful for `itensor_version="python"`, same restriction as `vev`/`correlator`):

```python
from dmrgpy.pyitensor import idmrg
from dmrgpy.pyitensor.index import Index
from dmrgpy.pyitensor.tensor import ITensor

sites_uc = ic._result.sites_uc
d = sites_uc.dim(1)
pauli_x = 2 * sites_uc.site_type(1).matrix("Sx")           # a chi_W=1 operator
link_l, link_r = Index(1, tags="Link"), Index(1, tags="Link")
s = sites_uc.si(1)
W = [ITensor((link_l, s, s.prime(1), link_r), pauli_x.reshape(1, d, d, 1))]

new_state = idmrg.apply_mpo(ic._result, W, cutoff=1e-12, maxdim=None)
idmrg.onsite_expectation(new_state, "Sz", 0)               # <Sz> of the new state
```

**Scope restriction — bounded operators only.** `W_bulk` must represent a *bounded* (non-extensive) periodic operator: the same tensor reused at every unit cell, with no unconditional "keep accumulating forever" self-loop — single-site products (as above), gates tiled once per unit cell (an SVD-split 2-site gate embedded with identity everywhere else), symmetry operators, and the like. `pyitensor/idmrg.py`'s own Hamiltonian automaton (built internally by `_build_periodic_mpo` for `gs_energy()`) is deliberately the *other* kind — its accumulator channel needs a genuine chain boundary to correctly represent an extensive sum, so feeding it into `apply_mpo`'s boundary-less periodic contraction does not compute "H|psi>" (see `pyitensor/idmrg.py`'s own "Applying a (bounded) MPO to the converged iMPS" section docstring for the specific failure mode). This makes `apply_mpo` the natural building block for e.g. an iTEBD-style real/imaginary-time evolution step (repeatedly apply a local Trotter gate + truncate) — not yet implemented as a public feature, but `apply_mpo` is exactly the primitive it would use.

`W_bulk`'s physical Indices must be the *same objects* as `result.sites_uc`'s own (`sites_uc.si(i+1)`), so a custom operator automaton must be built against `result.sites_uc` directly (as in the example above), not a freshly constructed `SiteX`. The truncation/regauging step's own numerical conditioning degrades the further the *raw* grown bond dimension (`chi_A * chi_W`, before any truncation) exceeds the state's real entanglement at that cut — keep `maxm`/`maxdim` modest relative to `chi_W` for the best-conditioned result.

**Overlap/fidelity between two converged infinite MPS (advanced).** `pyitensor.idmrg.imps_overlap(result_a, result_b, normalize=True)` computes the per-unit-cell overlap between two `IDMRGResult`/`PeriodicMPS` objects — the infinite-chain notion of a finite MPS inner product `<phi|psi>`. A literal `<phi|psi>` over an infinite chain is not, in general, a finite number (it scales as `eta**N` over `N` unit cells, `N -> infinity`, where `eta` is the dominant eigenvalue of the mixed transfer matrix built from the two states' own tensors), so by default `imps_overlap` returns the always-finite *per-site fidelity* instead — magnitude 1 iff the two iMPS represent the same physical state (any gauge or normalization convention on the raw tensors), magnitude `<1` otherwise:

```python
from dmrgpy.pyitensor import idmrg

idmrg.imps_overlap(ic._result, ic._result)          # 1 -- same state
idmrg.imps_overlap(ic._result, new_state)            # cross-check apply_mpo didn't change the physical state
idmrg.imps_overlap(ic._result, some_other_result)    # < 1 in magnitude for a genuinely different state
```

`result_a`/`result_b` must share the same `n_uc` and, at every sublattice position, the same local physical dimension — `imps_overlap` raises `ValueError` otherwise. The two states' bond dimensions need *not* match (e.g. comparing a ground state against an `apply_mpo` output truncated to a different `maxdim`). Pass `normalize=False` for the raw, un-normalized mixed-transfer eigenvalue instead — mainly a diagnostic, analogous to `apply_mpo`'s own returned `.eta`.

**Direct sum of two converged infinite MPS (advanced).** `pyitensor.idmrg.imps_sum(result_a, result_b, cutoff=..., maxdim=...)` is the periodic-chain analogue of the finite backends' `mpsalgebra.sum`: a block-diagonal-in-the-bond-space direct sum at every unit-cell cut (there is no open boundary to instead concatenate along the way a finite chain's two ends do), re-canonicalized/truncated via the same `apply_mpo`-style two-sided fixed-point procedure, returning a new `PeriodicMPS`:

```python
from dmrgpy.pyitensor import idmrg

result = idmrg.imps_sum(result_a, result_b, cutoff=1e-12, maxdim=None)
```

**Scope restriction — read before use.** Tiled to the thermodynamic limit, this "+" does *not* represent a literal Hilbert-space vector sum the way the finite-chain construction does. It is only well-posed (a single dominant branch surviving — the mathematically correct infinite-volume answer, not a truncation artifact) when `result_a`/`result_b` have a genuine per-site norm mismatch, i.e. different self-overlap transfer eigenvalues (`eta`, in the sense of `imps_overlap`'s own `eta_aa`/`eta_bb`). This mismatch never happens between two *ordinary* `IDMRGResult`s — every one is individually normalized to `eta=1` exactly by left-canonical SVD construction — so summing two separately-converged ground states (the most natural reason to want this operation, e.g. to combine two symmetry-related solutions of the same Hamiltonian) hits a genuine tie every time. That tied case is a "cat state" superposition of two macroscopically distinct branches, which this module's single-fixed-point machinery cannot represent as one canonical periodic MPS — `imps_sum` raises `RuntimeError` there (via `pyitensor.idmrg`'s internal degeneracy check on the combined transfer matrix's dominant eigenvalue) rather than silently collapsing to one arbitrary branch. Correctly evaluating observables on that common, physically meaningful case (e.g. via the thermodynamic-limit local-observable identity for two orthogonal, equally-weighted branches) needs correlator machinery this module does not have yet — a documented, deliberate scope limit, not something silently gotten wrong. See `pyitensor/idmrg.py`'s own "Summing two converged iMPS" section docstring for the full derivation and empirical confirmation.

**Direct sum of two converged VUMPS iMPS (advanced).** `pyitensor.vumps.imps_sum(result_a, result_b, cutoff=..., maxdim=...)` is the VUMPS-mixed-gauge analogue of `idmrg.imps_sum` immediately above — same construction in spirit (block-diagonal direct sum, re-canonicalized/truncated via `idmrg._canonicalize_periodic`), but working on `VUMPSResult`'s own grouped-supersite `AL` tensor (bond dimension `D_a+D_b`) rather than `idmrg.py`'s per-sublattice periodic `U_list`, and returning a new `pyitensor.vumps.UniformMPS` (accepted directly by `vumps.onsite_expectation`/`two_point_correlator`, same duck typing as `idmrg.PeriodicMPS`) whose full mixed gauge `{AL,AR,C,AC}` is then reconstructed via `vumps._complete_mixed_gauge` — the standard "bringing a uniform MPS to canonical form" procedure (Vanderstraeten, Haegeman, Verstraete, "Tangent-space methods for uniform matrix product states", arXiv:1810.07006, Eq.(9)-(17): factor the truncated `AL`'s own dominant right transfer-matrix fixed point `r = C C^dagger`, then `AR := C^-1 AL C` and `AC := AL C`):

```python
from dmrgpy.pyitensor import vumps

result = vumps.imps_sum(result_a, result_b, cutoff=1e-12, maxdim=None)
```

Same physical scope restriction as `idmrg.imps_sum` above, for the same reason: every converged `VUMPSResult` has `AL` exactly left-canonical *and* `AR` exactly right-canonical by construction of the mixed gauge, so both its left and right self-overlap transfer eigenvalues are exactly `eta=1` — summing two ordinary `VUMPSResult`s therefore always hits the same degenerate-dominant-eigenvalue tie `idmrg.imps_sum` does, and `imps_sum` raises `RuntimeError` there rather than silently returning one arbitrary branch (`arXiv:1810.07006`'s own Sec. 2.1 independently makes the same point: non-injective/"cat state" MPS tensors are exactly the ones with a degenerate dominant transfer-matrix eigenvalue). Only two states with a genuine per-site norm mismatch (e.g. one deliberately rescaled, exactly mirroring `idmrg.imps_sum`'s own worked example) have a well-posed sum — see `pyitensor/vumps.py`'s own "Summing two converged VUMPS iMPS" section docstring for the full derivation, and `examples/idmrg/vumps_imps_sum/main.py` for a worked example of both cases plus a cross-check of the surviving branch's `onsite_expectation`/`two_point_correlator` at genuinely entangled `D>1`.

Also ported to `itensor_version=3` (`Chain::vumps_imps_sum`, `mpscpp3/chain_session.h`), a dense-array translation reusing the same `vx_canonicalize_n1`/`vumps_complete_mixed_gauge` machinery `apply_mpo` below shares — same scope restriction and degeneracy behavior as the pyitensor version above. There is no `Infinite_Many_Body_Chain`-level wrapper on this backend either (matching pyitensor's own scope), so it is reached directly against a Chain's own converged VUMPS snapshot:

```python
D, d_g, AL, AR, C, AC, eta = ic._session3.vumps_imps_sum(D_b, AL_b, cutoff=1e-12, maxdim=0)
ic._session3.vumps_load_uniform_state(D, d_g, AL.flatten().tolist(),
                                       AR.flatten().tolist(), C.flatten().tolist())
```

`AL_b` is the second state's own flat `(D_b,d_g,D_b)` `AL` array (row-major) — e.g. obtained from another Chain's own snapshot via `Chain::vumps_get_snapshot()` (which returns `(D,d_g,AL,AR,C)` for `this` Chain's own current state); `maxdim<=0` means no cap (Python's `maxdim=None`). `vumps_load_uniform_state` writes the result back into a Chain's own snapshot so `vumps_onsite_expectation`/`vumps_two_point_correlator` (and, since `ic._session3_has_vumps` is untouched by any of this, `ic.vev`/`ic.correlator` too) see it. See `tests/test_vumps_imps_sum_v3.py` for the full cross-check against `itensor_version="python"`.

**Applying an operator/gate to a converged VUMPS iMPS (advanced).** `pyitensor.vumps.apply_mpo(result, W_bulk, cutoff=..., maxdim=...)` is the VUMPS-mixed-gauge analogue of `idmrg.apply_mpo` above: it groups `W_bulk` into a single grouped-supersite MPO tensor (the same `_group_automaton` routine that groups VUMPS's own Hamiltonian automaton), grows the converged `AL` by it via `idmrg.grow_by_mpo`, re-canonicalizes/truncates via `idmrg._canonicalize_periodic` (the identical two-sided fixed-point procedure `idmrg.apply_mpo` uses), and completes the resulting truncated left-canonical tensor to the full mixed gauge `{AL,AR,C,AC}` via `vumps._complete_mixed_gauge` (the same completion `imps_sum` immediately above uses). `W_bulk` takes *exactly* `idmrg.apply_mpo`'s own convention — a list of `n_uc` rank-4 `(Left, in, out, Right)` ITensors, one per unit-cell sublattice site — so the identical `W_bulk` list built for one backend can be fed to the other's own `apply_mpo` to cross-check both against each other on the same operator:

```python
from dmrgpy.pyitensor import vumps

new_state = vumps.apply_mpo(result, W, cutoff=1e-12, maxdim=None)
```

Returns a new `pyitensor.vumps.UniformMPS`, accepted directly by `vumps.onsite_expectation`/`two_point_correlator` exactly like a `VUMPSResult`. Same scope restriction as `idmrg.apply_mpo` above — `W_bulk` must represent a *bounded* (non-extensive) periodic operator; the Hamiltonian's own automaton (built internally by `vumps_ground_state` for `gs_energy()`) is out of scope, for the identical reason described there. There is no `Infinite_Many_Body_Chain`-level wrapper yet either, so `pyitensor.vumps.apply_mpo` is reached directly, working against `ic._vumps_result` (only meaningful for `itensor_version="python"` with `gs_method="vumps"`). See `examples/idmrg/vumps_apply_mpo/main.py` for a worked example: a `chi_W=1` unitary single-site flip cross-checked exactly against `idmrg.apply_mpo` at the exact `D=1` field-polarized point, the same flip's `<Sz>`/`<Sz(0)Sz(r)>`/`eta` invariants checked at a genuinely entangled `D>1` TFIM ground state, and a genuinely bond-growing `chi_W>1` two-site gate tiled once per an `n_uc=2` unit cell.

Also ported to `itensor_version=3` (`Chain::vumps_apply_mpo`, `mpscpp3/chain_session.h`): grows a Chain's own converged `AL` by a caller-supplied `W_bulk` (grouped via `vumps_group_automaton`, reused unmodified from the Hamiltonian-automaton path since the grouping contraction itself doesn't care what it's grouping), then shares `vumps_imps_sum`'s own `vx_canonicalize_n1`/`vumps_complete_mixed_gauge` machinery. Same bounded-operator scope restriction, and the same "no `Infinite_Many_Body_Chain`-level wrapper" caveat as `vumps_imps_sum` above — reached directly against a Chain's own converged VUMPS snapshot:

```python
W_bulk_flat = [W0.flatten().tolist(), W1.flatten().tolist()]  # one per unit-cell site
D, d_g, AL, AR, C, AC, eta = ic._session3.vumps_apply_mpo(
    W_bulk_flat, Dw_left, Dw_right, cutoff=1e-12, maxdim=0)
ic._session3.vumps_load_uniform_state(D, d_g, AL.flatten().tolist(),
                                       AR.flatten().tolist(), C.flatten().tolist())
```

`W_bulk_flat[p]` is a dense, row-major `(Left,in,out,Right)` array (size `Dw_left[p]*d_p*d_p*Dw_right[p]`) — the same convention as `idmrg.apply_mpo`'s own ITensor list, just flattened. See `tests/test_vumps_apply_mpo_v3.py` and `examples/idmrg/vumps_apply_mpo_v3_VS_python/main.py` for the full cross-check against `itensor_version="python"`, including the same three cases (`D=1` exact, `D>1` unitary invariants, `chi_W>1` bond growth at `n_uc=2`).

## 19. Performance: BLAS threads

DMRG spends its time in a great many *small* dense linear-algebra calls
rather than a few large ones — a two-site tensor at `maxm=30` on spin-1/2
sites is a 60×60 matrix, and that is what hits `svd`/`eigh`/`matmul` tens
of thousands of times in one ground-state solve. At that size a
multithreaded BLAS can lose more to thread barriers than it gains: measured
with MKL, going from one thread to two made a 60×60 complex `svd` 1.6x
slower and `eigh` 2.3x slower.

End to end that alone is minor (~1.13x for `itensor_version="python"`,
nothing measurable for `3`). It matters when the machine is
**oversubscribed** — a shared cluster node, or several dmrgpy runs launched
in parallel — because each process's BLAS assumes it owns every core. On a
14-core host with another job holding 10 of them, letting MKL use its
default thread count turned a 9.4s `"python"` solve into 28.2s, and a 1.6s
`itensor_version=3` solve into 26.7s.

The most reliable fix is to pin threads before numpy is imported:

```bash
MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 myscript.py
```

For finer control there is an opt-in context manager (needs `threadpoolctl`,
which is not a declared dependency):

```python
from dmrgpy import blasthreads

print(blasthreads.current_hint())     # how threads are configured now
with blasthreads.limit(1):
    e0 = sc.gs_energy()               # restored on exit, including on error
```

dmrgpy never changes thread counts on its own: threading does pay off once
the tensors get large enough, and a script may have chosen its setting
deliberately. Treat the numbers above as indicative — they come from a busy
shared host — and measure on your own machine before tuning around them.
See `dmrgpy/blasthreads.py` for the full measurements.
