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
`itensor_version=2` isn't yet, and `examples/mixed_spin_fermion_chain`
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
actually agree; see `examples/boson_maxnb_v3_VS_ED`.

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
NH-DMRG on all three backends against exact diagonalization and the
Arnoldi route on an interacting fermionic chain with a staggered
imaginary potential. The biorthogonal pair $|\psi_R\rangle,|\psi_L\rangle$
is also what feeds the non-Hermitian dynamical correlator,
`get_dynamical_correlator(submode="KPM")` for $H\neq H^\dagger$ — see §6.

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

which is then windowed with an exponential damping factor
$e^{-\delta t}$ (equivalent to a Lorentzian broadening of width $\delta$
in frequency) and Fourier transformed,

$$S_{AB}(\omega)=\frac1\pi\,\mathrm{Re}\!\int_0^{T}\!dt\;e^{i\omega t}\,C(t)\,e^{-\delta t}$$

The total simulated time $T$ (`damping_periods`/$\delta$) must be long
enough that the damping has suppressed truncation ringing by $t=T$.
Frequency resolution is set by $T$ (via the usual $\Delta\omega\sim
1/T$ time-frequency uncertainty), so this method is best when you want
fine resolution over a *narrow* frequency window, at the cost of a real
dynamical simulation (bond dimension grows with entanglement generated
during the evolution).

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
`"python"`, `tevol_method="TDVP"`, the paper's own setup), one-site TDVP
with global subspace expansion (`tevol_method="TDVP_GSE"`, see §7 — same
`itensor_version` support as `"TDVP"`), or falls back to the MPO-Taylor
propagator otherwise (`tevol_method="MPO"`, or `itensor_version=2`, which
has no TDVP) — the same TDVP-vs-Taylor choice `"TD"` already makes. Current
scope: only the "greater" branch of the correlator is computed (the same
simplification `"TD"` itself already makes), so this is best used the
same way as `"TD"`: high-resolution work in a narrow frequency window,
now reachable at a lower bond-dimension cost for a given simulated time.

**`submode="EX"` — exact diagonalization in a truncated DMRG subspace.**
Builds $A$, $B$, $H$ explicitly in the subspace spanned by the lowest
`nex` (orthogonalized) DMRG excited states, then evaluates the exact
Lehmann sum in that subspace,

$$G_{AB}(\omega)=\sum_{n=1}^{n_{ex}}\frac{\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle}{\omega-E_n+i\delta}$$

i.e.\ a small, explicit sum over poles at the computed excited-state
energies $E_n$, each with residue given by the transition matrix
elements. Cheap and exact *within* the truncated subspace; only as good
as how many/which excited states were computed.

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

**Choosing the propagator: `sc.tevol_method`.** Three options:

- `"TDVP"` (the default) — two-site TDVP, which grows the MPS bond
  dimension via SVD the same way ground-state DMRG does. Used whenever
  `itensor_version` is `3` or `"python"`; `itensor_version=2` falls back
  to `"MPO"` (below) even with this default.
- `"TDVP_GSE"` (`itensor_version` 3 or `"python"` only, same support as
  plain `"TDVP"` above) — one-site TDVP preceded, for the first
  `sc.tdvp_gse_sweeps` steps (default 3), by a *global subspace
  expansion* step: a Krylov subspace $\{\psi,H\psi,H^2\psi,\dots\}$ of
  dimension `sc.tdvp_gse_krylov_order` (default 3) is used to enlarge the
  MPS's bond dimension *without changing the state it represents*, using
  a cutoff `sc.tdvp_gse_cutoff` (default $10^{-8}$) — the scheme of Yang
  & White, [arXiv:2005.06104](https://arxiv.org/abs/2005.06104)/Phys.
  Rev. B 102, 094315 (2020). Most useful when the starting state's bond
  dimension is small (e.g.\ a product-state quench) and one-site TDVP
  alone (which conserves bond dimension exactly) wouldn't be able to grow
  into the entanglement the subsequent evolution generates.
- `"MPO"` — a hand-rolled 2nd-order Taylor expansion of $e^{-iH\,dt}$
  applied as an MPO each step; the only option on `itensor_version=2`
  (which has no TDVP at all), and available (if slower/less accurate for
  a given bond dimension) everywhere else too.

```python
sc.tevol_method = "TDVP_GSE"
sc.tdvp_gse_sweeps = 3
sc.tdvp_gse_krylov_order = 3
sc.tdvp_gse_cutoff = 1e-8
```

## 8. Density of states

**Many-body density of states.** The full many-body spectral density

$$\rho(E)=\sum_n\delta(E-E_n)$$

(sum over *all* eigenstates of $H$, not resolved by any operator) is
obtained via the same KPM machinery as §6 applied directly to $H$
itself:

```python
from dmrgpy import dos
es, rho = dos.get_dos(sc)
```

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

`T=0` recovers ordinary ground-state DMRG. For small systems ($n\lesssim
14$), `get_correlation_matrix(T=...)` instead computes the exact thermal
average directly by brute-force ED, $\rho=\sum_nZ^{-1}e^{-E_n/T}|n\rangle\langle n|$ — a useful cross-check of the purification approach.

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

$$\frac{\partial I}{\partial V}(eV)\propto\sum_{i,f}p_i\Big[\tfrac12|\langle f|S_-|i\rangle|^2+\tfrac12|\langle f|S_+|i\rangle|^2+|\langle f|S_z|i\rangle|^2\Big]\,\Theta(eV-\epsilon_{if})+U^2$$

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

**Scope and known limitations**, worth reading before trusting specific
numbers:

- Only a single chain site couples to the tip (`site=`); the paper's own
  model allows several sites with independent tip couplings, not
  implemented here.
- The third-order terms (`order=3`) are $\partial I^{t\to s}/\partial V$
  only, not the full bidirectional net-current derivative that
  `order=2` provides. The paper never gives a general $t\leftrightarrow
  s$ formula for them — its own worked $S=1/2$ example shows the two
  directions are related by more than a plain $eV\to-eV$ mirror, so no
  such generalization is attempted here. This mainly affects the
  *symmetry* of the returned third-order spectrum, not its second-order
  part.
- The paper's own closed-form equations for two numerical building
  blocks — the temperature-broadened step $\Theta(x)$ and the
  temperature-broadened Kondo log function $F(\epsilon,T)$ — do not
  reproduce the physics the paper itself describes for them (checked
  directly: the printed $\Theta(x)$ diverges rather than saturating, and
  the printed $F$ closed form drops its own temperature broadening).
  `kondospectrumtk/stepfunctions.py` uses corrected forms instead,
  re-derived from the paper's own unambiguous defining integrals and
  verified against digitized values from the paper's own figures (see
  that module's docstring, and `examples/kondo_spectrum_VS_paper/`).
- The potential-interference term's general-spin closed form (`U!=0` in
  `order=3`) is an extrapolation from the paper's own worked $S=1/2$
  example (only that special case is spelled out in closed form in the
  paper) and carries lower confidence than the other terms; treat it as
  provisional.

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
  that function's docstring).

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
