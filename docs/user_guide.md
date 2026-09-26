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
19. [Running the pure-Python backend on a GPU](#19-running-the-pure-python-backend-on-a-gpu)
20. [Performance: BLAS threads](#20-performance-blas-threads)
21. [What raises, and what changed in the 2026-08 audit](#21-what-raises-and-what-changed-in-the-2026-08-audit)

## 1. Physical models and Hilbert spaces

Every model is a chain of $n$ local Hilbert spaces $\mathcal H=\bigotimes_{i=1}^n \mathcal H_i$, on which a Hamiltonian and observables are built out of local operators.

| Chain class | Local Hilbert space | Key operators |
|---|---|---|
| `spinchain.Spin_Chain` | spin-$S$, $S\in\{\tfrac12,1,\tfrac32,2,\tfrac52,3\}$ per site | $S^x_i,S^y_i,S^z_i$ |
| `fermionchain.Fermionic_Chain` | spinless fermion (occupied/empty) | $c_i,c_i^\dagger,n_i=c_i^\dagger c_i$, Jordan-Wigner string $F_i$ |
| `fermionchain.Majorana_Chain` | Majorana fermion | Majorana operators built from `Fermionic_Chain` |
| `fermionchain.Spinful_Fermionic_Chain` | spin-$\tfrac12$ fermion (4 states: $0,\uparrow,\downarrow,\uparrow\downarrow$), built from two interleaved spinless sites per physical site | $c_{i\sigma},c^\dagger_{i\sigma},n_{i\sigma}$, plus derived $S^x_i,S^y_i,S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$, onsite pairing $\Delta_i=\tfrac12 c_{i\uparrow}c_{i\downarrow}$ |
| `fermionchain.Spinful_Fermionic_Chain_Native` | same physics as `Spinful_Fermionic_Chain`, but on a genuinely 4-dimensional local space (one tensor-network site per physical site; `itensor_version=3` and `"python"` only) | identical operator lists/formulas as `Spinful_Fermionic_Chain` |
| `bosonchain.Bosonic_Chain` | truncated boson Fock space, $n_i\in\{0,\ldots,\text{maxnb}_i-1\}$, per-site dimension `maxnb` (default 4, i.e. up to 3 bosons/site) settable via `Bosonic_Chain(n, maxnb=[...])` | $a_i,a_i^\dagger,n_i$, occupation projectors $\hat n_i^{(k)}=\lvert k\rangle\langle k\rvert$ for $k=0,\ldots,\text{maxnb}_i-1$ (`bc.D[i][k]`, plus `bc.D0`..`bc.D3` when every site has `maxnb`$\,\ge 4$) |
| `parafermionchain.Parafermionic_Chain` | $\mathbb Z_N$ parafermion (clock model), $N\in\{2,3,4\}$ | clock/shift operators $\sigma_i,\tau_i$ and composite parafermion operators $\chi_i,\psi_i$ built as $\tau$-string $\times\sigma_i$ |
| `bosonchain.SpinBoson_Chain` | mixes truncated-boson sites and genuine spin-$S$ sites *in the same chain*, one entry per location (the constructor takes the list of site labels, e.g. `SpinBoson_Chain(["B","S=1/2",...])`) | at a boson location: $a_i,a_i^\dagger,n_i$ and the occupation projectors `sb.D[i][k]` (plus `D0`..`D3` when every boson location has $\ge4$ levels); at a spin location: $S^x_i,S^y_i,S^z_i$. As in `Mixed_Spin_Fermion_Chain`, the operators that do not apply at a given location read as the integer `0` |
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
DMRG tolerance under `itensor_version=3` and `itensor_version="python"`
(the two DMRG backends this class wires up; `itensor_version=2` and
`"julia_live"` have no native spinful site). Despite halving the site count, it is *not* generally
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
~30% win). Leave `ctmode` at its default for this class — the resolver
already picks the fastest available method (`"fold"` on
`itensor_version=3`, `"batched"` under `"python"`), both of which are
faster than the `"full"` this comparison was made against. Otherwise prefer
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
factory) — `itensor_version=2` understands only the single fixed 4-level
boson site (ITensor's `BosonFourSite`), and now says so: any `maxnb`
entry other than 4 raises `ValueError` naming
`itensor_version=3`/`"python"` as the alternatives (with `mode="ED"` on
either of them for exact diagonalization, since `mode="ED"` next to
`itensor_version=2` does not avoid the check: a constructor keyword is
applied after the session is built), *before* a session is built. It used to abort the whole interpreter with SIGABRT
from inside ITensor instead, which no Python code could catch. So run a
non-default `maxnb` under `itensor_version=3` (the default when the
compiled C++ extension is available, `"python"` otherwise) or
`"python"`; see `examples/boson_models/boson_maxnb_v3_VS_ED`. On a
3-site chain at `maxnb=[3,5,3]` the `"python"` DMRG and ED ground-state
energies agree to machine precision (measured 8.9e-15). `itensor_version=3`
agrees on the same physics but neither to that tolerance nor
reproducibly: it starts from a random MPS with QNs off, so its agreement
varies run to run — three runs of that example here disagreed by
1.05e-03, 1.78e-03 and 1.88e-03, each one over its own `tol = 1e-3`.
Expect a loose agreement from v3 here, not a tight one, and re-run before
concluding anything from a single number. The same
restriction is *believed* to hold for `itensor_version="julia_live"`,
but that has not been verified here and is not guarded — treat it as
unchecked rather than as a documented limitation.

`SpinBoson_Chain(sitesin, maxnb=None)` takes one label per location: a
spin label (`"S=1/2"`, `"1/2"`, `"S=1"`, ... as in `Spin_Chain`) or a
boson label, which is `"B"` for the default 4-level site or `"B<k>"`
(`"B6"`) to name the local dimension — `"B4"` and `"B"` are the same
site. `maxnb=` is the alternative spelling of the same thing and is
honored rather than dropped: one entry per location with `None` at the
spin positions, so `SpinBoson_Chain(["B","S=1/2"], maxnb=[6,None])` and
`SpinBoson_Chain(["B6","S=1/2"])` both build sites `[106, 2]`. A list of
the wrong length, an `n=` that contradicts `len(sitesin)`, a `maxnb` at
a spin position, and an unrecognized label all raise `ValueError` naming
what was wrong (an unknown label used to surface as `RuntimeError: No
active exception to reraise`). Per-site occupation projectors are
`sb.D[i][k]`, an empty list at a spin location, mirroring
`Bosonic_Chain`. This class has no ED backend: `mode="ED"` raises
`NotImplementedError` pointing at `mode="DMRG"` (`itensor_version=3` or
`"python"`), or at `Bosonic_Chain` for a purely bosonic model.

`itensor_version=` is now accepted by *every* chain constructor.
`Bosonic_Chain`, `SpinBoson_Chain`, `Parafermionic_Chain` and
`spinfermionchain.Spin_Fermion_Hamiltonian` used to raise `TypeError:
unexpected keyword argument` for it — which is awkward precisely where
the paragraph above tells a boson user to pick a backend — and had to be
constructed first and switched afterwards with `setup_python()` /
`setup_cpp(version=3)`. Both spellings work now.

Every other constructor keyword names a setting of the chain and takes
effect as if it had been assigned right after construction, meaning that
`Spin_Chain(sites, maxm=60, nsweeps=20)` is `sc = Spin_Chain(sites)`
followed by `sc.maxm, sc.nsweeps = 60, 20`, through the same checks
(`maxm=0` raises `ValueError`). This holds for `mode=` too:
`Spin_Chain(sites, mode="ED")` builds its DMRG session all the same, so
`sc.mode = None` later goes back to DMRG, and
`Thermal_Spin_Chain(sites, T, mode="ED")` sets the wrapper's own `tc.mode`,
which its `get_gs()` hands to `tc.MBChain`. A chain's `mode="ED"`
overrides the `mode=` of every call made on it, while `mode="DMRG"` on the
chain means what no mode means, so an explicit `mode="ED"` call is still
answered by ED; until the 2026-09-25b pass a chain set to `"DMRG"`
answered that cross-check by DMRG, and `Thermal_Spin_Chain` now defaults
to `mode=None`. The settings are an explicit list, `sites.SETTINGS`, and a
keyword that is not one of them (a misspelling, a private name, a method,
the chain's state such as `hamiltonian`, `conserved_sector`, `wf0`, `e0`
or `computed_gs`, the Hamiltonian pieces `hopping`, `hubbard`, `pairing`
and `exchange` that the `set_*` methods accumulate, or an operator list
the model class builds, as in `Fermionic_Chain(4, N=5)`) raises
`TypeError` naming every offending key, `Thermal_Spin_Chain` included.
Until the 2026-09-25 fixes every such keyword was dropped, so a bond
dimension asked for at construction silently ran at `maxm=30` (§21), and
until the 2026-09-25b pass nine public attributes that are not settings
were admitted, so `Fermionic_Chain(4, hubbard=2.0)` followed by
`set_hoppings()` added the constant $2\,\mathrm{Id}$ to the free-fermion
Hamiltonian ($-0.236068$ against $-2.236068$). Note that `cvm_maxm` is set
from the default `maxm` when the chain is built, so it stays 30 unless it
is given too, exactly as when `maxm` is assigned afterwards.

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

**Hamiltonian shortcuts on fermionic chains.** Rather than writing every
term out, the fermionic classes offer builders that take a function of
the site indices and add the whole family of terms at once:

| Call | Adds |
|---|---|
| `fc.set_hoppings(fun)` | $\sum_{ij}f(i,j)\,c_i^\dagger c_j$ (spinless) |
| `fc.set_hubbard(fun)` | $\sum_{ij}f(i,j)\,n_in_j$; on a spinful chain $n_i=n_{i\uparrow}+n_{i\downarrow}$ |
| `sfc.set_hoppings_spinful(fun)` | the same hopping, spin-diagonal, with `fun` indexed by *physical* site rather than by interleaved spinless site |

```python
fc = fermionchain.Fermionic_Chain(4)
fc.set_hoppings(lambda i,j: 1.0 if abs(i-j)==1 else 0.0)
fc.set_hubbard(lambda i,j: 2.0 if abs(i-j)==1 else 0.0)
```

**Algebra on already-built operators.** `sc.toMPO(h)` compiles a
`MultiOperator` into an already-built matrix product operator, on every
backend: a `StaticOperator` for `itensor_version` 2, 3 and `"python"`,
an `EDOperator` under `mode="ED"`, and `mpsjulialive`'s own `MPO` for
`itensor_version="julia_live"`. The `+`, `-`, unary `-` and scalar
`*`/`/` algebra below is `StaticOperator`-only — the `julia_live` `MPO`
supports `*` alone — so it is available for `itensor_version` 2, 3 and
`"python"`, not yet for `"julia_live"`. Two `StaticOperator`s can be
combined directly, without going back through the symbolic
`MultiOperator` form:

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

**MPS and operator algebra.** Every chain also exposes the primitives
that act on wavefunctions and operators directly, which the worked
examples later in this guide use (e.g. §3's `promote_to_dense` example
builds a photoemission weight out of `applyoperator` and `overlap`).
Each takes the same `mode=`/`**kwargs` as the rest of the API, with one
qualifier: `applyoperator`, `summps`, `applyinverse`, `scale_mps` and
`exponential` take the backend from the *wavefunction's own type*
(`mps.MPS` → DMRG, `edtk.edchain.State` → ED), because for these the
wavefunction is the backend. An explicit `mode=` there is a consistency
check, not a switch: naming the other backend raises `TypeError` telling
you to rebuild the wavefunction with `get_gs(mode=...)` rather than
silently running the solver the state did not come from. `overlap`,
`aMb`, `trace`, `inverse_trace`, `operator_norm` and `is_zero_operator`
route on `mode=` directly — note `overlap`/`aMb` default to
`mode="DMRG"`, so handing them ED states without saying so still fails
inside the session rather than with a message.

| Call | Returns |
|---|---|
| `sc.overlap(a, b)` | $\langle a|b\rangle$ |
| `sc.aMb(a, M, b)` | $\langle a|M|b\rangle$ |
| `sc.applyoperator(A, wf)` | $A|\psi\rangle$ |
| `sc.applyinverse(A, wf)` | $A^{-1}|\psi\rangle$ (approximate — see §4) |
| `sc.exponential(h, wf)` | $e^{h}|\psi\rangle$ |
| `sc.summps(a, b)` | $|a\rangle+|b\rangle$ |
| `sc.scale_mps(x, wf)` | $x|\psi\rangle$ (a single-tensor rescale, not an MPO sweep) |
| `sc.trace(A)` | $\mathrm{Tr}\,A$, and `sc.inverse_trace(A)` for $\mathrm{Tr}\,A^{-1}$ |
| `sc.operator_norm(A)` | an estimate of $\lVert A\rVert$; `sc.is_zero_operator(A)` tests $A/c_{\max}$ against $10^{-20}$, $c_{\max}$ its largest coefficient (see below) |

`exponential(h, wf)` really is $e^{+h}$ on the DMRG backends now, as it
always was under `mode="ED"`. It was neither before: for a *multi-site*
`h` — every Hamiltonian-shaped operator, since the symbolic Hermiticity
test reported `False` for $S^x_iS^x_j+S^y_iS^y_j+S^z_iS^z_j$ back then,
which it no longer does — it fell
through to an uncontrolled two-term Taylor truncation (1.5%/16%/247%
relative error at $z=0.25/0.5/1$ on a 4-site Heisenberg chain under
`itensor_version="python"`, unbounded in $z$ and compounding if used in
a loop), and for a single-site `h`,
where the Hermitian branch did run, it computed $e^{-h}$. Anything built
on an anti-Hermitian argument, i.e. real-time evolution written as
`exponential(1j*dt*h, wf)`, was therefore evolving *backwards*. An
operator that is neither Hermitian nor anti-Hermitian now raises
`NotImplementedError` on the DMRG backends instead of printing a warning
and returning a number; use `mode="ED"` on a small chain for that case.

On the operator itself, `A.simplify()` returns it in canonical form, the
same operator with every term's factors sorted by site (carrying the sign
of the fermionic reordering) and equal terms collected, and
`A.is_hermitian()`, `A.is_antihermitian()` and `A.is_zero()` read the
answer off that form. Those three are one-sided: `True` is a proof and
needs no calculation at all, while `False` only means the canonical form
did not collapse, which a Hermitian operator can survive whenever its
Hermiticity rests on something the operator names do not carry, a
same-site identity ($S^xS^x=1/4$), two factors on one site commuting
without both being diagonal, or a term that names both a
pre-Jordan-Wigner fermion (`C`, `Cdag`, `Cup`, ...) and a bare
post-transform ladder operator (`A`, `Adag`, `Aup`, ...), which is left
spelled as written because the sign of exchanging the two depends on
which sits on the lower site. That last rule is new in the 2026-09-24
audit: such terms used to be sorted with the wrong sign, so
`simplify()` changed their value and an exactly anti-Hermitian operator
could be proven Hermitian (`docs/audit_2026_09_24_hole_hunt.md`,
finding 1); no chain class builds one, so no ordinary Hamiltonian lost
its proof. When a chain is at hand,
`sc.is_hermitian(A)` answers the same question without that caveat,
taking the proof when it lands and probing numerically when it does not,
which is what `gs_energy()` and `exponential()` gate on. Both halves decide
on $A/\max|c|$, $c$ running over the coefficients of $A$, so the answer is
a property of the operator and not of its units, and the probe resolves an
anti-Hermitian part down to about $10^{-10}$ of the largest coefficient; on
`mode="ED"` the same question is a relative test on the matrix,
$\lVert h-h^\dagger\rVert_F\le10^{-10}\lVert h\rVert_F$. Until the third
2026-09-24 pass the probe compared an unnormalized norm against an absolute
$10^{-4}$, so a weak loss term on an $O(1)$ Hamiltonian was called Hermitian
and its decay rate dropped (§21). `sc.is_zero_operator(A)` is relative in
the same way since the 2026-09-25b pass: an operator with no terms, only
zero coefficients or an empty canonical form is zero outright, and
otherwise $A/c_{\max}$, $c_{\max}$ its largest coefficient, is probed
against $10^{-20}$ on the mean squared norm, and on `mode="ED"` the exact
matrix is tested as
$\lVert A\rVert_F\le10^{-10}c_{\max}\lVert\mathbb 1\rVert_F$. It used to
threshold a squared norm at an absolute $10^{-4}$, so a nonzero operator
written at $10^{-2}$ read as zero; `cvm_solver="variational"` asks it
whether $A^\dagger=B$ when `canonical.is_dagger_pair` cannot prove it,
before choosing its solver (§6).

A term counts as absent only when its coefficient is at or below
$10^{-12}$ of the largest coefficient of the operator, and in the canonical
form at or below $10^{-12}$ of the larger of that and the magnitudes summed
into it, so an operator means the same thing in any unit of energy while the
rounding dust of $H-H^\dagger$ still cancels. Until the 2026-09-25 fixes it
was an absolute $10^{-8}$, so $10^{-9}S^z_0$ had no terms, its expectation
value and every correlator of it were exactly 0 on every backend, ED
included, and a Hamiltonian written below that scale was the zero operator.
The price of the relative rule is a hierarchy of more than twelve decades
inside one operator: in $10^6+10^{-7}S^z_0$ the field is dropped, which the
absolute rule kept.

On `itensor_version=2` and `3` three numbers inside ITensor are absolute,
calibrated for an operator of order one: `toMPO` discards squared singular
values that sum to its cutoff of $10^{-13}$ and skips a coefficient below
$10^{-14}$, and its Davidson solver replaces a Krylov direction by a random
vector once its residual is below $10^{-10}$. Since the 2026-09-25 fixes a
Hamiltonian whose largest coefficient is below 1 is handed to all three
multiplied by the power of two that brings that coefficient into $[1,2)$,
and the factor is divided back out exactly, so the ground state, the excited
states, the band edges and every KPM spectrum come out the same in any unit
of energy. Before, a Heisenberg chain written at $J\lesssim4\times10^{-7}$
(v3) or $10^{-7}$ (v2) lost its exchange channels and returned the Néel
energy, $-1.25J$ against $-2.4936J$ on 6 sites, and from about $J=10^{-6}$
down every energy lost digits. The largest coefficient is read after
ITensor has merged duplicate terms, the same reading the MPO is built
from, so $sH+S^z_0-S^z_0$ is solved exactly as $sH$; until the 2026-09-25b
pass the solver read the raw term list, ran that operator unscaled, and
returned $E_0/s$ 0.09 to 0.53 off on a 6-site Heisenberg chain at
$s=10^{-10}$. Since then `"python"` builds its operators at the same unit
scale too (a Heisenberg MPO at $s=10^{-7}$ was $1.1\times10^{-9}$ off
relative to itself on numpy with OpenBLAS, and is now at roundoff at every
scale), and NH-DMRG, on v2, v3 and `"python"`, solves $2^kH$ and divides
the energy back in the same way. What one scale cannot reach is a hierarchy
inside the Hamiltonian: `toMPO` truncates bond by bond, so a bond whose
strongest term is below about $4\times10^{-7}$ (v3) or $3\times10^{-7}$ (v2)
of the largest coefficient loses channels at any units, which is an $O(1)$
energy offset or one-site field next to exchange written at $10^{-7}$, or a
weak link $J'=10^{-7}$ in a $J=1$ chain, both read at a third of their
exchange, and `"python"` drops such a bond altogether. Three more routes
are not scale-free, and on those the Hamiltonian is best written in units
where its largest coefficient is of order 1: a Hamiltonian set as an
already-built MPO (`set_hamiltonian(toMPO(H))` on v3, $9.8\times10^{-6}$ off
at $J=10^{-8}$), real-time evolution ($6.4\times10^{-4}$ off with TDVP at
$J=10^{-8}$ on v3, the MPO-Taylor stepper of v2 diverging there, and v3's
Krylov exponentiator stopping on an absolute error; `"python"`'s has been
scale-free since 2026-09-26, see §21), and the Lanczos of `"python"`'s finite DMRG, which stops on an absolute
test ($1.2\times10^{-6}$ relative at $J=10^{-9}$, measured before its MPOs
were built at unit scale). NH-DMRG on v2 was a fourth until the 2026-09-25b
pass. `mode="ED"` does not meet ITensor's thresholds, but until that pass
it had absolute ones of its own: its iterative eigensolver above 2000
states, which now runs at unit infinity norm, `submode="ROOTN"`'s Lanczos
breakdown, now relative to the seed's own $\lVert Hq_0\rVert$ on ED and on
the MPS route alike, and the `dex` cutoff (§6) and `State.normalize()`'s
floor, which stay absolute and now warn. `"julia_live"` is unaffected by
the units.

`trace`, `operator_norm` and `is_zero_operator` take a `MultiOperator`,
not a compiled `StaticOperator`. On the wavefunction itself,
`wf.dot(other)` and `wf.norm()` give $\langle\psi|\phi\rangle$ and
$\lVert\psi\rVert$, `wf.normalize(tol=1e-8)` returns the normalized
state (or `None`, with a warning, if the norm is below `tol`), on every
backend, `mode="ED"` included since the 2026-09-25b pass, where
`normalize()` took no `tol` and `norm()` did not exist; `tol` is an
absolute floor on the norm, so for a state of another scale pass it
relative to that scale. `wf.get_conjugate()` gives the complex conjugate
state, and `wf.get_entropy(b)`
the entanglement entropy at bond `b` (§5). `wf * x` and `A * wf` are
shorthand for `scale_mps` and `applyoperator`. `wf.get_dm(inds=[...])`
returns the single-particle density matrix $\langle
c_i^\dagger c_j\rangle$ and is implemented for `Fermionic_Chain` only.

## 3. Ground-state properties

**Ground-state energy**

$$E_0=\langle\mathrm{GS}|H|\mathrm{GS}\rangle=\min_{|\psi\rangle}\frac{\langle\psi|H|\psi\rangle}{\langle\psi|\psi\rangle}$$

```python
e0 = sc.gs_energy()          # E0
wf = sc.get_gs()             # the |GS> wavefunction itself
```

`get_gs()` takes the same keywords as `gs_energy()` and returns the stored
state under the same condition, which is that the stored answer already
answers the call: a current state, no `wf0=`, and no keyword but
`reconverge=False`, `maxde=None` or `maxdepth=` (on `"julia_live"`, no
keyword at all). Every other call goes to the solver, as on a chain that
was never solved, so `reconverge=True` sweeps from the stored state,
`maxde=` refines it, and a misspelled keyword raises `TypeError`. Until the
2026-09-25b pass a solved chain returned the stored answer before reading
any of them: `gs_energy(maxde=1e-4)` gave the unrefined $-4.1432$ against
$-4.2580$ on a 10-site Heisenberg chain at `maxm=3`, and `wf=x` or
`reconverg=False` were swallowed. `get_gs(wf0=x)` sweeps from x, and
`get_gs(wf0=x, reconverge=False)` returns the unit vector
$x/\lVert x\rVert$ with `e0` = $\langle x|H|x\rangle/\langle x|x\rangle$,
on a solved chain as on a fresh one and on every route, those that resolve
to ED included, and every later reader measures that state. A state and its
multiples are one state: every setter normalizes what it is given, once,
and a state whose norm is zero or not finite raises `ValueError`. Until the
2026-09-25 fixes a solved chain returned its stored state without reading
x, and until the 2026-09-25b pass the readers met a state whose norm is not
one at two different normalizations, and the routes that resolve to ED
raised `TypeError` or dropped x.

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
re-entered with a wavefunction already in hand (`set_initial_wf_guess`,
or simply a previous `gs_energy()` call's own solution; `set_initial_wf`
instead takes the state unswept, normalized), the ramp's starting
bond dimension is floored at whatever that state already carries, so a
re-run can only improve the energy. On `itensor_version="julia_live"` the
setters follow the same contract since the 2026-09-25 fixes: a state set
with `set_gs()` or `set_initial_wf()` is taken unswept with its own energy
$\langle\psi|H|\psi\rangle$, `set_initial_wf_guess()` is swept from, and
`submode="KPM"` measures a set state from its own energy on the band-edge
window. Before, both setters were dropped and the next read solved from a
random start, `set_gs()` raised `AttributeError`, and `gs_energy(wf0=x)` on
a solved chain raised `TypeError`. The ramp applies to the ground-state
solve on all three DMRG backends (`itensor_version` 2, 3 and `"python"`);
`"julia_live"` keeps its own schedule, but since the 2026-09-25b pass a
solve there that has no state of its own starts from a random MPS of link
dimension `min(maxm, bond_ramp_start)` and puts `sc.noise` on the first
half of the sweeps, as the other backends do. It used to start from a
product state with no noise, so a Hamiltonian whose couplings skip a site
stopped at a classical product state with no warning: $-0.5$ against the
exact $-1.0$ for the Heisenberg chain on the even sites of 6, and
$\langle H\rangle=-0.5$ against $-1.0$ for a 3-site `Thermal_Spin_Chain` at
$T=0$, whose physical Hamiltonian lives on the even sites of the doubled
chain. Every `"julia_live"` solve now
varies from run to run within the solver's noise, and a warm 16-site solve
takes 2.1 to 2.35 s where it took 1.5 s. See
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
| `Fermionic_Chain` | `Nf` |
| `Spinful_Fermionic_Chain` (Jordan--Wigner) | `Nf`; also `Sz` under `mode="ED"` |
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
extracted by tuning a chemical potential or a field.

On `itensor_version=3` its cost relative to the unconstrained search
depends on how big the solve is: the quantum numbers restore the block
sparsity the default mode gives up, but they also add per-block
bookkeeping, and only the former scales. Measured on Heisenberg chains
(BLAS pinned to one thread), sector mode is **0.6x** — i.e. slower — at
$n=20$, `maxm=60`; **2.0x** faster at $n=40$, `maxm=100`; and **4.2x**
faster at $n=60$, `maxm=200`, always for the same energy. Expect a small
penalty on toy systems and a real speedup on the ones that actually cost
something. On `itensor_version="python"` the mechanism is different — that backend
stores its tensors densely and uses the quantum numbers only to *label*
basis states, so there is no block sparsity to gain — but a sector is not
a tax either. Measured on the same 40-site Heisenberg chain (`maxm=100`,
15 sweeps, one P-core, BLAS pinned to one thread, 3 repeats):

| | dense | sector |
|---|---|---|
| `itensor_version=3` | 2.85 s | 2.05 s |
| `itensor_version="python"` | 11.04 s | 8.91 s |

so the sector run is ~1.2x *faster* than the same backend's dense run, and
~4.3x slower than `itensor_version=3`'s sector run (about the same 4x
pure-Python-vs-compiled factor as everywhere else — the sector costs
nothing extra relative to it). The pure-Python speedup has a different
cause than v3's: confining the state lets it converge at a lower bond
dimension (chi=52 against 94 for the same energy), which more than pays
for the charge penalty's extra MPO bond dimension. That balance depends on
the sweep count — at 6 sweeps, before the sector run has settled to its
lower chi, it is the slower of the two — so treat "about the same, give or
take 20%" as the honest summary rather than a reliable speedup. What a
sector buys on this backend is the targeting itself: the sector's own
ground state, which a dense search cannot reach at all.

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

Sector targeting is implemented by three solvers: DMRG on
`itensor_version=3`, DMRG on `itensor_version="python"`, and `mode="ED"`
(see below). DMRG on `itensor_version=2` and on `"julia_live"` has no
quantum numbers at all, and a sector-mode chain deliberately raises rather
than letting one of those answer, since a solver without quantum numbers
would silently return the *global* ground state instead.

The two DMRG backends implement it differently, which matters in one place.
`itensor_version=3` confines the calculation structurally: its tensors are
block-sparse over quantum numbers, so an amplitude outside the sector has
nowhere to be stored. `itensor_version="python"` keeps dense storage and
confines the calculation instead by adding a charge penalty
$\lambda\sum_k(\hat Q_k-q_k)^2$ to the operator its *variational* solves
(ground state, excited states, band edges) minimize. That penalty is
identically zero on the target sector, so it changes no reported number —
`gs_energy()` reports $\langle H\rangle$ under the plain Hamiltonian
regardless — and it is not optional: a dense SVD mixes rows across charge
blocks at the $10^{-16}$ level on every truncation, and a variational
sweep amplifies that leak toward whichever sector is lower in energy. With
the penalty switched off, an $n=12$ chain asked for $N_f=2$ with an
attractive interaction converges to the *full* band instead; with it on,
it reproduces sector-restricted ED to $10^{-8}$. If a solve ever does end
up outside the requested sector, it raises rather than reporting the
energy of the wrong one. `chain.set_sector_penalty(lam)` overrides
the default strength (derived from the Hamiltonian's own coefficient
scale) on that backend.

One accuracy caveat, on that backend only: *excited* states inside a
sector are confined exactly (every one of them has the requested charge to
$10^{-8}$) but converge less accurately than on `itensor_version=3`. The
overlap-penalty excited-state solver is already the least accurate part of
the pure-Python backend without any sector, and the charge penalty widens
the local effective Hamiltonian's spectral range on top of that. Measured
on a 6-site Heisenberg chain at $S^z=0$, the second and third levels of
the sector came out $0.2$–$0.4$ above their exact values, where
`itensor_version=3` reproduces them exactly. Ground-state energies,
expectation values, correlators and time evolution are unaffected — those
match sector-restricted ED to $10^{-8}$ on both backends. See
`examples/backend_comparison/sector_v3_VS_python` for both backends' addition
spectra checked against sector-restricted ED, and for the charge-penalty
threshold measured directly.

What works under a sector: ground state and excited states, expectation
values and static correlators, entanglement entropies, real-time
evolution and dynamical correlators via the default
`tevol_method="TDVP"` (and `"TDVP_GSE"`) — all subject to the
conserving-operator rule above. What does not, and says so: METTS
(`metts_vev`, `metts_dynamical_correlator`), and the infinite-chain
algorithms (iDMRG, VUMPS). METTS averages over an ensemble sampled from
every sector at once, so confining the chain cannot restrict it; the
infinite-chain algorithms assemble tensors by hand in a way that is
dense-only under QN indices. They refuse rather than produce a wrong
number. `tevol_method="TEBD"` is the one entry that differs by backend: it
refuses on `itensor_version=3` (its gate assembly sums bond gates before
their fluxes agree, which aborts inside ITensor) and works normally on
`itensor_version="python"`, whose gates are dense and simply inherit the
Hamiltonian's charge conservation. See
`examples/groundstate/sector_targeted_groundstate` for the addition
spectrum and charge gap of a $t$--$V$ chain, checked against ED.

**Sector targeting under `mode="ED"`.** ED implements the same API, and by
the simplest mechanism of the three: a conserved charge is diagonal in
ED's product basis, so a sector is a *set of basis states*, and confining
a calculation to it means taking the corresponding submatrix of every
assembled operator. Same call, same names, same $2S^z$ units, same
answers -- checked sector by sector against `itensor_version=3` in
`tests/test_sector_conservation_ed.py`.

```python
fc = fermionchain.Fermionic_Chain(8)
fc.mode = "ED"                   # this chain is solved by ED
fc.set_hamiltonian(h)
fc.set_conserved_sector(Nf=3)
fc.gs_energy(mode="ED")          # the three-particle ground state
fc.vev(sum(fc.N), mode="ED")     # exactly 3.0
```

Three things about it differ in kind rather than in result:

- **It restricts, it does not shrink.** The full Hilbert space is still
  assembled before the mask is applied, so a sector buys a smaller
  eigenproblem, not a smaller construction. What limits ED is the full
  space either way.
- **An unreachable target can surface at the first calculation** rather
  than at `set_conserved_sector`, because the ED object is built lazily.
  When the chain's DMRG backend carries the same quantum numbers it still
  validates the request immediately, as before; the lazy case is the
  ED-only quantum number below.
- **`promote_to_dense()` is DMRG-only.** It exists because a sector
  re-solve is expensive there; in ED it is not, so the way to leave a
  sector and keep working is `set_conserved_sector()` and re-solve.

ED is also the only solver that can target a total-$S^z$ sector of
`Spinful_Fermionic_Chain`, the Jordan--Wigner spinful chain: its DMRG
representation is $2n$ *spinless* fermionic sites, which carry `Nf` and
know nothing about spin, while its ED representation has an explicit
up/down mode per orbital. Set `chain.mode = "ED"` before asking for it; a
DMRG call on that chain afterwards raises rather than answering with an
`Nf`-only sector.

Because ED has quantum numbers now, a sector also no longer forbids the
automatic DMRG-to-ED fallbacks: a chain that *explicitly* asked for
`itensor_version=2`/`3` on a machine where that extension was never
compiled, or an `itensor_version=3` chain too short for ITensor's
two-site DMRG ($n<3$), answers the sector correctly through ED instead of
refusing. (A chain that named no `itensor_version` at all no longer
reaches that fallback: it is resolved to `"python"` at construction time
when no extension is compiled, and the pure-Python backend implements
sectors itself.)

The conserving-operator rule is the same under ED, and enforced for the
same reason: restricting an operator to the sector is exact for a static
expectation value but identically *zero* for one that changes the charge,
so `vev(fc.C[0])` inside a fixed-$N_f$ sector raises instead of reporting
a clean, wrong zero. `examples/groundstate/ed_sector_addition_spectrum`
sweeps the addition spectrum and charge gap of a $t$--$V$ chain this way.

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
the expensive sweeps, not before. On `itensor_version="python"` the state
was stored densely all along and promotion only relabels its site
indices, so there is nothing to lose by promoting early — but nothing to
gain either. The Hamiltonian and the band-edge
caches are rebuilt on the next call that needs them, while the state the
chain holds is kept — the sector's ground state, or one set with
`set_gs()`/`set_initial_wf()`, taken unswept with its own energy
$\langle\psi|H|\psi\rangle$ on the next read, so a bare `gs_energy()`
afterwards returns that energy rather than re-solving; call `restart()` if
an unconstrained re-solve is what you want. Until the third 2026-09-24 pass
a state set by hand was dropped by the promotion, or re-swept on the dense
sites, which on `"python"` could leave the sector altogether (§21). A wavefunction Python is already holding is not reached by
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

`chain.get_sector_charge_operators()` returns the conserved quantities
this chain offers to `set_conserved_sector`, as a
`{name: MultiOperator}` dict measuring each one over the whole chain —
the same names and the same integer $2S_z$ units. They are ordinary
observables, so `chain.vev(ops["Nf"])` works directly. A chain that
conserves nothing (parafermions) raises.

`promote_to_dense()` is available on `itensor_version=3` and
`itensor_version="python"`, like `set_conserved_sector` itself, and does
nothing if no sector is set. Handing a chain a wavefunction built under a
different sector setting raises instead of quietly contracting the wrong
indices. See
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

`npow=` now means $\langle O^n\rangle$ on the ED route too. Under
`mode="ED"` it used to be accepted and ignored, returning $\langle
O\rangle$ for every $n$ — so any previously recorded ED value for
`npow>1` is a different quantity and is not comparable. `npow=0` returns
`1.0` on both routes. The two guards that were added are **ED-route
only**: on that route a negative `npow` now raises `ValueError` and, at
$T>0$ (`vev(..., T=...)`), an `npow` other than 1 raises
`NotImplementedError` telling you to pass the explicit product operator
instead — where both used to be accepted silently and answer with
$\langle O\rangle$. The DMRG routes are unchanged and still answer
rather than raise: under `itensor_version="python"`,
`vev(Sz[0], npow=-1)` returns $\approx0$ and
`vev(Sz[0], npow=2, T=1.0)` returns 0.25. This matters beyond an explicit `mode="ED"`: it also
reaches the *automatic* ED routes (no compiled C++ extension, or
`itensor_version=3` on a chain with fewer than 3 sites), which nobody
opts into.

The model classes wrap the most common of these profiles, so the loop
above rarely has to be written by hand. Each takes the same `mode=`/
`**kwargs` as `vev` and returns one value per site:

| Call | Returns | Available on |
|---|---|---|
| `.get_density()` | $\langle n_i\rangle$, summed over spin on spinful chains | fermionic chains |
| `.get_density_fluctuation()` | $\langle n_i^2\rangle-\langle n_i\rangle^2$ | fermionic chains |
| `.get_onsite_pairing()` | $\langle c_{i\uparrow}c_{i\downarrow}\rangle$ | spinful fermionic chains |
| `.get_magnetization()` | the $3\times n_s$ array of $\langle S^x_i\rangle,\langle S^y_i\rangle,\langle S^z_i\rangle$ | spin and spinful fermionic chains |

`get_magnetization()` also writes a `MAGNETIZATION.OUT` file into the
working directory as a side effect.

**Energy fluctuation** (a measure of how sharply the DMRG/ED state is an eigenstate, and physically the variance of $H$ in the prepared state):

$$\delta E=\sqrt{\langle H^2\rangle-\langle H\rangle^2}=\big\lVert(H-\langle H\rangle)|\psi\rangle\big\rVert$$

```python
de = sc.gs_energy_fluctuation()
```

It is computed as the right-hand side on every mode: $\langle H\rangle$ is
subtracted first and $(H-\langle H\rangle)|\psi\rangle$ is applied with no
bond-dimension cap, so the truncation, set only by `cutoff`, acts on the small
residual itself. Until the third 2026-09-24 pass it was the left-hand side,
with $H|\psi\rangle$ truncated to the chain's `maxm`, and that number was set
by the truncation rather than by the state: a state solved and measured at
the same `maxm` reported its fluctuation 10 to 51 times too low (6.7e-03
against 0.344 on an 8-site chain at `maxm=3`), and one wider than `maxm`
reported it at order one (§21). `get_gs(maxde=...)`/`gs_energy(maxde=...)`
reads the same quantity **per site**, $\delta E/n_s$, doubling `maxm` until it
drops below `maxde`, and returns the refined energy, the one it leaves on the
chain. It is therefore `gs_energy_fluctuation()/n_s`, not
`gs_energy_fluctuation()` itself, that is compared with `maxde`, and each
step of the loop prints it as `Energy fluctuation per site`.

`mode=` is now actually forwarded here: `gs_energy_fluctuation(mode="ED")`
runs ED, where it used to return the DMRG number byte for byte — so a
value recorded as an "ED energy fluctuation" from before this change is
a DMRG value. Together with the `npow` fix above, every ED route now
evaluates $\sqrt{|\langle H^2\rangle-\langle H\rangle^2|}$ rather than
$\sqrt{|\langle H\rangle-\langle H\rangle^2|}$: on a 4-site Heisenberg
chain `mode="ED"` reports 0.0 where it used to report 2.0561, and a
2-site chain on the *default* backend — which `mode.py` routes to ED by
itself — reports 1.05e-08 where it used to report 1.1456. Anything using
this number as a convergence criterion (`get_gs(maxde=...)`, or a script
that tightens `maxm` until it drops) therefore behaves differently on an
ED route. `npow=` is rejected with a `TypeError`: this function sets the
powers itself.

What the number measures below full bond dimension is the truncation of
the state, which is why it is worth watching while tuning `maxm`. On a
10-site Heisenberg chain at the stock `maxm=30`, one below the full bond
dimension 32, `itensor_version=3` reports 1.4e-05 and
`itensor_version="python"` 8.4e-06, and `mode="ED"` 4.2e-15. The two DMRG
numbers used to read 1.7e-07 and 6.3e-06, which this guide attributed to
backend-dependent floors (Lanczos accuracy on `"python"`, a
double-precision cancellation floor on v3); both were artefacts of
truncating $H|\psi\rangle$ at `maxm`.

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
Afterward, `fc.wf0` holds the eigenstate $|w_g\rangle$ of the
*generalized* problem, not a plain eigenstate of `fc.hamiltonian` alone,
and `fc.e0` that state's own energy
$\langle w_g|H|w_g\rangle/\langle w_g|w_g\rangle$ (on a non-Hermitian chain
the biorthogonal $\langle\psi_L|H|\psi_R\rangle/\langle\psi_L|\psi_R\rangle$
of the pair), while $\lambda$, still what the call returns, is kept as
`fc.lam_generalized` -- every other method that reads `fc.wf0` as an
ordinary ground state
(`get_excited_states()`, any dynamical/KPM correlator, ...) has no way to
tell the difference and will silently build on the wrong reference state
if called afterward. Since the 2026-09-25 fixes this is what every
correlator does, on a Hermitian and a non-Hermitian Hamiltonian alike and
whether or not the chain was solved before. It used to re-solve a plain
ground state over the generalized one at the first correlator on a
Hermitian chain never solved ($|\langle w_g|w_0\rangle|^2=0.60$ on a
6-site chain), and on a non-Hermitian chain whether it was solved first or
not (0.6877 on a 4-site chain), unless a correlator had run before the
generalized solve. A bare `gs_energy()` afterwards returns that energy,
since the stored state is current, and every submode measures the state's
lines from it, as TD and EX always did. Until the 2026-09-25b pass `e0` was
$\lambda$, so KPM, CVM, ROOTN and TDZ put every line $\lambda-E_{w_g}$ away
from where TD and EX put it: for $A=2\,\mathrm{Id}$, whose generalized state
is the ground state itself, the $(S^z_0,S^z_0)$ lines of a 4-site chain in
a field sat at $-0.15$ and $+0.56$ against the exact $+0.66$ and $+1.365$,
a line at negative frequency in a ground-state autocorrelator. Call
`gs_energy_generalized()` as the last step of a calculation, or
`restart()` and then `gs_energy()` if you need a genuine ground state for
one of those other methods.

Unlike plain `gs_energy()`, this method has no ED fallback for a chain
too short for the backend's two-site sweep, so on
`itensor_version="python"` it *raises* `RuntimeError` below 2 sites
rather than returning a number. That is deliberate: the outer
self-consistent iteration would otherwise still return the Rayleigh
quotient of a state no sweep ever touched — measured at $-0.3049$ for
an exact $-0.5$ — i.e. a silently wrong answer. Use `mode="ED"` for a
chain that short.

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

`get_correlation_matrix(dmmode=...)` picks how $C_{ij}$ is evaluated
(`"simple"`, `"fast"`, `"explicit"`, `"full"`). The default,
`dmmode=None`, is resolved from **the state being measured**, not from
the backend the chain was built with, and only after that state exists:
with no conserved sector it is `"fast"` (apply each $c_i$ to the state
and overlap the results — much the cheapest); with a sector, that step
would change the particle number, so it is `"full"` when the
wavefunction carries a live DMRG session handle that can answer the
matrix directly (`itensor_version=3` and `"python"`), and the
backend-agnostic `"explicit"` otherwise — `mode="ED"`, `"julia_live"`,
`basis="Nambu"`, or any state with no session handle. That last case is
why a conserved sector plus `mode="ED"` now answers
`get_correlation_matrix`, `get_correlation_eigenvalues`,
`get_correlation_entropy`, `get_correlation_entropy_density`,
`get_correlated_orbitals` and `get_correlated_density` instead of dying
with an `AttributeError`: the sector case used to hardcode the
session-only `"full"`, on the since-obsolete premise that only DMRG
could reach a sector at all.

A **misspelled** `dmmode`, `basis`, `ctmode` or `fpmode` now raises
`ValueError` naming the argument and listing the accepted values, where
it used to surface as `RuntimeError: No active exception to reraise` —
or, for `basis=`, as no error at all: the old `else` branch served both
`"electron"` and "you typo'd it", so a typo silently returned the
electron-basis matrix. A `dmmode` you pass explicitly is now also
rejected *before* the ground-state solve rather than after it. What has
not changed is that an explicitly-passed mode is still a hard request:
it raises rather than silently falling back when that method is not
available for the wavefunction at hand. The new error only distinguishes
"you misspelled this" from "this mode exists but not here".

`get_four_correlation_tensor(ctmode=...)` has five implementations:
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
single-threaded — see §20 on BLAS threads), i.e. roughly 8x and 19x,
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

`ctmode="batched"` (`itensor_version="python"`, spinless and native spinful
alike) computes the same tensor with every transfer contraction *batched*
over the tuples that share it, and is the default under `"python"` since it
is exact against both `"sweep"` and `"fold"` to machine precision and
15-28x faster. The reorganization is worth understanding, because it fixes
two separate problems at once. `"sweep"` and `"fold"` both evaluate one
tuple's worth of MPS transfer at a time -- a chain of $\chi\times\chi$
contractions -- and there are $O(n^4)$ of them, each far too small to keep
BLAS busy; and the *repeated*-index tuples (those whose four indices are
not pairwise distinct) get no environment reuse at all under `"sweep"`, so
despite being only $O(n^3)$ of them they measured at 61-65% of its total
runtime at $n=12..20$. `"batched"` writes each tuple in site-sorted order as
a sequence of at most four (local matrix, parity) steps; two tuples agreeing
on their first $r$ steps share a partial environment whatever sites those
steps landed on, so the whole $(N,N,N,N)$ tensor collapses onto a trie of a
few dozen matrix sequences, each holding one array of environments batched
over the site combinations that realize it. One site of the sweep is then
two large GEMMs per trie node instead of tens of thousands of small
contractions, and distinct and repeated tuples go through the same
machinery. Measured single-threaded on a spinless chain at `maxm=20`, full
tensor: $n=16$, 10.65\,s (`"sweep"`) / 6.69\,s (ITensor v3, C++) /
0.38\,s (`"batched"`); $n=30$, 117.3\,s / 72.4\,s / 7.4\,s. For native
spinful the comparison is against `"fold"`: 1.29\,s $\to$ 0.22\,s at 5
sites and 3.02\,s $\to$ 0.31\,s at 6, the crossover sitting at 4 sites
(below that the trie build is a fixed cost the tensor is too small to
amortize). The batch axis is also what makes this calculation worth putting
on a GPU at all, and it is the one place in this library where the device
wins *below* the usual chi ~ 120-160 crossover -- because the arithmetic per
array operation comes from the tuple batch rather than from bond dimension.
Measured on an H200 against one Xeon core at $n=30$ (warm, i.e. XLA's ~22\,s
of per-shape compilation already paid): 0.97/0.95/0.98\,s at `maxm` =
20/40/80 against 1.89/9.44/22.81\,s on the host, i.e. 2.0x, 9.9x and 23.3x
-- the device time does not move at all across the range, so every one of
those is a lower bound. Setting `backend.set_backend("jax")` immediately
before the call is enough; the ground state itself can stay on NumPy, since
the MPS is converted once. Full tables and the cold-run caveat in
`docs/gpu_cpu_performance.md`.

`ctmode=None` (the default) auto-selects the fastest method actually
available for the wavefunction's backend/chain type: `"batched"` whenever
it applies (`itensor_version="python"`, either flavour), else `"sweep"`
whenever it applies, else `"full"` whenever it applies (any `itensor_version` in
`2`/`3`/`"python"`, or `Spinful_Fermionic_Chain_Native` under
`itensor_version=3` via its own `four_correlation_tensor_spinful()`),
else the always-correct `"explicit"` fallback (e.g. for
`itensor_version="julia_live"`, which has no `"full"`/`"sweep"`
implementation). Passing a `ctmode` explicitly is still a hard request —
it raises rather than silently falling back if that method isn't
available for the wavefunction at hand. For `Spinful_Fermionic_Chain_Native` the only method that has no
native-spinful counterpart is `ctmode="sweep"`; the default resolver
picks `"batched"` under `itensor_version="python"`, `"fold"` under
`itensor_version=3`, and falls back to `"full"` and then `"explicit"`
if neither is available on the session at hand.

## 6. Dynamical (frequency-dependent) correlators

The quantity returned is the **complex Lehmann density** of the operator
pair $(A,B)$, Lorentzian-broadened by $\delta$:

$$S_{AB}(\omega)\;\equiv\;\sum_n M_n\,\frac{\delta}{\pi\big[(\omega-\Delta_n)^2+\delta^2\big]}\;\xrightarrow[\delta\to0]{}\;\sum_n M_n\,\delta(\omega-\Delta_n)\;=\;\frac{i}{2\pi}\Big[G^R_{AB}(\omega)-G^A_{AB}(\omega)\Big]$$

with $M_n=\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle$,
$\Delta_n=E_n-E_0$, and the retarded/advanced resolvents

$$G^{R/A}_{AB}(\omega)=\langle\mathrm{GS}|A\,\frac{1}{\omega-H+E_0\pm i\delta}\,B|\mathrm{GS}\rangle$$

$\delta$ is the small broadening. Here $|\mathrm{GS}\rangle$ is the chain's
own state, the solved ground state or one set with `set_gs()`, and $E_0$ is
its energy, the one `gs_energy()` reports, on every submode and both modes;
SECTOR, which measures its reference sector's ground state, is the exception
(see its entry). Until the third 2026-09-24 pass DMRG KPM measured a set state
from the solved ground-state energy and `mode="ED"` from three different
origins (§21). $M_n$ is complex in general, so the
returned array is complex. The discriminant that matters is
$\mathrm{Im}\,M_n=0$: **whenever every $M_n$ is real**, the density is
real and coincides exactly with the equally common convention
$-\frac1\pi\mathrm{Im}\,G^R_{AB}$. That is a weaker condition than the
**Hermitian pair** $A=B^\dagger$, which only adds
non-negativity ($M_n=|\langle n|B|\mathrm{GS}\rangle|^2\ge0$). A real
Hamiltonian with real $A$ and $B$ has real $M_n$ whether or not
$A=B^\dagger$: measured on a 6-site Heisenberg chain with $A=S^z_0$,
$B=S^z_2$, $\max_n|\mathrm{Im}\,M_n|=0$ and the two conventions agree to
1.1e-16 — while the density does dip negative there (min $-0.0582$
against a peak $0.1515$), which a Hermitian pair never does.
Every `name=(A,B)` example in this guide **on a Hermitian $H$** has real
$M_n$ — they are Hermitian pairs such as $(S^z_0,S^z_0)$ and
$(A,A^\dagger)$, plus the $S(q,\omega)$ sweep below, whose
$(S^z_i,S^z_j)$ at $i\neq j$ is *not* a Hermitian pair but is still real
— so if that is the case you work in, nothing below changes anything for
you. The qualifier is load-bearing: under a non-Hermitian $H$ (the
`submode="KPM"` example in §9) the $\Delta_n$ are themselves complex and
this whole construction — real $\Delta_n$, orthonormal $|n\rangle$ — does
not apply; that section defines its own quantity. The criterion is not vacuous:
mixing $S^y$ with $S^z$ on the same real Hamiltonian gives a purely
*imaginary* $M_n$ ($\max_n|\mathrm{Re}\,M_n|=0$).

The two conventions differ by the dispersive term:
$-\frac1\pi\mathrm{Im}\,G^R_{AB}=\sum_n[\mathrm{Re}(M_n)\delta-\mathrm{Im}(M_n)(\omega-\Delta_n)]/\pi[(\omega-\Delta_n)^2+\delta^2]$.
The operational reason dmrgpy picks the density is the
kernel-independent sum rule

$$\int d\omega\,S_{AB}(\omega)=\sum_n M_n=\langle\mathrm{GS}|A\,B|\mathrm{GS}\rangle$$

which any submode can be checked against, and which $-\frac1\pi
\mathrm{Im}\,G^R$ does *not* satisfy: its dispersive term has
principal-value tails that leak arbitrarily far outside any finite
frequency window. (Measured on a 4-site chain with complex hoppings and
$A=c^\dagger_0$, $B=c_2$, over `es=np.linspace(-12,18,1500)` at
$\delta=0.05$: the exact $\langle AB\rangle$ is
$0.11052176-0.27135291i$; `submode="KPM"` integrates to that within
2.6e-08 under `mode="ED"` and 7.5e-06 under `mode="DMRG"` on
`itensor_version="python"`, while the
$-\frac1\pi\mathrm{Im}\,G^R$ curve integrates to a real 0.13145.)

Choosing $A=B=S^z_i$ at the same
site gives the local dynamical spin structure factor $S^{zz}_{ii}(\omega)$
(what a local probe like NMR/ESR couples to); choosing $A=S^z_i$,
$B=S^z_j$ at different sites and Fourier-transforming over $i-j$ gives
the momentum-resolved dynamical structure factor $S(q,\omega)$ measured
in inelastic neutron scattering. All of the submodes below compute
$S_{AB}(\omega)$ as defined above — `submode="CVM"` under `mode="DMRG"`
being the one that only recently joined them, so that it now agrees with
its own `mode="ED"` counterpart (measured 1.9e-09 apart on a 6-site
Heisenberg chain — comfortably below the solver's own per-frequency
residual, which prints at ~1e-5 on that run; see §21), and the real-time
pair `submode="TD"`/`"TDZ"` being the one that joined them on
2026-09-22, having returned the complex one-sided Fourier transform
until then (see that section and §21). Sharing the convention is not the
same as agreeing pointwise, since two submodes can still differ in
kernel and in resolution, but for the default `submode="KPM"` they no
longer differ in resolution either: `delta` now sets the same width on
`mode="ED"` as on `mode="DMRG"`, so the two routes sit a fraction of a
per cent of the peak apart instead of a whole peak apart. On a 6-site
Heisenberg chain with $A=B=S^z_0$ at $\delta=0.2$, over
`es=np.linspace(0.01,4,120)`, the two curves sit 1.65e-04 apart on
`itensor_version="python"` against a resolvent (`submode="INV"`) peak of
0.194, and both integrate to 0.24999 against the exact
$\langle AB\rangle=0.25$. The audit that fixed this measured the same
quantity on `itensor_version=3` across three chains and found 4 to 7 per
cent of the resolvent peak, down from 95 to 124 per cent
(`docs/audit_2026_09_hole_hunt.md`, open item O2), and most of that
remainder turned out to be the DMRG routes reconstructing from two more
moments than the calibration asks for. Since they were cut to the
calibrated count (`docs/audit_2026_09_24_hole_hunt.md`, finding 4) the
same three chains sit 0.07 to 0.2 per cent of the peak apart on
`itensor_version="python"`, while `itensor_version=3` agrees as closely
in its median run and reaches about 1.3 per cent (worst of 25 runs on the
staggered 4-site chain) only in the runs where its
band-edge estimate, which moves by up to 2e-2 from run to run, lands
below the true upper edge. What remains between
KPM and a resolvent submode is the kernel itself, a near-Gaussian
Jackson line against a Lorentzian, which the `submode="KPM"` entry below
quantifies.
Otherwise the submodes differ only in
*how*, and therefore in what energy range/resolution/cost trade-off they
offer:

```python
(x, y) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM",
                                      name=(sc.Sz[0], sc.Sz[0]))
```

**What `name=` accepts.** Either a documented correlator string together
with the site indices it applies to (`name="ZZ", i=0, j=1`, and likewise
`"XX"`/`"YY"`/`"+-"`/`"cdc"`/`"densitydensity"`/...), or an explicit pair
`(A,B)` of operators. The pair may be `MultiOperator`s, as above, or the
*already-built* operators `toMPO()` returns (§2): a `StaticOperator`
under `mode="DMRG"`, an `EDOperator` under `mode="ED"`. Products of
those are accepted too, which is the natural way to write a composite
channel such as $c^\dagger_{0\uparrow}c^{\phantom\dagger}_{1\uparrow}$:

```python
Cdagup = [fc.toMPO(A, mode=mode) for A in fc.Cdagup]
Cup    = [fc.toMPO(A, mode=mode) for A in fc.Cup]
A01 = Cdagup[0]*Cup[1]
(x, y) = fc.get_dynamical_correlator(name=(A01, A01.get_dagger()), mode=mode)
```

Handing over an already-built operator saves rebuilding it on every
call, which is what makes it worth doing across a frequency sweep or a
parameter scan. It works for `mode="ED"` (every submode) and for
`submode="EX"`/`"CVM"`/`"ROOTN"` under `mode="DMRG"` — the paths that
consume an operator by *applying* it. The remaining submodes (`"KPM"`,
`"TD"`, `"TDZ"`, `"SECTOR"`, the METTS finite-temperature correlator of
§9, and `cvm_solver="variational"`) instead rebuild the operator inside
the backend from its symbolic term list, so they need the
`MultiOperator` itself and raise a `TypeError` saying so if given
`toMPO()` output.

**`submode="KPM"` — Kernel Polynomial Method.** $H$ is linearly rescaled
into $\tilde H=(H-b)/a$ so its full spectrum lies in $[-1,1]$
(`kpm_scale`, default 0.7, sets how much margin is left inside $[-1,1]$),
then $G_{AB}$'s spectral function is expanded in Chebyshev polynomials
$T_m(\tilde H)$:

$$\mu_m=\langle\mathrm{GS}|A\,T_m(\tilde H)\,B|\mathrm{GS}\rangle,\qquad S_{AB}(\tilde\omega)=\frac{1}{\pi\sqrt{1-\tilde\omega^2}}\left[\mu_0+2\sum_{m=1}^{N-1}g_m\,\mu_m\,T_m(\tilde\omega)\right]$$

with Jackson-kernel damping coefficients $g_m$ suppressing Gibbs
ringing from the truncation at $N$ moments. This is
the default, general-purpose method: robust, works across the whole
spectrum at once, cost grows only linearly with the number of moments.

*What `delta` buys you here.* The Jackson kernel turns an exact pole at
rescaled energy $x$ into a line of width $\sigma=\pi\sqrt{1-x^2}/N$
(Weisse, Wellein, Alvermann & Fehske, [RMP 78, 275
(2006)](https://doi.org/10.1103/RevModPhys.78.275), Sec. II-B), so $N$
and the requested broadening are the same knob read two ways. Every KPM
route in dmrgpy solves that relation for $N$
(`algebra/kpm.py::polynomials_for_broadening`), so the reconstructed
line comes out at FWHM $=2\delta$ at the centre of the rescaled band,
which is the width the resolvent submodes (`"CVM"`, `"INV"`, `"ROOTN"`, `"ED"`)
give for the same $\delta$. The relation holds only when the kernel is
handed exactly $N$ moments, and until the 2026-09-24 audit every DMRG
route reconstructed from $N+2$, its moment loop returning two more than
it was asked for, which made a DMRG line narrower and taller than the ED
one by about $2/N$. The DMRG routes are cut to $N$ now, meaning that
every DMRG KPM spectrum moved by roughly $2/N$ onto the ED one: on two
decoupled Heisenberg dimers the band-centre peak of
$\langle S^z_0;S^z_0\rangle$ went from 0.666761 to 0.620591 at
$\delta=0.2$, on `itensor_version=3` and `"python"` alike, and now
equals the ED peak (`docs/audit_2026_09_24_hole_hunt.md`, finding 4).
`kpm_n_scale` is a multiplier on that calibrated count, default 1, and
it has to be a positive integer (a Python or numpy integer, not a `bool`):
raise it to buy a line sharper than the $\delta$ you asked for, at
proportionally higher cost. Anything else raises a `TypeError` or
`ValueError` naming `kpm_n_scale` on every backend before a single
moment is computed, where `"python"`, `mode="ED"` and `"julia_live"`
used to round it down silently (1.5 gave exactly the calibrated count,
and anything below 1 the 16-moment floor at every $\delta$) while the
C++ backends raised; note that `kpm_n_scale=2.0`, which gave 2x on those
three routes, raises now too (finding 5 of the same record). Both the
calibration and that default are new in the 2026-09 audit's open item
O2 (§21). Before
it, `mode="DMRG"` counted `round((emax-emin)/delta)` moments and
multiplied by a `kpm_n_scale` that defaulted to 3, which put the line at
about $1.6\delta$, while `mode="ED"` used a rule of its own that read
neither of those and landed at about $0.93\delta$, on a Chebyshev window
three times as wide. Two residuals are properties of the kernel rather
than of this choice, and are documented rather than removed. The width
is exact at the band centre and tightens as $\sqrt{1-x^2}$ towards the
band edges, since no single moment count gives one width across a whole
band, and for a ground-state correlator that narrowing is the rule
rather than the exception. Every route rescales the same way, so with
$s$ the `kpm_scale` and $W$ the many-body bandwidth the ground-state
energy $E_0$ sits at $x_0=-1/(2s)$ on every chain ($-0.714$ at the
default $s=0.7$), the band centre is the energy $E_0+W/2$, and an
excitation at $\omega$ sits at $x(\omega)=x_0+\omega/(sW)$. As $W$ grows
with the chain every intensive excitation therefore tends to $x_0$,
where the line is $\sqrt{1-1/(4s^2)}=0.70$ of the requested width, and
on the ground-state-anchored window of `kpm_energy_truncate=True` the
factor at $E_0$ is 0.157. Measured on open Heisenberg chains with
$A=B=S^z_i$, the isolated lowest pole comes out at 0.755 of $2\delta$ on
12 sites and 0.725 on 20, and the delivered width averaged over the
spectrum is 0.87 to 0.90 of the requested one on 8 sites, 0.83 to 0.86
on 12 and 0.78 to 0.80 on 20. To get FWHM $=2\delta$ at a chosen
$\omega$, pass $\delta/\sqrt{1-x(\omega)^2}$ instead, which is the only
compensation available, since `kpm_n_scale` cannot ask for fewer moments
than the calibrated count (`docs/audit_2026_09_24_hole_hunt.md`,
finding 6, and `algebra/kpm.py::polynomials_for_broadening`'s docstring
for the rest of the measurements). And the line is near-Gaussian rather
than Lorentzian, so at equal FWHM and equal integrated weight its peak
stands about 1.5x higher than a resolvent submode's, 1.514 being the
Jackson kernel's own ratio on one pole at the band centre; the 0.315
against 0.195 once quoted here, on a 6-site Heisenberg chain at
$\delta=0.2$ with both integrating to the same 0.250000, is a ratio of
1.61 because its dominant pole sits at $x=-0.527$, where the line is
already 0.850 of $2\delta$. One guard rail sits underneath: a $\delta$ comparable to the
bandwidth itself calibrates to a handful of moments, which stops being a
spectrum at all, so the count is floored at 16 and the line there comes
out sharper than requested instead of unrepresentable.
The band edges entering the rescaling are obtained variationally: the
lower edge reuses the ground-state energy, and the upper edge runs a
deliberately reduced-effort DMRG on $-H$ (few sweeps at modest bond
dimension) — it is only a spectral *bound*, protected by the
`kpm_scale` margin, not a physical result, and a variational
underestimate only shrinks the number of moments. If the bound is ever
too tight for the chosen `kpm_scale`, or `kpm_scale` is below 1/2, where
the ground state itself sits at $x_0=-1/(2s)$, $s$ being `kpm_scale`,
outside $[-1,1]$, the moment recursion raises `RuntimeError` as soon as a
moment exceeds 1.5 times $\|v_i\|\|v_j\|$, the exact bound for a
spectrum inside $[-1,1]$, on every DMRG backend and on `mode="ED"` alike
(`julia_live` raises `juliacall.JuliaError`). Until the 2026-09-24
second-pass audit the DMRG threshold was $10^3(\|v_i\|\|v_j\|+1)$,
which let spectra up to 109 times the true peak through, and ED had no
check at all (findings 2 and 3). This catches every case where a pole
outside the window carries enough weight to show; in a sliver just below
1/2 (0.499 on a 4-site chain at $\delta=0.1$) the elastic line is
already distorted while the moments are still under the bound, so stay
above the ~0.5 floor unless `kpm_energy_truncate` is on. A harshly
truncated `kpmmaxm` (2 or 3 on a small chain) can also cross the bound,
and then its spectrum is tens of per cent wrong anyway.

**Reconstruction kernel (`kernel=`) — how the moments become a
spectrum.** The formula above is only one way to turn $\{\mu_m\}$ into
$S_{AB}$. `kernel="jackson"` (the default) is the classic choice;
`"lorentz"` and `"plain"` are the usual alternatives. Every one of them
damps the moments, $\mu_m\to g_m\mu_m$, which is equivalent to
convolving the exact spectrum with a positive kernel of width $\sim\pi/N$
— and *any* such kernel is only second-order accurate, so the error at a
smooth point decays as $O(N^{-2})$ no matter how many moments are spent.

`kernel="hodc"` selects the high-order delta-Chebyshev reconstruction of
Yi, Massatt, Horning, Luskin, Pixley & Kaye, [arXiv:2512.03149
(2025)](https://arxiv.org/abs/2512.03149), which breaks that
second-order barrier. Instead of damping the moments it replaces
$\delta$ itself by the order-$m$ rational regularization

$$K_\eta(E,x)=-\frac1\pi\sum_{l=1}^{m}\mathrm{Im}\left[\frac{w_l}{E-x+\eta z_l}\right],\qquad z_l=x_l+i,\qquad \sum_{l}w_l z_l^{\,j}=\delta_{j0}\ \ (j<m)$$

— $m$ complex-weighted Lorentzians of width $\eta$, whose weights are
fixed by a Vandermonde moment-matching system that annihilates the first
$m-1$ terms of the small-$\eta$ expansion, leaving $K_\eta\to\delta$ at
$O(\eta^m)$ rather than $O(\eta)$. Expanding $K_\eta(E,\cdot)$ in
Chebyshev polynomials gives energy-*dependent* coefficients
$\nu_m(E,\eta)$ and the reconstruction

$$S_{AB}(\tilde\omega)=\sum_{m=0}^{N-1}\nu_m(\tilde\omega,\eta)\,\mu_m$$

from the **same, undamped** moments. Nothing in the expensive half of the
calculation changes — no extra Chebyshev recursion, no backend change,
and `kernel=` is available on every backend that supports
`submode="KPM"`. `hodc_order` is $m$ (default 6, the paper's own choice;
1–8 are allowed, and $m=1$ degenerates to plain Lorentzian broadening),
and `hodc_eta` is $\eta$ in the same energy units as `delta`/`es`,
defaulting to `delta` itself — i.e. "you asked for resolution `delta`,
you get resolution `delta`, with an $O(\delta^m)$ error instead of an
$O(\delta^2)$ one". Two caveats, both inherited from the method rather
than from this implementation:

* The kernel is **not positive** for $m>2$, so the reconstructed
  spectrum can dip slightly negative next to a sharp feature. That is the
  price of the higher order.
* The gain is an asymptotic statement about *smooth* points. It is real
  where the spectrum is a continuum and disappears at a van Hove
  singularity or a band edge, and it disappears entirely for a spectrum
  made of isolated $\delta$-peaks, which is not smooth at any
  resolution. `examples/dynamical_correlator/hodc_VS_jackson_kernel`
  measures both regimes against an exactly solvable case.
* **It is more sensitive to MPS truncation noise than Jackson**, and this
  one is specific to using it inside DMRG rather than to the method. The
  Chebyshev recursion is only conditionally stable under truncation: the
  error in $\mu_m$ grows roughly exponentially in $m$, and a damping
  kernel's $g_m\to0$ as $m\to N$ happens to suppress exactly the worst
  moments. HODC takes them undamped — its only high-$m$ suppression is
  the $\nu_m\sim q^m\sim e^{-m\eta}$ decay of its own coefficients. So at
  modest `kpmmaxm` the useful moment count is capped by moment accuracy
  well before it is capped by cost, and if an HODC spectrum looks *worse*
  than the Jackson one the first thing to try is a larger `hodc_eta`
  (more damping), not more moments. Measured in that example, on an XX
  chain at `kpmmaxm=32` with $N=1200$ moments whose error has grown to
  $O(10^2)$ by the last of them: at the default $N\eta$ HODC is
  more than 10x *worse* than Jackson, and at $N\eta\approx14$ it is
  several times better (4-8x, drifting between DMRG runs).
  On exact moments the same sweep puts the optimum right at the default,
  which is where it belongs. That default $N\eta$ is whatever the moment
  calibration above implies for $\eta=\delta$, which is
  `JACKSON_FWHM_FACTOR`$/2=3.70$ at `kpm_n_scale=1`; the two
  ratios just quoted were measured before that calibration landed, when
  the same accounting gave 4.3, and they have not been re-measured
  against it. Both values sit in the flat minimum of the
  error-versus-$\eta$ curve, which is why the default was left where it
  is rather than chased. The same example measures the moment
  error directly against exact moments, which is the cleanest way to find
  where that horizon sits for a given chain.

Because the kernel only enters after the moments exist,
`sc.get_dynamical_correlator_moments(...)` returns the raw
$(\mu,E_\mathrm{min},E_\mathrm{max},a,N,\delta)$ and
`kpmdmrg.dynamical_correlator_from_moments(...)` reconstructs from them,
so several kernels can be compared without repeating the DMRG work.

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
Available for `itensor_version="python"` and `itensor_version=3`
(`mode="ED"` raises `NotImplementedError`, and so does `mode.py`'s
fallback onto it for a 2-site `itensor_version=3` chain; it used to
return the untruncated curve, finding 4 of the second-pass record); there
is no `itensor_version=2` port (mpscpp2 has no equivalent machinery to
build a per-site local effective Hamiltonian from, unlike mpscpp3's
`LocalMPO`/`diagHermitian`). Requesting it on `itensor_version=2` raises
`NotImplementedError` (it used to be silently ignored there, so the run
looked truncated and was not — see §21).
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
$b=-\eta B|\mathrm{GS}\rangle$. Second, the CG loop stops early once the
bond-dimension truncation (`cvm_maxm`) keeps it from improving, instead
of burning the full `cvm_nit` iteration budget. The CG starts from zero
($b=0$ returns 0 at once) and stops when its residual reaches `cvm_tol`
(default $10^{-5}$) relative to $\lVert b\rVert=\eta\lVert
B|\mathrm{GS}\rangle\rVert$. Its two early exits read the CG functional
$\phi(x)=\tfrac12\langle x|M|x\rangle-\mathrm{Re}\langle b|x\rangle$, $M$
the system matrix above, which exact CG lowers at every step: `cvm_patience`
(default 50) counts iterations without a new minimum of $\phi$, and
`cvm_blowup` (default 100) fires only when $\phi$ is above its minimum and
the running residual is that many times the best iterate's. The solver
returns the iterate that reached the tolerance or, failing that, the one
of lowest $\phi$. On an untruncated solve neither exit can fire, so
neither changes the answer there; on a truncated one they stop a
recurrence that has stopped improving, and the point warns. Each point
reports its CG iteration count and best residual; if that residual stalls
far above `cvm_tol`, the correction vector is not converged at this bond
dimension and the fix is a larger `cvm_maxm`, not more iterations.
(Warm-starting each point's CG from the neighboring point's correction
vector was tried and measured to *hurt* — truncated CG from a
nearby-but-wrong start can stagnate at a much worse residual than from
the cold start — so each point is solved independently.)

Until the 2026-09-25b pass `cvm_tol` was absolute, the CG started from $b$
itself, and both exits read the running residual, which exact CG does not
lower monotonically. Whenever $\eta\lVert B|\mathrm{GS}\rangle\rVert$ was
at or below $10^{-5}$ (an operator at scale $10^{-4}$, a Hamiltonian in
small units, or `get_kondo_spectrum`'s documented default `delta=2e-6`) the
start passed at iteration 0 and every frequency returned the flat
$\eta\langle AB\rangle/\pi$, and at full bond dimension the exits fired
on the stall that precedes CG's drop and returned the start: an on-line
point of a 6-site chain at $\eta=2\times10^{-3}$ gave
$1.5915\times10^{-4}$ against 13.2636, and 10 of the 121 points of a
10-site grid at $\eta=0.05$ were off by up to 0.989 of the peak. The
earlier claim that the exits never change the answer was therefore false
then; it holds now for untruncated solves. The cost moved with it: that
10-site grid takes 1214 s on v3 where it took 677 s, and a truncated
point about twice as long.

**`submode="CVM_explicit"` — the resolvent formed explicitly.** Computes
the same $S_{AB}(\omega)$ as `submode="CVM"`, but without the
positive-definite reformulation above: it applies
$(\omega+E_0\pm i\eta-H)^{-1}$ to $B|\mathrm{GS}\rangle$ directly, as two
`applyinverse` solves at $\pm i\eta$, and combines them. That is a more
literal transcription of the definition and needs no conjugate-gradient
machinery, but it inherits `applyinverse`'s accuracy rather than
`"CVM"`'s controlled CG residual, so `"CVM"` is the better default; use
this one to cross-check a suspicious `"CVM"` curve. `applyinverse`'s own
tolerance `delta=` is relative to the norm of the vector it inverts
against since the 2026-09-25b pass, so the solve is as accurate for a
small operator as for a large one. It takes any pair: it inverts against
$B|\mathrm{GS}\rangle$ and takes the overlap with
$A^\dagger|\mathrm{GS}\rangle$, which is
$\langle\mathrm{GS}|A(z-H)^{-1}B|\mathrm{GS}\rangle$, i.e. $C[A,B]$, and
the non-Hermitian CVM and INV do the same on DMRG and on ED. Until that
pass it assumed $A^\dagger=B$ behind an absolute test on a squared norm,
so it raised `NotImplementedError` for a non-adjoint pair at unit scale
(and, before that, a bare `RuntimeError` after a `print`) and answered one
at scale $10^{-2}$ with the adjoint pair's curve, 2.004 of the peak off
for $(10^{-2}S^x_0,10^{-2}S^y_1)$. On a non-Hermitian $H$ it is one of
the two submodes (with `"CVM"`) that has a genuine non-Hermitian
implementation.

"The same $S_{AB}$ as `"CVM"`" is a statement that only became true
recently, and it is the one number-changing item here that hits
*everyone*, Hermitian pairs included: `"CVM_explicit"` returned exactly
**twice** the correct value, at every frequency. This submode is
DMRG-only — `mode="ED"` raises `NotImplementedError: submode='CVM_explicit'
has no ED implementation` — but it was the same factor of 2 on every
DMRG backend, since the fix is one backend-agnostic line. Any spectrum,
peak height, integrated weight or figure produced with it before that
fix is a factor of 2 too large.
On the 6-site staggered Heisenberg chain the two now agree to 5.7e-07.
Its output is also plain complex now: on a non-Hermitian $H$ it used to
be passed through `np.abs()`, which discarded the phase and reported
$+|z|$ for values that are genuinely negative or complex.

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

**Both `mode=` values return the house convention**, the complex Lehmann
density $\frac{i}{2\pi}(G^R-G^A)$ defined at the top of this section, and
each gets there by running the fractional-resolvent recursion **twice**,
once at $+i\delta$ and once at $-i\delta$. That is the price of the
convention and it is close to a literal factor of two: `"ROOTN"` was
already the most expensive submode here ($N$ sequential Lanczos
subspaces of dimension `nkry` per frequency, each step under
`mode="DMRG"` a truncated MPO application over the whole chain), and it
now runs $2N$ of them.

There is no shortcut, and it is worth knowing why, because `"CVM"` does
have one. `cvm.py` solves a linear system whose $-\eta$ version is the
same system with the right-hand side negated, so both resolvents fall
out of one solve; root-$N$ applies a *function* of $H$, so the $-\eta$
pass genuinely re-seeds $N$ new Krylov subspaces. Nor does conjugation
help: $\overline{G^R_{A,B}}=G^A_{B^\dagger,A^\dagger}$ is the advanced
resolvent of a *different* operator pair, and it collapses to
$G^A_{A,B}$ only for $A=B^\dagger$ — a case in which the two conventions
already coincide, so the shortcut exists only where it is not needed.

Results produced before this change used $-\frac1\pi\mathrm{Im}\,G^R$
and are **not comparable when some $M_n$ is complex**, which for real
operators on a real Hamiltonian never happens (see the criterion at the
top of this section); where the two conventions coincide the curves
coincide too, up to the recursion's own tolerance.
`mode="ED"` moved first and `mode="DMRG"` followed; on the 4-site
complex-hopping chain with $A=c^\dagger_0$, $B=c_2$ at $N=6$,
`nkry=16`, the `mode="DMRG"` curve moves by up to 5.7e-01 against a
peak of 0.61, and the returned array is genuinely complex there where it
used to be a real float.

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

$$F_{AB}(\omega)=\frac1\pi\int_0^{T}\!dt\;e^{i\omega t}\,C(t)\,w_\delta(t)$$

A real-time run only ever produces $C(t)$ for $t\ge0$, and a one-sided
transform of it is a *resolvent*, not a spectral density: in the
$T\to\infty$ limit $F_{AB}$ is $-\frac{i}{\pi}G^A_{AB}(\omega)$, whose
real part is the density whenever every $M_n$ is real and whose
imaginary part is a dispersive term the density does not have. The
missing half of the transform is not a second simulation backwards in
time. For $t>0$ one has
$C(-t)=\langle\mathrm{GS}|B\,e^{i(H-E_0)t}A|\mathrm{GS}\rangle
=\overline{\langle\mathrm{GS}|A^\dagger e^{-i(H-E_0)t}B^\dagger|\mathrm{GS}\rangle}$,
so the backward half of the pair $(A,B)$ is the conjugate of the
*forward* half of the pair $(B^\dagger,A^\dagger)$, which the same
machinery computes with no change at all, and the two halves combine as

$$S_{AB}(\omega)=\tfrac12\Big[F_{AB}(\omega)+\overline{F_{B^\dagger A^\dagger}(\omega)}\Big]$$

which is what the returned array is. So `"TD"` returns the same density
every other submode returns, and there is nothing left to take the
`.real` of: on a Hermitian pair the imaginary part comes back exactly
zero. Measured on a 4-site Heisenberg chain with $A=B=S^z_0$
($\delta=0.3$, `dt=0.05`, `itensor_version="python"`, resolvent peak
0.1869), max $|\mathrm{Im}\,y|$ is 0.0 and the curve sits 2.97e-04 from
`submode="INV"`, which is the method's own discretization error, with
`mode="ED"` and `mode="DMRG"` agreeing to 1.6e-10 on the same chain.
When $A$ is *provably*
$B^\dagger$ the adjoint pair is the original pair and the combination
collapses to $\mathrm{Re}\,F$ exactly, so that case, which is every
`"TD"`/`"TDZ"` example in this guide, still costs one evolution;
anything else costs two. The test is
`multioperatortk.canonical.is_dagger_pair`, and
it refuses rather than guesses for an operator name with no known
adjoint (a parafermionic `Sig`, a name you invented yourself), which
buys a second evolution. That second evolution is built from
`get_dagger()`, so it is right exactly when `get_dagger()` knows the
name's adjoint, and it did not for the raw backend name `ISy`, which is
$iS^y$ and so anti-Hermitian, until the 2026-09-24 audit gave
`get_dagger()` a phase: before that, `ISy` in the first operator
returned exactly minus the correlator under `"KPM"`, `"EX"`, `"TD"` and
`"TDZ"` (`docs/audit_2026_09_24_hole_hunt.md`, finding 2). Note that the cost
discriminant is not the convention discriminant: what the combination
needs is a *proof* that $A=B^\dagger$, while what makes the imaginary
part vanish is $\mathrm{Im}\,M_n=0$, and the $S(q,\omega)$ sweep above
sits in the gap. Its $(S^z_i,S^z_j)$ at $i\neq j$ has real weights, so
the answer is real either way, and it is still not an adjoint pair, so
under `"TD"` it pays two evolutions per site pair. All of this is new in
the 2026-09 audit's open item O1 (§21), and it moves numbers: a
`"TD"`/`"TDZ"` spectrum from before it is not comparable, its imaginary
part most of all.

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
  spectral function than windowing alone where truncation of the time
  window, rather than the damping, sets the width, which is the standard
  fix in the DMRG-dynamics literature for finite-simulated-time resolution loss,
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
has the numbers). That comparison was made on the old frequency stage,
which evaluated the transform on an FFT grid of spacing
$2\pi/(n_t\,dt)$ and interpolated linearly onto `es`, and at the default
window ($\delta T=6$) most of the narrowing it saw was that grid, which
prediction's tenfold longer series made ten times finer. Since the
damped sum is evaluated at each requested frequency
(`docs/audit_2026_09_24_hole_hunt.md`, finding 7), a 4-site Heisenberg
chain at $\delta=0.05$ shows the line at the exact width with or without
prediction at the default window (FWHM 0.0975 on a 0.0025 grid, as the
exact density gives), and what prediction still
buys there is accuracy, $\max|y-y_{\rm exact}|$ going from 2.55e-03 to
3.3e-06 by carrying the series past the $e^{-6}$ cut; at `nt=200`
($\delta T=1$), where truncation sets the width, it narrows the line
from 0.2125 to 0.0975, the exact width. Note that neither knob changes the *asymptotic*
algebraic decay of the tail far from any resonance — `"KPM"`'s far tail
still decays much faster, since its Jackson-kernel-damped Chebyshev
reconstruction has no Lorentzian/algebraic tail to begin with; these
knobs sharpen `"TD"`'s peaks and suppress its near-tail ringing, not
close that specific gap.

Both `damping` and `predict`/`lp_*` are also accepted by `submode="TDZ"`
and by the internal `S(k,\omega)` reduction (`sxt_to_skomega`), since all
three share the same windowing stage, but their own defaults are
unchanged (`damping="exp"`, `predict=False`), since the empirical
comparison above was only run for `"TD"`. All three also share the last
step: they evaluate the damped sum directly at each requested frequency.
`sxt_to_skomega` stayed on the FFT grid for a while after `"TD"` and
`"TDZ"` left it, and at a converged window, $\delta T=6$, that put its
$S(k,\omega)$ 21 per cent of the peak off exact, with peak heights down
to 0.79 of the right ones; it moved onto the direct sum with the
2026-09-24 second-pass audit (`docs/audit_2026_09_24b_hole_hunt.md`,
finding 7, and §21).

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
scope: only the "greater" branch of the correlator is simulated along
the contour (the same simplification `"TD"` itself already makes), so
this is best used the
same way as `"TD"`: high-resolution work in a narrow frequency window,
now reachable at a lower bond-dimension cost for a given simulated time.

The returned array is the density of §6, assembled from the one-sided
transform exactly as in `"TD"` above and through the same shared code.
The contour does not obstruct that identity: the damping it puts on each
Lehmann term, $e^{-\Delta_n\alpha t}$, is real and is the same for a pair
and for its adjoint, so conjugating the adjoint run reproduces the
backward-time half envelope and all, and a pair that is provably its own
adjoint still costs one run rather than two. How close it lands is set
by the same finite time window as `"TD"`, not by the contour: on the
2026-09 audit's 4-site complex-hopping chain ($A=c^\dagger_0$, $B=c_2$,
$\delta=0.4$, `dt=0.1`, exact peak 0.2313, `itensor_version=3`) it lands
5.4e-04 from the exact Lehmann density, which is where `"TD"` lands at
`predict=False`, while `"TD"` at its default `predict=True` lands
3.2e-05. The contour and the Taylor-in-$\alpha_0$ reconstruction
together account for about 1e-06 of that, the distance between `"TDZ"`
and `"TD"` at `predict=False`, meaning that `alpha0` and `n_max` do not
move the floor. It used to read 3.31e-02, which the 2026-09 records
attributed to the contour and which was in fact the frequency stage's
linear interpolation of an FFT grid of spacing
$2\pi/(n_t\,dt)\approx1.05\delta$, shared with `"TD"` at `predict=False`
(`docs/audit_2026_09_24_hole_hunt.md`, finding 7). `"TDZ"` also takes
only the parameters it names since that audit, so a misspelled keyword
(`alpha=` for `alpha0=`, `nmax=` for `n_max=`) raises `TypeError`, where
it used to be accepted and ignored, returning the default-parameter
spectrum bit for bit, and a `toMPO()` operator raises the `TypeError`
naming `toMPO` again (finding 9 of the same record).

**`submode="EX"` — exact diagonalization in a truncated DMRG subspace.**
Builds $A$, $B$, $H$ explicitly in the subspace spanned by the lowest
`nex` DMRG excited states, then evaluates the exact Lehmann sum in that
subspace,

$$G_{AB}(\omega)=\sum_{n=1}^{n_{ex}}\frac{\langle\mathrm{GS}|A|n\rangle\langle n|B|\mathrm{GS}\rangle}{\omega-E_n+i\delta}$$

i.e.\ a small, explicit sum over poles at the computed excited-state
energies $E_n$, each with residue given by the transition matrix
elements. Cheap and exact *within* the truncated subspace; only as good
as how many/which excited states were computed. Here $|\mathrm{GS}\rangle$
is the chain's own state, the one `set_gs()` set if any, expressed in that
subspace, and $E_n$ is measured from its energy; a state outside the
subspace raises (until the 2026-09-24 second-pass audit, finding 15, EX
measured from the lowest vector of its own subspace whatever the chain
held). The `nex` excited-state
MPS are not assumed to be orthonormal (`dcex.py` builds their own overlap
matrix and solves a generalized eigenvalue problem $Hc=eSc$ rather than
assuming $S=1$), which is important in practice since Gram-Schmidt over
bond-truncated MPS is itself only approximate.

**`submode="SECTOR"` — the Lehmann sum over one quantum-number sector.**
The same explicit pole sum as `"EX"`, but the intermediate states are not
excited states of the whole Hilbert space: they are the lowest
eigenstates of the *one conserved-quantum-number sector that $B|\mathrm{GS}\rangle$
actually lands in*. That sector is not a choice. $B$ shifts every
conserved charge by a fixed amount, so $\langle n|B|\mathrm{GS}\rangle$
vanishes identically unless $|n\rangle$ carries the ground state's charges
plus $\mathrm{charge}(B)$:

| correlator | ground state in | intermediate states in |
|---|---|---|
| $(c_i,c_j^\dagger)$ | $N$ particles | $N+1$ |
| $(c_j^\dagger,c_i)$ | $N$ particles | $N-1$ |
| $(S^-_i,S^+_j)$ | $S_z$ | $S_z+1$ (i.e. `Sz`$+2$ in the 2$S_z$ units the sector API uses) |
| $(S^+_i,S^-_j)$ | $S_z$ | $S_z-1$ |
| $(S^z_i,S^z_j)$ | $S_z$ | $S_z$ |

```python
fc.set_conserved_sector(Nf=6)          # solve at exactly 6 particles
(x, y) = fc.get_dynamical_correlator(submode="SECTOR",
                                     name=(fc.C[0], fc.Cdag[0]), nex=20)
```

Each intermediate state is therefore a separately converged DMRG
eigenstate of a *smaller* Hilbert space. Symmetry — not the overlap
penalty `"EX"` has to rely on — is what keeps them out of the ground
state and out of each other's charge channels: they live in a *different*
sector, so the search cannot collapse back into $|\mathrm{GS}\rangle$.
Within the target sector they are still found by that same penalty
method, only in a much smaller and better-conditioned space, and the
usual per-state fluctuation warning still applies. There is no broadening/resolution tradeoff, no
expansion order, no time window and no Fourier artifact: the poles come
out at their converged energies with their exact weights, and
`return_poles=True` returns them as poles (energies and complex weights)
rather than only as a curve.

The single approximation is truncation of the sum at `nex` states, and it
is *measured* rather than assumed. The returned `info["captured"]` is

$$\text{captured}=\frac{\sum_{n\le n_{ex}}|\langle n|B|\mathrm{GS}\rangle|^2}{\langle \mathrm{GS}|B^\dagger B|\mathrm{GS}\rangle},$$

the fraction of the spectral weight of $B|\mathrm{GS}\rangle$ the returned
poles account for, and a warning fires below 0.9. `nex` is also capped
automatically at the target sector's exact dimension, so asking for more
states than exist is clipped rather than chased. Use this method for a
few sharp low-lying excitations; for a broad continuum the low-lying
states carry a small share of the weight and `"KPM"`/`"TD"` are the right
tools — exactly the opposite tradeoff.

Other keywords: `sector=` and `conserve=` override the reference sector
and which quantities are conserved (by default the chain's own
`set_conserved_sector`, or — with none set — the charges measured on the
unconstrained ground state, restricted to the quantities the Hamiltonian
actually conserves). Two operators whose charges do not cancel raise,
since the correlator then vanishes identically by symmetry.

The $|\mathrm{GS}\rangle$ here is always the ground state of the reference
sector, solved on an internal clone, never the chain's own state. After
`set_gs()`/`set_initial_wf()` the reference charges are therefore read from
the state that was set, and SECTOR proceeds only when that state is its
sector's ground state (by energy, to $10^{-6}$ relative); any other state
raises `NotImplementedError` naming the submodes that read the chain's
state. Until the third 2026-09-24 pass it returned the ground state's
spectrum instead, silently (§21).

An operator with *no* definite charge is handled rather than rejected:
$S_x$ raises **and** lowers $S_z$, so `name=(Sx,Sx)` has no single target
sector, but $S_x=(S^++S^-)/2$ splits it exactly into two pieces that do,
and the channels whose charges cancel are summed — an $S_z+1$ and an
$S_z-1$ contribution, which is also the physically right decomposition.
`info["target_sector"]` is then a list, one entry per channel. The
single-sector primitive `sectordc.sector_poles` still requires a definite
charge, since its `info` names one sector and carries that sector's
matrix elements; use it when you want those (to build $S(q,\omega)$ out
of them, say). `itensor_version=3`
and `itensor_version="python"` only, and never falls back to ED. The
reason is `promote_mps`: this submode solves two sectors on an internal
clone and contracts their states against each other, which needs both
rebased onto the chain's original dense indices. ED targets a sector by
a different mechanism (§3) and deliberately provides no
`promote_to_dense`/`promote_mps`, so there is nothing to fall back *to*
— and falling back to a plain ED solve would answer with the *global*
excited states, a different calculation.

**`get_spectral_function(i, j=None, spin=None, ...)` — the single-particle
spectral function.** The physics-facing wrapper, since "give me
$A(\omega)$" should not require assembling the hole part by hand:

$$A_{ij}(\omega)=\sum_n\langle \mathrm{GS}|c_i|n^{N+1}\rangle\langle n^{N+1}|c_j^\dagger|\mathrm{GS}\rangle\,L(\omega-\omega_n^+)+\sum_m\langle \mathrm{GS}|c_j^\dagger|m^{N-1}\rangle\langle m^{N-1}|c_i|\mathrm{GS}\rangle\,L(\omega-\omega_m^-)$$

with $\omega_n^+=E_n^{N+1}-E_0^N-\mu$ and $\omega_m^-=-(E_m^{N-1}-E_0^N)-\mu$.
The chemical potential $\mu=(E_0^{N+1}-E_0^{N-1})/2$ is subtracted by
default (`shift="mu"`, `shift=None` for the raw axis), which is the
convention that formula is normally written in: particle weight then sits
at $\omega>0$ and hole weight at $\omega<0$, separated by the gap. That
is not cosmetic — in the raw energies both families can land on the same
side of zero, for any Hamiltonian carrying an explicit chemical-potential
term or attractive interactions. `info` also returns `mu`, the charge gap
$E_0^{N+1}+E_0^{N-1}-2E_0^N$, and the two halves separately. On a spinful
chain pass `spin="up"`/`"dn"`; the target sectors are then $(N\pm1,S_z\pm1)$.
See `examples/dynamical_correlator/sector_spectral_function_hubbard`,
where the Mott gap opens with $U$ and $\mu=U/2$ by particle-hole
symmetry.

**`get_spin_spectral_function(i, j=None, ...)` — the dynamical spin
structure factor, resolved by $S_z$ channel.** $S^{+-}$ from the $S_z+1$
sector, $S^{-+}$ from $S_z-1$, $S^{zz}$ from the ground state's own
sector; returns $S^{zz}+(S^{+-}+S^{-+})/2$ (the same combination as
$S^{xx}+S^{yy}+S^{zz}$) and, with `return_poles=True`, the three channels
separately. Keeping them separate is half the point: which channel a
feature lives in is a selection rule no single broadened curve can show.
Note $S^{zz}$'s target sector *is* the reference sector, so its $n=0$
state is the ground state itself and the raw Lehmann sum carries an
elastic $\omega=0$ pole of weight $\langle S^z_i\rangle\langle S^z_j\rangle$
— that is the correct Lehmann sum (the exact ED reference contains it
too) rather than the connected structure factor; pass `connected=True` to
subtract it.

One property specific to this method makes momentum-resolved work cheap:
nothing in the expensive half depends on *which* operators are being
correlated, so one set of sector solves serves the entire matrix of
$(i,j)$ pairs. That is cached, so a full $S(q,\omega)$ or $A(k,\omega)$
map costs one pair of sector solves plus $L$ sets of matrix elements —
measured, a 12-site $S(q,\omega)$ map takes about as long as a single
site's correlator. See
`examples/dynamical_correlator/sector_spin_structure_factor`, which
reproduces the des Cloizeaux–Pearson lower edge of the two-spinon
continuum.

**The isotropic, degeneracy-averaged spin correlator.**
`sc.get_full_SS_correlator(mode="ED", i=0, j=None)` returns the
component-summed correlator

$$S_{ij}(\omega)=\tfrac13\sum_{\alpha\in\{x,y,z\}}S^{\alpha\alpha}_{ij}(\omega)$$

additionally averaged over the degenerate ground-state manifold
(`get_gs_manifold`), which is what an unpolarized INS or ESR measurement
of an isotropic magnet actually sees — neither the component sum nor the
manifold average is a one-liner over `get_dynamical_correlator`. `j`
defaults to `i`. It is **ED-only** and raises for any other mode;
remaining keyword arguments are forwarded to `get_dynamical_correlator`,
so `submode=` and the frequency grid still apply. See
`examples/dynamical_correlator/full_spin_correlator`.

**`submode="maxent"` — maximum-entropy reconstruction.** Reconstructs a
positive-definite spectral function from a finite set of moments
$\langle(H-E_0)^k\rangle$ using a maximum-entropy method
(`distribution.get_distribution_maxent`), rather than a Chebyshev
expansion — useful when positivity of the reconstructed $S(\omega)$
matters more than matching KPM's polynomial-expansion artifacts.
**Not available in a stock checkout:** the reconstruction itself lives
in the third-party `dmrgpy.maxenttk` (PyMaxEnt) module, which is not
distributed with this package, so `submode="maxent"` raises
`NotImplementedError` until you install it separately. The same is true
of `submode="CVMimag"`, whose Padé continuation needs `dmrgpy.padetk`.

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
spectral bound. The pair must come from an NH-DMRG solve: after `set_gs()`,
`set_initial_wf()` or `gs_energy(wf0=x, reconverge=False)` the chain holds a
right state with no left partner, and `submode="KPM"` raises
`RuntimeError`, where until the 2026-09-25 fixes it paired that state with
the left eigenvector of an earlier solve and absorbed
$\langle\psi_L|x\rangle$ (0.93 on a 4-site chain) into the normalization.
`gs_energy(wf0=x)` on a non-Hermitian chain raises `TypeError`, since
NH-DMRG takes no start state, and with `reconverge=False` takes x,
normalized, on every session backend. Implemented so far for the ED backend, `itensor_version=3`,
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
a factor of three of `dex` on either side, where the answer depends on
exactly where the cutoff was placed, and, since the 2026-09-25b pass, a
second one whenever the averaged manifold holds more than one state and is
wider than the broadening, `ex[nex-1] > delta`. The first alone cannot see
a cutoff that lies far above the levels it swallows, which is where the
answer depends on it most: a Hamiltonian written in small units puts its
whole spectrum below `dex/3`, and the call then averages over every state
with equal weight, i.e. returns the infinite-temperature spectrum, 0.738
of the peak off on a 6-site Heisenberg chain with $0.3S^z_0$ written at
$s=10^{-7}$, with no warning before. There the remedy is a `dex` in the
Hamiltonian's own units, `dex=s*1e-5`, which is exact and quiet. A
genuinely degenerate multiplet stays quiet under the second test, since
its width is roundoff.

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
combined with finite-temperature ED, see §9) *and* you have installed the
separate `maxenttk` module it needs. For a non-Hermitian $H$,
the dispatch is per submode, not wholesale: `"KPM"` runs a genuine
biorthogonal KPM (see above), `"CVM"`/`"CVM_explicit"` run the
non-Hermitian correction-vector resolvent (which *is* their
non-Hermitian implementation, not a substitution for it), and
`"EX"`/`"maxent"` are backend-agnostic enough to work as they are.
Every other submode — `"TD"`, `"TDZ"`, `"CVMimag"`, `"ROOTN"`,
`"SECTOR"` — raises `NotImplementedError`. Before the 2026-08 audit
these last ones silently returned the `CVM_explicit` resolvent instead,
so a caller's `submode=` was effectively ignored on a non-Hermitian
Hamiltonian (see §21).

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

`evolve_and_measure` and `evolution_ABA` take the same keywords on both
modes, `nt=1000`, `dt=1e-2` and `h=` (the Hamiltonian to evolve under, the
chain's own by default), and raise `TypeError` on one they do not read. On
`mode="ED"` the propagator is `scipy.sparse.linalg.expm_multiply`, exact to
rounding and unitary, so ED is the reference at any `dt` and whatever the
energy origin of $H$. Until the third 2026-09-24 pass a misspelled keyword
(`DT=0.2`) ran silently at the defaults on DMRG, `h=` raised on ED, the ED
default was `nt=100`, and ED integrated with RK45 at scipy's default
tolerances, whose error grew with `dt` times the absolute energy of the
state: a constant $+20$ added to $H$ moved an ED trajectory $5.6\times10^{-4}$
off exact at `dt=0.1` (§21).

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
  handles a bond-dimension-1 (product-state) start fine, and so does
  `itensor_version=3` since the third 2026-09-24 pass, where a start whose
  site 0 is one local basis vector (a product state, or any ladder
  operator or projector on site 0) used to lose the expansion at the left
  edge and leave site 0 frozen (§21).
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

and tracing out the ancillas gives the thermal (Gibbs) density matrix of
the physical chain,
$\rho(T)=\mathrm{Tr}_{\rm anc}|\Psi(\beta)\rangle\langle\Psi(\beta)|\propto e^{-\beta H}$,
up to the error of the imaginary-time stepper described below:

```python
from dmrgpy import thermal
tc = thermal.Thermal_Spin_Chain(spins, T=0.1)
tc.anneal_step = 2.0    # the default: largest dimensionless step dtau*W
tc.anneal_order = 8     # the default: Taylor order of each step
```

**The imaginary-time stepper.** $e^{-\beta H/2}$ is applied as $n_{\rm st}$
steps, each the order-$k$ Taylor polynomial of
$e^{-\Delta\tau(H-E_{\rm ref})}$, with $E_{\rm ref}$ the middle of the
spectrum $[E_{\min},E_{\max}]$, $\Delta\tau=\beta/(2n_{\rm st})$ and
$n_{\rm st}=\max(1,\lceil\beta W/(2a)\rceil)$, where $a$ is `anneal_step`
and $W=E_{\max}-E_{\min}$ comes from two band-edge solves on a clone of the
chain, on its own backend and mode. Every eigenvalue then enters a step as
$x=\Delta\tau(E-E_{\rm ref})\in[-a/2,a/2]$, so the error depends on the
dimensionless step $a$ and the order $k$ alone, in any units and at any
constant offset, and converges as $a^k$. Each polynomial is applied as its $k$ linear factors over its complex
roots, i.e. $k$ applications of $H$ per step, which runs on every backend.
At the defaults, `anneal_step = 2` and `anneal_order = 8` (class attributes
of `Thermal_Spin_Chain`, which an instance overrides as it sets `T`), the
thermal energy of the open $S=1/2$ Heisenberg chain is within
$6.8\times10^{-7}$ (relative) of the Boltzmann value from 3 to 12 sites and
$T$ from 0.05 to 10, for about $2\beta W$ applications of $H$. The
ground-state branch is taken at $T\le10^{-5}W$, and $T=\infty$ returns the
singlet purification without a step. `thermal.anneal(..., step=, order=)`
takes the same two directly, and its `dbeta=` is an optional absolute cap
on $\Delta\tau$.

Until the 2026-09-25b pass the stepper was first order,
$|\Psi\rangle\to(1-0.1H)|\Psi\rangle$ repeated `int(beta/2/0.1)` times on
the unshifted $H$, so its error was set by 0.1 times the extensive energy
($-2.777892$ against the exact $-3.166396$ on 10 sites at $T=0.5$), a
constant offset changed the answer, every $T>5$ took no step at all (a
`ZeroDivisionError` for a float $T$, the $T=\infty$ state for a numpy
scalar), and it returned early once one step changed the state by less
than $10^{-7}$, so a Hamiltonian in small units came back at the wrong
temperature ($E/s$ of 0.000000 at $s=10^{-3}$ against $-0.769460$ on 3
sites at $T=0.5s$). The ground-state switch was an absolute $T>10^{-5}$.
Finite-temperature results from before are not comparable (§21). The
cost now grows with $\beta W$: 8.6 s where it was 4.9 s on 12 sites at
$T=0.5$. See `examples/finite_temperature/thermal_purification_VS_exact`,
which compares 6 sites against the exact Boltzmann energy from $T=0.1$ to
10, with a constant offset, and plots the relative error.

`tc.MBChain` holds the annealed state as a state set by hand, so `vev()`
and the correlators on `tc.MBChain` measure $|\Psi(\beta)\rangle$ on every
mode, and on the DMRG backends `tc.MBChain.gs_energy()` is the purified
energy $\langle\Psi(\beta)|H|\Psi(\beta)\rangle$ (on `mode="ED"` it stays
the lowest eigenvalue of $H$). Until the 2026-09-25 fixes `gs_energy()`
returned the singlet Hamiltonian's energy ($-2.25$ against $-0.4164$ on a
3-site chain at $T=1$), a correlator re-solved the plain ground state, and
`vev()` and the KPM sum rule on `mode="ED"` read the singlet state
($\langle S^z_0S^z_1\rangle$ 0.0 against $-0.0694$). Those two reference
values are the first-order stepper's of the time; the present stepper
gives $-0.428230$ and $-0.071372$, both exact. Note that a dynamical
correlator of the purified state under $H$ on the physical sites alone is
not the thermal correlator, which needs $H$ minus its ancilla copy, so of
these only the static readers and the sum rule are thermal quantities.

`T=0`, like any $T\le10^{-5}W$, recovers ordinary ground-state DMRG. For small systems (Hilbert
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
than silently falling back to `njobs=1`. It is also unavailable while the
pure-Python backend is running on a device (`backend.set_backend("jax")`,
see the GPU section of `docs/documentation.md`), for a related reason:
the worker processes are started with `spawn`, so each one re-imports
`pyitensor` with the default NumPy backend rather than inheriting this
process's device selection, and the chains would silently run on the host
— so that combination raises too.

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
`topology.berry_phase_matrix(hkgen, nk=20)` computes this from a
callable `hkgen(k)` returning $H(k)$.

**Fermion parity.** For a fermionic chain, the total parity

$$P=\Big\langle\prod_i(1-2n_i)\Big\rangle$$

is $+1$ or $-1$ for a state of definite (even/odd) fermion number and is
the $\mathbb{Z}_2$ invariant distinguishing the two sectors a Majorana
chain's topological phase connects. It is evaluated on the wavefunction:

```python
wf = fc.get_gs()
p = wf.get_fermionic_parity()          # fpmode="full" (default)
p = wf.get_fermionic_parity(fpmode="iterative")
```

The two `fpmode` values are the same quantity by two routes — `"full"`
builds the whole parity string as one operator, `"iterative"` applies it
site by site — and any other value raises. See
`examples/topological/parity` for a sweep of $P$ against chemical
potential across the topological transition, and
`examples/topological/parity_modes` and `parity_long_chain` alongside it.

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

The couplings are read straight off the chain's own Hamiltonian — the
`MultiOperator` you passed to `set_hamiltonian()` — so any spin
Hamiltonian written out of one- and two-site `Sx`/`Sy`/`Sz` terms works,
including anisotropic exchange and bonds that are not
nearest-neighbour. `meanfield.decompose_spin_hamiltonian(sc)` exposes
that split directly, returning the on-site coefficients `b[i,a]`, the
exchange `J[i,j,a,b]` and any constant offset. A Hamiltonian this
decoupling is not defined for — a non-spin operator, a three-site term,
two factors on the same site — raises `ValueError` rather than being
silently approximated. Components that vanish are dropped below $10^{-10}$
of the largest coefficient of the Hamiltonian, not below an absolute
$10^{-10}$, so the decoupling reads the same in any unit of energy; before
the 2026-09-25 fixes a model written below that scale gave an empty
mean-field Hamiltonian and raised, and one at $10^{-9}$ was solved as the
zero operator.

The chain's own one-site terms are kept in full (they are an external
field, not something being decoupled), so at `p=1` the mean-field
Hamiltonian *is* the original model and the solver reproduces its exact
ground state. Keyword arguments are forwarded to
`gs_energy`/`get_magnetization`, so `mode="ED"` picks the ED solver; the
mixing scheme is `mixmode=` (`"default"` or `"broyden"`), since `mode=`
means the DMRG/ED solver here as everywhere else.

Simple mixing on an antiferromagnet oscillates if every site is pushed
the same way each iteration — the standard failure — so for an
antiferromagnetic model start from a staggered `m0` and/or lower `mix`.
The loop stops after `maxite` iterations with a warning rather than
spinning forever.

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
`sc.get_distribution_moments(...)` returns the raw $\langle X^k\rangle$
moments both reconstructions are built from — the distribution
counterpart of §6's `get_dynamical_correlator_moments`. It is a
DMRG/KPM quantity with no ED implementation (the ED path builds spectra
by explicit summation, not from moments), so it raises
`NotImplementedError` under `mode="ED"`. One thing it does *not* share
with §6: its `delta` is deliberately left off the broadening calibration
described there, and still only sets a polynomial count. This path
expands an arbitrary operator $X$ rather than the Hamiltonian, so there
is no rescaled band whose centre the FWHM $=2\delta$ relation could be
anchored to; read `delta` here as a resolution knob and compare
distributions computed at one value of it, not across values. What it
does share with every other distribution is the normalization,
$\int P(x)\,dx=1$, and under `mode="ED"` that held only since the
2026-09-24 audit: the ED route used to return a distribution whose total
weight was exactly $1/$`scale` (0.1000 at the default `scale=10`), and
with `xs=` given its imaginary part was a copy of its real part. An ED
distribution from before therefore moves by exactly the factor `scale`
and is now real on the requested points, a 2-site `itensor_version=3`
chain having reached that route without asking for it through the
automatic ED fallback, while nothing moved on any DMRG backend
(`docs/audit_2026_09_24_hole_hunt.md`, finding 3).

## 15. Post-processing tools

- **Analytic continuation** (`analyticcontinuation.py`): Padé
  continuation of a correlator known on the imaginary/complex-frequency
  axis (e.g.\ from a Matsubara-like or complex-shifted CVM calculation,
  §6's `submode="CVMimag"`) to the real frequency axis, where the
  physical spectral function lives. Requires the third-party
  `dmrgpy.padetk` module, which is not shipped with this package; without
  it both this tool and `submode="CVMimag"` raise `NotImplementedError`.
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
  `get_distribution_maxent`. Requires the third-party `dmrgpy.maxenttk`
  (PyMaxEnt) module, which is not shipped with this package; without it
  importing `reconstruct` raises `ModuleNotFoundError` and the two
  callers above raise `NotImplementedError`.

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
  chain's own `itensor_version`/mode setting. Works at any `T>=0`. The
  third-order sums run over the thermally occupied *initial* states only
  (those with $p_i>10^{-12}\max p$; one state at $T=0$, or the whole
  degenerate ground-state manifold with equal weights), so they cost
  $O(n_{\rm occ}\,{\rm dim}^2)$ per bias point rather than the
  ${\rm dim}^3$ they used to, and $F(\epsilon,T)$ is tabulated once per
  call (0.4 s) instead of being integrated per point: at 64 states and
  1 K a third-order term takes 0.05 s where it took 28 s and 4 GB. It
  reads only its named parameters, and since the 2026-09-24 audit any
  other keyword raises `TypeError` naming it, where it used to be
  accepted and ignored, so that `Jrho=` for `Jrho_s=` silently removed
  the whole third-order Kondo peak (`docs/audit_2026_09_24_hole_hunt.md`,
  finding 11); the moment count of its KPM correlator is set through the
  chain's `kpm_*` attributes, not through a keyword.
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
temperature-broadened logarithmic function $F(eV-\epsilon_m,T)$, one
for each of its two diagrams with different arguments (see below), that
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

**Direct and exchange diagrams.** Every third-order process comes in
two interaction orders (the paper's "normal" and "reversed"), and
reversing the order both reverses the electron-spin trace and turns
the intermediate electron from electron-like into hole-like, whose $F$
enters with the opposite sign. For the Kondo term the trace is the
antisymmetric Levi-Civita one, so the two sign flips cancel and both
orders add, $F(eV-\epsilon_{im})+F(eV-(\epsilon_f-\epsilon_m))$ (the
second argument is not the paper's, see "Where the exchange log sits"
below); for the potential-interference term the trace is the symmetric $\delta_{kj}$, so
only the hole-like flip survives and the two orders *subtract*,
$F(eV-\epsilon_{im})-F(eV+\epsilon_{im})$ — which is what the paper's
Fig. 7c shows (its 121u and 121uR curves are mirror images of opposite
sign). Two consequences: the elastic intermediate state $m=i$ drops out
of the potential term identically, and the whole term vanishes at $B=0$,
so the zero-field Kondo peak keeps its position and symmetry at $U\neq0$
(Fig. 7d at 0 T is exactly symmetric, peak 1.39 at $eV=0$).

> **Note (behaviour change, 2026-09-12).** Until then the potential
> term used the *summed* combination. That kept a spurious $m=i$ term,
> a $\mathrm{sign}(eV)\,F(eV)$ spike at zero bias that shifted the
> zero-field $U=0.25$ peak to $-0.2$ mV (1.47 instead of 1.39) and
> tilted its tails (1.156/1.052 instead of 1.145/1.135), and it gave
> the exchange diagram the wrong sign, so the 10 T step asymmetry came
> out $\sim0.27$ against the figure's $\sim0.05$ (now 0.047 against a
> digitized 0.054). Every `U!=0, order=3` result from before that date
> is not comparable; `U=0` and `order=2` are untouched.

**Where the exchange log sits, a departure from the paper.** The direct
diagram tunnels first: its intermediate state is the tunnelled electron
in the sample with the impurity in $m$, which goes on shell at the
sample's Fermi edge when $eV=\epsilon_m-\epsilon_i$, hence
$F(eV-\epsilon_{im})$. The exchange diagram scatters a sample electron
first, into the outgoing state, and leaves a hole behind; that goes on
shell when the hole reaches the Fermi edge, which fixes the *outgoing*
electron's energy at $\epsilon_i-\epsilon_m$, and the outgoing electron
has $eV-\epsilon_{if}$. So the exchange log is
$F(eV-(\epsilon_f-\epsilon_m))$. The paper's eq. 25, and the formula
under its Fig. 6, has $F(eV-\epsilon_{mi})=F(eV+\epsilon_{im})$ instead,
measured from the incoming energy as if the impurity absorbed nothing.
The two agree when $f=i$, so the potential term, whose $I_{fi}$ forces
$f=i$, is the same either way, and whenever every state the initial one
connects to is degenerate with it (a free $S=1/2$ at $B=0$, a Kramers
doublet). Any field, anisotropy or exchange splitting separates them.
`tests/test_kondo_spectrum_tmatrix.py` builds the second-order T-matrix
of the paper's own Hamiltonian as explicit fermion $\otimes$ impurity
operators, reading every intermediate state's energy denominator off
$H_0$, and the spectrum agrees with it to roundoff on anisotropic $S=1$
and exchange-coupled impurities in tilted fields, where the printed form
missed it by 5–20% of the third-order term.

> **Note (behaviour change, 2026-09-26).** Every third-order Kondo
> spectrum with an inelastic transition out of an occupied state moved,
> on `mode="ED"` and `mode="DMRG"` alike; a free $S=1/2$ at $B=0$, every
> zero-bias value, `order=2` and the potential term did not. On the
> paper's own Fig. 7d parameters at 10 T the step overshoots go from
> 1.226/1.183 to 1.248/1.203 and the $\pm4$ mV tails from 1.143/1.130 to
> 1.149/1.137, where the figure, drawn from eq. 25, reads 1.231/1.177 and
> 1.146/1.128; the step asymmetry, which is the potential term, does not
> move. At $T=1$ K on an anisotropic $S=1$ ($D=1$ meV, $E=0.3$ meV, 3 T)
> the largest move is 0.044 on a peak of 2.52. Results from before are
> not comparable.

**Scope and known limitations**, worth reading before trusting specific
numbers:

- Only a single chain site couples to the tip (`site=`); the paper's own
  model allows several sites with independent tip couplings, not
  implemented here.
- The paper's own closed-form equations for two numerical building
  blocks — the temperature-broadened step $\Theta(x)$ and the
  temperature-broadened Kondo log function $F(\epsilon,T)$ — are
  garbled as printed (checked directly: the printed $\Theta(x)$ diverges
  rather than saturating, and the printed $F$, its eq. 22, has a log that
  does not depend on the integration variable and an overall sign that
  makes it negative). `kondospectrumtk/stepfunctions.py` uses the
  evident intent of each instead: $\Theta$ is re-derived from the
  paper's current formula, and $F$ is the closed-form log
  $\ln[(\omega_0+|\epsilon|)/|\epsilon+i\Gamma_0|]$ convolved with the
  thermal kernel $\Theta'$, i.e. $F_0(\epsilon)=\ln(\omega_0+|\epsilon|)-\tfrac12\ln(\epsilon^2+\Gamma_0^2)$
  at $T=0$. Both are verified against digitized values from the paper's
  own figures: the six peak heights of its Fig. 5 (arXiv v1 numbering)
  and the $\pm4$ mV tails of its Fig. 7b (see that module's docstring,
  `examples/kondo/kondo_spectrum_VS_paper/` and
  `tests/test_kondo_spectrum_paper_fig7.py`).

  > **Note (behaviour change, 2026-09-12).** Until then $F$ was the
  > paper's *electron-like* defining integral (its eq. 20) evaluated
  > exactly, and used for the exchange diagram too, which the paper says
  > is the hole-like eq. 21. That form is not even in $\epsilon$ (its
  > band edge sits at $+\omega_0$ only, with a sharp-cutoff log
  > singularity there), and it put the Fig. 7b tails at 0.854 where the
  > figure reads 0.886. The thermal broadening is unchanged (eq. 20's
  > own double integral reduces to the same $\Theta'$ kernel: 7.459 vs
  > 7.460 at the Fig. 5 peak), so nothing changes at zero bias; every
  > third-order number away from it moves by $O(|eV|/\omega_0)$, up to
  > 3.5% of the total at $|eV|=\omega_0/5$.
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
  every step). $G$ carries the direct diagram; the exchange diagram, whose
  log sits at $eV-(\epsilon_f-\epsilon_m)$, needs that difference as a
  frequency, which $G$'s $t_2$ axis does not have, so it comes from the
  sheared function $G_x(s,\tau)=G(-s,\tau+s)
  =\langle\psi_l(-s-\tau)|S_k|\psi_j(-s)\rangle$, with
  $\psi_j(t)=e^{-i(H-E_0)t}S_j|\mathrm{GS}\rangle$ the $t_2$ trajectories
  already computed (one more short $\tau$ trajectory per $l$ per row,
  1.5x the cost). Both are then extracted via two closed-form time-domain
  kernels (derived by inverse-Fourier-transforming $\Theta_0$ and
  $F_0(eV-\cdot)$) rather
  than by evaluating those functions pointwise on a discrete frequency
  grid, which does not converge robustly for this construction — see
  `kondospectrumtk/twotime.py`'s module docstring for the full
  derivation and the numerical pitfalls it was built to avoid ($\Theta_0$'s
  kernel is a Cauchy principal value, computed via an FFT-based Hilbert
  transform for machine-precision accuracy; $F_0$'s, a sine/cosine-
  integral expression since 2026-09-12, is log-singular at $t_2=0$
  because $F_0$ decays only as $1/|\epsilon|$, and the grid points
  within eight cells of it are taken as cell averages). This is the expensive part:
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

What "the ground state" means at a degeneracy differs between the two
modes unless you say which one you want. `mode="ED"` takes the
$T\to0^+$ limit, the equal-weight average over the degenerate manifold,
while `mode="DMRG"` uses the single state its solver converged to, which
at an accidental crossing is whichever member, or superposition, the
random start reached, so the spectrum lands anywhere between those of
the members: on an $S=1$ impurity at the crossing of its $|0\rangle$ and
$|-1\rangle$ levels, between 1.0 and 2.0 against the 1.5 of the average.
Passing `n_gs=g` asks `mode="DMRG"` for that same average. Every term,
second order, third-order Kondo and potential, is then averaged with
equal weight over `get_excited_states(n=g, purify=True)`, which is
exactly the ED average, since every term is linear in the ground-state
density matrix, and it costs $g$ times the single-state run; on that
crossing `n_gs=2` gives 1.50002 to 1.50004 in every run. The default
`n_gs=1` is the unchanged single-state path, a `RuntimeWarning` fires
when a member lies more than `delta` from $E_0$, and `n_gs>1` runs on
`itensor_version` 2, 3 and `"python"`, measures each member exactly as
returned, with `gs_energy()` that member's own energy, restores the chain
afterwards, `gs_energy()` included, and accepts only the submodes that
read the chain's state (KPM, CVM, CVM_explicit, ROOTN, TD, TDZ, EX),
raising `NotImplementedError` for SECTOR, any other submode and
`"julia_live"` (`docs/audit_2026_09_24_hole_hunt.md`, finding 12, and
the 2026-09-24 second-pass record, findings 13 to 15: until then each
member was re-swept before it was measured, and on v2 and v3 an excited
member split below `delta` relaxed into the lower one in 16 of 24 runs).

Further `mode="DMRG"` limitations beyond the general ones above: the
`es` frequency grid has no safe default and must be supplied
explicitly, and what it has to cover is set by the chain's spectrum
rather than by the `eV` sweep range, and differs between the terms. The
second-order term needs only the transitions below $\max|eV|$, while
the potential-interference term needs every transition up to the top of
the tip-coupled site's $S_k$ spectrum, a few $\delta$ beyond it and
$\omega=0$ itself, since its kernel falls off only as about $2eV/\omega$
inside the band. One `es` is shared by both, so with `U!=0` it must meet
the second requirement, and the potential term issues a `RuntimeWarning`
when the sum rule $\sum_k\int S_{kk}(\omega)\,d\omega=S(S+1)$ says more
than 1e-2 of the weight is missing (see `second_order_dIdV_dc`'s and
`potentialdc`'s docstrings). The grid may be non-uniform and in any
order, since both terms sort it and then integrate with trapezoid weights
of the grid they are given (an unsorted grid, such as a refined block
appended after a coarse one, was silently wrong until the 2026-09-24
second-pass audit, 3.04 off a 6.97 peak on the second-order term,
finding 17); the
potential term weighted every point with the first spacing until the
2026-09-24 audit (finding 13), and the one example of it,
`examples/kondo/kondo_potential_term_dmrg_VS_ED`, violated the coverage
requirement and blamed the resulting 9.7 per cent DMRG-versus-ED gap on
KPM broadening, where with `es` spanning $\pm20$ meV it is 2.1 per cent,
measured before the moment-count cut of finding 4 (finding 14). The third-order term's `dt2`,
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
(`submode="KPM"`, $\delta=2\times10^{-5}$) agrees with the exact
excited-state sum to 0.2% of its maximum at every bias point on the
same chain. That figure predates the 2026-09-22 KPM broadening
calibration (§21) and has not been re-measured against it; what the
second-order term reads is a cumulative integral of $S(\omega)$, which
the sum rule pins independently of the moment count, so the
recalibration is not expected to move it. It was quoted as "a few tens of percent at thresholds"
until 2026-09-12; that error was the route's own cumulative sum, which
counted the whole frequency bin holding a threshold's delta-like peak
as lying below it (1.033 against an exact 0.808 at $eV=0$), not KPM --
it is a trapezoid rule now.

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

Every method above works on a *finite* chain of `n` sites. `Infinite_Many_Body_Chain`/`Infinite_Spin_Chain` (`infinitechain.py`) instead describe a translationally-invariant chain that repeats a fixed-size **unit cell** forever in both directions, and solve it with infinite DMRG (iDMRG, White's growing algorithm generalized to a multi-site unit cell) rather than sweeping a fixed-length system. A unit cell of any size is supported under `gs_method="vumps"` (the default), which routes `n_uc>2` to the sequential multi-site solver `pyitensor/vumps_ms.py`; `gs_method="idmrg"`, the growing algorithm, is limited to `n_uc<=2` and raises `NotImplementedError` above that — see `pyitensor/idmrg.py`'s module docstring for why.

`itensor_version="python"` (default) or `3` (`ic = infinitechain.Infinite_Spin_Chain(["1/2"], itensor_version=3)`) are both supported — `2` and `"julia_live"` have no iDMRG port and also raise `NotImplementedError`. `ic.gs_method` (`"vumps"` by default, on EITHER backend since 2026-08-08 — see below) picks the ground-state solver. Both solvers, on both backends, support `vev`/`correlator` (see "Static correlators" below): the v3 C++ backend's `gs_method="idmrg"` (`mpscpp3/chain_session.h`'s `Chain::idmrg_ground_state`) is a line-by-line port of `pyitensor/idmrg.py`'s own growing algorithm, and the three things the Python side had gained since that port — McCulloch's wavefunction prediction, the gauge-consistent unit-cell extraction described below, and the per-site energy baseline — are ported too, along with the static observables built on them (`Chain::idmrg_onsite_expectation`/`idmrg_two_point_correlator`/`idmrg_local_excitation_gap`). Relative speed depends on the model: v3's own `idmrg_ground_state` is ~2.6-2.7x slower than `"python"`'s for a gapless (critical) chain (the local 2-site solve needs close to its full Krylov dimension every macro-iteration), but can be *faster* for a gapped model, where that same solve converges quickly — see `docs/documentation.md`'s iDMRG section for benchmark numbers.

**Building the Hamiltonian: `L`/`C`/`R`-suffixed operators.** Instead of one operator list per absolute site index, an infinite chain exposes each operator three times per unit-cell site `i` (`i=0..n_uc-1`): `SxC[i]` (site `i` of the *central* cell — ordinary intra-cell use), `SxR[i]` (site `i` of the *next* cell), and `SxL[i]` (site `i` of the *previous* cell, provided purely so a coupling can be phrased in whichever direction reads most naturally). Terms of any finite range are accepted; each is canonicalized by translating it, as a whole, onto the cell its leftmost site lives in, so no bond is ever double-counted between a cell and its neighbour:

```python
from dmrgpy import infinitechain
ic = infinitechain.Infinite_Spin_Chain(["1/2"])       # n_uc=1, uniform chain
h = ic.SxC[0]*ic.SxR[0] + ic.SyC[0]*ic.SyR[0] + ic.SzC[0]*ic.SzR[0]  # NN Heisenberg
ic.set_hamiltonian(h)
```

A bare, single-site operator (e.g. `ic.SzC[0]`, no product) is also a valid term — a Zeeman-field-style onsite Hamiltonian, either on its own or added alongside bond terms.

**Couplings longer than one unit cell: `get_operator(..., group=c)`.** The three flat lists reach one cell either way. For anything further out, `ic.get_operator(name, i, group=c)` takes an *integer* cell offset — `c = -1, 0, 1` are exactly `L`, `C` and `R`, and `c = 2` is the cell after the next one — so a J1-J2 (next-nearest-neighbour) chain fits on a one-site cell:

```python
ic = infinitechain.Infinite_Spin_Chain(["1/2"])       # n_uc=1
J1, J2 = 1.0, 0.4
h = 0
for op in ("Sx", "Sy", "Sz"):
    o0 = ic.get_operator(op, 0)                        # site 0 of this cell
    h = h + J1*o0*ic.get_operator(op, 0, group=1)      # nearest neighbour
    h = h + J2*o0*ic.get_operator(op, 0, group=2)      # next-nearest neighbour
ic.set_hamiltonian(h)
print(ic.gs_energy())
```

Both ground-state solvers handle this, on both backends. `gs_method="vumps"` (the default) routes a Hamiltonian whose couplings exceed one unit cell to the *sequential* multi-site solver (`pyitensor/vumps_ms.py`, or `Chain::vms_ground_state` on `itensor_version=3`), whose channel-resolved environments carry one channel per site of a term's reach and so cost linearly in it — where the grouped path it otherwise uses would need the same chain rewritten on an `n_uc >= range` cell and folded into a `d**range` supersite, i.e. exponentially. `gs_method="idmrg"`'s growth loop consumes the same automaton and handles it too, needing no dispatch at all. `kpm_finite` follows along (it builds its own finite window).

One thing does not follow along. `excitation_energies`/`excitation_gap` raise `NotImplementedError` on either backend: the tangent-space ansatz is built on the reach-1 environment triple, so there is no sequential route for it to take — the ground state on that same chain still works, so rewrite it on a longer unit cell if you need the excitations. `vev`/`correlator`, by contrast, do follow along, on both backends and at any reach: `itensor_version=3` reads whichever snapshot the run left behind (`Chain::vms_onsite_expectation`/`vms_two_point_correlator` for a sequential answer, the grouped `vumps_*` pair otherwise), so nothing has to be rewritten to measure a long-range chain.

See `tests/test_infinite_long_range.py`, which pins a polarized chain carrying reach-2 and reach-3 `Sz`-`Sz` terms against its exact energy density (and its exact `vev`/`correlator`), and `examples/idmrg/long_range_infinite_chain`, which sweeps `J2` and checks the 1- and 2-site cells against each other.

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

**What these cost in memory** (`itensor_version="python"`). Every fixed point and every correlator above is a contraction along the chain's transfer matrix, and that chain used to be assembled as one `chi^2 x chi^2` array per unit-cell position, so the peak memory of a run grew as `chi^4` and was reached before any physics happened. Since 2026-09-22 the chain is carried lazily instead, as the `(chi_l, d, chi_r)` site tensors it is built from, and applied one site at a time, so nothing on that path allocates the rank-4 array at all. Measured on a 2-site Heisenberg cell at `maxm=64` under `gs_method="idmrg"`, threads pinned and seeded so the growth trajectory is identical either way: the tracemalloc peak of the first `vev` went from 517.8 MB to 5.8 MB and of the growth loop from 519.0 MB to 10.4 MB, peak RSS from 810.8 MB to 304.1 MB, and the growth loop itself from 23.92 s to 5.23 s, with `e0`, `vev` and an `r=1..7` correlator sweep identical to every printed digit. Nothing here changes a number, so read it as a budget you can spend on a larger `maxm` rather than as a result to recheck; what is left composing a full transfer matrix is the dense eigensolve small chains take anyway and the fallback for a non-converged ARPACK solve, so a small cell behaves exactly as it did. `gs_method="vumps"` goes through the same lazy chain once per iteration and is a wash in time at `D=16` (16.39 s to 16.11 s on a critical Heisenberg cell, against a run-to-run spread on that box several times larger), the win there being the memory alone.

**The default ground-state solver: VUMPS** (`ic.gs_method = "vumps"`, the default since 2026-08-08 — Variational Uniform Matrix Product States, Zauner-Stauber et al., arXiv:1701.07035; see `pyitensor/vumps.py`'s own module docstring for the algorithm) — instead of growing a finite window and truncating it down to `maxm` at every step (the `gs_method="idmrg"` growing algorithm above), VUMPS solves directly, in the thermodynamic limit, for the actual `maxm`-dimensional variational optimum (`ic.maxm` sets VUMPS's own target bond dimension `D` here too). Both `itensor_version="python"` (`pyitensor/vumps.py`) and `itensor_version=3` (`mpscpp3/chain_session.h`'s `Chain::vumps_ground_state`, a C++ port of the same algorithm — built from plain dense arrays closed over LAPACK rather than ITensor tensor-network objects, since the bond/physical dimensions this feature targets are always small; see that method's own doc comment) support `gs_method="vumps"`; the two are cross-checked directly against each other to ~1e-10 or tighter on TFIM/Heisenberg at `D=1,2,3` (`tests/test_vumps_v3.py`). Explicitly set `ic.gs_method = "idmrg"` instead for `local_excitation_gap`/`td_dynamical_correlator` (no VUMPS equivalent, see their own sections below) or if the growing algorithm's own more battle-tested behavior is preferred (VUMPS's former `D>1` convergence-robustness gap has since been traced to two bugs and fixed — see the reliability note below):

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

`vev`/`correlator` also work under `gs_method="vumps"`, on BOTH backends: `pyitensor.vumps.onsite_expectation`/`two_point_correlator` for `itensor_version="python"`, and `Chain::vumps_onsite_expectation`/`vumps_two_point_correlator` (a line-for-line C++ port of the same formula) for `itensor_version=3` — cross-checked directly against each other to ~1e-14 or tighter on TFIM at `D=2,3` (`tests/test_vumps_correlator_v3.py`). Both are computed directly from the converged mixed-gauge `{AC, AR}` rather than `pyitensor.idmrg`'s dominant-right-fixed-point eigenproblem: `AC` is already the exactly-normalized single-(super)site reduced state by construction of the mixed canonical gauge (Vanderstraeten, Haegeman, Verstraete, "Tangent-space methods for uniform matrix product states", arXiv:1810.07006, Eq.(34)), and `AR`'s exact right-orthonormality lets a two-point correlator spanning multiple unit cells close by a direct trace with no eigenproblem either (the mixed-gauge analogue of that same review's Eq.(37)-(39)). Unlike the growing-algorithm's own reconstructed-from-the-last-macro-iteration correlators above, these carry no `maxiter`/`<H_uc>`-self-consistency caveat — VUMPS solves directly at the target bond dimension in the thermodynamic limit, so once `.converged` is `True` the correlator is exact for that converged `{AL,AR,C}`, only limited by the bond dimension `D` itself (same caveat as the energy density, see below). `local_excitation_gap` is one of the two methods (with `td_dynamical_correlator`, below) that still require `gs_method="idmrg"` specifically (it re-diagonalizes the growing algorithm's own final 2-site effective Hamiltonian, which has no VUMPS equivalent), on either backend at `window=0` (`Chain::idmrg_local_excitation_gap` for `itensor_version=3`); its `window>0` variant is `itensor_version="python"`-only, being an explicit prototype rather than stable API; conversely `excitation_energies`/`excitation_gap` (below) *require* `gs_method="vumps"` (the default) — they need `VUMPSResult`'s own mixed-gauge `{AL,AR,C,GL,GR}`, which the growing algorithm's `IDMRGResult` (`gs_method="idmrg"`) has no equivalent of.

```python
ic.gs_method = "vumps"
ic.vev("Sz", 0)                    # <Sz> at site 0, from the converged VUMPSResult
ic.correlator("Sz", 0, "Sz", r)    # <Sz(0) Sz(0+r)>, same signature as gs_method="idmrg"
```

**A scoped reliability note.** `D>1` VUMPS was for a while genuinely unreliable here — independent calls at the same `D` could land on noticeably different energies, and a `D=4` TFIM run occasionally missed the exact answer by ~10%. That is no longer the state of the code: it was two identified bugs, both since fixed. The first was an environment fixed point closed against a missing conjugate, invisible at `D=1` and a real source of wrong `D>1` energies; with it fixed, the same `D=4` TFIM(g=1.5) case converges to ~1e-7 relative on 10/10 independent `nrestarts=6` calls. The second was the stopping criterion of the inner eigensolves: they stopped on the Ritz *value*, which leaves an eigen*vector* accurate only to the square root of its tolerance — fatal for VUMPS, whose convergence test compares `AC` and `C`, two independently-solved eigenvectors. That floored the gauge mismatch at ~1e-6, so `tol=1e-10` was unreachable at any number of iterations while the reported energy was perfectly good. Both backends now stop on the residual instead: `D=8` TFIM went from `converged` in 0/3 runs to 3/3 (and 32.6s to 6.7s at g=1.5, 38.4s to 4.6s at critical g=1.0), with the energies unchanged to every printed digit. What remains is ordinary restart-search difficulty rather than a known defect: the D-ramp, multi-restart and variational-principle safety-net machinery is load-bearing infrastructure, not vestigial, and for a `D>2` result on a harder or less-tested model it is still worth calling `gs_energy()` a few times independently and keeping the lowest reported density. Always check `.converged` — it is reported honestly and never silently assumed. See `pyitensor/vumps.py`'s own "Convergence robustness" docstring section for the full numerical account. See `examples/idmrg/vumps_TFIM/main.py` for a worked example sweeping `D` and comparing VUMPS against both `gs_method="idmrg"` and the transverse-field Ising model's own exact (free-fermion) energy density.

Entanglement/entropy are not implemented for infinite chains yet.

**Excited states: the tangent-space/quasiparticle excitation ansatz.** Unlike a finite chain, an infinite chain's excitations form a momentum-resolved band `E(k)`, not a single discrete state, so "the first excited state" here means the standard single-mode/quasiparticle ansatz (Haegeman et al.), built on top of a **VUMPS** ground state's own mixed-gauge `{AL,AR,C}` representation (`ic.gs_method = "vumps"` — required, see above): a tangent-space vector with one excitation tensor `B` (a function of a free matrix `X`) inserted at every unit-cell position, weighted by momentum `k`. `ic.excitation_energies(k, n=1)` returns the lowest `n` excitation energies (above the ground state) at momentum `k` (radians, per unit cell); `ic.excitation_gap(ks=None)` scans `k` (default `numpy.linspace(-pi, pi, 41)`) and returns the minimum — the scalar "gap", mirroring the finite-chain `get_gap()` naming (§4):

```python
ic.gs_method = "vumps"             # required for excitation_energies/excitation_gap
ic.excitation_energies(0.0, n=1)   # lowest excitation energy at k=0
ic.excitation_gap()                # min_k E(k), the scalar gap
```

**Scope.** Both `itensor_version="python"` (`pyitensor/idmrg_excitations.py`) and `itensor_version=3` (`mpscpp3/chain_session.h`'s `Chain::vumps_excitation_energies`, a C++ port of the same algorithm, same dense-array/LAPACK approach as `vumps_ground_state` above) are supported — both require `gs_method="vumps"` (`NotImplementedError` otherwise, or if a different `itensor_version` is used). Any converged bond dimension `D>=1` is supported — including a genuinely entangled ground state (`D>1`, e.g. the transverse-field Ising or a dimerized Heisenberg chain), which used to be an explicit, rejected scope limit here (see `pyitensor/idmrg_excitations.py`'s own module docstring, "History" section, for the eight-pass investigation that limit came from and how it was eventually resolved by rewriting the ansatz from scratch on top of VUMPS's own mixed-gauge state, mirroring MPSKit.jl's own architecture). For `D=1` the computed dispersion matches the exact free-fermion single-magnon dispersion of a field-polarized XX chain to ~14 digits across the whole Brillouin zone (`examples/idmrg/excitation_gap_xx/main.py`); for `D=2` it matches an independently-converged MPSKit.jl transverse-field Ising result to 6 significant figures (`examples/idmrg/excitation_gap_tfim/main.py`), and `H_eff(k)` is Hermitian to machine precision at every `D` tried. `itensor_version=3` is cross-checked directly against `itensor_version="python"` across a full momentum scan, matching to ~1e-10 or tighter on the gapped TFIM case and to a looser ~1e-4..1e-7 on the gapless/critical Heisenberg case (both backends' own non-convex VUMPS restart search can land on slightly different local optima there, not a discrepancy in the ported algorithm itself) — see `tests/test_vumps_excitations_v3.py` and `examples/idmrg/vumps_excitation_v3_VS_python/main.py`. A longer-range term (spanning more than 2 adjacent unit cells after `n_uc`-grouping) is rejected with `NotImplementedError` by `excitation_energies`/`excitation_gap` themselves: the ansatz's environments are the reach-1 `{GL, GR, bond channels}` triple. The ground state on that same chain is unaffected — `gs_energy()` routes it to the sequential multi-site solver instead of raising.

**Cost.** The eigenproblem `excitation_energies` solves has dimension `D*D*(d_g-1)`, so its cost grows quickly with the converged bond dimension. Two things bound that, and since 2026-09-22 they bound it on **both** backends rather than only on `itensor_version="python"`: the momentum-dependent channel resolvents are built once per momentum and cached with their LU factorization (they depend only on `k` and the momentum-independent environment, so rebuilding and refactorizing them inside every application of `H_eff(k)` was pure waste), and above a size threshold the eigenproblem is solved by Lanczos on that Hamiltonian's action rather than by assembling the matrix at all. Neither changes any returned number — the two solver paths are cross-checked against each other in `tests/test_infinite_chain.py` and `tests/test_vumps_excitations_v3.py`, agreeing to 4.6e-11 on `itensor_version=3` — and a momentum scan at large `D` is substantially cheaper than it was: over a 3-momentum scan on `itensor_version=3` with threads pinned to one core, a `D=16` TFIM chain went from 164.1 s to 10.3 s with the resolvent cache and to 1.8 s with both halves, and an `n_uc=2` Heisenberg chain at `D=10` from 15.6 s to 5.4 s to 1.2 s. The two backends keep different thresholds, and deliberately: `itensor_version=3` switches to Lanczos above dimension 64 rather than `"python"`'s 256, since one application of `H_eff(k)` there solves four channel resolvents, which puts the crossover below dimension 36. One design note if you ask for `n>1` on a cell whose dispersion is degenerate away from `k=0`, as the `n_uc=2` Heisenberg cell's is: a single Krylov space holds at most one direction out of a degenerate eigenspace, so both solvers run one deflated Lanczos per eigenvalue rather than one Krylov space for all of them. A plain single-vector Lanczos there returns distinct eigenvalues where two should coincide, every one of them a genuine eigenpair, so a residual check passes it. This guide used to say that `"python"` never had that exposure, and it did: above its dense threshold a single ARPACK call for all `n` from one constant start dropped one copy of a degenerate level and returned the next distinct one in its place, 7 of 16 calls wrong over four momenta on the `n_uc=2` critical Heisenberg cell at `maxm=10` (dimension 300), the earlier "not exposed" having been measured at dimension 12, where the Krylov basis is the whole space. Since the 2026-09-24 audit it runs one deflated Lanczos per value for `n>=2`, as the C++ solver does, and agrees with its dense path to 1.9e-15 there, while `n=1` is the unchanged single call (`docs/audit_2026_09_24_hole_hunt.md`, finding 15).

**Spectral weights: `S(k,w)` directly in the thermodynamic limit.** Each branch the ansatz returns is a genuine momentum eigenstate, so its contribution to a dynamical correlator is an exact δ-peak — not a broadened, windowed approximation to one. `ic.spectral_weights(opname, k, p=0, n=1)` returns that peak's position *and* its residue:

```python
ic.gs_method = "vumps"                                  # required, as above
e, w = ic.spectral_weights("Sx", 0.7, n=1)              # energies and weights
e, w, total = ic.spectral_weights("Sx", 0.7, n=1, return_total=True)
ks, es, S = ic.dynamical_structure_factor("Sx", delta=0.05)   # broadened map
```

with `w[a] = |<k,a|O(k)|Psi>|^2` for the normalized quasiparticle state `|k,a>` and `O(k) = N^{-1/2} sum_m e^{ikm} O_m`, i.e.

    S(k, w) = sum_a  w[a] * delta(w - e[a]).

`opname` is a local operator on sub-site `p` of the unit cell, and `k` is per unit cell. This is the same physical object `kpm_finite`/`td_dynamical_correlator` estimate, computed a completely different way: those embed a *finite* window in the infinite chain and therefore carry a window-size error plus a KPM/time-truncation broadening, while these peaks are exact in both momentum and energy. What they do not contain is multi-particle continuum weight — which is exactly what `return_total` measures, so the two families are complementary rather than one superseding the other. `ic.dynamical_structure_factor(opname, ks=None, energies=None, delta=0.05, p=0, n=1)` is the convenience wrapper that assembles those peaks into a Lorentzian-broadened `(k, w)` grid ready to plot as a heat map; the broadening there is cosmetic only.

**The sum rule, and what `return_total` is for.** A one-site operator applied to a uniform MPS lands *exactly* inside the same variational tangent space the excitations live in (`B = O.A` is itself a valid excitation tensor), so the total weight summed over **every** branch is exactly the per-site connected static structure factor `sum_r e^{ikr}(<O_0 O_r> - <O>^2)`. `return_total=True` returns that number — free, since it is the squared norm of the same source vector each weight is an inner product against — and `w.sum()/total` is then the fraction of the momentum-resolved response the branches actually returned account for. That is the practical test of whether a single-mode picture is adequate at that momentum, and it is a sharp one: on the paramagnetic transverse-field Ising chain, `sigma^x` (parity-odd, exactly one quasiparticle) puts >99.9% of its weight on the lowest branch, while `sigma^z` (parity-even) puts *exactly zero* there — measured at ~1e-21, an exact selection rule — and all of it on higher branches lying inside the exact two-particle continuum. See `examples/idmrg/dynamical_structure_factor_tfim/main.py`, which plots all three.

The connectedness is automatic: no `<O>` is subtracted anywhere, and none needs to be, because the disconnected piece lives entirely along the gauge direction the left gauge-fixing projector annihilates (see `pyitensor/idmrg_excitations.py`'s `_spectral_source_vector`/`_spectral_resolvent`). The `k=0`, `<sigma^z> = -1` product-state case — where the disconnected part is 1, not small — is pinned by a test.

**Validation and scope.** Four independent references, since a spectral weight has no single golden number. The sharpest is the **AKLT** point of the spin-1 chain (`H = S·S + (1/3)(S·S)²`), whose ground state is *exactly* a `D=2` MPS — so there is no variational error at all, and three closed forms are reproduced to machine precision (~5e-15): the static structure factor summed from `⟨S_0·S_r⟩ = 4(-1/3)^r`; Arovas–Auerbach–Haldane's single-mode dispersion `(5/27)(5+3cos k)` [PRL **60**, 531 (1988)], which the *first moment* must reproduce exactly since `S^z_k|Ψ⟩` lies entirely inside the tangent space; and the SU(2) content, which emerges rather than being imposed — the eight branches split into a magnon triplet and a quintuplet at every momentum, `S^z` reaches only the triplet (quintuplet weight ~1e-23), and the two multiplets cross near `k≈0.9`. On the pure **Haldane** chain the gap at `k=π` converges monotonically from below to the literature 0.4104789 (0.2132/0.4074/0.4094/0.4098 at `D=4/8/12/16`) while the sum rule holds identically at every `D` — including `D=4`, where the gap is 50% wrong, which is the point: it is an identity of the ansatz, not of the ground state's accuracy. The magnon triplet exhausts 97.3–97.6% of the `k=π` sum rule and falls to ~68% at small `k`. See `examples/idmrg/haldane_structure_factor/main.py`.

The other three, on spin-1/2: the exactly solvable `J=0` Ising product state (`sigma^x`/`sigma^y` weight exactly 1 at every momentum, `sigma^z` exactly 0); the static sum rule above, checked against a real-space sum of `ic.correlator` — machinery sharing no code with the excitation ansatz — matching to ~1e-13 at `D=2`; and the f-sum rule `sum_a e[a] w[a] = (1/2)<[O_k^dagger,[H,O_k]]>`, which brings the excitation *energies* in too and converges with the ground state's own bond dimension (2e-4 relative at `D=2`, 1e-7 at `D=4`). See `tests/test_infinite_chain_spectral.py`. Requires `itensor_version="python"` — unlike `excitation_energies`, the mpscpp3 port is energies-only, so it has neither the eigenvectors nor the mixed-transfer source vector a weight is built from, and that pair is now the whole of what is missing there: the solver work `excitation_energies` was waiting on (the cached, LU-factored channel resolvents and the Lanczos eigensolve) was ported on 2026-09-22, so a dispersion scan costs the same on both backends and only the weights are still one-sided (`docs/idmrg_improvement_plan.md`) — plus `gs_method="vumps"` and `n_uc <= 2` as above. A fermionic (parity-odd) `opname` is rejected: its Jordan-Wigner string would have to be closed at infinity, the same reason `correlator` rejects an odd-parity operator pair.

**Cost.** Independent of how many branches are asked for, and of `D`, `d_g` and the operator: exactly **two** linear solves per momentum, both cached on the chain across repeated calls at the same `k`. The naive way to build the weight — leaving the excitation tensor's legs open and solving once per basis element — would cost `D^2*d_g` solves per momentum instead; transposing each geometric series moves its solve to the operator end of the contraction instead. So the eigensolve for the branches themselves dominates, and adding weights to an existing dispersion scan is close to free.

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

`n_window` has no default — see the convergence caveat below. `window_chain_kwargs` is an optional dict of attribute overrides applied to the temporary finite `Many_Body_Chain` (`maxm`, `nsweeps`, `kpmmaxm`, `kpm_scale`, ...), independent of `ic`'s own `maxm`/etc.; a key the temporary chain does not hold as a setting (a misspelling, a method or a private name) raises `TypeError`, and so do `itensor_version` and `mode`, since the window always runs on `itensor_version="python"` (until the 2026-09-24 second-pass audit a misspelled key was stored where nothing reads it and the call returned the default spectrum, finding 6); remaining `**kwargs` (`delta`, `kernel`, `es`, ...) are forwarded to `kpmdmrg.get_dynamical_correlator` unchanged.

**Scope restriction — a finite-window approximation, read before use.** This is *not* an exact infinite-size method: results carry finite-size/open-boundary corrections that must be checked by convergence in `n_window`, exactly as a static `vev`/`correlator` caller would check `maxm`/`etol` convergence of the original iDMRG ground state. One Chebyshev moment corresponds to one application of the (nearest-neighbor) window Hamiltonian, so it can only move information by ~1 site per moment (a Lieb-Robinson-style bound) — but KPM's own moment count scales with the *window's own extensive bandwidth* divided by the requested `delta` (an ordinary finite chain's KPM already has this property, nothing new here), so a genuinely fine `delta` can require a moment count comparable to (or larger than) `n_window` itself, at which point open-boundary reflections contaminate the result regardless of how large `n_window` is. Prefer a coarser `delta`, or check that the correlator has visibly converged with growing `n_window`, for quantitative work (especially near a gapless point, where a fine `delta` is most tempting). Unlike `vev`/`correlator`, this does not need `ic._result` (no dependency on a previously converged `IDMRGResult`, or even on `ic.itensor_version`), so it works regardless of which backend `gs_energy()` itself used. See `examples/idmrg/dynamical_correlator_finite_window/main.py` for a worked example sweeping `n_window`.

**A genuinely infinite-chain dynamical correlator, via real-time TDVP (`td_dynamical_correlator`).** `kpm_finite`'s own open-boundary window has a real error source no amount of `n_window` alone can fix: an open chain's own ground state carries boundary artifacts (e.g. Friedel-oscillation-like features) that contaminate even the *central* region, not just the two edges. `ic.td_dynamical_correlator(opname_i, p_i, opname_j, n_window, dt=0.1, nt=200, x_values=None, maxdim=60, cutoff=1e-10, niter=50, connected=True, **kwargs)` fixes this by capping the window's two ends with the *converged* iDMRG growth environment (`idmrg_ground_state`'s own `HL`/`HR`, already computed during growth and exposed on `IDMRGResult`) instead of plain open boundaries — infinite boundary conditions (IBC), following Milsted/Vanderstraeten et al., "Infinite boundary conditions for response functions and limit cycles in iDMRG" (arXiv:1804.09163) — and evolves the perturbed window in real time via two-site TDVP rather than expanding in Chebyshev moments. Supports both `itensor_version="python"` and `itensor_version=3` (calls `gs_energy()` automatically if needed, like `vev`/`correlator`):

The native ITensor v3 backend (`Chain::td_dynamical_correlator_window`, `mpscpp3/chain_session.h`) reuses the vendored ITensorTDVP library's own boundary-tensor `tdvp(psi,H,t,LH,RH,sweeps,args)` overload directly against a tiled window MPS/MPO — unlike the `"python"` backend, which has to hand-roll its own window-aware TDVP sweep (pyitensor's generic TDVP infers a site's Link via a same-Index chain-neighbor lookup that cannot see a window's extra boundary legs). Two scope differences versus `"python"`: (1) `x_values` may not extend beyond the window's own explicit range (`center+x` must stay within the window, i.e. increase `n_window` instead of relying on padding) — the `"python"` backend pads beyond the window with extra unevolved unit-cell copies, not ported here; (2) as of this writing, `itensor_version=3`'s own `idmrg_ground_state` has a known, pre-existing convergence bug for Hamiltonians with an onsite ("field") term (energy diverging every macro-iteration) — unrelated to `td_dynamical_correlator` itself, but it means a v3 `td_dynamical_correlator` call inherits that limitation for such models; a purely bond-coupled Hamiltonian (e.g. plain Heisenberg) is unaffected.

```python
ks, es, Skw = ic.td_dynamical_correlator(
        "Sz", 0, "Sz", n_window=14, dt=0.05, nt=200,
        maxdim=60, cutoff=1e-10, x_values=range(-6, 7),
        ks=np.linspace(-np.pi, np.pi, 41), delta=0.1, window=[-1, 6])
```

`opname_j` is applied at sublattice position `p_i` and evolved forward in time under the window's own ground-state-energy-shifted Hamiltonian (`e^{-i(H-E_GS)t}`, matching this codebase's own established real-time correlator convention, e.g. `mpscpp3::quench_tdvp`'s `Hshift=H-EGS*Id`); `opname_i` is inserted at the shifted position (bra side, not itself evolved) — this is the paper's own headline efficiency result: every `x` (and, via the spatial Fourier transform below, every `k`) comes from this *one* window evolution, not one run per distance the way a naive real-time approach (or `kpm_finite`'s own one-run-per-`r` KPM calls) would need. Cross-checked against an exact non-interacting (free-fermion) reference for the XX chain, on systems far larger than many-body ED could reach — see `docs/documentation.md`'s own architecture-level notes on this method for the full derivation and what that check found. `connected=True` (default) subtracts the disconnected background `<opname_i><opname_j>` before the spatial Fourier transform — turning it off produces a spurious, dominant `k=0` contribution with no discernible dispersion (the raw correlator approaches `<opname_i><opname_j>`, not 0, at large separation). `Skw` (shape `(len(ks), len(es))`) is obtained via a spatial DFT (`S(k,t)=sum_x e^{-ikx}S(x,t)`), conjugated, and followed by the *same* damping and direct per-frequency sum (`delta` -> Lorentzian broadening) `submode="TD"` uses. The conjugation is the finite TD route's own conjugate-at-return: `S(x,t)` carries $e^{-i(H-E_{GS})t}$, and without it every line came out at $\omega=-D_n$ instead of $+D_n$ on both backends, mostly below the default `window=[-1,10]`, until the 2026-09-24 second-pass audit (finding 8). There is one deliberate difference at the end: what comes back is the raw one-sided transform, not the §6 density. The reduction that produces `S(k,t)` never names the operator pair, so the adjoint-pair identity `"TD"` completes its transform with (see that submode above) has nothing to act on here. Take the real part if you want the density, which is the right thing to do whenever the Lehmann weights are real, and note that this is the one dynamical-correlator route in the library where that instruction still stands.

**Fermionic operators, and the vacuum normalization** (both since 2026-08-29). A parity-odd pair (`"Cdag"`/`"C"`, either order, on any fermionic site type) computes the *physical* fermionic correlator: the Jordan-Wigner string is threaded across the window on the ket (before the evolution, so it is evolved along with the perturbation) and across the bra at measurement time, matching the convention `correlator` already used for the static case. The two are pinned to each other by an exact identity — `S(x, t=0)` equals `correlator(...)` to machine precision at every `x`, both signs, across the window's own padding — and against the exact free-fermion Green function `<c†_x(t) c_0> = Σ_l [e^{iht}]_{xl} P[l,0]` at `t>0` (`tests/test_idmrg_window_fermionic.py`, `examples/idmrg/fermionic_dynamical_correlator/main.py`). Note the anticommutation sign this implies: for two parity-odd operators at different sites `<A_x B_0> = -<B_0 A_x>`, so a fermionic `S(x,t)` and a static `correlator` written in the opposite site order differ by a minus sign. A pair with odd *total* parity (one fermionic operator against a parity-even one) raises rather than returning a number, exactly as `correlator` does: its string can never close. `connected=True` subtracts nothing for a parity-odd pair, whose disconnected background is zero by symmetry. Before this, the path applied a bare `C`/`Cdag` matrix with no string on either backend — a different number entirely, not a less converged one (measured: `+0.203` at `x=2` against an exact `-0.001`). Independently, every `S(x,t)` — spin included — is now divided by the vacuum amplitude `<ψ|ψ(t)>` measured through the identical contraction, which cancels a spurious global factor that had been inflating results on both backends (the shift `eshift` that keeps the evolution phase-stationary is measured with the window's boundary legs traced, while correlators are measured with them closed by the transfer-matrix fixed points, and the two see different energies). On the dimerized XX chain that alone took the residual against the exact free-fermion answer from ~0.07 to ~1e-5.

**Both backends are exact on the `S(x,0)` identity.** `S(x,t=0)` must equal the chain's own static `correlator(...)` exactly (no evolution has happened yet), and both `itensor_version="python"` and `itensor_version=3` satisfy it to ~1e-15. That was not always true: until 2026-08-29 the v3 window tiled the raw per-micro-step iDMRG factors rather than the gauge-consistent unit cell every other v3 static observable uses, and missed the identity by up to ~1.7e-1 on a plain Heisenberg chain whose energy density the two backends agreed on to 6.7e-11 — shape and finiteness looked fine throughout, which is how it went unnoticed. Fixing it also removed this backend's `exp(+i*eshift*t)` phase correction in favour of the exact vacuum normalization the Python backend already used, and a latent index-stride bug that only showed up when the converged unit cell is not square. See `examples/idmrg/td_dynamical_correlator_python_VS_v3/main.py`, which checks the identity at every `x` on both backends and then compares the two `S(x,t)` trajectories.

**Scope**: this simplifies the paper's own Eq. 7 to `t1=0` (the ground state is perturbed by `opname_j` and evolved forward only; `opname_i` is never itself time-evolved) rather than the full two-branch trick (evolving a *second*, independent window backward in time too, which doubles the accessible total time for the same TDVP cost) — a documented, straightforward follow-up. `n_window` and `x_values` are two *separate* convergence axes to check (not one): `n_window` controls how much environment margin surrounds the perturbation, `x_values` how far the spatial sum/Fourier transform reaches — growing `x_values` together with `n_window` can keep changing results (a slowly decaying connected-correlator tail keeps contributing as the range widens), so converge each independently. See `examples/idmrg/td_dynamical_correlator/main.py` for a worked `n_window`-convergence sweep and a cross-check against `kpm_finite`.

**Applying an operator/gate to the converged chain (advanced).** `pyitensor.idmrg.apply_mpo(result, W_bulk, cutoff=..., maxdim=...)` is the infinite-chain analogue of the finite backends' `applyMPO`: it contracts a periodic MPO onto every site of the converged unit cell and re-canonicalizes/truncates the grown bond dimension back down via the standard two-sided fixed-point infinite-MPS canonicalization procedure, returning a new `pyitensor.idmrg.PeriodicMPS` that `onsite_expectation`/`two_point_correlator` accept exactly like an `IDMRGResult`. There is no `Infinite_Many_Body_Chain`-level wrapper yet, so it is reached by working with `pyitensor.idmrg` directly against `ic._result` (`itensor_version="python"` only — unlike `vev`/`correlator`, which work on both backends):

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

## 19. Running the pure-Python backend on a GPU

`itensor_version="python"` can put its tensors on a GPU instead of in host
memory. It is one process-wide switch, set before building a chain:

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")          # "numpy" (default) puts them back
print(backend.device_info())        # "jax: cuda:0" if a device was found

from dmrgpy import spinchain
sc = spinchain.Spin_Chain(["S=1/2"]*30, itensor_version="python")
...                                  # everything else is unchanged
```

It needs `jax` with its CUDA plugin installed; nothing else changes about
your script, and results agree with the host run to ~1e-11 or better.

**Whether it is worth using depends on one number: your bond dimension.**

| bond dimension | what to expect |
|---|---|
| below ~120 | the GPU is *slower* (0.1-0.5x). Use `"numpy"`. |
| ~120-160 | break-even |
| 240 | ~5x (both ground states and KPM correlators) |
| 480 | ~20x on a ground state that needs it |

The reason is not the calculation type but the arithmetic-to-overhead
ratio: each device operation costs ~0.35 ms no matter how small, so many
small tensors lose and few large ones win. Note this also means longer
chains do not help by themselves (more operations, same size each) --
bond dimension does.

Ground states, static and KPM dynamical correlators, and **real-time
evolution** (TDVP, `timedependent.evolve_and_measure` and the `submode=
"TD"` correlators) all run on the device. Time evolution is the natural
fit: entanglement grows with time, so bond dimension climbs into the
paying range by construction, while a 1D ground state's does not (see the
warning at the end of this section).

So do excited states, conserved-sector ground states, entanglement
entropies, TEBD, TDVP-GSE, the four-point correlation tensor, and the
`"CVM"`/`"SECTOR"` correlators -- each checked against the NumPy result.
Three more **run and give the same answer but do not stay on the device**:
infinite chains (`gs_method="idmrg"` and `"vumps"`) and non-Hermitian DMRG
(`nhdmrg`) move data to the host on every iteration, so do not expect them
to get faster on a GPU (see `docs/gpu_cpu_performance.md`'s
device-compatibility section).

Three things to set when you use it:

* `backend.set_pad_bonds(K)` (with `K` your bond dimension) if the script
  does *one* calculation and exits: it makes every tensor shape identical
  so the GPU compiles each kernel once, worth 1.4-3.1x on such a run. Skip
  it for long sweeps inside one process, where it costs 6-31%. It pads the
  state's bonds, never the Hamiltonian's, whose bond dimension does not
  move. The padding is exact for the state, and it leaves every method's
  algorithm alone except one-site TDVP, where a padded zero direction is
  a live basis vector the integrator can populate: `tevol_method=
  "TDVP_GSE"` under padding used to run one-site TDVP on the manifold of
  bond dimension `K` instead of its Krylov expansion, 0.49 off the
  unpadded run on a 10-site quench at `K=maxm`. Since the 2026-09-24
  audit that route is exempt from the padding the way the Hamiltonian
  is, the padding being stripped once at trajectory entry from whatever
  state the evolution is handed (keyed on the state rather than on the
  flag since the 2026-09-24 second-pass audit, finding 18, so a state
  padded earlier and evolved after the flag is cleared is stripped too),
  so a padded `TDVP_GSE` run follows the unpadded one to roundoff
  (`docs/audit_2026_09_24_hole_hunt.md`, finding 16, and §21). The one
  case still started padded is `submode="TDZ"` with
  `tevol_method="TDVP_GSE"` at `tdvp_gse_sweeps=0`, since TDZ carries
  its wavefunction between calls; at `tdvp_gse_sweeps>0` its first
  expansion strips the padding. Under JAX with `set_jit("auto")` the
  exempt route retraces once per bond dimension it grows through, which
  is the cost padding is meant to remove, though that route never kept
  frozen shapes anyway. On a *consumer* card skip it above small bond dimension entirely:
  padding trades arithmetic for shape stability, and where FP64 is
  1/32-rate the arithmetic is the expensive half (measured on a GTX 1060:
  a win at maxm=30, 2.6x slower at maxm=120).
* `backend.set_jit(True)` fuses the engine's hot inner kernels into one
  compiled kernel each, which is what lowers the ~0.35 ms-per-operation
  floor. The default `"auto"` turns it on exactly when `set_pad_bonds` is
  set, because compiling only pays off once the tensor shapes stop
  changing -- and measured on an H200, the two together are worth
  **6.4-12.1x on a cold run** (a script that starts, computes and exits),
  where either alone is worth only 2.4-3.8x. Once everything is compiled
  the trade reverses: padding then costs arithmetic on known-zero blocks,
  so for a long sweep inside one process set `set_jit(True)` *without*
  padding, which measured 1.3-1.6x warm.
* for a KPM correlator, set `kpmmaxm` equal to `maxm`, so the ground-state
  solve and the moment recursion share one set of shapes.

`docs/gpu_cpu_performance.md` has the full measured comparison, including
which models are and are not meaningful GPU benchmarks -- a uniform
Heisenberg chain's ground state converges at chi ~ 60 and cannot show a
speedup at any `maxm`.

## 20. Performance: BLAS threads

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

## 21. What raises, and what changed in the 2026-08 audit

A cross-backend audit in August 2026 (`docs/audit_2026_08_hole_hunt.md`,
which records every reproduction) went looking for calls that silently did
something other than what was asked. Most of what it found is now either
fixed or loud. Two of the fixes **change numbers** and are worth knowing
about if you have results from before them:

- **`submode="TD"` and `submode="TDZ"` were a factor of $\pi$ too large.**
  The Fourier transform behind them omitted the $1/\pi$ of the convention
  $S_{AB}(\omega)=\frac1\pi\mathrm{Re}\int_0^T\!dt\,e^{i\omega t}C(t)w_\delta(t)$
  that this guide documents and that every other submode already followed,
  and gave the $t=0$ sample full rather than half weight. On a 4-site chain
  the exact sum rule $\int S(\omega)d\omega=\langle AB\rangle=0.25$ came
  out as 1.28. Peak positions and widths were never affected -- the error
  is uniform in $\omega$ -- so only absolute spectral weight changes.
  `sxt_to_skomega` and the infinite-chain $S(k,\omega)$ reduction share the
  same transform and inherit the fix.
- **`itensor_version=3` real-time evolution from a non-unit-norm state.**
  `evolution_ABA` (whose start state is $A|\mathrm{gs}\rangle$) and
  `evolve_and_measure(wf=...)` were rescaled by $1/\|\psi\|^2$ for every
  $t>0$ under `tevol_method="TDVP"`/`"TDVP_GSE"`, while $C(0)$ stayed
  correct. Against ED on a 5-site chain the error went from 0.176 to
  1.7e-9. `"TEBD"`, `"MPO"`, `itensor_version=2` and `"python"` were never
  affected.

Beyond those, the following now raise where they used to be silent:

| Call | Was | Now |
|---|---|---|
| `sc.maxm = m` with `m<1` | SIGABRT (v2/v3), wrong energy (`"python"`) | `ValueError` |
| changing `maxm`/`nsweeps` between `gs_energy()` calls | returned the first, less converged energy | re-solves |
| `sc.tevol_method = "tdvp"` (any unrecognized name) | ran the legacy MPO-Taylor integrator | `ValueError` |
| `submode=` on a non-Hermitian Hamiltonian | every submode returned the CVM resolvent | dispatches; `NotImplementedError` for the Hermitian-only ones |
| unknown keyword to `submode="KPM"` | silently discarded | `TypeError` |
| `name="ZZ"` with `submode="KPM"`/`"EX"`, or any `mode="ED"` | `RuntimeError: No active exception to reraise` | works |
| `submode="CVMimag"`/`"maxent"` | `exit()` -- terminated your process | `NotImplementedError` |
| `kpm_energy_truncate=True` on `itensor_version=2` | silently ignored | `NotImplementedError` |
| `kpm_energy_truncate=True` with far too small a `kpm_scale` | plausible spectrum at the wrong energy | `RuntimeError` in the worst regime (see `docs/known_issue_kpm_energy_truncation_window.md`) |
| an operator on an out-of-range site (`"python"`) | treated as the identity | `IndexError` |
| `get_excited_states(n=...)` beyond the Hilbert space | ground on for minutes | `ValueError` |
| `Thermal_Spin_Chain(..., T<0)` | treated as `T=0` | `ValueError` |
| `Thermal_Spin_Chain(..., itensor_version=...)` | dropped | honored |
| `get_distribution*(mode="ED")` | `AttributeError` deep inside | dispatches, or `NotImplementedError` |
| `vev(C[i]*C[i])` | `AttributeError` on every DMRG backend | `0` (as ED always answered) |
| `Many_Body_Chain.evolution()` | `AttributeError` (it called a function that does not exist) | removed |
| `gs_energy_generalized()` on a chain below 2 sites (`"python"`) | the Rayleigh quotient of an untouched random state (-0.3049 for an exact -0.5) | `RuntimeError` |

And these combinations now work where they used to fail:

- **conserved-sector mode + `get_rdm()`** on `itensor_version=3`, which
  returned uninitialized heap memory (non-Hermitian, different every run);
- **conserved-sector mode + a repeated-site operator product** such as the
  four-point correlator's `Cdag_0 C_1 Cdag_0 C_1`, which aborted the whole
  process inside ITensor and now evaluates to 0 as it should;
- **conserved-sector mode + site/pair entropy, mutual information, and the
  default `get_correlation_matrix()`**, all of which raised;
- **`get_rdm(i=ns-1)`** (the last site) on `itensor_version="python"`, and
  on `"julia_live"` since the 2026-09-25b pass, where it raised
  `BoundsError`.

One thing the audit's own `name=` fix broke, since fixed: centralizing
`name=` resolution made every submode go through the same normalizer,
whose explicit-pair branch admitted only `MultiOperator`s -- so
`name=(A,B)` with the already-built operators `toMPO()` returns started
raising `ValueError` on the paths that had always accepted them
(`mode="ED"`, and `submode="EX"`). They are accepted again, per element
and including products of them; see *What `name=` accepts* in §6 for
which submodes can and cannot take one.

**The ED Lehmann sum is now exact (2026-09-02).** `mode="ED",
submode="ED"` -- the reference every other submode is validated against --
truncated its final-state sum to `emu < max(es)`, comparing *absolute*
eigenvalues against a frequency measured from the ground state. Three
consequences, all reproduced: adding a constant to `H` changed the answer
(1.1e-3 at `H+3*Id` on a 6-site Heisenberg chain) although a constant can
move no pole; a spectrum sitting entirely above `max(es)` truncated itself
away and raised numba's `"zero-size array to reduction"` from deep inside
the sum; and the result depended on how wide an `es` window the caller
happened to ask for. `dynamical_correlator_finite_T` (the `T>0` route,
and the reference for the METTS correlator) carried the same line.

Both now sum over the full spectrum. **For almost every existing caller
this changes nothing measurable**: whenever the ground-state energy is
comfortably negative -- the ordinary case, and why this survived since
2025-11 -- the old threshold retained everything that mattered, and the
difference is float noise (1.4e-17 on the 5-orbital IETS atom of
`tests/test_atom_iets.py`, 8.3e-13 on a 6-site Heisenberg chain). It
differs only where the old code was wrong. Nor is it slower: the `T=0`
sum reads only `nex` rows of one matrix-element array and `nex` columns
of the other, so only those are now built, which is cheaper than the full
`n x n` products the truncated version was computing. The `T>0` route
*does* cost more, and unavoidably so: it sums over every initial state
with its own Boltzmann weight, so every row and column is genuinely read
and no such shortcut exists. That is O(dim^2) in the sum, the same order
as the dense diagonalization it already performs -- 0.47 s on a
1024-dimensional Hilbert space, against 0.004 s for a 6-site chain.

### The 2026-09 audit

A second cross-backend audit (`docs/audit_2026_09_hole_hunt.md`) recorded
36 confirmed findings, all of them now fixed. The last to land was the
iDMRG matrix-free item, which stood as PARTIAL for a while because its
compute half arrived without its memory half; the second pass closed on
2026-09-22, and what it buys an infinite-chain run is in §18. It also
left two *open items*, found while
re-measuring the record rather than during the hunt itself, and both of
those were closed on 2026-09-22: the real-time submodes sat off the
correlator convention (O1) and `delta` meant a different width under
`submode="KPM"` on each solver (O2). The list below is the
physics-facing half — what changes for a user who never reads the
architecture docs.

**Results that are not comparable across this change.** Read this the way
you read the two 2026-08 items above.

- **`submode="CVM_explicit"` returned exactly twice the correct
  spectrum**, at every frequency. The submode is DMRG-only (`mode="ED"`
  raises `NotImplementedError`), and the factor of 2 was the same on
  every DMRG backend, the fix being one backend-agnostic line.
  Verified directly on the 6-site staggered Heisenberg chain:
  the first four points came back as
  `[0.2169497, 0.0920095, 0.1731231, 0.1051759]` against `submode="CVM"`'s
  `[0.1084749, 0.0460046, 0.086561, 0.0525879]`, i.e. a ratio of 2 to six
  digits; they now agree to 5.7e-07. Any peak height, integrated weight
  or figure produced with `"CVM_explicit"` before this is a factor of 2
  too large. This is the one correlator change that hits Hermitian pairs
  as well.
- **Three dynamical-correlator routes moved onto the house convention**
  of §6 (the complex Lehmann density) from $-\frac1\pi\mathrm{Im}\,G^R$:
  `mode="DMRG" submode="CVM"`, `submode="ROOTN"` under both `mode=`
  values, and `mode="ED" submode="ED"` (at $T=0$ and $T>0$), which is
  the exact reference every other submode is validated against. **Where
  every $M_n$ is real nothing moves at all** — the two conventions are
  then identical (measured to 1.1e-16), and that covers every
  `name=(A,B)` example in this guide, Hermitian pairs and the
  $S(q,\omega)$ sweep's $(S^z_i,S^z_j)$ alike; §6 states the criterion
  and why $A=B^\dagger$ is not it. For a pair with complex $M_n$ the curves move by more than
  their own peak: on a 4-site complex-hopping chain with
  $A=c^\dagger_0$, $B=c_2$ at $\delta=0.15$ the two conventions differ by
  up to 5.720e-01 against a peak of 0.61. `submode="ED"` also returns a complex array
  now where it returned `float64` before (same real part, to the digit;
  the imaginary part it used to drop peaks at 0.5682 on that pair), so a
  caller that assumed a float array should take `.real`. `submode="KPM"`
  (the default), `"INV"`, `"CVM"` under `mode="ED"`, `"EX"`, `"TD"`,
  `"TDZ"` and `"SECTOR"` were untouched *by that change*; `"TD"`,
  `"TDZ"` and `"KPM"` moved later, on 2026-09-22, for the two reasons in
  the following two entries.
  `mode="ED" submode="ROOTN"` and `mode="DMRG" submode="ROOTN"` are also
  roughly **2x slower** now, for the reason given in §6.
- **`submode="TD"` and `submode="TDZ"` now return the density too**,
  where they used to return the complex one-sided transform it is
  assembled from (the audit's open item O1, closed 2026-09-22). §6's
  `"TD"` entry carries the identity and the cost rule; what matters for
  a result you already have is that the array moved, and that the
  imaginary part moved most. On the audit's 4-site complex-hopping chain
  ($A=c^\dagger_0$, $B=c_2$, $\delta=0.4$, `dt=0.1`, exact peak 0.2313)
  the distance to the exact Lehmann density went from 1.16e-01 to
  2.57e-04 under `mode="DMRG"`, from 2.76e-01 to 2.57e-04 under
  `mode="ED"`, and from 1.13e-01 to 3.31e-02 for `"TDZ"`, a residual
  this entry used to put down to the contour and which was in fact the
  frequency stage's FFT grid; since the 2026-09-24 audit (finding 7, in
  the subsection below) `"TDZ"` lands 5.4e-04 from exact and `"TD"`
  3.2e-05. On a Hermitian pair the
  imaginary part used to run at 60% of the peak (1.11e-01 against
  0.1869 on a 4-site Heisenberg chain at $\delta=0.3$) and is now
  exactly zero, with the curve itself landing 2.82e-04 from exact. Taking
  `.real` of the result, which this guide used to instruct, is now a
  no-op on a real-weight pair and *wrong* on a complex-weight one.
  `mode="ED"` moved further than `mode="DMRG"` because that route
  additionally read the operator pair in the opposite order, returning
  $S_{BA}$ where every other route returns $S_{AB}$; it was invisible
  whenever $A$ and $B$ are the same operator, which is every example
  here. A pair that is not provably its own adjoint now costs two
  evolutions rather than one.
- **`delta` under `submode="KPM"` is now the same broadening it is in
  the resolvent submodes**, FWHM $=2\delta$ at the band centre, on both
  solvers (the audit's open item O2, closed 2026-09-22), and
  `kpm_n_scale`'s default changes from 3 to 1 because the base moment
  count it multiplies is now the calibrated one. Neither solver realized
  the intended width before, and they missed it differently: `mode="ED"`
  rescaled onto a window three times the bandwidth and counted moments
  by its own rule, giving about $0.93\delta$, while `mode="DMRG"` gave
  about $1.6\delta$. Every KPM curve therefore changes width, and the
  two solvers stop disagreeing: measured with $A=B=S^z_0$ on three
  chains, $\max|\mathrm{ED}-\mathrm{DMRG}|$ against the resolvent peak
  went from 95 to 124 per cent down to 4 to 7 per cent, and to 0.07 to
  0.2 per cent on `itensor_version="python"` once the 2026-09-24 audit
  cut the DMRG routes to the calibrated moment count (finding 4, in the
  subsection below). Peak heights
  change with the width and the sum rule does not, which is exactly why
  nothing caught this: both routes integrated to 0.250000 against an
  exact 0.250000 before and after. A script that was tuning `delta` by
  eye against a KPM curve, or carrying `kpm_n_scale=3` explicitly, gets
  a different spectrum now. `get_distribution()`'s own KPM path is
  deliberately left off the calibration (§14).
- **`vev(op, npow=n)` on an ED route** returned $\langle O\rangle$ for
  every $n$; it now returns $\langle O^n\rangle$. Measured on a 4-site
  Heisenberg chain, `vev(Sz[0], npow=2, mode="ED")` went from 5.6e-17
  (i.e. $\langle S^z_0\rangle$) to 0.25.
- **`gs_energy_fluctuation` on an ED route** inherited that and also
  ignored `mode=`. On the same chain `sc.mode="ED"` reported 2.0561 and
  now reports 0.0; a 2-site chain on the *default* backend (which
  `mode.py` routes to ED by itself) reported 1.1456 and now reports
  1.05e-08. Anything using it as a convergence criterion behaves
  differently on those routes. DMRG-route values are bit-unchanged.
- **`exponential(h, wf)` on the DMRG backends** — see §2. Measured as
  $\langle\mathrm{GS}|e^{zH}|\mathrm{GS}\rangle$ on a 4-site Heisenberg
  chain under `itensor_version="python"`, against the ED value: 0.6776
  vs 0.6676 at $z=0.25$, 0.5184 vs 0.4457 at $z=0.5$, 0.6897 vs 0.1987
  at $z=1$ (1.5%, 16%, 247%). It now agrees with ED to ~1e-8 at all
  three.
- **Boson occupation projectors under `mode="ED"`** were wrong on any
  Hamiltonian that does not conserve the total boson number. On a
  3-site, 4-level chain with an $a_i+a^\dagger_i$ drive, site 0's
  $P(n)$ read `[0.0042, 0.0054, -0.048, 0.0706]`, summing to 0.032 —
  a negative probability included. It now reads
  `[0.0322, 0.2393, 0.4732, 0.2553]`, summing to 1.0000 and matching
  DMRG. A number-*conserving* Hamiltonian is unaffected to machine
  precision, which is why nothing caught this; $\langle n_i\rangle$,
  the energies and every DMRG backend were always right.
- **`SpinBoson_Chain(..., maxnb=[...])`** silently built the 4-level
  site regardless: `SpinBoson_Chain(["B","S=1/2"], maxnb=[6,None])` gave
  sites `[104, 2]` and now gives `[106, 2]`. Earlier results at a
  non-default `maxnb` are results at `maxnb=4`. With no `maxnb=` nothing
  changes.
- **`itensor_version="python"`**, several ways. Its `applyMPO` used to
  truncate in a gauge where discarding the smallest singular values is
  not the optimal truncation, so every repeated-application consumer
  moved at the 1e-5..1e-9 level: `applyoperator`, `vev(..., npow>1)`,
  `gs_energy_fluctuation` (5.06e-05 → 6.27e-06 on a 10-site Heisenberg
  chain at `maxm=30` — the old number was mostly the MPO-application
  error rather than the state's), the KPM Chebyshev recursion and
  therefore every KPM correlator and `get_distribution`, CVM,
  `applyinverse`, and the MPO-Taylor evolution paths. The new values are
  the ones that match `itensor_version=3` and ED. The product is now
  formed exactly and truncated once afterwards, which is more work per
  application; the cost is workload-dependent, so budget for it rather
  than expecting a fixed factor. Separately, its TDVP real-time
  evolution did not canonicalize the state before the first half-sweep,
  so a trajectory whose start state was not already in that gauge took
  its first step under an unprojected generator; the audit record
  (`docs/audit_2026_09_hole_hunt.md`, finding #1) carries the
  reproduction and the mechanism. Which start states escaped it is not
  recorded there and is not claimed here. And a chain object reused across
  several `set_hamiltonian()` calls could keep a start state belonging to
  the previous Hamiltonian. Both are fixed; no
  `itensor_version` other than `"python"` was affected by any of this.
- **`promote_to_dense()` now keeps the sector's energy on
  `itensor_version="python"`** — the guarantee §3 states ("a bare
  `gs_energy()` afterwards still returns the sector's energy rather than
  re-solving") was true on `itensor_version=3` and `mode="ED"` but not
  there: promotion invalidated the session's cached energy, and
  pyitensor's unconstrained re-solve then left the sector. On a 3-site
  Hubbard chain confined to `Nf=3` the audit recorded a post-promotion
  energy of $-2.3399130755$ with $\langle N\rangle=2$ against the
  sector's own $-1.9408140222$ with $\langle N\rangle=3$; it now returns
  the sector value, re-measured here. Any `"python"` result that called
  `gs_energy()`/`get_gs()` after `promote_to_dense()` — the documented
  workflow — is a different number than it was. The other backends are
  unchanged, and a chain where the sector and global ground states
  coincide never saw it.

**And these now raise, or now work, where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| `Bosonic_Chain(n, maxnb=[...])` with `maxnb != 4` on `itensor_version=2` | SIGABRT inside ITensor | `ValueError` naming the working backends |
| `Bosonic_Chain`/`SpinBoson_Chain`/`Parafermionic_Chain`/`Spin_Fermion_Hamiltonian(..., itensor_version=...)` | `TypeError: unexpected keyword argument` | honored |
| `SpinBoson_Chain(["boson", ...])`, or any unknown site label | `RuntimeError: No active exception to reraise` | `ValueError` listing the accepted labels |
| `SpinBoson_Chain(..., maxnb=[...])` of the wrong length, or `n=` contradicting `len(sitesin)` | silently dropped | `ValueError` |
| `SpinBoson_Chain` under `mode="ED"` | failure inside an unimplemented stub | `NotImplementedError` naming the DMRG alternatives |
| `vev(..., npow=n)` with `n<0`, or with `n!=1` at `T>0`, **on an ED route** | accepted silently, answering with $\langle O\rangle$ | `ValueError` / `NotImplementedError` (the DMRG routes still answer) |
| `gs_energy_fluctuation(npow=...)` | accepted and misused | `TypeError` |
| `exponential(h, wf)` with `h` neither Hermitian nor anti-Hermitian, on DMRG | a warning and a number | `NotImplementedError` |
| `applyoperator`/`summps`/`applyinverse`/`scale_mps`/`exponential` with a `mode=` naming the other backend | ran the wavefunction's own backend anyway (or a bare `raise`) | `TypeError` saying which backend the state belongs to |
| `get_rdm()` on any ED route | `AttributeError` several frames deep | `NotImplementedError` (on `"julia_live"` since the 2026-09-25b pass) |
| the bond entanglement entropy on any ED route | `AttributeError` several frames deep | `NotImplementedError` pointing at `get_site_entropy`/`get_pair_entropy` (on `"julia_live"` since the 2026-09-25b pass) |
| `get_bond_entropy(wf, i, j)` with a site out of range | the audit recorded a SIGABRT on the C++ backends | `IndexError` |
| `get_hamiltonian()` before `set_hamiltonian()` | `'int' object is not iterable`, or `takes 0 positional arguments` | `ValueError` naming the fix |
| a misspelled `ctmode`/`dmmode`/`fpmode` | `RuntimeError: No active exception to reraise` | `ValueError` listing the accepted values |
| a misspelled `basis=` in `get_correlation_matrix` | the electron-basis matrix, silently | `ValueError` |
| `submode="CVM_explicit"` with $A^\dagger\neq B$ | a `print` then `RuntimeError: No active exception to reraise` | `NotImplementedError` naming the submode; since the 2026-09-25b pass it returns $C[A,B]$ |
| `get_distribution()` with no `X=` | `RuntimeError: No active exception to reraise` | `TypeError` naming `X=` |
| `MultiOperator / "a"`, `*` by an ndarray, or `multioperator.obj2MO` on an object it does not recognize | `RuntimeError: No active exception to reraise` | `TypeError` naming what was passed |
| `sc.mode = "dmrg"` (any unrecognized spelling) | accepted, and silently decided by the automatic fallbacks | `ValueError` |
| a conserved sector under `mode="ED"` + the `get_correlation_*` family | `AttributeError` | works, via the backend-agnostic `dmmode="explicit"` |
| `get_dynamical_correlator()` on `Parafermionic_Chain` without an explicit `mode=` | the audit recorded a process abort | works, and obeys `self.mode` and the automatic ED fallbacks like every other chain class |

`Parafermionic_Chain` lost its own `get_dynamical_correlator` override in
the process — it was the base-class method minus the mode resolution and
minus the `name="..."` string resolution — so both the documented string
form of `name=` and `mode=`/`self.mode` work there now.
`Parafermionic_Chain.test(ntries=...)` is honored too, as
`Spin_Chain.test(ntries=...)` already was.

### The 2026-09-24 audit

A third hunt (`docs/audit_2026_09_24_hole_hunt.md`, five lenses over the
commits `1c2606d..765b537` that closed the 2026-09 audit's open items)
recorded 16 confirmed findings. Every code defect among them is fixed,
with one residual kept by design, `submode="TDZ"` under `set_pad_bonds`
at `tdvp_gse_sweeps=0` still being started from the padded state (§19),
and finding 6, which was documentation only, is closed by §6's
`"KPM"` entry.
Several of the findings corrected what the entries above say, the
`"TDZ"` residual and the KPM agreement between the two solvers among
them, and those entries now point here. Each item below names its
finding, and the record carries the reproduction that was actually run.

**Results that are not comparable across this change.**

- **`submode="TDZ"`, and `submode="TD"` at `predict=False`, were
  evaluated on an FFT grid** of spacing $2\pi/(n_t\,dt)$, about
  $1.05\delta$ at the default `damping_periods=6`, and interpolated
  linearly onto `es`, meaning that a Lorentzian of half-width $\delta$
  got about one sample per $\delta$. The damped sum is now evaluated at
  each requested frequency. On the 2026-09 audit's complex-hopping chain
  at $\delta=0.4$ (exact peak 0.2313) the distance to the exact density
  goes from 3.314e-02 to 5.408e-04 for `"TDZ"` and for `"TD"` at
  `predict=False`, and from 2.571e-04 to 3.169e-05 for `"TD"` at its
  default, and the returned curves move by up to 14 per cent of the peak
  for the first two (19 per cent at $\delta=0.25$) and by 2.89e-04 for
  the third. `sxt_to_skomega`, the infinite-chain $S(k,\omega)$
  reduction, stayed on the FFT stage and returned the same bits as
  before (§6, finding 7), until the second pass below.
- **Every DMRG KPM spectrum moves by roughly $2/N$ onto the ED one**,
  $N$ being the calibrated moment count, since the DMRG routes
  reconstructed from $N+2$ moments. On two decoupled Heisenberg dimers
  the band-centre peak goes from 0.666761 to 0.620591 at $\delta=0.2$ and
  from 1.266336 to 1.220236 at $\delta=0.1$, on `itensor_version=3` and
  `"python"` alike, both now equal to the ED peak, and the three
  residuals of the 2026-09 record's O2 table go from about 1.3e-02 to
  2.55e-04, 2.09e-04 and 3.50e-04 on `"python"`. `mode="ED"` does not
  move (§6, finding 4).
- **`get_distribution(mode="ED")` moves by exactly the factor `scale`**:
  on a 4-site open Heisenberg chain with $X=S^z_0$ at the default
  `scale=10` and $\delta=0.05$ the total weight goes from 0.1000 to
  1.0000 and the peak from 0.517316 to 5.173159, and with `xs=` the
  imaginary part, which was a copy of the real one, goes from
  $\max|\mathrm{Im}\,y|=4.96$ to exactly 0. Nothing moves on any DMRG
  backend (§14, finding 3).
- **`ISy` in the first operator of a correlator** returned exactly minus
  the correlator under `"KPM"`, `"EX"`, `"TD"` and `"TDZ"`, since
  `get_dagger()` treated this anti-Hermitian name as Hermitian. On a
  4-site Heisenberg chain in a field, the KPM curve of
  `(get_operator("ISy",0), Sz[1])` goes from exactly minus the curve of
  `(1j*Sy[0], Sz[1])` to identical to it, and its first moment from
  +0.079470 to -0.079470 (exact -0.079685); `sc.is_hermitian(ISy)` is
  now `False` and `sc.is_hermitian(1j*ISy)` `True`, where both were the
  other way round. Only operators that contain `ISy` move (§6,
  finding 2).
- **An operator that mixes pre-Jordan-Wigner `C`-type names with bare
  `A`-type ones** was canonicalized with the wrong sign, so `simplify()`
  changed its value and the Hermiticity, zero and adjoint-pair proofs
  could be false. Such terms are now left as written (§2). A `"TD"`
  pair of that kind on a 4-site complex-hopping chain goes from
  3.616e-02 off exact, with one evolution, to 3.063e-05 with two, on
  `"python"` and v3 alike, and no dmrgpy chain class builds such a term,
  so no ordinary Hamiltonian moves (finding 1).
- **`mode="ED"` `submode="TD"` on a degenerate ground state** built the
  two halves of the density on two different randomly chosen ground
  states: three identical calls on a degenerate 3-site Heisenberg chain
  spread by 2.03e-01 against a 1.768e-01 peak. Both halves now use the
  cached ground state, so the calls are identical and agree with
  `"KPM"`, `"INV"`, `"CVM"` and `"ROOTN"`, though still not with
  `submode="ED"`'s average over the manifold, a split between the two ED
  conventions that predates this fix. `evolution_ABA(mode="ED")` and
  `evolve_and_measure(mode="ED")` without `wf=` change the same way, and
  a non-degenerate ground state moves by at most 3e-11 (finding 8).
- **The Kondo potential-interference term on a non-uniform `es`**
  weighted every point with the grid's first spacing. On a single
  $S=1/2$ at 10 T the peak on a sinh-mapped grid goes from 68.2332 to
  0.6704 against an exact 0.6709, and a uniform grid moves by 2.4e-08
  (§17, finding 13).
- **`excitation_energies(k, n)` with `n>=2` on `itensor_version="python"`**,
  and `spectral_weights`/`dynamical_structure_factor` with it, could
  drop one copy of a degenerate level of $H_{\rm eff}(k)$ above the dense
  threshold (dimension 256) and return the next distinct level in its
  place. At $k=0.37$ on the `n_uc=2` critical Heisenberg cell at
  `maxm=10`, `[0.289896072 0.291597913 0.313428946]` becomes
  `[0.289896072 0.289896072 0.291597913]`, and the $S^x$ and $S^z$ weight
  fractions go from 0.058 to 0.184014 and from 0.339 to 0.041282, both
  now equal to the dense path's. `n=1`, `excitation_gap` and dimensions
  up to 256 are byte-identical (§18, finding 15).
- **`tevol_method="TDVP_GSE"` under `set_pad_bonds(K)`** on
  `itensor_version="python"` ran a different algorithm, one-site TDVP on
  the manifold of bond dimension $K$, its bond growth coming from a QR
  completion of the padded zero directions rather than from the Krylov
  expansion. Padded runs through `quench_tdvp_gse` and
  `evolve_and_measure_tdvp_gse` now follow the unpadded ones: on a
  10-site XXZ quench at $K=$`maxm`$=4$ the two went from 0.4928 apart to
  8.9e-16 at `tdvp_gse_sweeps=0` and from 4.1e-06 to 1.7e-07 at
  `tdvp_gse_sweeps=3`. Note the direction at `tdvp_gse_sweeps=0`, where
  the padded run used to sit 9.8e-05 from ED only because the larger
  manifold held the answer, and now sits 0.4929 from it, like the
  unpadded frozen product state, which is the correct one-site answer
  (§19, finding 16).

**And these now raise, or now work, where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| `kpm_n_scale` not a positive integer, on `"python"`, `mode="ED"` or `"julia_live"` | rounded down silently (1.5 gave 1x, anything below 1 the 16-moment floor) | `TypeError`/`ValueError` naming it, as v2/v3 already did; `kpm_n_scale=2.0` raises too |
| a misspelled keyword to `submode="TDZ"` (`alpha=`, `nmax=`) | ignored, the default-parameter spectrum bit for bit | `TypeError` |
| a `toMPO()` operator to `submode="TDZ"` | a crash several frames deep | `TypeError` naming `toMPO` |
| `get_dynamical_correlator_MB(name="ZZ", i=1, j=1)` under `"TD"`/`"TDZ"` | $C[S^z_0,S^z_0]$, the sites dropped | the sites asked for (the public route was never affected; `"ROOTN"` on this route followed on 2026-09-25) |
| `get_kondo_spectrum(mode="ED")` with a keyword it does not read | ignored, so `Jrho=` for `Jrho_s=` removed the Kondo peak | `TypeError` naming the unknown keywords |
| `get_kondo_spectrum(mode="DMRG")` at a degenerate ground state | one arbitrary member's spectrum | the same by default; `n_gs=g` gives the equal-weight average over the manifold that `mode="ED"` takes |
| the potential term with an `es` that misses spectral weight | a silently low spectrum | `RuntimeWarning` from the sum rule |

### The 2026-09-24 second pass

A fourth hunt (`docs/audit_2026_09_24b_hole_hunt.md`, five lenses over
the single commit `30200a4` that closed the third) recorded 18 confirmed
findings, twelve of them older than that commit and reached by probing
next to it, and every one is fixed. Each item below names its finding;
the record carries the reproduction that was actually run.

**Results that are not comparable across this change.**

- **Every infinite-chain $S(k,\omega)$**, on v3 and `"python"`, was
  mirrored in frequency and interpolated off an FFT grid. It is now
  conjugated after the spatial sum and evaluated directly at each
  frequency: on the transverse-field paramagnet the $k=\pi$ magnon moves
  from $\omega=-1.045$, with 0.928 of its weight below zero, to $+0.905$
  against $\varepsilon=0.900$, and on a single magnon at $\delta T=6$ the
  curve goes from 21 per cent of the peak off the closed form to 2e-14
  (§18, findings 7 and 8).
- **`evolve_and_measure` and `evolution_ABA` on `mode="ED"` ran backwards
  in time**, and on every DMRG backend returned the conjugate of
  $\langle O\rangle$. Larmor precession of a +x chain in a field now gives
  $\langle S^y_0\rangle(t=1)=+0.4207$ on ED, where it gave $-0.4207$, and
  $O=S^z_0+iS^x_0$ on an eigenstate now gives $+0.5i$ on DMRG, where it
  gave $-0.5i$. A real Hamiltonian, a real start and a Hermitian
  observable with real matrix elements, which is what every earlier test
  measured, do not move (§7, findings 9 and 10).
- **A state set by hand now reaches the solver.** After `set_gs()` every
  DMRG correlator measured the session's own solved state and then
  overwrote the chain's with it (TD on a pure doublet member goes from
  2.04e-01 to 3.3e-06 off exact), and `set_initial_wf`/
  `set_initial_wf_guess` never reached the session, so the transverse-field
  Ising example plotted $M_z/n\approx0.001$ in the ordered phase where the
  seeded branch gives 0.4998 (findings 11 and 12).
- **`get_kondo_spectrum(mode="DMRG", n_gs>1)`** re-swept each member before
  measuring it, so on v2 and v3 an excited member split below `delta`
  relaxed into the lower one (1 of 6 runs at the average 1.5, now 6 of 6),
  left `gs_energy()` at the last member's energy afterwards, and was inert
  under `submode="EX"`, which now measures from the chain's own state
  (1.1392 to 1.8711 across seeds, now 1.512780 on every one; §17,
  findings 13 to 15).
- **Smaller moves**: the Kondo terms on a non-increasing `es` (3.035 to
  0.03183 off exact on a refined block appended after a coarse one,
  finding 17); a padded ground state evolved with `set_pad_bonds` cleared
  under `TDVP_GSE` at `tdvp_gse_sweeps=0` (0.4929 to 1.2e-15 off the
  unpadded run, finding 18); `mpsalgebra.disentangle_manifold` on a
  Hermitian operator the canonical proof cannot see (a non-orthonormal
  basis, Gram error 0.32, now 1e-15, and proven operators move at the
  gauge level, finding 1).

**And these now raise where they used to be silent:**

| Call | Was | Now |
|---|---|---|
| KPM with a pole outside the rescaled window (`kpm_scale` below 1/2) | DMRG spectra up to 109 times the true peak, ED integrals up to 1e6 | `RuntimeError` on every backend, the guard being 1.5 times the exact moment bound (findings 2 and 3) |
| `kpm_energy_truncate=True` on `mode="ED"` | the untruncated curve | `NotImplementedError`, also through `mode.py`'s fallback (finding 4) |
| `kpm_n_scale=1.5` on a non-Hermitian `mode="ED"` KPM | `TypeError`, where `mode="DMRG"` accepted it | accepted on both, neither reading it (finding 5) |
| a misspelled key, a method or `itensor_version` in `kpm_finite`'s `window_chain_kwargs` | ignored, the default spectrum bit for bit | `TypeError` (finding 6) |
| `name=(A,B)` with `i=`/`j=` in `get_dynamical_correlator` and every wrapper of it | the sites silently dropped | `TypeError` (finding 16) |
| `n_gs>1` under `submode="SECTOR"` | the single-state value | `NotImplementedError` (finding 15) |

### The 2026-09-24 third pass

A fifth hunt (`docs/audit_2026_09_24c_hole_hunt.md`, four lenses over the
single commit `867e2b4` that closed the second pass) recorded 18 findings,
six of them from that commit and twelve older, and every one is fixed; two
needed a rebuild of the compiled extensions. Each item below names its
finding.

**Results that are not comparable across this change.**

- **`gs_energy_fluctuation()` below full bond dimension, on every DMRG
  backend**, was set by truncating $H|\psi\rangle$ to `maxm`. It is now
  $\lVert(H-\langle H\rangle)|\psi\rangle\rVert$ with an uncapped
  application: on an 8-site Heisenberg chain solved and measured at
  `maxm=3` it goes from 6.7e-03 to 0.344 on `"python"` and from 1.5e-02 to
  0.347 on v3, and an exact state measured at `maxm=3` from 1.64 (v3) to
  1e-11. `mode="ED"` subtracts $\langle H\rangle$ first too, so its value at
  an exact eigenstate goes from the 1e-07 roundoff floor to 1e-15. With it,
  `gs_energy(maxde=...)` stops later, since it used to stop on an
  under-reported number, and returns the refined energy it leaves on the
  chain, -3.374933 where it returned the unrefined -3.279373; a
  correlator afterwards measures that refined state instead of re-solving
  at the original `maxm` (§3, findings 6 and 7).
- **After `set_gs()` of a state off the ground manifold** every submode on
  both modes measures it from its own energy. DMRG KPM and ED KPM, CVM,
  INV, ROOTN and TD move by $E_x-E_0$ (0.3 on the 3-site chain in a field,
  where the elastic line sat at +0.3), and ED `submode="ED"`, which read no
  state at all, now gives the set state's density, 2.834 away from what it
  returned on a 2.835 peak (§6, findings 1 and 2).
- **`promote_to_dense()` keeps the state the chain holds.** A state set in
  a sector and promoted with no read in between reads -0.957107, its own
  energy, where it read the sector ground energy -1.616025, and a carried
  state that is not the sector ground state is kept (-1.131526) where
  `"python"` re-swept it out of the sector to the global ground state
  -1.857107 (§3, findings 3 and 4).
- **`set_hamiltonian(H2, restart=False)`** now drops the stored ground state,
  so `gs_energy()`, `vev()`, `get_excited()` and the direct KPM moments answer
  for H2: -1.780099 where they gave H1's -1.616025 on a 4-site chain (finding
  5).
- **Weakly non-Hermitian operators.** The Hermiticity probe is scale-free,
  so an anti-Hermitian part below about 1e-2 in absolute size is no longer
  called Hermitian: a weak loss term keeps its decay rate (Im $E_0$ =
  -0.002134, where v2 and v3 returned 0 and `"python"` a real part 3 to 18
  per cent off), a non-Hermitian chain written in small units is solved as
  one, and `disentangle_manifold` diagonalizes such an operator rather than
  its Hermitian part (§2, finding 12).
- **`"python"` operators whose coefficients are all below about 2e-7**
  were built as a different operator, a lone `1e-7*Sz0` as the zero MPO;
  `vev(1e-7*Sz0)` goes from 0 to -1.8936e-08 (finding 13).
- **v3 `TDVP_GSE` from a start whose site 0 is one local basis vector** (any
  ladder operator or projector on site 0, or a product state) left site 0
  frozen: `evolution_ABA(A=C_0)` goes from 4.7e-01 to 1e-07 off exact, the
  public (Cdag_0, C_0) TD spectrum's peak from 0.081 to 0.47, and a Néel
  start is exact in 8 of 8 runs where up to 5 of 8 failed (§7, finding 14).
- **Every `mode="ED"` real-time route** (`evolve_and_measure`,
  `evolution_ABA`, `submode="TD"`) uses an exact propagator, and moves by
  1e-08 to 1e-06 on the small chains the tests use and up to 1e-03 at
  `dt=0.2` on an 8-site chain with an extensive energy (§7, finding 17).

**Faster, with no number changed.** The first read after `set_gs()` no
longer solves the upper band edge and an energy fluctuation it discarded
(7.1 s on a 24-site `"python"` chain for a 0.04 s result), and `"python"`'s
upper band edge caps its local Krylov dimension at 20, costing about one
ground-state solve where it cost 5 to 9 on a Heisenberg chain (findings 9
and 10). A malformed KPM call raises before any ground-state work again,
and a SECTOR call makes no solve on the caller's chain (finding 11).

**And these now raise, or now work, where they used to do something else:**

| Call | Was | Now |
|---|---|---|
| a misspelled keyword to `evolve_and_measure`/`evolution_ABA` on DMRG | ignored, the run at the defaults | `TypeError`, as on ED (finding 15) |
| `h=` to `evolve_and_measure`/`evolution_ABA` on ED | `TypeError`, "multiple values for argument 'h'" | evolves under `h`, as on DMRG (finding 16) |
| `submode="SECTOR"` after `set_gs()` of a state that is not its sector's ground state | the ground state's spectrum | `NotImplementedError` (finding 8) |
| `disentangle_manifold` on an ED manifold | `AttributeError` for every operator | the eigenbasis (finding 18) |
| an unknown keyword to `gs_energy_fluctuation()` on DMRG | forwarded and dropped | `TypeError` (finding 7) |

### The 2026-09-25 open items

A sixth pass (`docs/audit_2026_09_25_open_items.md`) took ten of the items
the five hunts had left open or recorded as unreviewed leads, reproduced
each on `8dd2198`, fixed it and handed every fix to a reviewer briefed to
refute it; one needed a rebuild of both compiled extensions. Each item
below names its entry in that record.

**Results that are not comparable across this change.**

- **A setting passed to a chain constructor** now takes effect: a 12-site
  Heisenberg chain built with `maxm=4, nsweeps=10` on v3 goes from
  -5.1420906326, the `maxm=30` answer, to -5.1323602278, the number the
  same request made by assignment returns, and a constructor `mode="ED"`
  answers every read that does not pass `mode=` itself by ED (the KPM ZZ
  peak of a 4-site Heisenberg chain in a field of 0.2, 0.205144 on
  `"python"` and 0.205151 on v3, both DMRG, becomes 0.205190). The same
  holds for `Thermal_Spin_Chain(..., mode="ED")`, which ran DMRG (§1,
  init-kwargs).
- **v2 and v3 Hamiltonians in small units** (largest coefficient below
  about $10^{-6}$): the ground state, the excited states, both band edges,
  the KPM spectra and `gs_energy_generalized` now match the unit-scale
  calculation, where a 6-site Heisenberg chain written as $sH$ returned the
  Néel $-1.25$ for $E_0/s$ against $-2.4936$ from $s\approx4\times10^{-7}$
  (v3) and $10^{-7}$ (v2) down. Every v2 or v3 calculation whose largest
  coefficient is below 1 now runs on $2^kH$, which moved nothing beyond the
  run-to-run noise at ordinary scales (§2, small-units).
- **Operators with coefficients below $10^{-8}$** keep their terms on every
  backend: `vev(eps*Sz0)/eps` at `eps=1e-9` goes from exactly 0 to
  -0.189361 on a 6-site chain in a field, and $E_0/s$ of a Hamiltonian at
  $s=10^{-8}$ from 0 to the right answer on ED and `"python"`. A term at or
  below $10^{-12}$ of its operator's largest coefficient is now dropped
  where the absolute rule kept it when it was above $10^{-8}$, and an
  anti-Hermitian part between about $10^{-10}$ and $10^{-8}$ of the largest
  coefficient now makes `gs_energy()` take the non-Hermitian route
  (§2, clean-threshold).
- **`"julia_live"` warm starts**: `set_initial_wf_guess()` of an exact
  doublet member goes from $|\langle t|\mathrm{gs}\rangle|^2=0.0008$ (0.0008
  to 0.68 over runs) to 1.0000, and `set_initial_wf(x)` followed by
  `gs_energy()` from the solved $-1.0$ to $\langle x|H|x\rangle$ (§3,
  julia-warm-start).
- **`get_gs(wf0=x, reconverge=False)` on a solved chain** returns x with
  `e0` = $\langle x|H|x\rangle$ (0.0388284989 where it returned the stored
  $-2.5231886435$ on a 6-site `"python"` chain), and every later reader
  measures x; `get_gs(wf0=x)` sweeps from x (§3, get-gs-wf0). Since the
  2026-09-25b pass the state returned and measured is $x/\lVert x\rVert$,
  and the call works on the routes that resolve to ED, where it raised
  `TypeError`.
- **Readers after a correlator that followed `gs_energy_generalized()`**
  measure the generalized state: on a Hermitian chain never solved before,
  `e0` goes from $-2.493577$ to $\lambda=-3.597994$ and $\langle
  S^z_0\rangle$ from 0 to $-0.4801$ (6 sites, $A=1+0.8S^z_0$); on a
  non-Hermitian chain with no correlator before the generalized solve,
  solved first or not, `e0` goes from $-1.596396$ to
  $\lambda=-2.17581-0.213083i$ (§3, generalized-cache). Since the
  2026-09-25b pass `e0` is the generalized state's own energy instead of
  $\lambda$, which is kept as `lam_generalized` (next subsection).
- **`tc.MBChain` after `Thermal_Spin_Chain.get_gs()` at $T>10^{-5}$**:
  on a 3-site chain at $T=1$, `gs_energy()` goes from $-2.25$ to $-0.416388$
  and the KPM sum rule from $-0.166370$ to $-0.069291$ on the DMRG
  backends, and on `mode="ED"` `vev(Sz0*Sz1)` from 0 to $-0.069398$ and
  `gs_energy()` from $-2.25$ to $-1.0$, the lowest eigenvalue (§9,
  thermal-bypass). The 2026-09-25b stepper moved these again, to
  $-0.428230$ for `gs_energy()`, $-0.071372$ for `vev(Sz0*Sz1)` (both
  exact) and $-0.071331$ (DMRG) and $-0.071188$ (ED) for the sum rule, and
  the switch is now at $T>10^{-5}W$, $W$ the spectral width.
- **`gs_energy(wf0=x, reconverge=False)` on a non-Hermitian chain** is
  $\langle x|H|x\rangle$ on `"python"`, v2 and v3, where it returned the
  NH-DMRG $-1.596396$ (§6, nh-injected).
- **`submode="ROOTN"` on the lower-level route** (`get_dynamical_correlator_MB`
  with a string name) honours `i=`/`j=`: at (1,1) on a 4-site `"python"`
  chain it goes from the (0,0) curve, $4.07\times10^{-2}$ off ED's (1,1) on
  a 0.103 peak, to $4\times10^{-16}$ (rootn-ij).

**And these now raise, where they used to do something else:**

| Call | Was | Now |
|---|---|---|
| a chain constructor keyword that is not a setting (a misspelling, the chain's state, a model-built operator list) | ignored | `TypeError` naming every offending key (init-kwargs) |
| `maxm=0` or an unknown `mode=` at construction | ignored | `ValueError`, as on assignment (init-kwargs) |
| an unknown keyword to `submode="ROOTN"` | ignored | `TypeError` (rootn-ij) |
| `get_gs(best=True, wf0=x)` on a solved chain | the stored state | `TypeError`, as on a fresh chain (get-gs-wf0) |
| non-Hermitian `submode="KPM"` after `set_gs()` | a spectrum built on the left state of an earlier solve | `RuntimeError` naming the missing left state (nh-injected) |
| `gs_energy(wf0=x)` on a non-Hermitian chain | the NH-DMRG energy, x dropped | `TypeError`; `reconverge=False` takes x (nh-injected) |

**And these now work:** on `"julia_live"`, `set_gs()` and
`gs_energy(wf0=x)` on a solved chain, which raised `AttributeError` and
`TypeError`; `gs_energy_generalized()` on a v3 MPO Hamiltonian, which raised
`AttributeError` (-3.59799449, the `MultiOperator` route's value); and the
first correlator after a plain NH solve no longer repeats the solve.

### The 2026-09-25b hole hunt

A sixth hunt (`docs/audit_2026_09_25b_hole_hunt.md`, four lenses over the
single commit `e7b1196` that fixed the open items above) recorded 29
findings and two unreviewed leads, and its fix pass found one more
(finding 30). Three of the 29 came from that commit and the rest are
older, reached by probing next to it, and seventeen are one family: a
threshold in absolute energy units, which the relative term cutoff of
`e7b1196` let a Hamiltonian or an operator written in small units reach
for the first time. Every one is fixed, and both compiled extensions were
rebuilt. Each item below names its finding; the record carries the
reproduction that was actually run, and every number is from it.

**Results that are not comparable across this change.**

- **A state whose norm is not one**, set with `set_gs()`,
  `set_initial_wf()` or `set_initial_wf_guess()`, or passed as `wf0=` with
  `reconverge=False`, is now normalized once, where it becomes the chain's,
  on every backend. Before, `e0`, the session's `vev()` and the energy
  fluctuation divided by $\langle x|x\rangle$ while KPM, CVM, TD, TDZ,
  ROOTN, `evolve_and_measure()` and, on `mode="ED"` and `"julia_live"`,
  `vev()` did not, so each of those now drops by exactly
  $1/\langle x|x\rangle$. On a 4-site Heisenberg chain at $x=2s$, $s$ its
  ground state: the KPM sum rule of $(S^z_0,S^z_0)$ goes from 0.999995 to
  0.249999, `evolve_and_measure(Sz0*Sz1)` at $t=0$ and ED's `vev(Sz0*Sz1)`
  from $-0.910684$ to $-0.227671$, `gs_energy_fluctuation()` on ED and
  `"julia_live"` from 9.696152 to 0 ($8\times10^{-9}$ on `"julia_live"`),
  and `"julia_live"`'s
  `get_excited_states(n=2)` from $[-6.464102, -0.957107]$ to
  $[-1.616025, -0.957107]$. `get_gs()` hands back $x/\lVert x\rVert$, and
  `get_rdm()` is the density matrix of the ray on all four backends, where
  it was $\rho/c^2$ for a state of norm $c$ (trace 0.25 at $c=2$). An
  explicit `wf=` to `evolve_and_measure`/`evolution_ABA` is left as given
  (findings 1 and 11).
- **Readers after `gs_energy(wf0=x, reconverge=False)` on a route that
  resolves to ED** (the chain's `mode="ED"`, a `mode="ED"` call, or v3's own
  fallback below three sites) measure x, where the ED ground state stayed
  on the chain: on a 4-site Heisenberg chain with $0.4S^z_0+0.3S^x_3$,
  `vev(Sz0)` goes from $-0.2286763037$ to $\langle x|S^z_0|x\rangle$,
  $-0.1089088536$ and $0.0390101221$ for the random states drawn on
  `"python"` and v3. The returned energy stays the lowest eigenvalue
  (finding 2).
- **A chain whose `mode` is `"DMRG"` no longer overrides a `mode="ED"`
  call**, so the ED cross-check is ED again: `gs_energy(mode="ED")` on an
  8-site Heisenberg chain in a staggered field at `maxm=2` goes from
  $-3.6734578613$ to the exact $-3.7040879103$, and
  `get_correlation_matrix(T=1)`, which is ED only, from occupations up to
  0.149 off the Fermi function to $5.6\times10^{-16}$, and from 4.2 s to
  0.1 s (finding 3).
- **`gs_energy(reconverge=True)` and `gs_energy(maxde=...)` on a solved
  chain** reach the solver, with `get_gs()` alike: the first sweeps from the
  stored state ($-4.2548001917$ to $-4.2548333556$ on a 10-site Heisenberg
  chain at `maxm=6` and one sweep, `"python"`), and the second refines
  ($-4.1431954920$ to $-4.2580352072$ at `maxm=3`, exact $-4.2580352073$)
  (finding 5, and the 2026-09-25 record's `maxde=` lead).
- **`"julia_live"` solves on a Hamiltonian whose couplings skip a site**
  leave the product state they were trapped in (§3): the ground energy goes
  from $-0.5$ to $-1.0$ on the Heisenberg chain on the even sites of 6 and
  from $-1.5$ to $-3.232051$ on an 8-site chain with next-nearest-neighbour
  exchange only, `Thermal_Spin_Chain` at $T=0$ from
  $\langle H\rangle=-0.5$ to $-1.0$, `gs_energy_generalized` from
  $-0.833333$ to $-1.496331$, the non-Hermitian solve from
  $-1.4997+0.0994i$ to $-3.219401$, and the lowest excited state from
  $-3.190943$ to $-3.232051$, all now exact. Every other `"julia_live"`
  solve moves within its run-to-run noise, and a warm 16-site solve costs
  about 40 per cent more (finding 7).
- **`gs_energy()` and `e0` after `gs_energy_generalized()`** are the
  generalized state's own energy (§4), where they were $\lambda$: on a
  4-site Heisenberg chain in a field of 0.3, from $-2.162930$ to $-1.385986$
  for $A=1+0.8S^z_0$, from $-0.808013$ to $-1.616025$ for $A=2\,\mathrm{Id}$
  and from $-1.115761$ to $-1.565921$ for $A=1.5+0.4S^z_0$, on `"python"`,
  v3 and `"julia_live"`, and on a non-Hermitian chain the biorthogonal
  energy of the pair. $\lambda$ is still the return value and is kept as
  `sc.lam_generalized`. Every KPM, CVM, ROOTN and TDZ spectrum afterwards
  moves rigidly by $\lambda-E_{w_g}$ onto the lines TD and EX already gave
  ($-0.777$ for $A=1+0.8S^z_0$, $+0.808$ for $A=2\,\mathrm{Id}$), and on
  `"julia_live"` the KPM window takes its lower edge from a solve of $H$,
  so the calls that raised "KPM moments diverging" return, and the
  $A=1+0.8S^z_0$ line's height goes from 1.2248 to 1.3037
  (findings 8 and 10).
- **`submode="CVM"` whenever $\eta\lVert B|\mathrm{GS}\rangle\rVert$ is
  small, and on exact solves at small $\eta$** (§6): on a 6-site Heisenberg
  chain with $0.3S^z_0$, every frequency of $(S^z_0,S^z_0)$ at an operator
  scale of $10^{-4}$ was the flat $\eta\langle AB\rangle/\pi$ (0.0159 after
  dividing by $\epsilon^2$) and is now the ED curve, and the same for
  $sH$ at $s=10^{-4}$; `get_kondo_spectrum(mode="DMRG", submode="CVM",
  order=2, T=0)` at its default `delta=2e-6` on three spins at
  $J=10^{-3}$ goes from about $10^{-8}$, flat, to the exact curve on the
  same grid (4.6981 at the ends, 0.6763 in the middle); the on-line point
  $\omega_0=0.52710639$ of that 6-site chain goes from
  $1.5915\times10^{-4}$ to 13.2636 at $\eta=2\times10^{-3}$ and from
  $1.5915\times10^{-5}$ to 132.6333 at $\eta=2\times10^{-4}$, as ED; and on a
  10-site chain on a 121-point grid, v3, 10 wrong points (up to 0.989 of
  the peak) become none at $\eta=0.05$, and 45 become none at $\eta=0.02$.
  At unit scale the tolerance tightens by $1/\lVert b\rVert$ (tenfold at
  $\eta=0.2$), so ordinary curves move within the old tolerance, the
  6-site one from $6.8\times10^{-7}$ to $6.7\times10^{-9}$ of the peak on
  `"python"`. A truncated solve still stops early and warns, and is not
  converged either way: a 20-site chain at `maxm=cvm_maxm=30` gives
  0.30333897 at $\omega=0.3$ where it gave the flat 0.01193662. The
  121-point grid takes 1214 s where it took 677 s, and a truncated point
  about twice as long. `applyinverse` stops relative to the norm of the
  vector it inverts against, so `CVM_explicit` on
  $(\epsilon S^z_0,\epsilon S^z_0)$ is $1.4\times10^{-6}$ of the peak off at
  every $\epsilon$ on v3, where it was 0.37 off at $\epsilon=10^{-4}$
  (findings 13 and 14).
- **Non-adjoint pairs under `CVM_explicit`, the non-Hermitian CVM and INV,
  and `cvm_solver="variational"`** return their own $C[A,B]$ (§6), where an
  absolute adjoint gate answered small pairs with another pair's density:
  `CVM_explicit` on $(10^{-2}S^x_0,10^{-2}S^y_1)$ from 2.004 of the peak off
  to $3.5\times10^{-5}$, the variational solver on
  $(10^{-2}S^z_0,10^{-2}S^z_3)$ from 2.124 to $5.9\times10^{-8}$, and ED's
  non-Hermitian CVM/INV on $(10^{-5}S^x_0,10^{-5}S^y_1)$ from 2.064 to
  $1.4\times10^{-15}$ (finding 16).
- **Every finite-temperature `Thermal_Spin_Chain` result** (§9). At
  $T=0.5$, on 10 sites ($\langle H\rangle$, `"python"` DMRG at `maxm=64`)
  from $-2.777892$ to $-3.166391$ (exact $-3.166396$), on 12 sites from
  $-3.287097$ to $-3.849203$ (exact $-3.849213$), on 6 sites on `mode="ED"`
  from $-1.674494$ to $-1.800863$, and on 3 sites from $-0.754223$ to
  $-0.769459$ (exact $-0.769460$), now the same at every scale and under a
  constant offset, where an offset of 20 gave 0.393848. On that 3-site
  chain in a field $B$ at $T=B$, $\langle S^z_{\rm tot}\rangle$ goes from $-0.089480$,
  $-0.013179$ and $-0.002545$ at $B=10^{-2}$, $5\times10^{-3}$ and
  $10^{-3}$ to $-0.231057$ at all three (exact $-0.231059$). A float $T$
  above 5 returns ($-0.071119$ at $T=5.5$) where it raised
  `ZeroDivisionError`, and a numpy scalar above 5 gets its temperature
  ($-0.055408$ at `np.float64(7.0)`) where it got $T=\infty$. The
  2026-09-25 thermal-bypass numbers move with it, `tc.MBChain.gs_energy()`
  at $T=1$ on 3 sites from $-0.416388$ to $-0.428230$ and `vev(Sz0*Sz1)`
  from $-0.069398$ to $-0.071372$. On `"julia_live"` a 4-site chain at
  $T=0.5$ comes out $2.3\times10^{-5}$ off where the other backends are
  $2\times10^{-8}$ off, since its band-edge solves returned a narrower
  width (not examined further). The cost grows with $\beta W$: 8.6 s where
  it was 4.9 s on 12 sites at $T=0.5$, and 19.8 s where it was 0.7 s for
  the field chain at $B=10^{-3}$ (findings 17 to 19).
- **`get_excited(n, mode="ED")` above 2000 states in small units** keeps
  every copy of a degenerate level: on a 12-site ferromagnet at
  $s=2\times10^{-7}$ the six-fold ground level came back with a magnon at
  $-2.715926s$ in place of a copy of $-2.75s$ in 8 of 8 calls, and the
  calls that stalled at $s=10^{-5}$ and $10^{-6}$ return. Below an infinity
  norm of 1 the levels also move at roundoff (finding 20).
- **`submode="ROOTN"` in small units** keeps its Lanczos basis, where it
  collapsed to the seed alone and returned one Lorentzian: the peak of
  $sC_s$ on a 6-site chain goes from 0.393371 to the unit-scale 0.157205 at
  $s=10^{-12}$ on ED and $s=10^{-10}$ on v3 (finding 21).
- **A KPM cross-correlator of two small operators** is no longer read as
  an autocorrelator: $C[\epsilon S^z_0,\epsilon S^z_3]/\epsilon^2$ at
  $\epsilon=1.2\times10^{-10}$ and $10^{-11}$ was $C[S^z_3,S^z_3]$, 2.070 of
  the peak off, and is now $C[S^z_0,S^z_3]$ to $10^{-14}$ on v2, v3 and
  `"python"` and to $7\times10^{-5}$, its run-to-run noise, on
  `"julia_live"`; `get_distribution` likewise, from 1.935 to $10^{-14}$
  (finding 23).
- **VUMPS in small units or next to a large constant**, on v3 and
  `"python"`, reports `converged=True` where it could not: TFIM at $D=8$
  from $s=0.1$ down to $10^{-8}$ and at constant offsets of $+10$, $+100$
  and $-100$ on v3; below about $s=10^{-6}$ the energy moves too, $e_0/s$
  from $-0.440126707917$ to $-0.440127030549$ at $s=10^{-8}$
  (finding 24).
- **v3 iDMRG in small units**: the energy density of TFIM is
  $4.8\times10^{-12}$ off the closed form at every $s$ down to $10^{-14}$,
  where it was 0.12 off at $s=10^{-13}$ and of the wrong sign at $10^{-14}$
  (finding 27).
- **NH-DMRG in small units or next to a large constant**, on v2, v3 and
  `"python"`, lands on the lowest level where it returned an excited one:
  $E_0/s$ on a 6-site Heisenberg chain with $0.3iS^z_0$ at $s=10^{-7}$ goes
  from $-0.9904-0.0568i$ to $-2.4610410218$ on v2, $E_0-c$ at
  $c=5\times10^5$ from $-1.9919-0.1109i$ to $-2.4610374186$ on `"python"`,
  and v3 at $s=10^{-13}$ from 0.570 off to $7\times10^{-15}$. At unit scale
  the tie window between Ritz values is about three times wider than it
  was, which moved nothing measured (findings 25 and 27).
- **v2 and v3 Hamiltonians whose duplicate terms cancel** are solved at the
  scale of the MPO they build (§2): $sH+S^z_0-S^z_0$ at $s=10^{-10}$ from
  0.09 to 0.53 off in $E_0/s$ on a 6-site Heisenberg chain to
  $5\times10^{-15}$. On v3 the XX chain at $J=1$, whose realified
  coefficients are 0.5, now also solves at scale 2, moving within its
  run-to-run noise (finding 28).
- **`"python"` operators whose largest coefficient is below 1** are built
  at unit scale (§2): a Heisenberg MPO at $s=10^{-7}$ goes from
  $1.1\times10^{-9}$ off relative to itself ($2.5\times10^{-4}$ at
  $s=2^{-40}$) to $10^{-15}$ at every scale, on numpy with OpenBLAS, and
  every `"python"` `MultiOperator*MPS` in small units moves with it,
  $\lVert(sH)|u\rangle\rVert/s$ reading 0.897242774928 at every $s$ down to
  $10^{-16}$, where it read 1.432 at $10^{-16}$ (finding 30).
- **`mpsalgebra.lowest_energy_arnoldi` in small units**, a non-default
  route: $E_0/s$ on a 4-site Heisenberg chain at $s=10^{-6}$ goes from
  $-1.6150192976$ to $-1.6160253791$, and at $s=10^{-7}$ with
  `delta=s*1e-3` from $+0.3510023868$ to $-1.6160254038$ (lead
  `scale-arnolditk-absolute-stop`).

**New attributes.** `sc.lam_generalized` holds $\lambda$ after
`gs_energy_generalized()`; `tc.anneal_step` and `tc.anneal_order` set the
thermal stepper (§9); `sites.SETTINGS` is the explicit list of the 36
settings a chain constructor accepts (§1); and on `mode="ED"` the state
gains `normalize(tol=)` and `norm()`, as on the MPS backends (§2).

**Quieter changes.** With `verbose` on, v2 and v3 print one line before
each solve that runs at $2^kH$, saying that the energies logged below it
are $2^k$ times those of the operator (finding 29); the returned energies
never were. NH-DMRG's convergence certificate is relative to the
Hamiltonian's units and constant, so a run it used to accept silently in
small units or next to an offset now retries up to `ntries` times and
warns with its "best relative residual" (finding 26). `State.normalize()`
on `mode="ED"` prints the MPS backends' warning where it was silent below
its floor (finding 22).

**And these now raise, or now work, where they used to do something
else:**

| Call | Was | Now |
|---|---|---|
| a constructor keyword naming one of nine attributes that are not settings (`hopping`, `hubbard`, `pairing`, `exchange`, `fields`, `resorder`, `resordered_indexes`, `hubbard_matrix`, `fit_td`) | accepted; `hubbard=2.0` added $2\,\mathrm{Id}$ to the Hamiltonian | `TypeError` naming the key, `Thermal_Spin_Chain` included; the last five are gone from the chain (finding 4) |
| a misspelled keyword (`wf=`, `reconverg=`) to `gs_energy()`, `get_gs()` or `get_excited(n=1)` on a solved chain | swallowed, the stored answer returned | `TypeError`, as on a fresh chain (finding 5) |
| a reader on a chain with no Hamiltonian | an error not naming the cause | `ValueError` saying to call `set_hamiltonian()` (finding 6) |
| `gs_energy(H=H2)` on a non-Hermitian chain | H2's pair stored as the chain's state | `TypeError` naming `nhdmrg(H=...)` (finding 9) |
| `get_gs(wf0=x, reconverge=False)` on a route that resolves to ED | `TypeError` | $x/\lVert x\rVert$ (finding 2) |
| an MPS as `wf0=`, or a misspelled keyword, to `gs_energy()` on a route that resolves to ED | ignored | `TypeError` (finding 2) |
| a set state whose norm is zero or not finite | `None` stored as the state | `ValueError` (finding 1) |
| an unknown single-site name on an ED spin chain (`ISy`, `Adag`) | `RuntimeError: No active exception to reraise` | `ValueError` naming it (finding 6) |
| `mode="ED"`, `submode="ED"` with the averaged manifold wider than `delta` | silent | `RuntimeWarning` (finding 15) |
| `Thermal_Spin_Chain` with a negative or NaN `T` set after construction | the ground state | `ValueError` (finding 19) |
| `get_rdm(i=ns-1)` on `"julia_live"` | `BoundsError` | the matrix (finding 12) |
| `get_rdm()` or the bond entropy under an ED mode on `"julia_live"` | the DMRG matrix, or `AttributeError` | `NotImplementedError`, as on the other backends (lead `session-julia-ed-guard-carveout-rdm-bond-entropy`) |
| `CVM_explicit`, or the non-Hermitian CVM/INV, on a non-adjoint pair | `NotImplementedError` | $C[A,B]$ (finding 16) |

Results from before this pass are not comparable wherever a number above
moved; everywhere else the calculation is the one it was.

### The 2026-09-26 Kondo exchange diagram

Not an audit: a check of the Kondo spectrum against the paper found that
its eq. 25 puts the third-order Kondo term's exchange-diagram log at
$eV+\epsilon_{im}$ where the second-order T-matrix of its own
Hamiltonian puts it at $eV-(\epsilon_f-\epsilon_m)$ (§17, "Where the
exchange log sits"). `get_kondo_spectrum` now uses the latter on both
`mode="ED"` and `mode="DMRG"`. Every third-order spectrum with an
inelastic transition out of an occupied state moves, by up to 5–20% of
the third-order term at $T=0$ and a few per cent of the total at 1 K
(the Fig. 7d 10 T overshoots from 1.226/1.183 to 1.248/1.203); zero-bias
values, a free $S=1/2$ at $B=0$, `order=2` and the potential term do not.
Results from before are not comparable. The lower-level
`kondospectrumtk.twotime.kondo_term_from_two_time` now takes
`(t, G, Gx)` triples and raises on the old `(t, G)` pairs, and its
kernel `K_W` is replaced by the one-sided `K_F`.

### The 2026-09-26 review of `"python"` real-time evolution

Not an audit: a check of `itensor_version="python"`'s real-time evolution
against exact propagation. The integrators themselves were right. Two-site
and one-site TDVP reproduce $e^{-iH\,dt}$ to $10^{-11}$ at full bond
dimension on spin-1/2 (long-range and Dzyaloshinskii-Moriya), spin-1,
Jordan-Wigner fermion and boson chains, at real and complex $dt$, and below
full bond dimension they track ED as closely as v3 does. Three things
around them changed numbers:

- **The Krylov exponentiator** behind every `"python"` TDVP route (TD,
  TDZ, `evolve_and_measure`, `evolution_ABA`, METTS, the infinite-chain
  window, the Kondo two-time construction) now stops when
  $|z|\,\beta_k\,|[e^{zT_k}]_{k1}|<10^{-10}$ for $e^{zH}$, an error
  relative to the evolved vector and free of units. It used to stop on
  $\lVert v\rVert\,\beta_k\,|[e^{zT_k}]_{k1}|$, which carries the units of
  $H$ and the norm of the state, and it returned a step its 50-vector
  budget could not take unconverged (0.31 off at $|z|$ times the spectral
  width equal to 100); such a step is now split. This closes the "absolute
  Krylov error goal" the 2026-09-25 record left open, on `"python"`:
  $C[\epsilon S^z_0,\epsilon S^z_3]/\epsilon^2$ under `submode="TD"` was
  $8.8\times10^{-2}$ off $C[S^z_0,S^z_3]$ at $\epsilon=10^{-8}$ and $0.24$
  at $10^{-9}$, and is $10^{-11}$ at both, and a Néel quench under $sH$ at
  $dt/s$ was $0.48$ off ED at $s=10^{-12}$ and is at the $2.9\times10^{-6}$
  of $s=1$. At ordinary units results move below $10^{-10}$. v3's
  `applyExp` is vendored ITensor and keeps its absolute goal.
- **`evolve_and_measure` and `evolution_ABA`** on `"python"` (TDVP,
  TDVP_GSE, TEBD) and on v3 TEBD now put the state back to its initial
  norm after every step, as v3's TDVP and every `"julia_live"` loop
  already did. Real-time evolution conserves the norm, so what the
  truncation removes from it is error, and it used to scale every later
  $\langle\psi(t)|O|\psi(t)\rangle$: on a 12-site Néel quench at
  `maxm=8`, $\langle\psi|\psi\rangle$ fell to 0.982 by $t=5$ under TDVP and
  to 0.970 under TEBD, and $\langle H\rangle$ drifted by
  $4.9\times10^{-2}$ where v3 drifted by $1.8\times10^{-4}$ (now the
  same $1.8\times10^{-4}$ on both). Nothing moves where nothing truncates,
  and the `quench` loops behind `submode="TD"` always renormalized.
- **`tevol_method="TDVP_GSE"` on `"python"`** now normalizes each Krylov
  vector $H^k|\psi\rangle$ before choosing the new directions, as v3's
  expansion does. The choice used to depend on $\lVert H\rVert$: that
  8-site Néel quench was $1.0\times10^{-4}$ off ED at $s=1$ and
  $1.3\times10^{-3}$ at every $s\le10^{-4}$, and is $2.7\times10^{-5}$ at
  all of them.

Results from before are not comparable where these apply. The
regressions are in `tests/test_pyitensor_time_evolution_review.py`. One
thing was measured and left as it is, on every backend: under
`TDVP_GSE` the bond dimension stops growing after `tdvp_gse_sweeps` steps
(3 by default), so a quench whose entanglement keeps growing saturates.
On the 12-site Néel quench the largest error up to $t=5$ stays near
$2\times10^{-3}$ from `maxm=16` to 64, where `"TDVP"` reaches $9\times10^{-7}$ at
`maxm=64`. Raise `tdvp_gse_sweeps` to the number of steps, or use
`"TDVP"`, for such a run.
