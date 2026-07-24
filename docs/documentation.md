# DMRGPY Documentation

## 1. Overview

DMRGPY is a Python library for quasi-one-dimensional quantum lattice
models — spin chains, spinless and spinful fermionic chains, Kondo
chains, parafermion chains, and bosonic chains — solved with matrix
product state (MPS) methods, principally the density matrix
renormalization group (DMRG). Most calculations can also be run with
exact diagonalization (ED) on small systems, so DMRG results can be
cross-checked against an exact reference.

The library exposes a single, backend-agnostic Python API
(`src/dmrgpy/`). The same physics — Hamiltonian, observables, sites — is
written once, and DMRGPY dispatches the actual numerical work to
whichever solver is requested or available:

- **DMRG**, run in-process through a compiled C++ extension
  (`itensor_version=2` or `3`, built against a vendored copy of the
  [ITensor](https://itensor.org/) library);
- **DMRG**, run through a live, in-process Julia session using
  ITensors.jl (`itensor_version="julia_live"`);
- **DMRG**, run in a from-scratch, pure-Python/NumPy/SciPy
  reimplementation of the ITensor v3 API subset DMRGPY needs
  (`itensor_version="python"`) — no compiler or external language runtime
  required at all;
- **ED**, a pure-Python fallback for small systems (dense/sparse matrices
  via NumPy/SciPy), used automatically whenever a requested compiled
  backend is not available, and useful on its own as a correctness
  reference for DMRG.

## 2. Installation

There are two independent install paths, and they cover different
backends:

- **`pip install dmrgpy`** (via the repo's `pyproject.toml`) installs the
  pure-Python part of the package only — the ED backend and the
  pure-Python `itensor_version="python"` DMRG/TDVP backend both work
  immediately, no compiler needed. The compiled C++ backends
  (`itensor_version=2`/`3`) are *not* on PyPI: they're built against a
  ~400MB vendored copy of ITensor, and `pyproject.toml` excludes
  `mpscpp2`/`mpscpp3` from the distributed package entirely (see its own
  comments). Without them, DMRGPY still works end to end — it just falls
  back to ED/`"python"` (see mode.py) — but if you want the compiled
  backends you need the git-clone path below instead.
- **A git clone + `python install.py`** is required for the compiled C++
  backends (and is how this repository is developed). DMRGPY is used
  directly from `src/`, added to `PYTHONPATH` (or symlinked into
  `site-packages`) by the installer.

```bash
python install.py                        # compiles ITensor v3 (mpscpp3, the default) + its pybind11 extension
python install.py --itensor-version=2     # compiles ITensor v2 (mpscpp2) instead
python install.py --itensor-version=both  # compiles both v2 and v3 backends
python install.py --gpp=g++-6             # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install.py --doctor                # check build requirements only, skip the build
python install_julia.py                   # alternative/additional: Julia backend
```

`install.py` first verifies every build requirement by actually trial
compiling and trial linking it (a C++ compiler, LAPACK/BLAS, `pybind11`,
`make`) rather than just checking version strings, auto-detecting a
working compiler (preferring a conda-provided one when running under
conda Python, since the extension shares a process with conda's own
numpy/scipy). Only once every requirement is confirmed does it compile
the vendored ITensor static library and the pybind11 extension for the
requested backend version(s). If nothing is compiled at all (the
extension build fails, or the user only ran `install_julia.py`), DMRGPY
still works: `itensor_version="python"` requires no compiled extension,
and every call transparently falls back to ED when a requested compiled
backend isn't present.

No compiler or pybind11 is needed at all to use `itensor_version="python"`
— only NumPy/SciPy/SymPy/`numba` (all four are core `pip install dmrgpy`
dependencies; `numba` is required rather than optional because
`algebra/kpmextrapolate.py` imports it at module level, and that module
is reached unconditionally from `manybodychain.py` via
`dynamics.py`/`kpmdmrg.py` — this is separate from `pyitensor`'s own
opt-in JAX/numba-accelerated matvec kernel, which *is* genuinely optional
and off by default, see §5.4). `jax` is optional either way.

## 3. Quick start

```python
from dmrgpy import spinchain

spins = ["S=1/2" for i in range(30)]      # 30-site spin-1/2 chain
sc = spinchain.Spin_Chain(spins)          # create the chain object

h = 0                                     # build the Hamiltonian symbolically
for i in range(len(spins) - 1):
    h = h + sc.Sx[i] * sc.Sx[i + 1]
    h = h + sc.Sy[i] * sc.Sy[i + 1]
    h = h + sc.Sz[i] * sc.Sz[i + 1]
sc.set_hamiltonian(h)

print("Ground state energy:", sc.gs_energy())
```

Every calculation follows this pattern: build a chain object for the
model of interest, build a Hamiltonian (and any observables) out of the
chain's own operators, call `set_hamiltonian`, then call one of
`Many_Body_Chain`'s methods (`gs_energy`, `get_gs`, `vev`, `get_excited`,
`get_dynamical_correlator`, time-evolution methods, ...). Most methods
accept a `mode="DMRG"|"ED"` keyword so a result can be cross-checked
against exact diagonalization on small systems, e.g.:

```python
print("Energy with DMRG:", sc.gs_energy(mode="DMRG"))
print("Energy with ED:",   sc.gs_energy(mode="ED"))
```

The `examples/` directory contains 100+ self-contained scripts, one per
physical model or feature, and doubles as the project's de facto
regression suite (there is no separate unit-test suite for the Python
code). See `examples/readme_examples/` for the snippets shown in
`README.md`, and `examples/v2_VS_v3_*` /
`examples/backend_timing_gs_energy/` for scripts that directly compare
backends against each other on correctness and timing respectively.

## 4. Architecture

### 4.1 Entry-point classes

Each physical model is a thin subclass of `Many_Body_Chain`
(`manybodychain.py`), which holds all shared state — bond dimension,
sweep count, cutoffs, mode, temp folder, etc. — and implements the
generic operations: `set_hamiltonian`, `gs_energy`, `vev`, `get_gs`,
correlators, time evolution, entanglement entropy, and more. Model
modules just define the local Hilbert space and convenience operators
for that statistics:

| Module | Class | Notes |
|---|---|---|
| `spinchain.py` | `Spin_Chain` | alias of `Many_Body_Chain` itself |
| `fermionchain.py` | `Fermionic_Chain` | spinless fermions; `C`/`Cdag`/`N`, Jordan-Wigner string `F` |
| `spinfermionchain.py` | `Spin_Fermionic_Chain` | spinful fermions |
| `bosonchain.py` | `Bosonic_Chain` | bosonic sites |
| `parafermionchain.py` | `Parafermion_Chain` | parafermion sites |
| `mixedchain.py` | `Mixed_Spin_Fermion_Chain` | mixes spin sites and spinful-fermion locations in one chain |

All of the above build `self.sites`, the list of per-site type codes
passed to `Many_Body_Chain.__init__`, from a **uniform** repetition of a
single type code (e.g. `[0]*n` for `Fermionic_Chain`). Nothing in
`sites.py`/the C++ or pyitensor `Chain` constructors actually requires
this — the site-construction layer (`mpscppN/get_sites.h`'s
`SpinX(std::vector<int> const&)`, `pyitensor/sites/siteset.py::SiteX`)
already accepts an arbitrary, heterogeneous per-site type-code list.
`bosonchain.py::SpinBoson_Chain` and `mixedchain.py::Mixed_Spin_Fermion_Chain`
are the two model classes that build such a heterogeneous list directly
(spin+boson and spin+spinful-fermion respectively), each following the
same pattern: build every operator that could be meaningful at *any*
site, then mask the invalid combinations to the literal int `0` per
site. `Mixed_Spin_Fermion_Chain` represents a spinful-fermion location as
a pair of physical spinless-fermion (type-code 0) sites, exactly like
`fermionchain.Spinful_Fermionic_Chain`, so it reuses the existing
`C`/`Cdag`/`A`/`Adag`/`F` operator vocabulary and Jordan-Wigner code
unchanged. It currently only supports `itensor_version` `3` and
`"python"`: a Jordan-Wigner string that has to cross a non-fermionic
(spin) site relies on `SiteSet::op()`'s "unrecognized `F` request
resolves to identity" fallback, present in ITensor v3
(`mpscpp3/ITensor/itensor/mps/siteset.h`) and in pyitensor's own
`SiteType.matrix()` (`pyitensor/sites/base.py`), but not in ITensor v2
(`mpscpp2/ITensor/itensor/mps/siteset.h`), which would hard-abort
instead.

`fermionchain.Spinful_Fermionic_Chain_Native` uses site-type code `1`
instead (`HubbardSite`/`ElectronSite`, ITensor v3's own stock 4-state
site, `mpscpp3/ITensor/itensor/mps/sites/electron.h`) -- one
tensor-network site per physical site instead of `Spinful_Fermionic_Chain`'s
two interleaved type-0 sites. `mpscppN/get_sites.h`'s site-type dispatch
already had a branch for code `1` (`SpinX`'s `HubbardSite` case, both
constructors) before this class existed; only `mpscpp3` (the C++ side)
and ED (`pyfermion.mbfermion.MBFermion.get_operator`, extended to
understand `Cup`/`Cdn`/`Cdagup`/`Cdagdn`/`Nup`/`Ndn`/`Ntot` at a physical
site index) actually wire it up today -- `mpscpp2` never registered a
Hubbard/Electron site, and neither pyitensor nor `mpsjulialive`
implement site-type code 1. See `fermionchain.py`'s class docstring for
why this alternative isn't a general-purpose speedup over
`Spinful_Fermionic_Chain` despite the halved site count (a directly
measured result, not a theoretical claim).

`bosonchain.Bosonic_Chain` uses a *range* of type codes rather than a
single one: `100+dim` per site, `dim` being the requested local Hilbert
space dimension (`Bosonic_Chain(n, maxnb=[...])`, default `dim=4`, i.e.
code `104`). This lets `mpscpp3/get_sites.h`'s `BosonFourSite`
(`extra/bosonfour.h`, generalized from a hardcoded 4-level site to an
`Args("MaxOcc",...)`-driven arbitrary occupation cutoff) build a site of
the actual requested dimension instead of always the fixed 4-level one
regardless of `maxnb` — previously a real, silent bug: the DMRG session
only ever saw type code `104` no matter what `maxnb` was, so DMRG and ED
results silently diverged for any non-default value. `itensor_version=3`
and `itensor_version="python"` both understand the extended `102..199`
range today — pyitensor's `siteset.py::_site_class()` routes it to
`sites/boson.py::get_boson_site(dim)`, a small factory that builds (and
caches) a `SiteType` subclass per requested dimension the same way
`BosonFourSite` is generalized on the C++ side, since pyitensor's
`SiteType` machinery is otherwise entirely class-level (no per-instance
state) and has no other way to parameterize a site's dimension.
`mpscpp2` and `mpsjulialive` still only handle the single code `104`, so
a non-default `maxnb` should be run under `itensor_version=3` (the
default) or `"python"` for cross-backend agreement. ED
(`pyboson/boson.py::BosonChain`) always builds the exact requested
per-site dimension regardless of backend.

### 4.2 Operator representation: `MultiOperator`

Hamiltonians and observables are built as `MultiOperator` objects
(`multioperator.py`) — sums of products of named single-site operators,
e.g. `Sx[i]*Sx[j]`. This is a backend-agnostic intermediate
representation: the *same* `MultiOperator` is later either

- converted into a plain list of `(coefficient, [(opname, site), ...])`
  terms and written directly into an ITensor `AutoMPO`/`HTerm` by the
  C++/pyitensor DMRG backends, or
- converted into a sparse matrix via `multioperator.MO2matrix` for the ED
  backend (`edtk/edchain.py`).

`multioperatortk/` holds the supporting machinery: Jordan-Wigner string
threading for fermionic operators, static/long-range operator
construction, and sympy-based symbolic building. Spinless fermionic
operators (`C`/`Cdag`, one fermionic mode per site, `Fermionic_Chain`/
`Spinful_Fermionic_Chain`'s interleaved sites) are threaded by
`jordanwigner.py`; flavor-resolved operators on a native spinful site
(`Cup`/`Cdn`/`Cdagup`/`Cdagdn`, two fermionic modes per site,
`Spinful_Fermionic_Chain_Native`) are threaded by the separate
`jordanwigner_spinful.py` instead, dispatched by operator name in
`multioperator.jordan_wigner()`. The two never mix within one term (the
operator name sets are disjoint), and `jordanwigner_spinful.py` always
uses the generic per-factor dressing recipe (no separate optimized
2-/4-point path the way `jordanwigner.py` has for plain `C`/`Cdag`),
since that recipe is already exact for an arbitrary product/order of
factors — see its module docstring.

### 4.3 Backend dispatch

`mode.py` decides, per call, whether a calculation runs via DMRG or ED
(the `get_mode`/`run` functions): DMRG unless `self.mode` forces ED, or
unless `self.itensor_version` is `2` or `3` and the corresponding
compiled pybind11 extension isn't available (`cppext.available(version)`)
— in which case DMRGPY silently falls back to ED.
`itensor_version="python"` never falls back this way, since it has no
compiled-extension precondition at all. Most public `Many_Body_Chain`
methods accept `mode="DMRG"|"ED"` so results can be cross-validated
between the two solver families.

| `itensor_version` | Engine | Requires | Fallback |
|---|---|---|---|
| `2` | ITensor v2, in-process C++ (`mpscpp2`) | compiled pybind11 extension | ED |
| `3` (default) | ITensor v3, in-process C++ (`mpscpp3`) | compiled pybind11 extension | ED |
| `"python"` | pure-Python `pyitensor/` | NumPy/SciPy only | none |
| `"julia_live"` | live in-process Julia session (`mpsjulialive/`), ITensors.jl | `juliacall`/PythonCall.jl (self-provisions Julia) | none (feature-by-feature; missing methods simply aren't implemented) |

Regardless of `itensor_version`, if `self.mode` is forced to `"ED"`, or
the requested backend isn't available, calculations run through
`edtk/edchain.py` (`EDchain`): dense/sparse operators built directly in
Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`, `pyzn/` provide
per-statistics many-body operator construction; `edtk/one2many.py`
promotes single-site operators to the full Hilbert space), diagonalized
with `scipy.sparse.linalg`.

### 4.4 The in-process C++ backend (ITensor v2 and v3)

`mpscpp2/bindings.cc` and `mpscpp3/bindings.cc` each compile a pybind11
extension (`mpscpp2/_dmrgcpp*.so` / `mpscpp3/_dmrgcpp*.so`) that runs DMRG
entirely in-process against its own vendored ITensor copy — no task
files, no operator files, no subprocess, no temp directory. `mpscpp3` is
a line-by-line port of `mpscpp2` to the ITensor v3 API (which merges
`IQTensor`/`IQIndex` into a single `ITensor`/`Index` type carrying QN
blocks on the Index), exposing an identical `Chain` class and Python
surface. `cppext.py` lazily imports whichever compiled extension a chain
needs (`get_backend(version)`, `available(version)`).

`mpscppN/chain_session.h`'s `Chain` class is the stateful, in-process
session — one instance per `Many_Body_Chain` (held as `self._session`) —
with a method per DMRG task: `gs_energy`, `vev`, `apply_operator`,
KPM/CVM dynamical correlators, `correlation_matrix`, `reduced_dm`,
`excited_states`, time evolution, and more.

The session caches its ground-state energy and (for KPM) both band
edges, but `Chain::set_hamiltonian` unconditionally invalidates those
caches — so the Python side (`groundstate.py::gs_energy_single`) only
re-sends the Hamiltonian when its `to_terms()` output, any solver
parameter a re-run would pick up (`maxm`, `nsweeps`, `cutoff`, `noise`,
the effective MPO bond dimension), or the session object itself actually
changed since the last send (`_session_ham_cache`). Repeated calculations on an
unchanged Hamiltonian (e.g. successive `get_dynamical_correlator` calls,
each of which re-verifies the ground state) then hit the session's
caches instead of re-running warm DMRG sweeps and band-edge solves. Code
paths that *want* a fresh solve of the same Hamiltonian either pass
through `restart()`/`set_hamiltonian` (which force DMRG via
`skip_dmrg_gs=False`) or clear `_session_ham_cache` explicitly (see
`groundstate.py`'s best-of-`n` loops).

Notable, deliberate implementation details (not bugs to "fix"):

- `mpscpp3` builds every site with `ConserveQNs=false` and starts DMRG
  from an actual `randomMPS`, not a plain product state — matching v2's
  real (unconstrained-search) behavior rather than ITensor v3's stricter,
  QN-conserving-from-the-start convention, at the cost of losing the QN
  block-sparsity speedup for `itensor_version=3`.
- `mpscpp3`'s (and `pyitensor`'s) real-time MPS evolution defaults to a
  proper two-site TDVP integrator (vendored under `mpscpp3/TDVP/`;
  `pyitensor/tdvp.py` for the pure-Python backend), selectable via
  `Many_Body_Chain.tevol_method` (default `"TDVP"`, `itensor_version` 3 or
  `"python"` only); `mpscpp2` has no TDVP and always uses a hand-rolled
  2nd-order Taylor expansion of `exp(-i dt H)` as an MPO instead.
- `tevol_method="TDVP_GSE"` (`itensor_version` 3 or `"python"` only,
  same support as plain `"TDVP"`) runs one-site TDVP with Krylov *global
  subspace expansion* (GSE) beforehand for the first `tdvp_gse_sweeps`
  steps (default 3) — the scheme of Yang & White, arXiv:2005.06104/Phys.
  Rev. B 102, 094315 (2020): a Krylov subspace `{psi, H*psi, H^2*psi,
  ...}` of dimension `tdvp_gse_krylov_order` (built via repeated MPO
  application) enlarges the MPS's local bond bases *without changing the
  represented state*, giving one-site TDVP (which conserves bond
  dimension exactly on its own, unlike two-site) room to grow into the
  entanglement the subsequent evolution generates. `itensor_version=3`
  (`Chain::global_subspace_expand()`/`Chain::tdvp_step(...,num_center=1)`,
  `mpscpp3/TDVP/basisextension.h`, vendored unmodified from upstream) and
  `"python"` (`pyitensor/gse.py`) implement this. A v2-API port
  (`mpscpp2/TDVP/`) was attempted and briefly landed: numerically correct
  (verified against ED and against v3/`"python"`, exact agreement on a
  6-site cross-check), built on a from-scratch Lanczos `applyExp` (v2's
  own `itensor/iterativesolvers.h` never had one) since v2's stock
  `LocalMPO` also turned out to have no working one-site local
  Hamiltonian at all (confirmed in `LocalMPO::position()`, which
  unconditionally builds a two-site block regardless of `numCenter()`),
  paired with two-site TDVP instead as a result. It was reverted after a
  severe, unresolved performance regression at `n≳10` sites (the
  dynamical-correlator step didn't finish in 25 minutes at `n=12`, versus
  under a second for the same computation on `itensor_version=3`) that
  couldn't be root-caused in the time available — see git history around
  `mpscpp2/TDVP/` if picking this up again.
- All three session backends implement a real non-Hermitian DMRG
  (`Chain::nhdmrg`, driven by `nhdmrg.py`, exposed as
  `Many_Body_Chain.nhdmrg()`): a port of ITensorNHDMRG.jl's
  "onesided"+"fidelity" algorithm that optimizes a biorthogonal
  left/right eigenpair of a non-Hermitian `H`, targeting the eigenvalue
  with smallest real part. `mpscpp3/chain_session.h`'s copy is the
  annotated original (including two deliberate deviations from the
  reference for Re-degenerate spectra); `mpscpp2`'s is its v2-API
  back-port, and `pyitensor/nhdmrg.py` is the pure-Python port (built on
  `dmrg.py`'s environment/matvec machinery; unlike `dmrg.py` it *does*
  implement the noise term, because for NH-DMRG it measurably matters).
  `groundstate.py`'s non-Hermitian `gs_energy` branch routes to it for
  `itensor_version` 2, 3 and `"python"`; the MPS Arnoldi route
  (`algebra/arnolditk.py`) remains as the fallback for other backends.
  The adjoint MPO is built from `MultiOperator.get_dagger()`'s terms on
  the Python side, and since the non-Hermitian energy is not variational,
  `nhdmrg.py` certifies each run by the eigen-residual and redraws a
  fresh random start when a run stalls.
- A few pre-existing bugs in the original (pre-refactor) file-based
  backend are deliberately reproduced rather than silently fixed —
  see the call-site comments in `chain_session.h` for both versions.
- `itensor_version` 2, 3 (and `"python"`) expose MPO algebra directly on
  already-built operators: `StaticOperator` supports `+`, `-`, unary `-`,
  and scalar `*`/`/` between two `StaticOperator`s (or a scalar), on top
  of the pre-existing `*` for `StaticOperator*MPS`/`StaticOperator`
  contraction. `A + B` (`Chain::sum_operators`, a public wrapper each
  backend already had internally as a private `sum_mpo()` helper used by
  `custom_exp()`/`evoloperator()`) is a compressed direct sum at the
  tensor-network level — algorithmically the same construction as
  ITensorMPS.jl's `+(::MPO, ::MPO)` (`abstractmps.jl`'s default
  `"densitymatrix"` algorithm) — not a MultiOperator-level symbolic sum
  (which already existed, see §4.2, and remains the preferred way to
  combine Hamiltonians before ever building an MPO). It exists for
  combining operators that only exist as already-built `StaticOperator`s,
  e.g. two independently constructed products or exponentials.
  `"julia_live"` doesn't implement this yet (`StaticOperator.__add__`
  raises `NotImplementedError` rather than silently doing something
  backend-specific).

### 4.5 The pure-Python backend (`pyitensor/`)

`pyitensor/` is a from-scratch, pure-Python/NumPy/SciPy reimplementation
of exactly the ITensor v3 API subset `mpscpp3/chain_session.h` uses: an
`Index`/`ITensor` tensor core, the same ten site types, an `AutoMPO` term
compiler with its own Jordan-Wigner threading, MPS/MPO algebra, two-site
Lanczos-based DMRG (ground and excited states via the overlap-penalty
method), two-site Krylov-based TDVP, and a `chain.Chain` facade exposing
the identical method surface as the compiled backends' `self._session`.
Because of that shared surface, every call site elsewhere in DMRGPY that
accepts `itensor_version in (2, 3)` treats `"python"` as a third,
always-available option with no separate code path.

It exists so DMRG/TDVP work with zero compiler/pybind11 dependency, at
the cost of being slower than compiled ITensor by default (no
block-sparsity, no JIT) — see §5 for how much slower in practice, and how
`numba`/`jax` narrow that gap.

### 4.6 The Julia backend (`mpsjulialive/`)

`itensor_version="julia_live"` drives a live, in-process Julia session
(via [`juliacall`/PythonCall.jl](https://github.com/JuliaPy/PythonCall.jl),
`mpsjulialive/juliasession.py`) with its own set of modules mirroring the
top-level ones (`mpsjulialive/groundstate.py`, `vev.py`, `mps.py`,
`mpo.py`, `dynamics.py`, ...) that talk to Julia's ITensors.jl instead of
the C++ pybind11 extension. This is an independent implementation, not a
shared protocol with the C++ path — a feature missing on one side simply
isn't implemented there yet.

`juliacall` replaced an earlier PyJulia-based bridge: PyJulia requires the
Julia build's PyCall to be linked against a matching libpython (the same
class of ABI landmine documented in §2 for the C++ extension) and needs a
`compiled_modules=False` workaround for slow/broken precompilation.
`juliacall` avoids both, and self-provisions Julia plus the required
packages (`ITensors`, `ITensorMPS`, `ITensorNHDMRG`, pinned in
`src/dmrgpy/juliapkg.json`) into its own managed project the first time
`mpsjulialive` is imported — no manual Julia download step. One porting
detail worth knowing: unlike PyJulia, `juliacall` does *not* implicitly
convert a Python `list` of `str` into a Julia `Vector{String}` (it wraps
it as a lazy `PyList{Any}`, which fails to dispatch against
strictly-typed Julia functions like `sites.jl`'s `get_sites`); every call
site that used to rely on that now goes through
`juliasession.to_julia_strvec()`, which forces the conversion explicitly
via `juliacall.convert`.

The KPM dynamical correlator's Chebyshev-moment recursion runs natively
in Julia (`mpsjulialive/kpm.jl`'s `kpm_moments_full`/
`kpm_moments_accelerated`, driven through `kpm.jl`'s own `apply_op`/
`mpsalgebra.jl`'s `summps`), one call per correlator rather than one
Python↔Julia round trip per Chebyshev step. `mpsjulialive/dynamics.py`
only builds the scaled-Hamiltonian MPO and the two seed wavefunctions in
Python (mirroring `pyitensor/chain.py`'s own pure-Python KPM algorithm,
since `julia_live` has no single `Chain` session object comparable to
`self._session` on the C++/pure-Python backends for `kpmdmrg.py` to
dispatch to directly), then hands the whole moment loop to Julia in one
call. Measured directly on a 30-site Heisenberg chain: eliminating the
per-step round trip only closed part of the gap to the compiled ITensor
v3 backend (v3: 23.4s, Julia native loop warm: 34.0s, vs. 37.7s for the
earlier per-step Python loop) — most of the remaining time is genuine
Julia-side tensor-contraction/truncation cost (`alg="densitymatrix"`
SVD-based `add`/`contract`), not Python↔Julia marshaling overhead.

**Follow-up optimization**: `kpm_moments_full`/`kpm_moments_accelerated`
originally called `mpsalgebra.jl`'s `applyoperator()` once per Chebyshev
step, which does *two* truncated MPO-MPS contractions per call
(`contract(A,psi)` plus a second contraction against an identity MPO —
the same "fixes a bug" workaround documented above for `tdvp.jl`'s
`apply_clean`, which switched to `ITensorMPS`'s own higher-level
`apply()` for the same reason, on the TDVP side, to fix a *different*
bug). `mpsalgebra.jl`'s `apply_op()` does the same swap, cutting one of
the two contractions per moment (not a clean 2×, since the identity-MPO
contraction is against a trivial bond-dimension-1 operator, much cheaper
than the contraction against the real, bond-growing scaled Hamiltonian).
Measured directly on the same 30-site chain: 34.5s → 28.4s warm (~18%
faster), narrowing the gap to v3 from ~1.47× to ~1.21× slower. Peak
positions confirmed unchanged against the existing ED cross-check after
the change. `apply_op()` is shared between `kpm.jl`'s recursion and
`tdvp.jl`'s `apply_clean()` (which layers its own `orthogonalize!()` on
top, needed for `tdvp()`/`dmrg()` inputs but not for KPM's) — it lives in
`mpsalgebra.jl`, loaded before both, rather than as two independently
maintained near-duplicates. `mpsalgebra.jl`'s `summps()` (the `add()`
wrapper the Chebyshev recursion sums consecutive Chebyshev vectors with)
now also takes an explicit `cutoff` argument: `add()`'s own default
(`1e-15`) is three orders tighter than dmrgpy's usual configured cutoffs,
so every moment step used to silently keep far more Schmidt values than
`self.kpmcutoff` called for, up to `kpmmaxm`. The KPM seed vectors
(`psi1`/`psi2` in `mpsjulialive/dynamics.py`) are now also built via
`apply_op()` rather than the older `applyoperator()`, closing the one
spot the original optimization pass left on the old primitive.

Real-time TDVP evolution (`timedependent.py`'s `evolve_and_measure_dmrg`/
`evolution_dmrg_DC`, submode `"TD"`) is implemented the same way —
`mpsjulialive/tdvp.jl`'s `evolve_and_measure_tdvp`/`quench_tdvp` run the
whole nt-step trajectory in one Julia call, driving `ITensorMPS.jl`'s own
`tdvp()` (already exported by the `ITensorMPS` dependency; no separate
`ITensorTDVP.jl` package needed). One real-time step is `exp(-i·dt·H)`,
i.e. `tdvp()`'s complex evolution parameter is `-im*dt`. Two real,
non-obvious `ITensorMPS`/`ITensors` bugs had to be worked around to get
this correct, both confirmed by direct inspection of the intermediate
tensors' index structure rather than guessed at from the stack trace
alone:

- `mpsalgebra.jl`'s `applyoperator()` (`contract(A,psi)` plus a
  second contraction against an identity MPO, "this fixes a bug" — a
  *Site*-index priming issue) leaves the *Link* indices of its result
  carrying a stale nonzero prime level. This is invisible to KPM (which
  only ever feeds `applyoperator()`'s output back into more
  `applyoperator`/`summps`/`inner` calls), but `tdvp()`'s internal
  environment bookkeeping (`ProjMPO`) cannot handle it: it aborts deep
  inside `KrylovKit.expintegrator` with "... but the indices are not
  permutations of each other" on the affected Link index. `tdvp.jl`'s
  `apply_clean()` uses `ITensorMPS`'s higher-level `apply()` instead
  (`contract()` is the lower-level primitive both `applyoperator()` and
  `apply()` are built on, but only `apply()` avoids the stale prime).
- Whenever `tdvp()`'s input MPS did not come straight from `dmrg()`
  (e.g. `evolution_ABA()` first applies an operator via
  `mpsjulialive/mpo.py`'s `MPO.__mul__`, which still goes through
  `applyoperator()` for other callers' sake), the same stale Link prime
  reappears — and, confirmed directly, `orthogonalize!()` alone does
  *not* clear it. `evolve_and_measure_tdvp()` therefore explicitly
  `noprime(copy(wf),"Link")`s (never mutating the caller's own MPS) before
  `orthogonalize!()`-ing and starting the TDVP loop.

Both paths were validated directly against exact diagonalization
(quench from a Néel-favoring field to a Heisenberg chain, mirroring
`examples/tdvp_VS_ED_time_evolution`): `evolve_and_measure` agreed with
ED to ~8e-7 on a 4-site chain, and `evolution_DC`'s quench correlator
agreed to ~9e-11 on a 6-site chain.

Excited states (`get_excited_states`/`get_excited`, `n>1`, Hermitian) and
the single-site reduced density matrix (`get_rdm`) are implemented the
same way: `mpsjulialive/excited.jl`'s `excited_states_dmrg` runs the
whole *n*-state loop in one Julia call (one new random-start warm state
per additional excited state, deflated against every wavefunction found
so far via `ITensorMPS.jl`'s own orthogonality-penalty
`dmrg(H, wfs, psi0, sweeps; weight)`, mirroring
`mpscpp3/chain_session.h`'s `Chain::excited_states` and
`pyitensor/chain.py`'s `excited_states`); `mpsjulialive/densitymatrix.jl`'s
`reduced_dm` is a direct Julia port of `pyitensor/chain.py`'s
`reduced_dm` (itself a port of `Chain::reduced_dm`), including its
"divide by the norm squared, not its square root" quirk, preserved for
cross-backend consistency (in practice a no-op, since the ground state
is essentially always already unit-norm). Validated directly: excited
energies match the golden regression values in
`tests/test_excited_states.py` to ~3e-15 on a 4-site chain, and
`get_rdm` matches the existing backend-agnostic `reduced_dm_projective`
to ~1e-15. `n==1` and non-Hermitian excited states already worked before
this — they route through `gs_energy()`/the generic Arnoldi method
(`mpsalgebra.mpsarnoldi`), neither of which is backend-specific.

Bond entanglement entropy (`MPS.get_bond_entropy`, used e.g. by
`entropy.py`'s `central_charge`) was missing entirely on the Julia
backend (`mpsjulialive/mps.py`'s `MPS` class had no `get_bond_entropy`/
`get_site_entropy` methods at all — a plain `AttributeError`, not a
crash inside a dispatch branch). `mpsjulialive/entropy.jl`'s
`bond_entropy` computes it the standard MPS way, via SVD of the two-site
tensor at that bond (mirrors `pyitensor/chain.py`'s `bond_entropy`,
itself a port of `Chain::bond_entropy`) — much cheaper than the generic
`get_correlation_entropy`/`"explicit"`-dmmode fallback the Julia backend
already had (still the only option for multi-site correlation entropy,
just not for a single bond). Validated directly: a 2-site Heisenberg
singlet gives exactly ln 2 (the analytically known Bell-pair
entanglement), and on a 6-site chain the bond entropy at the first bond
exactly matches the existing generic `get_site_entropy` (which computes
the same quantity a completely different way, through
`reduced_dm_projective`).

The CVM dynamical-correlator submode (`get_dynamical_correlator(...,
submode="CVM")`, `cvm.py`) needed almost no Julia-specific code at all:
its correction-vector CG solve (`cvm.py::cvm_correction_vector`) is
already implemented purely with backend-agnostic MPS/MPO algebra
(`self.toMPO()`, MPS `+`/`-`/scalar-`*`, `.dot()`), which already works
against `julia_live`'s `mpsjulialive.mpo.MPO`/`mpsjulialive.mps.MPS`
classes without any changes. The only thing standing in the way was that
`cvm.py` unconditionally called `self._session.set_sweep_params(...)` to
point every MPO application at `cvm_maxm` rather than the DMRG `maxm`
(no `_session` object exists for `julia_live` to call this on);
`cvm.py::_set_cvm_sweep_params` now temporarily overrides `self.maxm`
instead for that backend, since that's what
`mpsjulialive/mpo.py`/`mps.py` actually read. Validated against ED with
the same tolerance `examples/dynamical_correlator_VS_ED/main.py` uses
for the C++/pure-Python backends (`1e-6`) — matched to `~2e-15` here, in
3 CG iterations per frequency point.

Wiring this up surfaced a real bug in
`mpsjulialive/dynamics.py::get_dynamical_correlator`: its first version
declared `name=`/`delta=`/`es=` as its own named parameters (with
KPM-flavored defaults) purely so it could compute the KPM moments
directly in the same function, which silently captured a caller's
`es=`/`name=`/`delta=` kwargs out of `**kwargs` before they could reach
`cvm.dynamical_correlator` — confirmed directly, `submode="CVM"` ran on
a completely different, wrong frequency grid as a result (CVM's own
`np.linspace(0.,10.0,100)` default instead of the caller's requested
window). Fixed by splitting the KPM implementation out into its own
`_kpm_dynamical_correlator` and making the public
`get_dynamical_correlator(self,submode="KPM",**kwargs)` take no named
parameters of its own — mirroring the top-level
`dynamics.py::get_dynamical_correlator`'s own signature, which was never
vulnerable to this because it never declared submode-specific parameters
either.

The TDZ dynamical-correlator submode (complex-time evolution +
perturbative real-axis reconstruction, `tdz.py`, arXiv:2311.10909) needed
only one new backend primitive, same as CVM: the whole algorithm is
otherwise generic MPS/MPO algebra (`self.toMPO()`, `.dot()`) already
working on `julia_live`. That one primitive is a single complex-time-step
propagator (`tdz.py::_advance_complex_time_step`, formerly gated to
`itensor_version in (3,"python")`) — `mpsjulialive/tdvp.jl`'s `tdvp_step`
already generalizes to it with **no changes**, since `-im*dz` is exactly
the right formula whether `dz` is real (the `evolve_and_measure_tdvp`/
`quench_tdvp` case) or genuinely complex (TDZ's per-step contour
increment); only a thin Python wrapper
(`mpsjulialive/timedependent.py::advance_complex_time_step`) was needed.

Wiring TDZ up caught a second occurrence of the same stale-Link-prime bug
documented above for TDVP: `tdz.py`'s first evolved state,
`self.toMPO(A)*wf_g`, goes through the same `applyoperator()` path
`evolution_ABA()` does, and `tdvp_step` had only been sanitized against
that inside `evolve_and_measure_tdvp`'s own wrapper loop, not inside
`tdvp_step` itself — so any *other* direct caller (TDZ's driver loop in
`_complex_time_correlator`) hit the identical failure. Fixed by moving
the `noprime(copy(psi),"Link")` + `orthogonalize!()` sanitization into
`tdvp_step` itself, so it runs on every step regardless of caller,
instead of relying on each call site to remember to pre-clean its input
— a small, one-time cost next to the sweep itself. Validated against the
same golden test `tests/test_dynamical_correlator.py` uses for the other
backends: peak within `0.03` of the exact 4-site Heisenberg gap
(`0.658919`) — landed at `0.68`.

The last two dynamical-correlator submodes needed no new Julia code at
all, only dispatch routing: `dcex.py` (submode `"EX"`, correlator via
exact diagonalization in the excited-state subspace) only calls
`self.get_excited_states()` (already backend-agnostic, see above) and
generic `MultiOperator`/MPS algebra (`.dot()`, `A*wf`) before dropping
into plain NumPy/SciPy (`eigh`); `distribution.py`'s maxent path
(submode `"maxent"`) is built the same way on top of
`vev.py::power_vev`. Validated `"EX"` the same way
`tests/test_dynamical_correlator.py` already does for the other
backends — cross-checked directly against `"KPM"` and `"CVM"` on the
same chain/operator/frequency grid, landing on the exact same peak
(diff `0.0`) rather than just the analytic gap. `"maxent"` is wired up
for parity but not actually functional on *any* backend right now —
`distribution.py::get_distribution_maxent` imports
`dmrgpy.maxenttk.pymaxent`, which doesn't exist in this checkout
(confirmed directly: `ModuleNotFoundError`), so it hits the same
`except: print("Not functional yet"); exit()` regardless of
`itensor_version` — a pre-existing, cross-backend gap, not something
introduced or left broken by this work.

The fermionic 4-point correlator tensor (`MPS.get_four_correlation_tensor`,
`<Cdag_i C_j Cdag_k C_l>`) was missing the same way bond entropy was — no
`get_four_correlation_tensor` method at all on `mpsjulialive/mps.py`'s
`MPS` class — and needed the same fix: just add the method, delegating
to `entropytk/correlationentropy.py`'s default `ctmode="explicit"` (a
Python loop of `MultiOperator` products + `.aMb()`/`.dot()`, already
generic; `ctmode="full"`, a native per-element AutoMPO build, mirrors
`pyitensor/chain.py`'s own `four_correlation_tensor` but isn't ported to
Julia). Validated directly against ED on the same well-gapped free-fermion
model `tests/test_four_point_correlator.py` uses (agreement to `~2e-15`).

`get_four_correlation_tensor_explicit`'s Python loop
(`entropytk/correlationentropy.py`) didn't exploit the tensor's own
Hermitian symmetry
(`<Cdag_i C_j Cdag_k C_l>^dagger = Cdag_l C_k Cdag_j C_i`, so
`ct[i,j,k,l]` and `ct[l,k,j,i]` are complex conjugates) the way
`ctmode="full"`'s C++/pyitensor implementations already did via their
own `accelerate` flag — added the same `accelerate=True` default there
too (skip one representative of each `(current,conjugate)` pair,
fill in its mirror from the other), an exact, not approximate,
speedup. This mattered for `fermionchain.Spinful_Fermionic_Chain_Native`
(added alongside the interleaved `Spinful_Fermionic_Chain`, see the
model table above): its self.C/self.Cdag are flat, single-flavor-per-
entry lists (`Cup`/`Cdn` interleaved as mode `2*i`/`2*i+1`, matching
`Spinful_Fermionic_Chain`'s own indexing exactly, so the two classes'
tensors compare index for index) added purely so this already-generic
`ctmode="explicit"` path works unchanged for it.

`ctmode="full"` (the C++-accelerated path) was then added for this
class too, via a dedicated `Chain::four_correlation_tensor_spinful()`
(`mpscpp3/chain_session.h`/`bindings.cc`) — the existing
`Chain::four_correlation_tensor()` couldn't be reused as-is, since it
hardcodes the literal `"Cdag"`/`"C"` operator names, undefined on
ITensor's `ElectronSite` (only `Cup`/`Cdn`/`Cdagup`/`Cdagdn` are). The
new method instead hands ITensor's own `AutoMPO` the flavor-resolved
names directly and relies on its built-in automatic fermionic-sign
insertion (`autompo.cc`'s `isFermionic()`/`fermionicTerm()`, triggered
by any operator name starting with `'C'`) to do the Jordan-Wigner
threading — a deliberate, one-off exception to this codebase's usual
rule of threading Jordan-Wigner strings explicitly at the Python level
for backend-agnosticism (see `multioperatortk/jordanwigner_spinful.py`);
safe here only because this calculation always builds and discards its
own fresh, self-contained `AutoMPO`, not a change to the general
Hamiltonian/MPO pipeline. `entropytk/correlationentropy.py::
get_four_correlation_tensor_cpp` dispatches to whichever C++ method
matches `type(wf.MBO)`.

Measured (n=3,4,5,6,12 orbitals): `Spinful_Fermionic_Chain_Native`'s
`ctmode="full"` is the fastest of all four combinations (native/
interleaved × explicit/full) tried at every size, including n=12 (24
flat modes: ~620s vs ~890s for `Spinful_Fermionic_Chain`'s own
`ctmode="full"`, a ~30% win — see
`examples/four_correlation_tensor_spinful_native`). Before this
C++-accelerated path existed, native's only option
(`ctmode="explicit"`) still beat interleaved's `ctmode="full"` at
n=3..6, by a margin that grew with n there but did not continue to
n=12 (the two `ctmode="explicit"` numbers alone came back to
essentially tied, ~700s each) — adding `ctmode="full"` for this class
is what restores and extends the win at n=12. This remains the one
calculation checked so far where the native-site class wins at all,
since it is a static overlap computation rather than an iterative
two-site search, so it never pays the two-site combined-local-
dimension penalty documented in `Spinful_Fermionic_Chain_Native`'s own
class docstring.

Investigated whether `get_correlation_matrix`'s default `dmmode="fast"`
(`entropytk/correlationentropy.py::correlation_matrix_fast`, `n` MPO
applications total rather than `explicit`'s `n(n+1)/2` two-operator
products) could also work on `julia_live` — `mpsjulialive/mps.py`'s
`get_correlation_entropy`/`get_correlation_entropy_density` already
force `dmmode="explicit"`, and this turned out to be necessary, not an
arbitrary choice: `dmmode="fast"` builds an MPO for a single bare
fermionic operator (`Cdag[i]`, odd particle-number parity) before
combining it with another, and `ITensorMPS.jl`'s OpSum-to-MPO compiler
rejects that outright ("Parity-odd fermionic terms not yet supported by
OpSum to MPO conversion", confirmed directly) — `explicit` avoids the
problem by always building the two-operator *product* (`Cdag_i C_j`,
even parity) before ever converting to an MPO. Not pursued further: it's
a real `ITensorMPS.jl` limitation, not a bug on this side, and
`explicit` is already fast in absolute terms (0.29s for a 10-site
chain's full correlation matrix) — a from-scratch fermionic-JW-string
MPO builder to work around it would be a disproportionate amount of new
code for a performance-only win on an already-cheap operation.

A final sweep of every remaining `self._session.` call site reachable
from `julia_live` (`grep -rl "self\._session\." *.py`, checked against
what each call site's own top-level dispatcher actually routes for that
`itensor_version`) turned up two more real gaps, deliberately left
unfixed rather than pursued further, since both fall outside the
"already-generic code, just needs a dispatch branch or a small native
Julia primitive" pattern every fix above followed:

- **`vev.py`'s `npow=` kwarg** (`sc.vev(op, npow=2)`, computes
  `<X^2>`/`<X^n>` via repeated MPO application instead of building the
  `O(n^2)`-term squared operator directly) — `mpsjulialive/vev.py`'s
  `vev(MBO,MO)` takes no `**kwargs` at all, so this raises a plain
  `TypeError` for `julia_live` rather than silently misbehaving. The one
  example exercising it, `examples/power_vev/main.py`, sets
  `sc.itensor_version = "julia"` directly (bypassing
  `setup_julia()`/`_reset_dmrg_state()` entirely) — the legacy,
  already-inert subprocess-based backend documented below, not
  `julia_live` — so it wouldn't validate this fix even if made.
- **`get_distribution`/`kpmdmrg.py::general_kpm`** (spectral distribution
  of an arbitrary operator, not just the Hamiltonian; has example
  coverage in four `examples/magnetization_distribution*` scripts and
  `examples/dynamical_correlator/dynamical_correlator_shift`, but no
  `pytest` coverage) — unlike everything above, this is not a
  wire-up-already-generic-code fix: `kpmdmrg.py::general_kpm_moments`
  depends on `scale_operator()`, which calls `self.lowest_eigenvalue()`/
  `self.bandwidth()` (`manybodychain.py`), both of which route through
  `self.clone()` → `deepcopy(self)` — and `julia_live`'s live Julia
  session handles (`self.jlsites`, every `MPS.jlmps`/`MPO.jlmpo`) are not
  established to be deepcopy-safe (this is exactly why
  `dynamics.py::_max_energy_bound`/`excited.py`'s energy-bound helpers
  were written as their own mutate-and-restore functions instead of
  calling `bandwidth()`/`lowest_eigenvalue()` directly, see their
  docstrings). Doing this properly means generalizing that same
  mutate-and-restore pattern from "the Hamiltonian, negated" to "an
  arbitrary caller-supplied operator" — real, but new design work rather
  than a small addition, so left as a documented follow-up rather than
  implemented here.

**Fixes from an independent multi-angle code review of the whole Julia
backend** (8 finder passes + per-finding verification, run once the above
was otherwise feature-complete): six real, confirmed issues, all fixed
directly rather than left as follow-ups —

- `cvm.py`'s `_set_cvm_sweep_params` (the `julia_live` equivalent of
  `self._session.set_sweep_params`) is now a context manager
  (`_cvm_sweep_params`) instead of a plain set-and-return-the-old-value
  function: `self.maxm` is restored via `finally`, so a CG blow-up mid-loop
  (or a standalone call to the public `cvm_correction_vector`, which
  previously never restored it at all) can no longer leave `self.maxm`
  permanently stuck at `cvm_maxm`.
- `dynamics.py`'s `julia_live` branch now checks `self.is_hermitian(...)`
  before dispatching, the same way the `(2,3,"python")` branch already
  does — previously a non-Hermitian Hamiltonian silently ran the
  Hermitian-only KPM/CVM/TDZ math and returned numerically wrong output.
  There's no non-Hermitian route to fall back to for `julia_live` yet
  (that needs `applyinverse()`, C++-only today), so this raises a clear
  `NotImplementedError` instead.
- `mpsjulialive/excited.py`'s `get_excited_states_dmrg` now honors
  `self.excited_gram_schmidt` (previously silently dropped — `excited.jl`
  has no Gram-Schmidt step of its own), applying the same generic,
  backend-agnostic `algebra.arnolditk.gram_smith` the top-level
  `get_excited_states`'s own `purify=True` path already uses, and
  recomputing energies against the orthogonalized states to match.
- `mpsjulialive/dynamics.py`'s `_max_energy_bound` now restores
  `self.maxm`/`self.nsweeps`/`self.hamiltonian` in a `finally` block —
  previously a `self.gs_energy()` failure partway through (a real
  possibility for a live juliacall session) would leave the chain
  object with a negated Hamiltonian and clamped bond dimension/sweep
  count for every later call.
- `mpsjulialive/tdvp.jl`'s `evolve_and_measure_tdvp` now explicitly
  renormalizes every step (matching `quench_tdvp`'s own existing pattern,
  and `mpscpp3/chain_session.h`'s `tdvp_step`, which passes
  `"DoNormalize",true`) — `tdvp_step`'s own `normalize=false` stays
  unchanged, since `tdz.py`'s `advance_complex_time_step` also calls
  `tdvp_step` directly and its physically-meaningful norm decay (the
  whole point of TDZ's damping mechanism, arXiv:2311.10909) must not be
  renormalized away.
- `mpsjulialive/tdvp.jl`'s `tdvp_step` no longer forces `mindim=2`: that
  floor was a guessed fix for a bug later found to be the stale-Link-prime
  issue (fixed separately, see above) and was never confirmed necessary:
  neither `mpscpp3`'s own `tdvp_step` nor `pyitensor/tdvp.py`'s `svd`
  calls force a floor above ITensor's own default of 1, and forcing 2
  could pad in a spurious singular value at a bond whose true Schmidt
  rank genuinely is 1 (e.g. a near-product state on a small chain).

One more finding was fixed as a cleanup rather than a correctness bug:
`mpsjulialive/excited.jl`'s `excited_state_dmrg` had copy-pasted the
`Sweeps`-construction + stdout-silencing boilerplate a third time
(`get_gs.jl` already had two near-identical copies); both are now built
on shared `make_sweeps`/`run_quiet` helpers in `get_gs.jl`.

A seventh candidate (`densitymatrix.jl`/`entropy.jl` lacking a bounds
check when `site`/`b` is the chain's last site) was investigated and
found to be a pre-existing, deliberately-preserved limitation shared
identically by `pyitensor/chain.py` and `mpscpp3/chain_session.h` (see
`reduced_dm`'s own docstring above) — not a new risk, so left as-is.

(A separate, older subprocess-based Julia path, `itensor_version="julia"`
via `juliarun.py`, is not reachable through the normal public API and
should be treated as legacy/inert.)

**A second review round** (run against the KPM optimization + the six
fixes above together) found two of those fixes were themselves
incorrect, plus several smaller cleanup items — all fixed, and this time
backed by persisted regression coverage (`tests/test_julia_live.py`,
`julia_live`'s first, since the first round's validation was ad hoc and
non-committed):

- `evolve_and_measure_tdvp`'s new renormalization (previous round, above)
  forced the state to *unit* norm every step instead of preserving its
  own starting norm — diverging from `quench_tdvp`'s existing
  `norm0`-target pattern in the same file, and silently rescaling every
  result by `1/‖wf‖²` whenever the input wasn't already unit-norm (e.g.
  `evolution_ABA()`'s `wfA = A*wf`, for any non-unitary `A`). Now
  captures `norm0 = ‖wf‖` once at entry and renormalizes to that every
  step, matching `quench_tdvp` and the pure-Python reference's physical
  behavior.
- The non-Hermitian guard added to `dynamics.py`'s `julia_live` branch
  (previous round, above) ran *before* submode dispatch, so it
  unconditionally blocked every submode — including `submode="EX"`,
  which already has a working, itensor_version-agnostic non-Hermitian
  path (`dcex.py` → `excited.py`'s `excited_states_non_hermitian`, not
  gated on `itensor_version` at all) and `submode="maxent"`. The guard
  now only fires for `submode in ("KPM","CVM","TDZ")`, matching what its
  own error message already claimed.
- Both renormalization sites now use `norm(psi)` instead of
  `sqrt(abs(inner(psi,psi)))`: `tdvp_step()`'s own `orthogonalize!(psi,1)`
  guarantees a well-defined orthogonality center on its output, so
  `ITensorMPS`'s `norm()` short-circuits to a cheap `O(D²)` single-tensor
  norm instead of a full `O(N·D³)` chain contraction — confirmed via a
  live benchmark (`N=40`, `D=60`: ~4ms vs ~1.18s per call).
- `mpsjulialive/dynamics.py`'s `_same_mps` had an unused `jlsites`
  parameter left over from the same cleanup that removed it from its
  sibling `_kpm_moments_full`/`_kpm_moments_accelerated` — dropped, and
  `same_mps`/`summps` on the Julia side gained the `cutoff` argument
  described above.

The commit that introduced this round's fixes claimed "new targeted
tests" without actually adding any (`tests/`/`examples/` were untouched)
— `tests/test_julia_live.py` now covers exactly the four behaviors that
commit described validating: CVM's `maxm` restoration, the non-Hermitian
dispatch split, `evolution_ABA`'s norm preservation, and the KPM peak
position after the `apply_op`/`summps` changes. The whole file is skipped
if `juliacall`/Julia isn't available, the same way `tests/` already skips
itensor_version 2/3 when the corresponding compiled extension isn't
present.

### 4.7 Supporting `*tk` packages

Functionality is generally split into a top-level module (the public
method) and a `<topic>tk/` package with the implementation, e.g.
`correlator.py` / `fermionchaintk/staticcorrelator.py`, `dynamics.py` /
`dynamicstk/`, `entropy.py` / `entropytk/`, `mpsalgebra.py` /
`mpsalgebratk/`. When changing behavior for a specific feature, check
both the thin dispatch module and its `tk` counterpart.

### 4.8 KPM / dynamical correlators

Dynamical correlators and generic operator distributions are computed
with the Kernel Polynomial Method (`kpmdmrg.py`;
`Chain::kpm_dynamical_correlator` / `Chain::general_kpm` in
`chain_session.h` on the C++ side; `algebra/kpm.py`-style
moment recursion on the ED side) rather than by direct spectral
decomposition, since exact diagonalization of the full spectrum is
infeasible for large chains.

### 4.8b Non-Hermitian KPM (NH-KPM)

For a non-Hermitian Hamiltonian, `dynamics.py`/`edtk/dynamics.py`'s
`submode="KPM"` gate routes to `nonhermitian/kpm.py` instead of
`kpmdmrg.py`/`algebra/kpm.py`'s Hermitian moment machinery (previously
every submode fell through unconditionally to the correction-vector
fallback `nonhermitian/dynamics.dynamical_correlator_non_hermitian`,
which assumes `A^\dagger=B` and never used the left ground state at all).
NH-KPM is a port of the biorthogonal Chebyshev-moment algorithm in
[NHKPM.jl](https://github.com/GUANGZECHEN/NHKPM.jl) (Phys. Rev. Lett.
130, 100401): rather than a single self-dual ground state, it uses the
biorthogonal pair `(psi_L,psi_R)` already produced by NH-DMRG (§4.4,
`nhdmrg()`/`gs_energy_nhdmrg`) on the DMRG side, or
`algebra.biorthogonal_ground_state` (a small addition to `algebra.py`
diagonalizing both `h` and `dagger(h)` and matching/normalizing the pair
so `<psi_L|psi_R>=1`) on the ED side.

The algorithmic difference from Hermitian KPM: a non-Hermitian rescaled
operator `hs=(z*Id-H)/E_max` has no real spectrum, so the ordinary
single-operator Chebyshev recursion (valid once `H` is rescaled onto
`[-1,1]`) doesn't apply. NH-KPM instead builds a *coupled*
forward/adjoint recursion from both `hs` and `hs_dag` (`get_mu_n_nh`/
`spec_from_moments_nh` in `algebra/kpm.py`, `Chain::nhkpm_moments` in
`mpscpp3/chain_session.h`, ported line-for-line from the reference's
`get_vn_NH`/`get_mu_n_NH`/`get_spec_kpm_NH`). The practical consequence:
unlike Hermitian KPM (moments computed once, reused at every frequency
via cheap Chebyshev polynomial evaluation), NH-KPM's expansion operator
depends on `z` itself, so `nonhermitian/kpm.py` rebuilds the MPO/matrix
and recomputes the full moment recursion at *every* requested frequency
point — mirroring the reference algorithm's own cost profile, not a
missed optimization. `E_max` must be supplied by the caller: there is no
automatic spectral-radius estimator yet for a non-Hermitian operator
(`chain_session.h`'s `maximum_energy()`, used for the Hermitian case,
runs ordinary variational DMRG on `-H` and is meaningless once `H` isn't
Hermitian).

On the DMRG side, the `(z*Id-H)/E_max` operator and its dagger are built
once per frequency as ordinary `MultiOperator`s in Python (`z*identity()
- self.hamiltonian`, then `.get_dagger()`), converted via the existing
`to_terms()`/`build_mpo` machinery — no new MPO arithmetic was needed in
C++, only the new coupled-recursion primitive `Chain::nhkpm_moments`
(public, alongside `general_kpm`) and its `bindings.cc` wrapper. One
non-obvious bookkeeping fix was needed to reuse NH-DMRG's own state:
`gs_energy_nhdmrg` renormalizes `wf0=psir.normalize()` to unit norm but
leaves `nh_left_wf` at NH-DMRG's own `<psi_L|psi_R>=1` scale, so the two
attributes no longer satisfy that biorthogonal relation together;
`nonhermitian.kpm.dynamical_correlator_nhkpm` restores it explicitly
(`psil = psil/conj(psil.dot(psir))`) before using the pair, rather than
assuming the convention still holds after `gs_energy()`.

Implemented so far for the ED backend, `itensor_version=3`, and
`itensor_version="python"` (`itensor_version=2` raises
`NotImplementedError` from the DMRG-side driver, matching the same
scope limitation `nhdmrg()` itself doesn't have but this feature's C++
side does — no `Chain::nhkpm_moments` was ported to `mpscpp2`). The
pure-Python backend's own `Chain.nhkpm_moments`
(`pyitensor/chain.py`) is a line-for-line transcription of
`mpscpp3/chain_session.h`'s version, built on the same MPO/MPS algebra
(`applyMPO`/`sum`/`inner`) `general_kpm` already used there; no new
pyitensor primitives were needed. All three backends agree to machine
precision (`examples/non_hermitian/nhkpm_v3_VS_ED` for ED vs v3,
`examples/non_hermitian/nhkpm_python_VS_v3_timing` for v3 vs
`"python"` — the latter on a non-uniform-hopping, non-uniform
imaginary-onsite-energy interacting fermionic chain, chosen to avoid
the Re-degenerate-ground-state pitfall a uniform staggered pattern can
trigger; `"python"` measured ~1.8-2x slower than v3 for this workload,
consistent with NH-KPM's per-frequency moment recursion being far more
matvec-heavy than the Hermitian KPM path, which amortizes its moments
over the whole spectrum instead of recomputing per frequency).
Validated directly against the live upstream Julia implementation, not
just by hand: the reference's own shipped `.OUT` output files did not
reproduce under the parameters recorded in the currently-checked-in
example scripts (plausibly stale, given the upstream repo ships
multiple undated `_backup`/`_old` source variants of the same
algorithm), so instead the exact `get_vn_NH`/`get_mu_n_NH`/
`get_spec_kpm_NH`/`dos_kpm_NH` function bodies were run directly (Julia
was available; only the `mode="matrix"` branches were needed, which
have no `ITensors.jl`/`Arpack.jl` dependency) on the same single-particle
Hatano-Nelson-type chain used in the reference's own `Hatano.jl`
example (asymmetric-hopping non-Hermitian SSH chain, N=20, all-real
spectrum). Comparing the resulting tDOS(E) (E scanned along the real
axis at fixed small imaginary broadening) against this dmrgpy port on
the identical matrix and identical E_max/N/kernel/energy grid: all 22
detected peaks land at bit-for-bit identical energies and heights
(max relative difference ~3e-15, i.e. floating-point round-off, not an
algorithmic discrepancy), and every peak sits at the real part of a
genuine eigenvalue of the underlying matrix, to within the energy
grid's resolution. (Earlier, weaker checks are kept for context: the
coupled recursion also reproduces an independently-derived scalar
recursion exactly for a 1x1 test case, and the reconstructed density
correctly sharpens with more moments.)

### 4.8c STM/Kondo tunneling spectra (`kondospectrumtk/`)

`Spin_Chain.get_kondo_spectrum` (physics documented in the user guide,
§17) is a different kind of ED calculation from the dynamical-correlator
machinery above: it needs the *full* spectrum (every eigenstate, as a
possible virtual intermediate state in the third-order perturbation
theory), not just the ground state or a resolvent at a target frequency.
It reuses the ED backend's existing `EDchain.get_diagonalized_hamiltonian`
(dense `eigh`, already used elsewhere for small systems) rather than
adding a new diagonalization path, and is deliberately independent of a
chain's own `itensor_version`/mode setting -- DMRG cannot supply the
full spectrum this method needs regardless of which backend a chain
would otherwise use, so `KondoSpectrum` (`kondospectrumtk/edkondo.py`)
always goes through `Spin_Chain.get_ED_obj()` directly.

Two of the source paper's own closed-form equations for supporting
numerical functions (the temperature-broadened step and a
temperature-broadened logarithmic Kondo function) do not reproduce the
behavior the paper itself describes and plots for them; `kondospectrumtk/
stepfunctions.py` re-derives both from the paper's own unambiguous
defining integrals instead, each independently verified against
digitized values from the paper's own figures. See that module's
docstring and `docs/user_guide.md`'s §17 for the full derivation notes
and the feature's other scope limitations (single tip-coupled site,
third-order terms are single-direction only, the potential-scattering
interference term's general-spin form is an extrapolation from the
paper's own worked example).

**`mode="DMRG"`** (`T=0` only) is a second, independent route through
this same feature that *does* stay within `itensor_version=3`, never
diagonalizing beyond the ground state -- viable because at `T=0` both
terms simplify enough to avoid needing the full spectrum `mode="ED"`
requires:

- The second-order term (`secondorder_dc.py`) becomes a `Theta0`-weighted
  cumulative integral of the ordinary `T=0` dynamical structure factor,
  so it reuses `get_dynamical_correlator` (`submode="KPM"`/`"CVM"`)
  as-is.
- The third-order Kondo term (`twotime.py`, `dmrgtwotime.py`) is a
  genuinely new construction: a Heisenberg three-point function
  `G(t2,tau)=<GS|Sl(t2+tau)Sk(t2)Sj(0)|GS>` built via real TDVP time
  evolution (a "checkpoint-and-branch" pattern: evolve `Sj|GS>` forward
  and backward through `t2`, apply `Sk` at each checkpoint, evolve each
  branch further through `tau`, overlap with a fixed `Sl|GS>` reference
  at every step -- reusing the single-step `tdvp_step` primitive
  `tdz.py` already drives manually the same way), then converted to a
  frequency-domain quantity via two closed-form time-domain kernels
  (`Theta0`'s Cauchy-principal-value kernel, computed via an FFT-based
  Hilbert transform; `F0`'s kernel, smooth and closed-form) rather than
  by evaluating those functions pointwise on a discrete frequency grid --
  an initial FFT-grid attempt at the same physics did not converge
  robustly (a discontinuous step function and a log singularity both
  have slowly-decaying spectral content that a finite discrete grid
  resolves poorly), which is why this construction exists in the first
  place rather than a simpler direct approach. `twotime.py`'s functions
  are backend-agnostic (they only ever consume a discretely-sampled
  `G(t2,tau)` array or batches thereof) -- `edtwotimeref.py` supplies
  the same interface via exact ED eigenbasis time evolution, used both
  as the development-time reference this construction was validated
  against and as a fast, exact ground truth for `dmrgtwotime.py`'s own
  test.

`dmrgtwotime.py` was written against this codebase's existing, verified
DMRG API and validated once a compiled ITensor v3 backend became
available: `G(t2,tau)` matches the ED reference to ~1e-9-1e-10
pointwise, and the swept third-order Kondo term matches a
grid-consistent ED reference to ~1e-10
(`tests/test_kondo_spectrum_dmrgtwotime.py`, skipped automatically when
no compiled `itensor_version=3` backend is available). Getting there
surfaced three real bugs, none of which showed up in the ED-only testing
this module was originally written against:

- `tdvp_step` silently renormalizes its output to unit norm on *every*
  call. `Sj|GS>` (the state each two-time trajectory starts from) is
  generally not unit-norm (`Sj` isn't norm-preserving), so letting that
  renormalization stand discarded the true amplitude at every step --
  confirmed directly, it produced overlaps exactly a factor of 2 too
  large end to end for a state of norm 0.5. Fixed by evolving a
  normalized copy and rescaling every returned checkpoint by the true
  original norm.
- The forward/backward time-stepping loop shared one running `times`
  list, so the "backward" branch kept counting down from the forward
  branch's endpoint instead of from `t=0` -- it never actually reached
  negative times at all. Fixed by walking forward and backward from
  `t=0` independently and merging afterward.
- `twotime.py`'s `t2`-integral used `np.trapz` per chunk, which is
  exactly 0 for a single-point chunk (no width to integrate over) -- and
  DMRG chunks are always single-point (each `t2` checkpoint is its own
  real trajectory), so this silently zeroed the entire third-order Kondo
  term end to end. Fixed with a uniform Riemann-sum accumulation
  (trapezoidal endpoint correction only at the true global first/last
  `t2` points) that is well defined for any chunk size.

A 1-site test chain also hit an unrelated internal ITensor v3 error
(building the Hamiltonian MPO, distinct from the two-site-`dmrg()`
short-chain crash `mode.py`'s ED fallback already guards against) --
chains need at least 3 sites for this feature. The second-order
term (`submode="KPM"`) was spot-checked too, agreeing to within a few
tens of percent at thresholds, consistent with the expected
delta-broadening/moment-truncation error on top of what the ED path
already has.

### 4.9 TDZ / complex-time-evolution dynamical correlator

`tdz.py` implements `submode="TDZ"` (Cao, Lu, Stoudenmire & Parcollet,
"Dynamical correlation functions from complex time evolution",
arXiv:2311.10909): instead of evolving along the real-time axis
(`submode="TD"`, `timedependent.py`), it evolves along a complex-time
contour `z(t,alpha0)`, whose negative imaginary part damps high-energy
content as the simulation proceeds, keeping the MPS bond dimension
needed for a given accuracy much smaller than under real-time evolution
alone. The true real-time correlator is then recovered order by order
via a perturbative Taylor expansion in `alpha0` around the simulated
contour (`_reconstruct_real_axis`, hardcoding the paper's own explicit
n<=4 Appendix-B formulas). Each order needs only a fixed set of
precomputed overlap targets (`{H^n(B|GS>)}`, built once via repeated
`toMPO`/`StaticOperator` application, independent of t) plus pure scalar
contour integrals -- no new tensor-network machinery beyond a single new
per-step propagator primitive:

- **`itensor_version` 3 or `"python"`, `tevol_method="TDVP"`** (the
  paper's own setup): a single two-site-TDVP step with a *complex* time
  argument -- `pyitensor/tdvp.py`'s Krylov-exponentiation core
  (`_lanczos_expm_multiply`) was already fully generic to complex
  coefficients with no changes needed; `mpscpp3/chain_session.h`'s
  `Chain::tdvp_step` (moved from a private helper to a public method and
  widened from `double dt` to `Cplx dt`) simply forwards to the vendored
  `TDVP/tdvp.h`, which documents its own time argument as natively "real,
  imaginary, or complex".
- **`itensor_version=2`** (mpscpp2 has no TDVP at all): a single
  MPO-Taylor step, `Chain::evolve_taylor_step` in both
  `mpscpp2/chain_session.h` and `mpscpp3/chain_session.h` (the latter as
  a cross-check / non-TDVP alternative), built on the existing
  `evoloperator()` Taylor-expansion of `exp(z*H)` (also mpscpp2-only,
  now widened from a real `dt` to a complex `z` throughout, including
  `pyitensor/chain.py`'s own `_evoloperator`) -- the pre-existing
  deliberately-reproduced "z^3/6 multiplies H2 not H3" quirk (see
  CLAUDE.md) is unaffected by this widening.

Current scope: only the "greater" branch of the correlator is computed
(the same simplification `submode="TD"` already makes), fed into the
same windowing/FFT tail as `"TD"` (factored out into
`timedependent._fourier_transform_correlator` so both submodes share it).

## 5. Backend performance: v3 vs the pure-Python backend

The pure-Python backend (`itensor_version="python"`) trades raw speed for
zero build dependencies. Numbers below were measured directly on one
machine, with `numba` installed *and* explicitly opted in
(`pyitensor.kernels.USE_NUMBA = True`) for its accelerated
matvec/effective-Hamiltonian kernel (see `pyitensor/kernels.py`). This
kernel is **not** used automatically just by having `numba` installed --
`USE_NUMBA` defaults to `False`, because the one-time JIT compile cost is
a net regression for this library's typical one-shot-script usage
pattern (see `kernels.py`'s own docstring for measured numbers). Set
`kernels.USE_NUMBA = True` yourself before running DMRG if you want this
path. Absolute times will vary by machine and load, but the qualitative
trends should hold.

### 5.1 Ground state energy (Heisenberg spin-1/2 chain)

| n  | v3 (s) | python (s) | ratio (python / v3) |
|----|--------|------------|----------------------|
| 8  | 0.10   | 0.16       | 1.6x  |
| 12 | 0.38   | 0.51       | 1.3x  |
| 16 | 0.66   | 1.12       | 1.7x  |
| 20 | 1.09   | 2.44       | 2.2x  |
| 24 | 1.58   | 4.06       | 2.6x  |
| 28 | 1.84   | 6.93       | 3.8x  |
| 30 | 2.57   | 10.84      | 4.2x  |
| 32 | 2.80   | 14.20      | 5.1x  |

The ratio grows with system size: v3 benefits from real ITensor's
QN-block-sparse tensors, which `pyitensor`'s dense NumPy tensors don't
have even with numba JIT acceleration on the hot matvec loop. For quick,
small-system exploratory work the pure-Python backend is quite
competitive (1.3–2x); for production-size chains (n ≳ 30) expect it to be
several times slower.

### 5.2 KPM dynamical correlator (`⟨Sz₀(t)Sz₀(0)⟩`, Heisenberg chain)

| n  | v3 (s) | python (s) | ratio |
|----|--------|------------|-------|
| 6  | 0.19   | 0.35       | 1.8x  |
| 8  | 0.39   | 1.58       | 4.0x  |
| 10 | 1.11   | 4.54       | 4.1x  |
| 12 | 3.63   | 8.17       | 2.2x  |

Dynamical correlators are generally worse for the Python backend than
plain ground-state DMRG (2–4x vs 1.3–2.6x at comparable n): KPM performs
many moment-recursion sweeps at a fixed, comparatively large bond
dimension (`kpmmaxm=50` by default, vs `maxm=30` for ground states), so
the missing block-sparsity cost is paid on every sweep rather than
amortized. Systems around n≈30 were not benchmarked to completion here —
extrapolating from this table, expect a KPM correlator run at that size
on the pure-Python backend to take on the order of minutes or more.

### 5.3 Other calculation types (single-run timings, n=6–8)

| Calculation | v3 | python | ratio |
|---|---|---|---|
| Excited-state gap (TFIM, n=8) | 0.41s | 2.35s | 5.7x |
| TDVP quench evolution (n=6, 30 steps) | 0.11s | 0.97s | 8.9x |
| Hubbard ground state (3 sites) | 0.04s | 0.11s | 3.1x |

These single-run numbers include a one-time numba JIT compilation cost
(~0.1–0.4s) for each distinct kernel the first time it runs in a process
— a script that repeats the same calculation type many times (e.g. a
parameter sweep) amortizes this and approaches the steady-state ratios
shown in §5.1/§5.2, not these first-call numbers.

### 5.4 Practical guidance

- Use `itensor_version="python"` for quick prototyping, CI/environments
  without a C++ toolchain, or small systems (n ≲ 20) where the 1.3–2x
  overhead doesn't matter.
- `numba` is already installed (it's a core dependency, see §2); also
  install `jax` and set `pyitensor.kernels.USE_NUMBA = True` (off by
  default, see §5's note above) if you want the accelerated matvec
  kernel — without that opt-in, `pyitensor` falls back to plain NumPy
  loops and is substantially slower than the numbers above (see
  `pyitensor/__init__.py`'s own docstring).
- For production-size chains, dynamical correlators, or anything
  performance-sensitive, compile the C++ extension (`python install.py`)
  and use `itensor_version=2` or `3`.
- `examples/backend_timing_gs_energy/main.py` reproduces the n=8..20
  rows of the §5.1 table directly against the current codebase and
  current machine (n=24..32 aren't included, to keep the example fast to
  run) — re-run it rather than trusting these numbers verbatim on a
  different setup.

## 6. Directory reference

| Path | Contents |
|---|---|
| `src/dmrgpy/` | the Python library |
| `src/dmrgpy/manybodychain.py` | `Many_Body_Chain`, the shared base class |
| `src/dmrgpy/multioperator.py`, `multioperatortk/` | backend-agnostic operator representation |
| `src/dmrgpy/mode.py`, `cppext.py` | DMRG/ED and backend dispatch |
| `src/dmrgpy/mpscpp2/`, `mpscpp3/` | vendored ITensor C++ (v2, v3) + pybind11 bindings |
| `src/dmrgpy/pyitensor/` | pure-Python ITensor-v3-subset reimplementation |
| `src/dmrgpy/mpsjulialive/`, `mpsjulia/` | Julia/ITensors.jl backend modules |
| `src/dmrgpy/edtk/`, `pyfermion/`, `pyspin/`, `pyboson/`, `pyzn/` | exact-diagonalization backend |
| `src/dmrgpy/kpmdmrg.py` | Kernel Polynomial Method dynamical correlators |
| `src/dmrgpy/nonhermitian/kpm.py` | non-Hermitian KPM dynamical correlator (NHKPM.jl port, ED + itensor_version=3) |
| `src/dmrgpy/tdz.py` | complex-time-evolution dynamical correlator ("TDZ", arXiv:2311.10909) |
| `examples/` | 100+ self-contained example scripts (also the regression suite) |
| `examples/v2_VS_v3_*` | backend-vs-backend correctness comparisons |
| `examples/backend_timing_gs_energy/` | backend-vs-backend timing comparison |
| `installtk/` | build-requirement checking and ITensor/extension compilation |
| `docs/` | this document, and its LaTeX counterpart |

## 7. Further reading

- `README.md` — installation quick-start and further usage examples
  (static correlators, CFT central charge, bond-dimension convergence).
- `CLAUDE.md` — a more detailed, implementation-level architecture
  reference (written for AI coding assistants working in this
  repository, but equally useful to a human maintainer diving into a
  specific subsystem).
- Tutorials linked from `README.md` cover many-body quantum magnetism,
  correlated fermionic systems, and tensor-network methods more broadly.
