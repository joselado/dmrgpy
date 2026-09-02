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

### 2.1 Releasing a new version to PyPI (maintainers)

`release_pypi.py`, run from the repo root, automates cutting a release of
the pure-Python distribution described above:

```bash
python release_pypi.py                     # build + validate + upload to TestPyPI
python release_pypi.py --repository pypi    # the real thing
python release_pypi.py --no-upload          # build + validate only, upload nothing
python release_pypi.py --skip-tests         # skip `pytest tests` (faster iteration)
```

It runs the full `pytest tests` suite, rebuilds `dist/*.whl`/`dist/*.tar.gz`
from a clean `dist/`/`build/`, runs `twine check` on both artifacts, then
installs the built wheel into a throwaway virtualenv (not from `src/`, so
missing package-data or import bugs actually surface) and cross-checks the
pure-Python DMRG backend against ED on a small Heisenberg chain before
asking for interactive confirmation (`--yes` skips the prompt) and
uploading. See CLAUDE.md's "Packaging / PyPI" section for the specific
`pyproject.toml` pitfalls this guards against (the `numba`-driven
`requires-python` floor, the `juliapkg.json` package-data entry, etc.). A
given version can never be re-uploaded to (Test)PyPI once published — bump
`pyproject.toml`'s `[project].version` for the next attempt, then tag it:
`git tag vX.Y.Z && git push origin vX.Y.Z`.

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
physical model or feature, organized into thematic subfolders
(`groundstate/`, `staticcorrelators/`, `dynamical_correlator/`,
`spin_models/`, `fermion_models/`, `boson_models/`, `time_evolution/`,
`excited_states/`, `entanglement/`, `topological/`, `kondo/`,
`non_hermitian/`, `finite_temperature/`, `magnetization/`,
`backend_comparison/`, `utilities/`, `algebra/`, ...), and doubles as the
project's de facto regression suite (there is no separate unit-test
suite for the Python code). See `examples/readme_examples/` for the
snippets shown in `README.md`, and `examples/backend_comparison/` plus
the `v2_VS_v3_*`-named examples living inside each theme folder (e.g.
`examples/spin_models/v2_VS_v3_heisenberg_spin_half`) for scripts that
directly compare backends against each other on correctness, and
`examples/groundstate/backend_timing_gs_energy/` for timing.

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

### 4.1b Infinite chains (`infinitechain.py`, `pyitensor/idmrg.py`)

`Infinite_Many_Body_Chain`/`Infinite_Spin_Chain` (`infinitechain.py`) are
a deliberately **independent** object, not a `Many_Body_Chain` subclass:
`Many_Body_Chain`'s whole design (mode dispatch between DMRG/ED, a fixed
`self.ns`-length site list, KPM/dynamics/entanglement/excited-states
machinery) assumes a *finite* chain throughout, and none of that has any
meaning for a translationally-invariant, infinite unit-cell chain.
`itensor_version="python"` (default) or `3` are supported (constructing
anything else raises `NotImplementedError`) — `2` and `"julia_live"` have
no iDMRG port. The public surface is deliberately narrow: `set_hamiltonian`,
`gs_energy` (energy *per site*, the physically meaningful quantity for an
infinite system), static `vev`/`correlator`, and (via a finite-window
reduction to the existing finite-chain KPM stack, see below)
`kpm_finite` — no excited states, no entanglement/entropy in v1.

**Backend dispatch, unlike the finite-chain `mode.py`, is explicit rather
than automatic**: `gs_energy` branches directly on `self.itensor_version`
rather than falling back silently (there's no ED oracle for an infinite
system to fall back *to*). `itensor_version="python"` runs
`pyitensor/idmrg.py`'s own growing algorithm (see below) and keeps the
resulting `IDMRGResult` in `self._result` for `vev`/`correlator` to reuse.
`itensor_version=3` instead calls the compiled `mpscpp3` backend's
`Chain::idmrg_ground_state` directly — a line-by-line C++ port of the same
algorithm (built from the start with fresh Index objects on every
automaton tensor to sidestep the duplicate-Index bug class described
below, rather than needing the same retrofit). The two implementations had
diverged for a while — the Python side gained McCulloch's wavefunction
prediction, the theta-cell extraction and the energy baseline described
below, and none of the three was ported — but **all three are ported
now** (`idmrg_wavefunction_prediction`, `idmrg_theta_cell` +
`ic_canonicalize_cell`, `idmrg_subtract_energy_baseline`), and with the
unit cell in hand the static observables built on it follow:
`Chain::idmrg_onsite_expectation`, `Chain::idmrg_two_point_correlator` and
`Chain::idmrg_local_excitation_gap`. So `vev`/`correlator` work under
`gs_method="idmrg"` on this backend too; `self._result` is still left
`None` on the v3 path (it is a `"python"`-backend-only cache), because
everything the v3 observables need lives inside the C++ `Chain` as private
snapshots. `test_itensor_version3_matches_python_backend` now compares the
converged *state* (`vev` plus a short `correlator` sweep) rather than only
the energy density, which was never the part that differed, so the
divergence is guarded rather than merely documented; `tests/
test_idmrg_correlator_v3.py` covers the new surface directly. What stays
`"python"`-only: `local_excitation_gap`'s `window>0` variant (an explicit
prototype rather than stable API) and the iMPS algebra over `IDMRGResult`
(`imps_overlap`/`apply_mpo`/`imps_sum`), which `infinitechain.py` does not
expose for `gs_method="idmrg"` anyway and which the VUMPS path already
covers on this backend.
Benchmarked (S=1/2 uniform Heisenberg chain, `maxm=30`, `maxiter=60`,
`etol=1e-9`): v3 is ~2.6-2.7x slower than `"python"` at both `n_uc=1` and
`n_uc=2`, with energy-density agreement between the two backends to
~1e-8. This is the *gapless* (critical) case, the hardest for any
Krylov-based local solve since the correlation length diverges and the
local 2-site Arnoldi solve genuinely needs close to its full
`krylovdim` to converge at every macro-iteration (confirmed by direct
instrumentation) — for a *gapped* model the same local solve typically
converges in far fewer Krylov vectors, and v3 can end up *faster* than
`"python"` (confirmed directly on a dimerized `n_uc=2` chain: v3 0.08s
vs `"python"` 0.23s at `maxm=30`, `maxiter=60`). This asymmetry is
`arnoldi_smallest_real`'s own `early_tol` early-exit at work (see its
own doc comment): every pre-existing caller of that shared Arnoldi
routine (`nhdmrg_one_sweep`) opts out (`early_tol<0`, the default) and
is completely unaffected; only `idmrg_local_solve`'s own call opts in,
checking the same residual-tolerance criterion the between-restart
check already trusted, just checked incrementally (every 3rd
Krylov-vector addition, from `m>=8` on) instead of only once a full
`krylovdim`-size subspace has been built — mathematically the same
"stop once already converged" logic the pre-existing "happy breakdown"
branch already relied on (an exactly-zero residual), just generalized
to a tolerance. `TagSet` string-tag parsing (`"Site"`/`"Link"`/`"a"`)
was also confirmed by profiling to be a real, multi-percent cost on its
own (repeated on every Krylov iteration and every micro-step) and is
now cached in function-local statics instead of re-parsed per call.
`Chain::idmrg_ground_state` is reached by constructing a `Chain`
directly from the unit cell's own `site_types` (`cppext.get_backend(3).
Chain(self.site_types)`), the same low-level pattern `nhdmrg.py` uses for
its own compiled-backend calls (`H.to_terms()` fed straight into
`self._session.<method>(...)`) — terms are shifted from `infinitechain.py`'s
own 0-based `h_intra`/`h_inter` site convention to the codebase's usual
1-based one via the same `MultiOperator.to_terms()` every other compiled
backend call already uses -- but, uniquely among this codebase's compiled
backend calls, with **`jordan_wigner_transform=False`**: the infinite path
hands the C++ side the raw, untransformed operator names.

That is not an optimization but a correctness requirement, and getting it
wrong is what left `itensor_version=3` unable to run a fermionic infinite
chain at all until 2026-08-22. `multioperator.jordan_wigner()` is the
*finite*-chain transform: it makes every fermionic operator at site `s`
carry an explicit `F` on sites `1..s-1`, a string anchored at site 1 of
the chain. An infinite chain has no site 1, and the strings also break the
documented "at most 2 distinct sites per term" contract, since each `F`
factor lands on a site of its own -- so `Chain::idmrg_ground_state`
rejected any term whose endpoints were not adjacent (the adjacent case
happens to survive, because the global string then collapses to exactly
the correct endpoint-`F` composition with nothing in between, which is why
the spin-only test suite never noticed). Both backends now thread the
string themselves, locally and translation-invariantly between each term's
own two endpoints: `pyitensor/idmrg.py`'s `_term_site_matrices` (which
always did) and `mpscpp3/chain_session.h`'s `idmrg_classify_terms`, a port
of it -- fermionic reordering sign, the endpoint `F` factor,
`IdmrgBond::carry_ferm` (consumed by `idmrg_build_row`, which propagates
`F` rather than `Id` along that bond's pending channel), and rejection of
odd total fermion parity. `idmrg_classify_terms` serves
`Chain::vumps_ground_state` too, so one implementation covers both
`gs_method`s. Odd-parity terms are additionally caught in
`infinitechain.py`'s own `_canonicalize_hamiltonian`, on the raw terms
every backend shares, because the C++ rejection is an ITensor `Error()`
that aborts the process rather than raising into Python. Cross-site
correlators of fermionic operators are a separate, still-unthreaded code
path and are refused outright -- see `Infinite_Many_Body_Chain.
correlator`'s own guard.

**A third `gs_energy` dispatch option: VUMPS.** `pyitensor/vumps.py`
implements Variational Uniform Matrix Product States (Zauner-Stauber,
Vanderstraeten, Fishman, Verstraete, Haegeman, arXiv:1701.07035) as a
completely independent ground-state solver, reached via
`self.gs_method="vumps"` (the default since 2026-08-08, on EITHER
`itensor_version` -- `itensor_version=3` calls `Chain::vumps_ground_state`
instead, see below) rather than `pyitensor/idmrg.py`'s growing algorithm
(`self.gs_method="idmrg"`). Unlike that growing algorithm
— which extends a finite window indefinitely and truncates it back down
to `maxm` at every micro-step, so its own converged unit cell is only an
approximation to the true `maxm`-dimensional optimum — VUMPS solves
directly, in the thermodynamic limit, for the mixed-gauge `{AL, AR, C}`
fixed point at a FIXED target bond dimension (`self.maxm` doubles as
VUMPS's own `D`). `vumps.py` reuses `idmrg.py`'s automaton builder
(`_build_automaton`) and low-level transfer-matrix primitives
(`_apply_transfer`/`_apply_transfer_from_left`/
`_dominant_right_fixed_point`/`_dominant_left_fixed_point`/
`_canonicalize_periodic`) and `idmrg_excitations.py`'s own dense-
regularized-environment machinery (`_op_transfer_matrix`/`_cap_left`/
`_cap_right`/`_apply_op_ket`/`_dense_linear_map`/`_pending_channels`/
`_onsite_matrix`) rather than re-deriving any of that — the only
genuinely new derivations are the one-sided (`GL` from `AL` alone, `GR`
from `AR` alone) regularized-environment solves and the `H_AC`/`H_C`
effective-Hamiltonian actions themselves (see `vumps.py`'s own module
docstring for the full algorithm and diagram list). `self._vumps_result`
(a `pyitensor.vumps.VUMPSResult`, distinct from `self._result`'s
`IDMRGResult`) holds the converged state. Despite being a grouped-
supersite representation with no per-sublattice tensor list, `vumps.py`
itself now also implements `onsite_expectation`/`two_point_correlator`
directly on a `VUMPSResult` (mirroring `idmrg.py`'s own same-named
functions' public signature): unlike `idmrg.py`, which needs left- and
right-fixed-point eigenproblems because its own converged cell is only
approximately (left-)isometric, VUMPS's converged `AC` is already the
exactly-
normalized single-(super)site reduced state (Vanderstraeten, Haegeman,
Verstraete, arXiv:1810.07006, Eq.(34)) and `AR` is exactly right-
orthonormal, so both the left closure (via `AC`) and the right closure
(via `AR`, for a correlator spanning multiple unit cells) are exact
contractions with no eigensolve at all — a physical sub-site operator is
embedded into the grouped `d_g x d_g` space via `_embed_group_operator`
(`np.kron` in the same sequential order `_group_automaton` itself groups
physical indices in), and multi-cell separations propagate the open
right-bond object through plain `AR` transfer tensors
(`idmrg._apply_transfer_from_left`, reused as-is) before a final trace
closure. `Infinite_Many_Body_Chain.vev`/`correlator` dispatch on
`self.gs_method` to call either `idmrg.py`'s or `vumps.py`'s version
directly. `local_excitation_gap` remains `gs_method="idmrg"`-only (it
re-diagonalizes the growing algorithm's own final 2-site effective
Hamiltonian, which has no VUMPS equivalent; its `window>0` variant is
additionally `itensor_version="python"`-only, being an explicit prototype
rather than stable API). `itensor_version=3` mirrors all of this: both
`gs_method` values support `vev`/`correlator` there
(`Chain::vumps_onsite_expectation`/`vumps_two_point_correlator` and
`Chain::idmrg_onsite_expectation`/`idmrg_two_point_correlator`), and
`local_excitation_gap` at `window=0` maps onto
`Chain::idmrg_local_excitation_gap`. `excitation_energies`/`excitation_gap`
are
the opposite of `vev`/`correlator`: they require `gs_method="vumps"` (see
below), since the tangent-space excitation ansatz needs exactly the
mixed-gauge `{AL,AR,C,GL,GR}` `VUMPSResult` carries, which `IDMRGResult`
has no equivalent of.

**A confirmed, documented convergence-robustness gap, not fully closed.**
Plain single-attempt VUMPS from a random initial tensor was confirmed
directly, during development, to reliably land on self-consistent
(gauge-mismatch-converged) fixed points at an energy *worse* than the
same model's own smaller-`D` result — impossible for the true variational
optimum, so a genuine failure to find it, not merely slow convergence.
`vumps_ground_state` mitigates this with (1) a `D`-ramp (warm-start each
bond dimension from the previous, smaller one's own best solution,
`_grow_initial_state`) plus multiple restarts at every step, and (2) an
explicit variational-principle safety net that spends a bounded extra
attempt budget, warm-started from the same known-good smaller-`D`
solution, whenever a given `D`'s own search still lands worse than that
reference. Both measurably help (eliminating the worst, qualitatively
wrong failures observed with (1) alone, e.g. a `D=4` transverse-field
Ising run landing at a *positive* energy) but do not fully close the gap:
repeated calls on the same simple, gapped model at `D=4` were still
confirmed to occasionally (roughly 1 in 10 calls) land noticeably further
from the exact answer than the typical run. See `pyitensor/vumps.py`'s
own "Convergence robustness" module-docstring section for the full,
numerically-confirmed account, and `docs/user_guide.md`'s own iDMRG
section for the user-facing version of this caveat.

**How the `D`-ramp is walked, and where the two VUMPS backends now
diverge in cost.** The ramp visits `1, 2, 4, 8, ..., D` (doubling,
always ending exactly at `D`: `vumps.py`'s `_d_ramp`, mirrored in
`chain_session.h`'s own `d_ramp` vector). It used to step one integer at
a time, which made the *ramp*, not the solve at the target bond
dimension, dominate the whole driver: every step is a full multi-restart
VUMPS solve costing roughly `O(D^3)`, so a unit ramp is `O(D^4)` overall
and needs ~100 complete solves to reach `D=30` where doubling needs 6.
Warm-starting is unaffected -- `_grow_initial_state`/`vumps_grow_init`
embed the previous solution in the new tensors' leading block plus noise
and assume nothing about the size of the jump.

The remaining `O(D^6)` costs inside a single solve have been removed on
the `"python"` side only, and this is a real divergence between the two
ports rather than a difference of style. `pyitensor` computes the
transfer matrix's dominant fixed point with ARPACK on a matvec that
applies each site's own transfer tensor in turn (`idmrg.py`'s
`_dominant_fixed_point`), and solves the two environment systems with
GMRES on the same kind of matvec (`idmrg_excitations.py`'s
`_solve_linear_map`), instead of composing/materializing a
`chi^2 x chi^2` matrix and calling `np.linalg.eig`/`np.linalg.solve` on
it. Both keep the dense route as an exact fallback below a size
threshold (`_DENSE_EIG_MAX`, `_DENSE_SOLVE_MAX`) and on non-convergence,
so neither is less robust than what it replaced; `tests/
test_infinite_chain.py` pins the thresholds to force each route and
checks they agree. `Chain::vumps_ground_state` still uses the dense
`zgeev` fixed point (no ARPACK is available there, and a hand-rolled
Arnoldi is not worth the risk), so **the C++ VUMPS scales worse in `D`
than the Python one** even though both share the geometric ramp. Prefer
`itensor_version="python"` with `gs_method="vumps"`, or
`itensor_version=3` with `gs_method="idmrg"`, for a large-`D` infinite
chain.

**VUMPS and the tangent-space excitation ansatz, ported to
`itensor_version=3`.** `mpscpp3/chain_session.h`'s `Chain::
vumps_ground_state`/`Chain::vumps_excitation_energies` are C++ ports of
`pyitensor/vumps.py`'s `vumps_ground_state` and
`pyitensor/idmrg_excitations.py`'s `excitation_energies`, reached the
same way `idmrg_ground_state` already is (`infinitechain.py`'s
`gs_energy`/`excitation_energies` branch on `self.gs_method` for
`itensor_version=3` too now, not just `"python"`). Architecturally these
two new methods are a deliberate departure from every other `Chain`
method in this file: `idmrg_ground_state` (and everything else in
`chain_session.h`) is built from genuine ITensor `Index`/`ITensor`
tensor-network objects, but `vumps_ground_state`/
`vumps_excitation_energies` instead work entirely with plain dense
row-major `std::vector<Cplx>` arrays, closed over four LAPACK routines
already vendored with ITensor v3's own build (`itensor::zgeev_wrapper`
for the transfer-matrix dominant-fixed-point eigenproblem,
`zheev_wrapper` for the `H_AC`/`H_C` ground-state eigenproblem in place
of pyitensor's own truncated Lanczos, `zgesv_wrapper` for the
regularized `GL`/`GR`/`GBL(k)`/`GBR(k)` environment linear solves, and
`zgesvd_wrapper` for both the `AL`/`AR` orthogonal-Procrustes update and
the null-space isometry `V_L` the excitation ansatz needs) rather than
ITensor's own tensor contraction/decomposition machinery. This is a
deliberate simplification, not a shortcut taken for lack of a better
option: `D` and `d_g` (the grouped supersite's own physical dimension)
are always small in this feature's scope (`n_uc<=2`, reach-1 bonds
only, the same scope `idmrg_ground_state` itself has), so exact,
dense-matrix linear algebra is both simpler to get right and
dramatically lower-risk than re-deriving VUMPS's and the excitation
ansatz's own already extremely subtle fixed-point/gauge bookkeeping
(the eight independent investigation passes documented in
`pyitensor/idmrg_excitations.py`'s own "History" section, and the
missing-conjugate sign bug found along the way — see above) a second
time against ITensor v3's stricter `Index`/`IndexSet` API. Concretely,
this means the C++ port is close to a direct, line-for-line translation
of the numpy-array-based Python algorithm (same `(D,d_g,D)`/`(D,D)`
shape conventions, same automaton channel layout) rather than an
independent re-derivation in tensor-network form — the translation risk
this trades away is "did I copy the formula correctly", not "did I
re-derive VUMPS/the excitation ansatz correctly", which is exactly the
risk category the eight-pass investigation showed to be the dangerous
one. One initialization detail is intentionally NOT ported faithfully:
pyitensor's own `_random_initial_state`/`_grow_initial_state` canonicalize
a raw random tensor via `idmrg.py`'s fixed-point-based
`_canonicalize_periodic` (needed there because that function is also
used for genuine bond truncation elsewhere in `idmrg.py`); the C++ port
instead uses plain Gram-Schmidt orthonormalization of a random matrix's
columns (`vx_gram_schmidt_columns`) to build an exactly isometric
starting `AL0`/`AR0` directly — mathematically equally valid for VUMPS
initialization purposes (any exactly isometric starting tensor is a
legitimate starting point; VUMPS's own iteration, not the initialization
method, is what drives the state to a self-consistent, Hamiltonian-
optimal fixed point), and considerably simpler since the C++ port never
needs `_canonicalize_periodic`'s own bond-truncation machinery at all.
The D-ramp/multi-restart robustness strategy (`vumps_ground_state`'s own
outer loop) IS ported; the "spend a bounded extra attempt budget to beat
an already-known smaller-`D` energy" safety net documented above is not
(see `Chain::vumps_ground_state`'s own comment) — a caller doing
quantitative work at a larger `D` should still call it several times
independently and keep the lowest reported `e0`, the same discipline
pyitensor's own docstring already recommends.

Validated directly against `itensor_version="python"`'s own VUMPS
(cross-checked, not just internally self-consistent): ground-state
energy density matches to ~1e-10 or tighter at `D=1,2,3` on the
transverse-field Ising model (`n_uc=1`) and the uniform Heisenberg chain
(`n_uc=2`, exercising `vumps_group_automaton`'s two-sublattice grouping
branch); the excitation dispersion matches across a full momentum scan
to the same tolerance on the gapped TFIM case, and to ~1e-4..1e-7 on the
gapless/critical Heisenberg case (attributed to the two backends'
independent, non-convex VUMPS restart searches landing on slightly
different local optima there, not a discrepancy in the ported algorithm
itself — see `tests/test_vumps_v3.py`/`tests/test_vumps_excitations_v3.py`).

**Static correlators, ported to `itensor_version=3` too.**
`Chain::vumps_onsite_expectation`/`Chain::vumps_two_point_correlator` are
the same kind of line-for-line dense-array C++ port as `vumps_ground_state`/
`vumps_excitation_energies` above, this time of `pyitensor/vumps.py`'s own
`onsite_expectation`/`two_point_correlator` (see this section's own earlier
paragraph for that algorithm — `AC`'s exact normalization and `AR`'s exact
right-orthonormality, both direct consequences of the converged mixed
gauge, mean no eigensolve is needed at either the left or the right
closure). `Chain::vumps_embed_group_operator` (a new private helper) is the
C++ analogue of `_embed_group_operator`'s `np.kron`-based embedding, built
from the SAME per-sublattice dense operator matrices `idmrg_op_dense`
already reads off `SiteSet::op()` for `idmrg_ground_state`'s own automaton
construction, combined in the identical sequential (site 0 slowest, site
`n_uc-1` fastest) order `vumps_group_automaton` itself uses for the grouped
physical index — no separate re-derivation of that ordering convention.
`Infinite_Many_Body_Chain.vev`/`correlator` route to it whenever
`itensor_version=3` and `self.gs_method=="vumps"`, calling `gs_energy()`
first (or re-running it, via the same `self._session3_has_vumps`
bookkeeping `excitation_energies` already relies on) if the cached session
does not already hold a converged VUMPS snapshot. Validated directly
against `itensor_version="python"`'s own VUMPS correlators (not just the
same exact-D=1 closed-form cases both backends already had — a
field-polarized product state and a decoupled Heisenberg singlet dimer):
matches to ~1e-14 or tighter across the ground-state energy, `vev`, and a
multi-`r` correlator scan at `D=2,3` on the transverse-field Ising model —
see `tests/test_vumps_correlator_v3.py`.

**A real, previously-undetected bug this same audit found and fixed in
the already-shipped `idmrg_build_row`** (shared by `idmrg_ground_state`
and the new `vumps_ground_state`, since both build their automaton from
it): a site's own onsite/field-operator term was being added onto the
automaton's `F,F` self-loop transition instead of the direct `S,F`
transition — exactly the bug `pyitensor/idmrg.py`'s own
`_build_periodic_mpo` docstring documents finding and fixing on the
Python side (putting it on `F,F` instead of `S,F` either silently drops
every onsite term for a bond-term-free Hamiltonian, since `W` then stays
block-diagonal between `{S}` and `{F,pending...}`, or causes an
exponential-in-macro-iteration blow-up once at least one bond term
activates `F`, since every further-absorbed site re-adds the onsite
content on top of an already-summed total). Never caught by the existing
test suite because no `itensor_version=3` `idmrg_ground_state` test
exercises an onsite/field term (every `tests/test_infinite_chain.py` v3
case is a pure bond-only Heisenberg/XXX-style Hamiltonian) — found while
building this feature's own VUMPS test cases, since a TFIM Hamiltonian
(which has a field term) is the natural first thing to validate a new
VUMPS port against.

The Hamiltonian is built from `SxL[i]`/`SxC[i]`/`SxR[i]`-suffixed
operators (previous/central/next unit cell, `i=0..n_uc-1`) rather than
absolute site indices; `_canonicalize_hamiltonian` validates every term
touches at most two adjacent cells and shifts away the `L` suffix (`L→C`,
`C→R`) before storing, so the stored Hamiltonian is uniformly "intra-cell
`C` terms" + "`C`-to-`R` inter-cell terms", each physical bond attributed
to exactly one cell.

`pyitensor/idmrg.py` implements the actual growing/infinite-size
algorithm (White 1992, with an `n_uc`-site unit cell): `_build_automaton`
builds the periodic per-sublattice MPO tensor directly as a hand-rolled
finite-state automaton (not extracted from a compressed finite reference
system — an earlier design that was tried and abandoned, see the module's
own docstring for the two independent bugs that approach hit), then
`idmrg_ground_state` grows two environments `HL`/`HR` from `HL=HR=None`,
one unit cell at a time, reusing `dmrg.py`'s own `_lanczos_ground_state`
and `kernels.make_matvec` for each micro-step's local 2-site solve (no
separate Lanczos implementation).

Three ingredients beyond the plain growth loop are what make the *state*
usable, not just the energy. All three were missing originally, and the
failure was silent — the energy was correct throughout — so the history is
worth keeping:

- **The wavefunction prediction** (`_wavefunction_prediction`, McCulloch
  2008 and the reference implementation vendored at
  `mpscpp2/ITensor/itensor/mps/idmrg.h`). Each micro-step's Lanczos starts
  from the previous step's converged state translated into the new bond
  bases, `theta = lambda_k . B_k . lambda_{k-1}^{-1} . A_k . lambda_k` —
  `idmrg.h`'s `lastV = PseudoInvert(dag(D))` together with its
  `swapUnitCells`. Substituting the canonical forms collapses it to
  `lambda.Gamma.lambda.Gamma.lambda`, i.e. it is not a heuristic warm start
  but the exact translation of the converged state into the enlarged bases.
  Note the swap (the previous step's right-canonical `V` supplies the new
  *left* physical slot): that is what makes the sublattice labels line up,
  and it only does so for `n_uc<=2`, one more reason `n_uc>=3` is rejected.
- **A gauge-consistent extraction** (`_theta_cell`). The prediction fixes
  the trajectory but not the readout. Chaining one left-canonical factor
  per micro-step, as `U_list` does, still compares bond bases minted by
  *different* micro-steps, so tiling it inserts a spurious unitary at every
  unit-cell boundary. The converged cell is instead taken from a single
  micro-step's own `theta`, whose two outer legs are the same bond by
  construction, then re-gauged (`_canonical_theta_cell`, reusing
  `_canonicalize_periodic`) and normalized to a transfer eigenvalue of
  exactly 1. This is `IDMRGResult.cell_list`, and it is what
  `_correlator_cell` hands to every static observable, `apply_mpo`,
  `imps_overlap` and `imps_sum`. `IDMRGResult.cell_raw` is the same cell
  *before* re-gauging — needed by `pyitensor/idmrg_window.py`, since only
  there are the outer bonds still literally in the `env_HL_ket`/
  `env_HR_ket` bases its environment caps are expressed in. `U_list` is
  retained but is not a valid periodic MPS on its own; nothing should tile
  it.
- **An energy baseline** (`_subtract_energy_baseline`, the equivalent of
  `idmrg.h`'s `HL += -energy*IL`). No separate identity environment has to
  be accumulated for it: the automaton's channel space already carries one,
  since `chans[p]` is `[S, F] + pending` and a growing environment's `S`
  and `F` components *are* its identity and energy accumulations. The whole
  operation is therefore a slice update — `HL[F] -= e*HL[S]` on the left
  and the mirrored `HR[S] -= e*HR[F]` on the right. Without it the
  superblock eigenvalue grew extensively while `_lanczos_ground_state`
  stopped on a *relative* criterion, so the absolute noise in the
  finite-difference density grew linearly with the iteration count:
  measured on a gapped `n_uc=2` chain, `|E| = 603` after 400
  macro-iterations with the density jittering over ~9e-11 long after it had
  physically converged, against a default `etol` of 1e-10. It is now
  `|E| = 3.8` and ~4e-15, so `etol` is meaningful down to machine
  precision. The finite-difference estimator is kept and simply adds the
  shift back, which is exact as long as the shift is constant across the
  two macro-iterations being differenced — hence it is latched once and
  never revised.

Static correlators are reconstructed after convergence from `cell_list`
via the standard infinite-MPS transfer-matrix formalism. `_expectation`
closes the operator string against the transfer operator's own left fixed
point *at the string's own starting cell position* and divides by the same
contraction with the operators dropped, rather than assuming the left fixed
point is the identity. That assumption holds only for exactly
left-canonical tensors; using the cell's position-0 fixed point instead left
`<Sz(0)Sz(0)>` at 0.2499994784 rather than the exactly-0.25 it must be for
spin-1/2.

**How the static correlators got fixed, and how they were diagnosed.**
Before the prediction and the theta-cell extraction above, the reported
*energy* was correct while the extracted *state* was not — a failure mode
worth recognizing, since every energy-based check passes straight through
it. The diagnostic that exposes it is the model-agnostic
`<H_uc>`-equals-`n_uc*density` identity: exact for any translationally
invariant state, and checkable with the module's own correlator machinery
without needing `e0` in closed form.

An intermediate fix threaded each macro-iteration's own local ground vector
back in as the next one's Lanczos warm start, reasoning that a
gapless/SU(2)-symmetric model has a (near-)degenerate local ground manifold
and that a fresh random restart lets Lanczos land on an arbitrary member of
it each time. That story is real, but it was not the cause: the same
failure persisted at 4.5e-2 on a *gapped, non-degenerate* TFIM whose energy
was exact to 1e-12, where there is no degenerate manifold to wander in. The
raw flat vector being reused carried no basis information, so it could not
reconcile two macro-iterations' bond bases — which is what actually
corrupted `U_list`. The prediction subsumes the continuity argument and is
additionally a genuine basis change; the extraction fixes the readout.

Measured across gapped TFIM and the exactly solvable XX chain, with the
energy correct to ~1e-14 throughout: `<H_uc> - n_uc*e0` went from missing by
up to 1.2e-1 to 1e-15..1e-9. On the XX chain at `n_uc=1`, `<Sz(0)Sz(1)>`
went from +0.028 against an exact -0.101 — the wrong *sign* — to
reproducing it within the bond dimension's own truncation error, and `<Sz>`
from -0.13 to 1e-13 against an exact 0. The error was also not monotone in
`maxm`/`maxiter` beforehand (4.3e-5 at `maxiter=60` versus 1.2e-1 at 240),
the signature of an unfixed gauge rather than an unconverged one.
`IDMRGResult.state_overlap` — now the overlap between the solve's result
and the prediction it started from, i.e. the direct analogue of `idmrg.h`'s
own "Overlap of initial and final psi" — reaches 1-1e-13 instead of
plateauing around 0.5-0.65.

**Currently limited to `n_uc<=2`** (enforced at `Infinite_Many_Body_Chain`
construction) — the micro-step loop's sublattice pairing
(`p_L=mstep`, `p_R=n_uc-1-mstep`) only produces two genuinely *adjacent*
active sites for `n_uc` in `{1, 2}`; for `n_uc>=3` an intermediate
micro-step's two active sites are several real, not-yet-inserted sites
apart, and `_build_automaton`'s own per-sublattice pending-channel
content differs between them (not just its count), so there is no way to
even wire them together correctly without a genuine redesign of the
growth scheme — flagged as a known follow-up rather than attempted.

**A confirmed, fixed bug worth knowing about if this module is touched
again**: `idmrg_ground_state` re-mints a fresh mpo-axis Index for `HL`'s
and `HR`'s own environment tensor at the end of *every* micro-step,
rather than reusing whatever Index `_build_automaton` happened to label
that channel space with. This is required because `_build_automaton`'s
`boundary_idx` pool has only `n_uc` distinct Index objects total, reused
verbatim every time a given sublattice position comes up again (every
macro-iteration past the first) — left un-broken, `HL`'s and `HR`'s own
mpo axes can end up being the *same* Index object (guaranteed for
`n_uc=1`, where a single sublattice's left and right automaton legs are
literally identical), corrupting the effective Hamiltonian from the
second macro-iteration onward once both environments are non-trivial and
used together in the same local 2-site solve. Confirmed directly by
extracting the dense effective Hamiltonian at macro-iteration 2 of a
uniform `n_uc=1` Heisenberg chain and finding its full eigenvalue
spectrum did not match exact 6-site ED at all, while macro-iteration 1's
spectrum matched to machine precision — the same self-referencing-Index
hazard already seen once in this module (`_project_channel`'s own
docstring), here manifesting across the growth loop's iterations rather
than within a single tensor. `idmrg.py`'s module docstring and the inline
comments in `idmrg_ground_state`/`_relabel_pos` have the full derivation.

A related symptom, now gone: `U_list[0]`'s own left bond and
`U_list[n_uc-1]`'s own right bond are, in a truly periodic MPS, the *same*
wraparound bond, but each was independently truncated by its own
micro-step's SVD, so nothing guaranteed they even came out the same
*dimension*. `_dominant_right_fixed_point` raises a clear `RuntimeError`
rather than a cryptic `reshape` crash when that happens, and its comment
once described the case as intermittent — measured, it fired in 10 of 20
`(n_uc, maxm, maxiter)` combinations on a *gapped* model, i.e. it was the
common case, and it is the same root cause surfacing as a hard failure
instead of a wrong number. Since the static path tiles `cell_list`, whose
two ends are one bond by construction, none of those 20 raise any more. The
check is kept as defense in depth for anything that still hands
`_transfer_matrices` a raw `U_list`.

**Applying an MPO to the converged iMPS**: `idmrg.py` also implements
`apply_mpo(result, W_bulk, cutoff, maxdim)` — the infinite-chain analogue
of `mpsalgebra.applyMPO`, and the primitive a future iTEBD-style
real/imaginary-time evolution feature would build on. `grow_by_mpo`
Kronecker-merges `W_bulk`'s own per-site bonds with the converged cell's (mirroring
`mpsalgebra._apply_chain`'s per-site body, generalized to a periodic
index range including the wraparound cut), giving raw, non-canonical
tensors at bond dimension `chi_A*chi_W`; `_canonicalize_periodic` then
re-canonicalizes/truncates them back down via the standard two-sided
fixed-point infinite-MPS procedure (Orús & Vidal, PRB 78, 155117 (2008)):
factor the dominant left/right transfer-matrix fixed points at every cut
into square-root factors `X_p`/`Y_p` (`_psd_sqrt_factor`), SVD
`M_p=X_p@Y_p` and truncate (reusing `svd.py`'s own `_truncate`), then
build an *asymmetric* gauge (`G_left[p]=U_p^H X_p`, no `S` factor at all;
`G_right[p]=Y_p V_p^H S_p^{-1}`, the *full* inverse) — asymmetric, rather
than a textbook symmetric `S^{-1/2}` split on both sides, specifically so
the result comes out plain left-canonical (matching `IDMRGResult.cell_list`'s
own convention) with no separate bond-weight layer to thread through,
confirmed numerically before committing to the formula (a first,
symmetric-split derivation looked plausible by analogy to Vidal's
canonical form but reproduced Identity only to ~1, not a small residual).
`_all_left_fixed_points` (the new, symmetric counterpart to the existing
`_all_right_fixed_points`) needed one more fix beyond the obvious
transpose-the-eigenproblem mirroring: its own per-cut renormalization
(needed to avoid blow-up across repeated propagation, exactly like the
existing right-fixed-point code) discards the *relative* scale between
different cuts that `_canonicalize_periodic`'s left-gauge construction
turns out to depend on (unlike the right-gauge one, which is provably
scale-invariant) — tracked via an explicit accumulated-scale return value
and undone before use. A second, independent bug (an index-order/
transpose mismatch between `_apply_transfer_from_left`'s own (ket,bra)
convention and the (bra,ket) convention the isometry derivation needs)
was only found by deriving the exact identity the gauge must satisfy by
hand and checking it numerically term-by-term — a genuine case where
tests alone (which caught that *something* was wrong, via a hand-crafted
internal left-canonicality check) weren't enough to *localize* the bug,
only real algebra was. **Scope restriction**: `W_bulk` must be a
*bounded* (non-extensive) periodic operator — the Hamiltonian's own
`_build_periodic_mpo`-built automaton is deliberately the *other* kind
(an unconditional accumulator self-loop that needs a genuine chain
boundary to mean "sum", not "product") and is not a valid `apply_mpo`
input; see `idmrg.py`'s own "Applying a (bounded) MPO to the converged
iMPS" section docstring for the specific failure mode this would hit.

**Overlap between two converged iMPS**: `idmrg.py` also implements
`imps_overlap(result_a, result_b, normalize=True)` — the thermodynamic-
limit analogue of a finite MPS inner product. Since a literal
`<phi|psi>` over an infinite chain scales as `eta**N` over `N` unit
cells (not a finite number unless `|eta|==1` exactly), the default
return value is instead the *per-site fidelity*
`eta_ab/sqrt(eta_aa*eta_bb)`, where `eta_ab` is the dominant eigenvalue
of the *mixed* transfer matrix built from the two states' own tensors
(`_transfer_matrices` generalized with a new optional `bra_list`
parameter, defaulting to the self-overlap case every existing caller
still uses) and `eta_aa`/`eta_bb` are each state's own self-overlap
(via the existing `_dominant_right_fixed_point`, reused verbatim — valid
since a self-overlap transfer tensor is always square). The mixed case
needed a dedicated eigen-solve, `_dominant_eigenvalue_mixed`, rather
than reusing `_dominant_right_fixed_point` directly: two independently
converged/truncated iMPS need not share a bond dimension, so the
composed mixed transfer tensor is `(chi_a,chi_b,chi_a,chi_b)`, square
only in the `chi_a==chi_b` self-overlap special case
`_dominant_right_fixed_point` itself assumes.

**A second, real bug found while implementing `imps_overlap`** (needed a
Hamiltonian with an onsite/field term to exercise a meaningfully
non-trivial "orthogonal states" test case, which surfaced two bugs in
`_build_periodic_mpo`'s onsite-term handling, both pre-dating and
unrelated to `imps_overlap` itself): `_onsite_matrix`'s own unpacking
loop expected 3-tuples (`rel_site, coef, mat`) but `_build_periodic_mpo`
had already consumed `rel_site` as its `onsite_by_p` dict key, storing
only `(coef, mat)` pairs — a `ValueError: not enough values to unpack`
for *any* Hamiltonian with an onsite term (dead code before this, since
every existing test used bond-only Hamiltonians). Fixing just that
crash would have been unsafe: the underlying automaton wiring was also
wrong. Onsite content was added on the accumulator's own `F,F`
self-loop (`Id + onsite_mat`) rather than a direct `S,F` transition — no
existing branch ever populates `S,F` at all, so (1) for a bond-less
Hamiltonian, `W` stays block-diagonal between `{S}` and `{F,...}`
forever, and the boundary-extracted `<S|W^N|F>` (exactly what
`_project_channel`'s boundary tensors compute) is identically 0 for
every `N` — confirmed directly: a pure-field `n_uc=1` chain converged to
energy density 0 and `<Sz>~0` regardless of field strength, instead of
the exact, closed-form `-|B|/2`/`sign(B)*0.5`; and (2) once a bond term
*did* activate `F`, every further site absorbed while `F` was already
hot re-added the onsite content on top of an already-summed running
total (an exponential, multiplicative blow-up, not the intended linear
sum) — confirmed directly: energy density reached `-1e23` by `B=0.3` and
`-1e69` by `B=1.0` on an otherwise-normal dimerized `n_uc=2` chain, with
`.converged` staying `False` throughout instead of a modest,
field-dependent shift. The fix moves the onsite content onto the direct
`S,F` transition (leaving `F,F` as plain `Id`) — the standard textbook
automaton-MPO construction for summing an onsite term into a running
total (the same upper-triangular-in-channel-index structure any
finite-range Hamiltonian MPO derivation uses), verified both by a
by-hand matrix-product induction argument and numerically (the same
field tests above now reproduce the exact closed-form field-only result
and a finite, converging energy density with a bond term present).

**Direct sum of two converged iMPS**: `idmrg.py` also implements
`imps_sum(result_a, result_b, cutoff, maxdim)` — the periodic-chain
analogue of `mpsalgebra.sum`. `_periodic_direct_sum` block-diagonally
concatenates `result_a`'s and `result_b`'s own bond spaces at *every*
unit-cell cut, including the wraparound (unlike `mpsalgebra.sum`'s finite
chain, which has two dimension-1 boundaries to instead concatenate along
— a periodic chain has no such boundary, so the sum is block-diagonal
everywhere); the raw result is then re-canonicalized/truncated by reusing
`apply_mpo`'s own `_canonicalize_periodic` verbatim, no new gauge-fixing
algorithm needed.

**A real, load-bearing correctness issue found and fixed while
implementing this** (not merely a theoretical concern — reproduced
directly before the fix): `_canonicalize_periodic`'s two-sided fixed-point
procedure assumes the combined transfer matrix has a *unique* dominant
eigenvalue, picked via `np.argmax(np.abs(w))` in
`_dominant_right_fixed_point`. Every `IDMRGResult` is individually
normalized to self-overlap eigenvalue `eta=1` exactly (left-canonical SVD
construction), so summing two *different*, ordinarily-converged
`IDMRGResult`s — the natural use case — produces a combined transfer
matrix with an *exactly degenerate* leading eigenvalue (two blocks, each
eigenvalue 1). Before a guard was added, `np.argmax` silently picked
*one* of the two tied eigenvectors (confirmed directly: summing two
independently-converged, oppositely Sz-polarized `n_uc=1` product states
reliably collapsed to just one branch, with *which* branch survived
sensitive to floating-point tie-breaking noise inside `np.linalg.eig` —
not even reproducible run to run), silently discarding the other branch
entirely while still returning a plausible-looking, correctly-normalized
`PeriodicMPS`. This is the classic "looks right, isn't" failure mode this
codebase's own testing culture exists to catch. The fix: a shared helper,
`_check_dominant_eigenvalue_nondegenerate`, checks the top two eigenvalue
magnitudes of the composed transfer matrix and raises `RuntimeError` when
their ratio exceeds `1 - 1e-9`. This threshold is set far tighter than
the loosest value that would have passed every legitimate case tried —
the genuine tie this check exists to catch lands at a `~1e-15` relative
gap (exact to machine precision), so `1e-9` still rejects it with 6
orders of magnitude to spare, while leaving a wide safety margin below:
a gapped dimerized `n_uc=2` chain's top two magnitudes were `(1.0,
0.093)`; a gapless/critical uniform `n_uc=1` chain showed a clear
maxm-dependent narrowing (top-two gap `0.376` at `maxm=40`, `0.116` at
`maxm=60`, `0.160` at `maxm=80` — not perfectly monotonic, since `n_uc=1`
has its own known convergence limitation, see above, but never close to
tight) — a threshold merely "loose enough for today's test suite" would
risk a *future* false positive on an ordinary, non-degenerate correlator
call at a larger `maxm` a user might reasonably pick for better accuracy
near criticality, so `1e-9` was chosen deliberately over a looser value
like the originally-tried `1e-6`. Physically, the tied case is a "cat
state" superposition of two macroscopically distinct branches;
representing it as a single injective/canonical periodic MPS is out of
scope for this module's existing (single-fixed-point) correlator
machinery — `imps_sum` therefore raises clearly there rather than
silently returning one arbitrary branch, mirroring `apply_mpo`'s own
"bounded operators only" scope restriction in spirit. The shared helper
is used by `_dominant_right_fixed_point`, `_dominant_left_fixed_point`,
and `_dominant_eigenvalue_mixed` (the mixed-transfer-matrix eigenvalue
`imps_overlap`'s own cross term uses) alike — all three previously used
an unguarded `np.argmax`, the identical latent bug, just not yet exercised
by an existing caller for the left/mixed cases — so this same guard now
protects `onsite_expectation`/`two_point_correlator`/`imps_overlap` for
free too; confirmed not to trip on any pre-existing test or on
`imps_overlap`'s own realistic (orthogonal and gauge-comparison) cases
(all previously-passing `tests/test_infinite_chain.py` cases still pass).

**Correction: that "wide safety margin" did not hold, and the guard was
rejecting physics.** The margin argument above was measured on Heisenberg
chains, where the top-two magnitude gap stayed at `0.376`/`0.116`/`0.160`.
It does not generalize. A *half-filled free-fermion* chain is critical with
`2k_F = pi`, so its transfer matrix carries a peripheral eigenvalue at
phase `pi` whose magnitude approaches the leading one as the correlation
length diverges. Measured on exactly that model at `D=16` (the `n_uc=1`
fermionic case in `tests/test_infinite_chain.py`), *every* observed firing
had the identical signature

```
|lambda| = (1.0, 0.999999999)      arg = (0, +-pi)
```

— a `1e-9` gap, i.e. sitting precisely on the threshold this section
argued was 6 orders of magnitude clear of any legitimate case. The cost was
not subtle: **~42% of all individual VUMPS attempts and ~16% of whole
ground-state solves** on that model died here, surfacing as an
intermittent `RuntimeError: vumps_ground_state: every attempt at D=16
failed (degenerate transfer-matrix spectrum)`. The test had already been
given `vumps_nrestarts = 10` to paper over it.

The error was conceptual rather than a bad constant. `+rho` and `-rho` are
**distinct eigenvalues with distinct eigenvectors** — nothing is
ill-defined about them. A transfer matrix's peripheral spectrum is
`rho * exp(2 pi i k / p)` for a period-`p` state, and only the `k=0`
member carries the positive(-semidefinite) eigenvector that is the density
matrix every caller wants. So a magnitude tie there is *well-posed and
must be resolved*, not rejected; only a repeated **eigenvalue** (the cat
state, where both copies sit at `+1`) is genuinely ambiguous. Worse, since
`np.argsort` orders magnitude-tied entries arbitrarily, the old code could
also have returned the `-rho` eigenvector as a "density matrix" — not
positive, and with a trace near zero that callers divide by.

`_check_dominant_eigenvalue_nondegenerate` therefore now takes a `perron`
flag selecting between two genuinely different questions:

- `perron=True` (`_dominant_fixed_point`): degenerate means the top two
  entries are close **as complex numbers**, and magnitude ties within
  `_PERRON_TIE_RTOL` (`1e-6`, deliberately wider than `_DEGENERACY_RTOL`
  so the selection is not decided by rounding at the observed `1e-9`) are
  broken toward the largest real part — the Perron root.
- `perron=False` (`_dominant_eigenvalue_mixed`): unchanged magnitude test.
  `imps_overlap` extracts a per-unit-cell factor `|eta|^N`, which
  oscillates rather than converging when two eigenvalues share a
  magnitude, so there a magnitude tie really is ambiguous.

The cat-state case the guard exists for is unaffected
(`tests/test_vumps_imps_sum.py` still passes unchanged), and
`tests/test_transfer_degeneracy_guard.py` now pins both halves directly:
period-2 and period-3 peripheral spectra accepted with the Perron root
selected, genuine repeated eigenvalues still rejected, the mixed caller
still rejecting magnitude ties, and the returned fixed point verified
positive semidefinite.

**The same bug existed in `mpscpp3`, where it was worse.** `chain_session.h`
carried the identical magnitude-only test at three sites (`vx_dominant_eig`,
`ic_dominant_fixed_point`, `vms_environments`). There it did not merely
flake: the grouped (`n_uc<=2`) v3 VUMPS path could not complete a `D=8`
Heisenberg ground state at all, failing with *"every attempt at D=4 failed
(degenerate transfer-matrix spectrum)"*, and `n_uc=1` failed at every bond
dimension tried. All three now share `vx_perron_reorder` /
`vx_check_perron_nondegenerate`, the C++ counterparts of the Python fix.

One extra step is needed in the matrix-free path that the dense one does not
require. `ic_arnoldi_dominant` reconstructs and returns exactly ONE
eigenvector, so a caller handed the magnitude-dominant Ritz pair cannot
recover the Perron one if `-rho` edges ahead — the selection has to happen
inside. Selecting on the Ritz *eigenvalue* alone is also not enough, because
Ritz values are Krylov approximations and which of a `+rho`/`-rho` pair comes
out larger in magnitude is partly numerical accident: measured on a critical
chain at `D=16`, it returned the `-rho` branch and the caller then died on
*"dominant left fixed point has ~zero trace"*. So the peripheral candidates
are separated by the property the caller actually needs — a period-2 state's
`-rho` eigenvector carries opposite signs on the two sublattices and its
trace cancels to ~0, while the Perron one is positive with an O(1) trace, and
every caller of this function divides by that trace. `ic_arnoldi_tie_rtol_`
(`1e-3`) is correspondingly much looser than `_PERRON_TIE_RTOL`, since it is
applied to approximate Ritz magnitudes rather than exact eigenvalues.

That intermediate state — Perron reorder in place but selection still driven
by the eigenvalue — was caught by `test_idmrg_correlator_v3.py`'s
`test_local_gap_survives_a_positive_stored_superblock_energy`, which flipped
the sign of its starved run's stored superblock energy. The test was right;
it went green again once the selection used the trace.

**Performance, which is what the guard fix unblocked.** With the path usable
again, the grouped `n_uc<=2` route was still far slower than the pure-Python
one it ports, for the reason the multi-site (`vms_*`) path had already fixed:
a dense `zgeev` on the whole `(D^2,D^2)` transfer matrix (O(D^6), run twice
per iteration, once per side), and a `vx_build_linear_map` + dense LU in
`vx_regularized_solve` — that helper applies the action `D^2` times purely to
MATERIALIZE a matrix which then gets another O(D^6) factorization. Both now
have matrix-free branches above their thresholds (`vx_dominant_pair` using
the same `ic_arnoldi_dominant`, and `vx_bicgstab` respectively), reusing the
multi-site path's own helpers rather than adding second implementations, and
both keep the dense route below threshold so the small cases every v3 VUMPS
test was validated against stay on the path they were validated on.

Measured on a uniform Heisenberg chain, `nrestarts=2`, against the
pure-Python backend as reference:

| case | pyitensor | v3 before | v3 after |
|---|---|---|---|
| `n_uc=1`, D=8/16/24 | 6.8 / 36.2 / 171.5 s | *failed at every D* | 3.3 / 16.5 / 83.2 s |
| `n_uc=2`, D=8 | 5.3 s | 4.3 s | 4.0 s (dense, unchanged) |
| `n_uc=2`, D=16 | 22.9 s | 133.7 s | 29.3 s |
| `n_uc=2`, D=24 | 116.8 s | 1483.4 s | 45.8 s |

i.e. `n_uc=2` at D=24 went from 12.3x SLOWER than pure Python to faster than
it, and the `n_uc=1` path from unusable to working, with `e0` agreeing with
pyitensor to all printed digits at D=24 (-0.443136013 on both). Note the
pure-Python column is itself noisy run to run (random restarts), so treat the
ratios as order-of-magnitude.

**Direct sum of two converged VUMPS iMPS (`pyitensor/vumps.py`'s
`imps_sum`)**: the VUMPS-mixed-gauge port of `idmrg.py`'s own `imps_sum`
above, following the same literature check this project's own standing
practice requires before new physics work (arXiv:1810.07006, "Tangent-space
methods for uniform matrix product states", Sec. 2.1's Eq.(9)-(17)
"Algorithm for finding canonical forms", plus its own remark that
non-injective/"cat state" MPS tensors — exactly what a degenerate dominant
transfer-matrix eigenvalue signals, per that section's own footnote 4 —
appear with "measure zero" and are excluded throughout that reference too,
independently confirming the same scope limit `idmrg.imps_sum` already
documents). Since VUMPS already works at a single grouped-supersite level
(`n_uc` sites folded into one `d_g`-dimensional site via
`_group_automaton`), there is only ever one cut to sum across — no
per-sublattice list construction like `idmrg._periodic_direct_sum` is
needed, just a single block-diagonal `(D_a+D_b, d_g, D_a+D_b)` tensor built
from `result_a.AL`/`result_b.AL` directly. Both are already individually
left-canonical (isometric) by construction of any converged `VUMPSResult`,
so this raw direct sum is *already* exactly left-canonical too — a
block-diagonal direct sum of two isometries is itself an isometry
(confirmed algebraically: `sum_p AL_sum_p^dagger AL_sum_p =
block_diag(sum_p AL_a_p^dagger AL_a_p, sum_p AL_b_p^dagger AL_b_p) =
block_diag(I,I) = I`). `idmrg._canonicalize_periodic` is still reused on
this raw tensor (as a trivial `n_uc=1` "periodic chain") rather than
skipped, for two reasons: it applies the caller's requested
`(cutoff, maxdim)` truncation via a genuine two-sided fixed-point SVD (not
a naive per-block truncation), and — crucially — it is where the exact
same `idmrg._dominant_right_fixed_point` degeneracy check `idmrg.imps_sum`
relies on fires automatically for the "two ordinary/tied-eta branches"
case, with no new degeneracy-detection code needed. A new helper,
`_complete_mixed_gauge`, then completes the resulting truncated
left-canonical `AL` to the full mixed gauge `{AL,AR,C,AC}`: factor `AL`'s
own dominant right transfer-matrix fixed point `r = C C^dagger` (Hermitian
PSD square root via `eigh`, mirroring `idmrg._psd_sqrt_factor`'s own
Hermitizing-before-square-root reasoning), then `AR := C^-1 AL C` and
`AC := AL C` (`== C @ AR` as an exact algebraic identity, not merely
approximate, since `AR` is defined in terms of `C` in the first place).
This is deliberately *not* a further VUMPS energy-minimization step
(`imps_sum`'s output is generally not an eigenstate of anything) — purely
gauge bookkeeping, unlike `_random_initial_state`'s own two-direction
canonicalization trick (which leaves `C` at the identity, relying on a
follow-up VUMPS iteration to fix it up; `imps_sum` has no such follow-up,
so `C` is solved for properly here instead). The result is a new
lightweight `UniformMPS` class (mirrors `idmrg.PeriodicMPS`'s relationship
to `IDMRGResult`), accepted directly by `onsite_expectation`/
`two_point_correlator` with no changes to either (both only ever read
`.sites_uc`/`.n_uc`/`.AC`/`.AR`). Verified end-to-end against the
untouched dominant branch's own `onsite_expectation`/`two_point_correlator`
at genuinely entangled `D=2,3` (TFIM, gapped) to ~1e-6, not just the
trivial `D=1` scalar case — see `tests/test_vumps_imps_sum.py` and
`examples/idmrg/vumps_imps_sum/main.py`.

**Applying an MPO to a converged VUMPS iMPS (`pyitensor/vumps.py`'s
`apply_mpo`)**: the VUMPS-mixed-gauge port of `idmrg.py`'s own `apply_mpo`
above, reusing that function's construction and literature reference
(Orus & Vidal, PRB 78, 155117 (2008)) unmodified rather than re-deriving
it — no new physics beyond what `idmrg.apply_mpo` and `vumps.imps_sum`
(immediately above) already establish independently. `W_bulk` takes
*exactly* `idmrg.apply_mpo`'s own convention (a list of `n_uc` rank-4
`(Left, in, out, Right)` ITensors, one per unit-cell sublattice site), so
the identical `W_bulk` list can be fed to either function — reused
directly in `tests/test_vumps_apply_mpo.py`/
`examples/idmrg/vumps_apply_mpo/main.py` to cross-check both backends'
results against each other on the same operator at the exact `D=1`
field-polarized point. `W_bulk` is first grouped into a single
`Dw x Dw x d_g x d_g` grouped-supersite MPO tensor via `_group_automaton`
(the same routine `vumps_ground_state` already uses to group the
Hamiltonian automaton — grouping site-local MPO tensors into one is a
purely structural contraction, independent of whether the result feeds a
Hamiltonian or an arbitrary bounded operator), then `idmrg.grow_by_mpo`
grows `AL` by it on the trivial `n_uc=1` "periodic chain" VUMPS already
works at (the same reduction `imps_sum` uses for its own
`idmrg._canonicalize_periodic` call), which re-canonicalizes/truncates the
grown tensor exactly as `apply_mpo`/`imps_sum` already do. The resulting
truncated left-canonical tensor is completed to the full mixed gauge via
`_complete_mixed_gauge` — the identical bookkeeping step `imps_sum`'s own
construction needs, for the identical reason: neither operation is itself
a VUMPS ground-state solve, so there is no outer iteration left to fix up
an inconsistent gauge the way `vumps_ground_state`'s own iteration does.
Returns a `UniformMPS`, same as `imps_sum`. Scope restriction is identical
to `idmrg.apply_mpo`'s own (bounded/periodic operators only, the
Hamiltonian's own unbounded automaton is out of scope) — see that
function's own docstring. Verified against exact `D=1` closed-form cases,
a genuinely entangled `D=2,3` TFIM ground state (a unitary single-site
flip must preserve `eta`, flip `<Sz>`, and leave `<Sz(0)Sz(r)>` unchanged),
and a genuinely bond-growing `chi_W>1` two-site gate tiled once per an
`n_uc=2` unit cell — see `tests/test_vumps_apply_mpo.py` and
`examples/idmrg/vumps_apply_mpo/main.py`.

**Excited states: the tangent-space/quasiparticle excitation ansatz
(`pyitensor/idmrg_excitations.py`)**: an infinite chain's excitations
form a momentum-resolved band `E(k)`, not a single discrete state a
finite chain has, so `Infinite_Many_Body_Chain.excitation_energies
(k, n=1)`/`excitation_gap(ks=None)` implement the standard single-mode/
quasiparticle ansatz (Haegeman et al.), built on top of a **VUMPS**
ground state's own mixed-gauge `{AL, AR, C, GL, GR}` (requires
`gs_method="vumps"`, the default -- the opposite requirement from
`local_excitation_gap`/`td_dynamical_correlator`, which need
`gs_method="idmrg"` specifically; `vev`/`correlator` work under either
— see above): a tangent-space vector
`|Phi_k(X)> = sum_n e^{ikn} |...A_L A_L B_n(X) A_R A_R...>` with one
excitation tensor `B(X)` inserted at every unit-cell position, weighted
by momentum `k`. Multi-site unit cells are handled entirely inside
`vumps.py` (which groups the `n_uc` sites into one effective supersite
before ever building `AL`/`AR`/`C`), so this module only ever sees an
already-grouped, single-supersite `VUMPSResult`; every bond term must
have reach exactly 1 after grouping (checked once, inside
`vumps.vumps_ground_state` itself, `idmrg_excitations._check_reach_one`).
This mirrors MPSKit.jl's own architecture (`QuasiparticleAnsatz`,
`src/algorithms/excitation/quasiparticleexcitation.jl` /
`src/environments/qp_envs.jl`), not just its published equations: a
null-space gauge-fixing tensor `V_L` (`B(X) = reshape(V_L @ X)`,
`V_L^dagger @ AL = 0`); the ordinary, momentum-independent
channel-resolved background environments `GL_full`/`GR_full` (built from
VUMPS's own `GL`/`GR` plus the one-more-site bond content
`vumps._precompute_bond_environments` already computes); the
momentum-dependent, channel-resolved `GBL(k)`/`GBR(k)` ("the excitation
has already happened somewhere to the left/right", MPSKit's own
`lBs`/`rBs`) — each solved as a single self-consistent linear fixed-point
system per momentum (`_build_GBL`/`_build_GBR`), rather than a
hand-summed geometric series; and a 3-term `H_eff(k)` (`_h_eff_action`)
built by reusing one generic contraction helper
(`_full_channel_contraction`) with three different environment triples.
No metric reduction is needed for the resulting eigenproblem (unlike an
earlier, uniform-gauge version of this module, see below) — the
mixed-gauge tangent-space norm is already the trivial Euclidean one.

**How that eigenproblem is actually solved.** `H_eff(k)` has dimension
`dim = D*D*(d_g-1)`, and `excitation_energies` picks between two solvers
by size (`_DENSE_EIG_MAX`): assemble the matrix one basis vector at a time
and call `eigh` (`_lowest_dense`), or run Lanczos directly on
`_h_eff_action` (`_lowest_iterative`, `scipy.sparse.linalg.eigsh` over a
`LinearOperator`). The dense route costs `dim` applications of
`_h_eff_action`; the iterative one costs a few hundred regardless of
`dim` (measured: ~155-225 per momentum for `n=2` on a D=12 chain, not the
"few tens" a first estimate suggests -- which is why the threshold is where
it is rather than at 0).
Three things keep the dense route load-bearing rather than vestigial:
`dim` can legitimately be `1` (a `D=1` spin-1/2 chain), where ARPACK
cannot be used at all since it needs `nev < dim`; it is genuinely faster
below the threshold; and it is what the iterative route falls back to when
its eigenpairs fail a residual check, the same "never less robust than the
path it replaces" contract `_solve_linear_map` already has. The iterative
route uses a fixed, deterministic start vector rather than ARPACK's random
default, so a near-degenerate or gapless dispersion does not come out
differently run to run.

The bigger cost this sits on top of is `_channel_resolvent`, which
`_build_GBL`/`_build_GBR` need once each per momentum and direction.
Nothing in it depends on the excitation tensor `B` — only on `k` and the
momentum-independent environment — yet it used to be rebuilt inside
*every* `_h_eff_action` call: a dense `(D*D, D*D)` build plus a fresh
`O(D^6)` solve, twice per application, `dim` applications per momentum.
It is now built once per momentum and cached on the environment
(`_resolvents_for`, `ExcitationEnvironment.resolvent_cache`), LU-factored
once (`scipy.linalg.lu_factor`) so every later solve is two triangular
solves. Because the factorization is now amortized over many solves rather
than paid per solve, this resolvent keeps its own dense-vs-iterative
threshold (`_RESOLVENT_DENSE_MAX`) set by *memory* rather than by flops,
far above `_solve_linear_map`'s own `_DENSE_SOLVE_MAX` — the two solve
genuinely different problems, one-shot versus reusable.

`excitation_energies(env, k, n, return_vectors=True)` additionally returns
the tangent-space parameters `X` (the excitation tensor itself being
`B = (V_L @ X).reshape(D, d_g, D)`), which the energies alone previously
discarded. `spectral_weights` is what consumes them. This is
`pyitensor`-only: the
`itensor_version=3` port (`Chain::vumps_build_h_eff_dense`) still assembles
`H_eff(k)` densely and rebuilds its resolvents per application, so the two
backends now differ in cost (not in results — the cross-checks in
`tests/test_vumps_excitations_v3.py` are unchanged and still pass).

**Spectral weights (`spectral_weights`, `_spectral_source_vector`)**:
with the eigenvectors in hand, each branch's exact delta-peak residue
`|<k,a|O(k)|Psi>|^2` follows from one linear functional
`X -> <Phi_k(V_L X)|O(k)|Psi> = trace(X^dagger @ v)`, and
`_spectral_source_vector` builds that `(Dx, D)` matrix `v` once per
(momentum, operator). `Infinite_Many_Body_Chain.spectral_weights` /
`dynamical_structure_factor` are the public surface, gated exactly like
`excitation_energies` except that `itensor_version=3` is excluded too
(`Chain::vumps_excitation_energies` returns energies only).

Putting the ket in mixed canonical form with its center where the bra's
excitation tensor sits leaves only three regions — operator on the
excitation, operator somewhere to its right through an `(AR,AR)` chain,
operator somewhere to its left through an `(AL,AL)` chain — with both
semi-infinite backgrounds closing to the identity by canonicality. The
two chains are geometric series in `e^{ik} E_RR` / `e^{-ik} E_LL`.

Two implementation decisions are worth knowing, because both are about
cost or about the class of bug this module's history is made of:

- **The solves are transposed.** Written forwards, each series takes a
  `B`-dependent source, so keeping `B`'s legs open (which is what
  producing `v` rather than one scalar means) would need `D^2*d_g` solves
  per momentum. `idmrg._apply_transfer` and
  `idmrg._apply_transfer_from_left` are exact transposes of each other
  under the bilinear pairing with the *same* `E`, so transposing each
  solve moves it to the operator end of the contraction: **two** solves
  per momentum, independent of `D`, `d_g`, the operator, and the number
  of branches requested. They are cached per momentum
  (`_spectral_resolvents_for`) alongside the H_eff ones.
- **`_channel_resolvent` is deliberately not reused.** Its deflation
  projector hardcodes the `conj(r_mixed)` dual convention validated for
  the *mixed* `E_RL`/`E_LR` transfers; the transfers here are the
  ordinary `(AR,AR)`/`(AL,AL)` ones, whose duals are exactly the identity
  and whose non-trivial fixed points are `C^dagger C` / `C C^dagger` in
  closed form (no eigensolve). `_spectral_resolvent` writes that
  projector explicitly in the bilinear pairing, and — unlike
  `_channel_resolvent`, which can only deflate at the singular momentum —
  applies it at *every* momentum: what deflation adds to the solution is
  a multiple of the identity, which contracts into a multiple of `AC`,
  which `V_L^dagger` annihilates exactly. The returned weight is
  unchanged and the conditioning is bounded uniformly in `k`.

That same annihilation is why no `<O>` is subtracted anywhere: the
disconnected part of the correlator lives in exactly that direction, so
the `k=0` answer comes out connected on its own. It also makes the
branch-complete total (`return_total=True`) exactly the connected static
structure factor — a one-site operator applied to a uniform MPS lands
exactly inside the tangent space, so nothing of `O(k)|Psi>` is outside
the span the branches diagonalize — which is the sum rule
`tests/test_infinite_chain_spectral.py` checks against an independent
real-space `correlator` sum, alongside the f-sum rule and the exact
product-state limit.

**History.** An earlier version of this module built the ansatz on top
of a single uniform-gauge tensor from `idmrg.py`'s growing algorithm
instead, with a hand-derived, closed-form-resummed diagram list; it
matched the exact free-fermion dispersion of a field-polarized XX chain
(`D=1`) to $\sim$14 digits, but was wrong (anomalously flat) for a
genuinely entangled (`D>1`) ground state despite eight independent
investigation passes (three real bugs found and fixed along the way,
none of which closed the gap; see `git log` on this file and
`docs/idmrg_excitation_mpskit_port_plan.md` for the full account). That
investigation's own conclusion — rewrite from scratch mirroring
MPSKit.jl's actual architecture on top of a genuine mixed-gauge
`{AL,AR,C}` ground state (VUMPS), rather than patching the diagram list
further — is what the current version implements. Getting it working at
`D>1` additionally required fixing a real, pre-existing sign bug in
`vumps.py` itself (`_solve_left_environment`/`_solve_right_environment`/
`_energy_density_and_source_from_left`/`_right` were all missing a
conjugate on the dominant fixed point used to close an
`apply_transfer_from_left`/`apply_transfer` output — invisible at `D=1`,
where that fixed point is always real; see `vumps.py`'s own docstrings)
and a subtler point specific to the excitation ansatz itself: the
constant subtracted from `H_eff(k)`'s raw eigenvalues must be `H_AC`'s
own Rayleigh quotient on the converged `AC=AL@C` (`lambda_AC`), not the
ground state's physical energy density `e_cell` (`=lambda_AC-lambda_C`,
the difference of `H_AC`'s and `H_C`'s own eigenvalues — the standard
VUMPS energy-density identity) — see `ExcitationEnvironment`'s own
docstring. With both fixes, the `D=2` transverse-field Ising dispersion
matches an independently-converged MPSKit.jl result to 6 significant
figures at every momentum tested, and `H_eff(k)` is Hermitian to machine
precision at `D=1`, `2` and `3`. See
`examples/idmrg/excitation_gap_xx/main.py` (`D=1`) and
`examples/idmrg/excitation_gap_tfim/main.py` (`D=2`).

**A cheaper, cruder alternative: the local superblock gap
(`pyitensor/idmrg.py::local_excitation_gap`)**: rather than a fresh
tangent-space construction, this reuses the growing algorithm's own final
micro-step exactly as it already runs — `idmrg_ground_state` now stashes
the last micro-step's local-solve ingredients (`HL`/`HR` environments, the
two MPO tensors, the two physical Indices, the converged local ground
vector/energy) on `IDMRGResult.local_superblock` — and re-diagonalizes
that *same* 2-site effective Hamiltonian (rebuilding its matvec via
`kernels.make_matvec` on the stored pieces) for its second-lowest
eigenvalue instead of only its ground state, returning the difference.
Framed explicitly as the infinite-chain analogue of finite DMRG's own
Lagrange-multiplier/overlap-penalty excited-state method (`dmrg.py`'s
`dmrg_excited`/`Chain.excited_states`): that method needs a penalty term
threaded through a whole re-sweep because it enforces orthogonality
against a *separate*, externally-held ground-state MPS one local update at
a time; here the "state to stay orthogonal to" is just the local ground
vector already found, in the very same local Hilbert space, by the very
same solve, so the constraint is enforced *exactly* via deflation
(`P = I - |psi0><psi0|`; `deflated_matvec(v) = P(matvec(P(v)))`, Hermitian
since both `P` and the underlying local Hamiltonian are) rather than
approximately via a finite penalty weight — mathematically, a Hermitian
operator's constrained stationary points (extremize `<psi|H|psi>` subject
to `<psi|psi>=1`, `<psi|psi_0>=0`) are exactly its *other* eigenvectors, so
this is what a penalty method converges to anyway as its weight -> infinity,
obtained directly with nothing to tune. Implementation mirrors
`_local_two_site_solve`'s own small-dimension fallback (dense `eigh` for
`dim<=3`) and otherwise reuses `dmrg._lanczos_ground_state` (imported
already) on the deflated matvec, rather than introducing a second Lanczos
implementation or a new `scipy.sparse.linalg.eigsh` dependency.

Deliberately positioned as a *cruder, opt-in* alternative, not a
replacement: it has no momentum label at all (a single scalar, not a
dispersion), does not require `D=1`, and — critically — never lets `HL`/
`HR` relax for whatever the second local eigenstate actually represents,
since they are exactly the ground state's own converged environments.
Measured directly against the two cases the tangent-space ansatz was
itself validated against: for the field-polarized XX chain (`D=1`, exact
answer known), it comes out $\sim$10% too high (5.5 vs.\ the exact 5.0)
— confirming it is a genuinely cruder approximation, not a bug; for a
genuinely entangled (`D>1`) dimerized Heisenberg chain — exactly the case
`excitation_gap` cannot handle — it lands within $\sim$0.5% of a large
finite open chain's own ED gap (`Spin_Chain.get_gap(mode="ED")` at
`n_sites=12,14,16`, the same finite-size-extrapolation cross-check style
`test_n_uc2_dimerized_chain_matches_finite_size_extrapolation` already
uses for the ground-state energy density), a reassuring result precisely
where the tangent-space ansatz offers nothing at all. See
`examples/idmrg/local_excitation_gap/main.py` for the worked comparison.

**Sequential multi-site VUMPS** (`pyitensor/vumps_ms.py`). The unit cell may
have any number of sites. `vumps.py`'s original support for a multi-site
cell was `_group_automaton`: fold the `n_uc` sites into one supersite of
dimension `prod(d_p)` and run the single-site algorithm, which is exact but
exponential in the cell size and is why it rejected `n_uc > 2`. The
literature names that as the trap rather than an implementation detail
(arXiv:2003.01142: "the cost of a naive application of the VUMPS algorithm
would scale exponentially with the size of the unit cell ... a computational
effort that scales linearly rather than exponentially"), so `n_uc > 2` now
routes to the sequential algorithm instead: per-site `AL[n]`/`AR[n]`/`C[n]`,
one `H_AC[n]` eigenproblem per site and one `H_C[n]` per bond, swept across
the cell.

Two things generalize, not one. Beyond the site count, the environments
become fully channel-resolved -- `GL[n]`, `GR[n]` are `(Dw, D, D)`, one
matrix per MPO channel -- because `vumps.py`'s reach-1 specialization (one
accumulated `GL`, one `GR`, a list of one-site-away bond channels, guarded
by `_check_reach_one`) is exactly what grouping bought: any coupling inside
the cell is on-site once the cell is one supersite, and ungrouped it is not.
The channel-resolved form subsumes the reach-1 one; reading `vumps.py`'s own
`_h_c_action` in that language shows its `_cap_right(C,GR)` is the `S`
channel, `_cap_left(GL,C)` is `F`, and its bond terms are the pending
channels.

`n_uc <= 2` deliberately stays on the grouped path so its validated values
do not move, and because `idmrg_excitations`' tangent-space ansatz consumes
the grouped mixed gauge (`excitation_energies` raises for `n_uc > 2` rather
than reading a list of per-site tensors as one tensor). Cross-checks in
`tests/test_vumps_multisite.py`: machine-precision agreement with the
grouped path at `n_uc = 1, 2`; exact TFIM energy and magnetization; a
trimerized free-fermion chain against its own 3-band integral with genuinely
non-uniform per-site occupations; and cell-size invariance, including the
same 3-periodic chain written on a 6-site cell.

**Product-state traps in the growing algorithm, and the noise that breaks
them** (`pyitensor/idmrg.py::_noise_perturbed_split`,
`Chain::idmrg_noisy_isometry`). A particle-number-conserving Hamiltonian has
the vacuum as an exact eigenstate, hence an absorbing fixed point of the
growth loop: rank-1 state, warm-started solve returns it unchanged, no
mechanism to regrow `chi`, and `etol` trips immediately because consecutive
densities are identically equal. Both backends collapsed identically on a
dimerized spinless-fermion chain at uniform `mu=+1.25` (6 runs of 6),
returning `e = n = 0` where the band integral gives `e/site = -0.01648450`
at `n/site = 0.166092` — a shared algorithmic gap, not a porting bug.

The cure has two halves and the second is the load-bearing one here: White's
density-matrix perturbation built from the individual MPO *channels* (using
`H|theta>` would be useless — it is parallel to `|theta>` at an eigenstate),
plus a random admixture of weight `O(noise)` on the local solve's start
vector, because a product state stays an exact eigenvector of the local
effective Hamiltonian in any basis however enriched, and a Krylov solve
started on an exact eigenvector breaks down immediately. Measured: the
density-matrix term alone grew `chi` off 1 and left the energy identically 0
for every macro-iteration.

The schedule is demand-driven rather than a fixed iteration count, for two
measured reasons: escape depends on duration not magnitude (`1e-4` over 12
iterations escapes, `0.1` over 10 does not), and a fixed schedule perturbed
healthy runs enough to regress a correlator test by 1.3e-6 against a 1e-6
tolerance. Noise therefore arms only while the noise-free reduced state has
purity 1 (i.e. is a product state), plus a short tail; an entangled model
never arms it. `noise_iters` caps the total so a legitimately-product ground
state stops re-arming, and `etol` is suppressed while noise is active so a
perturbed state is never reported as converged.

**The deflated excited solve must shift, not merely project**
(`pyitensor/idmrg.py::_lowest_two_eigenvalues`, and its C++ port inside
`Chain::idmrg_local_excitation_gap`). The obvious implementation — build
`P = I - |psi0><psi0|` and hand `P H P` to a smallest-eigenvalue solver —
is wrong. `P H P` agrees with `H` on `psi0`'s orthogonal complement, but it
also carries `psi0` itself at eigenvalue *exactly zero*. That extra
eigenvalue is invisible while the rest of the spectrum lies below zero, and
fatal the moment the stored operator's own ground eigenvalue is positive:
zero is then `P H P`'s smallest eigenvalue and precisely what the solver
returns, so the reported gap is `0 - e0` — minus the stored superblock
energy. Reported from an external Kondo-chain session at `maxm=16` as
−358.003 meV against a stored energy of +0.358002587681, digit for digit,
and reproduced in isolation on a random Hermitian matrix shifted to a
positive ground energy (`tests/test_idmrg_local_gap_deflation.py`), where
the deflated solve returns 0.0 with `|<psi0|v>|**2 = 1`, i.e. `psi0` itself.
Deflation excludes `psi0` from the Krylov space in exact arithmetic only;
rounding regrows the component every iteration and a restarted solver locks
onto it. Which models this hits is not under anyone's control — the per-site
baseline (`_subtract_energy_baseline`) leaves the stored energy as a small
residual boundary term of either sign, so the same model was correct at
`maxm=8` (energy −1.05) and wrong at `maxm=16` (+0.358). Both backends now
diagonalize `P H P + sigma |psi0><psi0|` with `sigma` far above the bottom
of the spectrum. Two safety nets sit on top, both no-ops normally: if the
solve still returns something sitting on `psi0`, `sigma` is raised and the
solve retried; if it returns something *below* the accepted ground
eigenvalue, that lower state is promoted and the deflation redone.

**Both eigenvalues come from the stored superblock, at the same solver
strength.** The ground one is not read back from the growing algorithm even
though its stored value is by construction the Ritz value of the stored
local ground vector: that loop warm-starts every micro-step from McCulloch's
prediction and stops on a *relative* Ritz-value criterion, which bounds the
distance to *some* eigenvalue rather than to the lowest one, so mixing its
number with a strong random-started solve is not a difference of two
eigenvalues of one operator. The detail entry point
(`Chain::idmrg_local_excitation_gap_detail`, `detail=True` on the Python
side) hands back both ground eigenvalues so `infinitechain.py`'s
`_warn_if_growth_missed_local_ground_state` can raise a `RuntimeWarning`
when they disagree.

**The local superblock gap is not a fixed-charge quantity.** Neither backend
conserves QNs (`ConserveQNs=false` in `mpscpp3/get_sites.h`, and `pyitensor`
likewise), so the 2-site local spectrum contains particle-number-changing
excitations against the frozen environment. Measured exactly, by dense
diagonalization of the stored superblock on a gapped SSH spinless-fermion
chain at `maxm=16`: at `mu=0` the second and third eigenvalues are exactly
degenerate (the ±1-electron pair) and the gap is 0.631321; adding `mu` to
every on-site energy splits them and the gap becomes 0.431325 / 0.431324 at
`mu=±0.2` — exactly `|mu|` lower — while the converged state, its energy
density per the shift, and its correlators to 8 decimals are unchanged, and
the stored ground vector remains the exact ground state (overlap 1.0000)
throughout. So "add a constant to every on-site term and assert the gap is
unchanged" is *not* a valid invariance test for this estimator, however
reasonable it looks: a correct implementation fails it. It is a real
limitation of the quantity, not of the code.

**Tightening the local superblock gap further: `window=`
(`pyitensor/idmrg.py::local_excitation_gap_windowed`)**: rather than adding
Krylov vectors to the same frozen 2-site block — `local_excitation_gap`'s
own `dim>3` branch already diagonalizes that to Krylov convergence, so more
iterations there do not change the answer — this grows the local
block itself by `window` extra *free* physical sites on each side, using
fresh copies of the periodic per-sublattice MPO tensor (`_build_automaton`)
threaded onto `HL`/`HR`'s own mpo axis via the same `_relabel_pos`/
`_fresh_physical_copy` machinery `idmrg_ground_state`'s own growing loop
uses (recovering `HL`/`HR`'s mpo-axis Index by elimination from their own
3 legs, since `local_superblock` does not store it directly). Both the
ground state and the deflated first excited state are then re-solved fresh
within this larger block via `_lanczos_ground_state`, rather than reusing
the growing algorithm's own ground vector — `window=0` reduces to exactly
the same effective Hamiltonian `local_excitation_gap` diagonalizes (used as
an internal consistency check on the construction, confirmed to agree to
Lanczos precision). Costs `d**(2*window)` more local Hilbert space
dimension per extra site pair (`d` = the physical dimension), and is
supported only for `n_uc=1` — widening needs to know which sublattice
position each extra inserted site takes, which is not tracked for `n_uc=2`.

Measured directly, this is a real, converging refinement, not just noise:
on the field-polarized XX chain (`D=1`, exact answer known) the error
drops monotonically from 10\% at `window=0` to 3.8\%/2.0\%/0.81\%/0.44\% at
`window=1/2/4/6`; on a gapped, genuinely entangled (`D>1`)
transverse-field Ising chain, `window=3` (an 8-site local block, gap
2.047) matches an 18-site open finite chain's own ED gap (2.049) to
$<$1\%, converging at least as fast as growing the finite chain itself
does (`n`=10/12/14/16/18 give ED gaps 2.133/2.099/2.076/2.060/2.049).
Reproducible to $\sim$1e-9 across repeated calls despite the randomized
Lanczos start, exactly like `local_excitation_gap` itself. The convergence
*rate*, however, depends on the model's correlation length, not just its
gap: on the `S=1` Heisenberg chain (the Haldane gap, correlation length
$\sim$6 sites), `window=0/1/2` only move the estimate from
$\sim$29\%/24\%/21\% too high — still improving, but far more slowly,
since a handful of extra sites barely dents a correlation length that
long, and the larger physical dimension (`d=3` vs.\ `d=2`) makes each
extra site pair 9$\times$ more expensive rather than 4$\times$, making
large `window` impractical there. See
`examples/idmrg/local_excitation_gap/main.py` for the worked comparison
across models.

**Dynamical correlators, via a finite-window reduction to the existing
finite-chain KPM stack (`infinitechain.py`, not `pyitensor/idmrg.py`)**:
unlike the static-correlator machinery above, `Infinite_Many_Body_Chain.
kpm_finite` deliberately does **not** add any new Chebyshev-recursion
code to `idmrg.py`. Named `kpm_finite` (not `get_dynamical_correlator`,
the finite-chain method it wraps) so the name itself flags that this is
a finite-window approximation, not the exact infinite-size method
described further below. `H` is extensive/unbounded in the
thermodynamic limit, so there is no infinite-chain analogue of the
transfer-matrix formalism `vev`/`correlator` use — a literal KPM
expansion of the full infinite `H` has no meaning (the same reason
`apply_mpo` restricts `W_bulk` to *bounded* periodic operators, see
above). Instead, `kpm_finite` builds an ordinary finite, open-boundary
`Many_Body_Chain` (`itensor_version="python"`, `n_window` repeats of
the unit cell) with a Hamiltonian tiled from `self._h_intra`/
`self._h_inter` (`_window_hamiltonian`/`_shift_terms`, new module-level
helpers in `infinitechain.py` — pure term-list arithmetic, reusing the
same `MultiOperator.op` 0-based-site-index format `_canonicalize_
hamiltonian` already produces, just concatenated across cells with an
increasing site-index shift), places the two operators at the window's
central cell, and calls `kpmdmrg.get_dynamical_correlator` on that
temporary chain directly — the exact same code path (and, transitively,
`pyitensor/chain.py`'s `Chain.kpm_dynamical_correlator`/`_kpm_moments_
full`/`_scaled_hamiltonian`) an ordinary finite `Spin_Chain`/
`Fermionic_Chain` already uses, completely unmodified. `kpmdmrg.
get_dynamical_correlator` is called directly rather than through the
usual public `Many_Body_Chain.get_dynamical_correlator`/`dynamics.
get_dynamical_correlator` dispatch, to sidestep `dynamics.py`'s own
`is_hermitian()` gate — confirmed (this is `set_hamiltonian`'s own
already-documented gap, not a new finding) that `is_hermitian()`'s
`simplify()` step false-rejects an ordinary cross-site
`Sx[i]*Sx[j]+Sy[i]*Sy[j]+Sz[i]*Sz[j]`-style term as non-Hermitian, which
would otherwise misroute a perfectly ordinary Heisenberg-type window
Hamiltonian to the non-Hermitian KPM/dynamics path.

Correctness of the term-tiling itself (as opposed to the pre-existing,
independently-tested KPM math it feeds into) is checked in `tests/
test_infinite_chain.py` by an *exact* ED ground-state-energy comparison
between `_window_hamiltonian`'s output and an independently, by-hand-built
equivalent finite chain — the same finite-size-extrapolation pattern
`test_n_uc2_dimerized_chain_matches_finite_size_extrapolation` already
uses for `gs_energy`, just applied to the Hamiltonian construction
directly instead of to iDMRG's own converged energy density.

This is a genuine finite-size **approximation**, not an exact
infinite-size dynamical-correlator method — see `kpm_finite`'s own
docstring (and `docs/user_guide.md`'s corresponding section) for the
light-cone/moment-count argument for why a fixed `delta` needs
`n_window` chosen large enough (and, past a point, `delta` itself
coarsened) to avoid open-boundary reflections contaminating the result,
and `examples/idmrg/dynamical_correlator_finite_window/main.py` for a
worked `n_window`-convergence sweep.

**A genuinely infinite-chain KPM is possible in principle, but not
implemented here.** Two separate ingredients would be needed, both
reusing machinery this module already has, not a small patch on top of
`kpm_finite`:

1. Cap the window's two ends with the *converged environment* — the
   `HL`/`HR` accumulated Hamiltonian blocks `idmrg_ground_state` already
   builds during growth, or an equivalent construction from the
   fixed-point machinery `apply_mpo`/`imps_overlap` already use —
   instead of plain open boundaries: infinite boundary conditions (IBC,
   Phien/McCulloch/Vidal, PRB 86, 245107 (2012)). This removes a real
   error source `n_window` alone cannot fix no matter how large: an open
   chain's own ground state carries boundary artifacts (e.g.
   Friedel-oscillation-like features) that contaminate even the
   *central* region, not just the two edges, because the window's own
   DMRG ground-state search has no way to know it should match the true
   infinite bulk rather than terminate at a physical edge.
2. Independently, let the *active* window grow by roughly one site per
   Chebyshev moment as the recursion proceeds (reusing the converged
   unit-cell tensors to extend the state on demand, similar in spirit to
   a "growing window" real-time TEBD/TDVP simulation on an infinite
   chain), rather than fixing a static `n_window` up front. This removes
   the `n_window`-vs-moment-count constraint entirely, bounding cost
   only by the usual bond-dimension/computational budget of any DMRG
   calculation, not by a fixed geometric window guessed in advance.

Assembling these into a working IBC-window KPM implementation is a new,
nontrivial feature in its own right — flagged as a known, deliberate
follow-up rather than attempted here (see `kpm_finite`'s own docstring
for the same note, addressed to a caller rather than a maintainer).

**Infinite-boundary-condition (IBC) window construction and real-time
dynamical correlator (`pyitensor/idmrg_window.py`)**: following
Milsted/Vanderstraeten et al., "Infinite boundary conditions for response
functions and limit cycles in iDMRG" (arXiv:1804.09163), Sec. V.1 — the
item-1 ingredient flagged just above as unimplemented future work. The
paper's method time-evolves a finite window of tensors
around a local perturbation, with everything outside the window frozen
at the uniform ground state and capped by "the Hamiltonian terms outside
the window, projected onto the reduced D-dimensional Hilbert space of the
left/right block" — which turns out to be *exactly* what
`idmrg_ground_state`'s own growth loop already computes as `HL`/`HR`,
just discarded once converged. `idmrg.py` was extended (`IDMRGResult`'s
own `env_HL`/`env_HL_bra`/`env_HL_ket`/`env_HL_mpo` and mirrored `env_HR_*`
fields, plus `W_bulk`) to snapshot these — specifically, `HL`/`HR` as
they stood *entering* the last executed macro-iteration (captured via a
per-macro-iteration `env_window_boundary` snapshot, overwritten every
iteration so whatever remains at loop exit is the one consistent with the
cell the window tiles: the snapshot is taken *entering* the very micro-step
whose `theta` `_theta_cell` is built from, so that `theta`'s own two outer
legs literally are this snapshot's `HL_ket` and `HR_ket`) — no new
environment solve needed, only exposing what already exists.

`idmrg_window.py`'s `build_window(result, n_window)` tiles periodic repeats
of `IDMRGResult.cell_raw` (the window's MPS) and `W_bulk` (its MPO)
into an ordinary `mpscontainer.MPS`/`MPO` pair (`_tile_periodic`, fresh
Link Indices at every copy boundary to avoid the identity-collision
hazard `_relabel_pos`/`_project_channel` already warn about elsewhere in
this module), with the very first/last tensors' outer legs left as `env_
HL_ket`/`env_HL_mpo`/`env_HR_ket`/`env_HR_mpo` so they attach directly to
the environment caps. Environment extension *through* the window
deliberately reuses `idmrg.py`'s own explicit-Index `_extend_HL`/
`_extend_HR` (`_extend_through_left`/`_extend_through_right`), not
`dmrg.py`'s generic `_Chain`/`_link_at`-based `_extend_left`/
`_extend_right`: the latter infers a site's left/right Link purely by
looking for a same-Index neighbor *within the chain itself*, correct for
an ordinary finite chain (whose boundary sites truly have no further leg,
see `mpscontainer.py`'s own docstring) but wrong here, since the window's
own edge sites *do* carry one more leg each (connecting to `env_HL`/
`env_HR`, tensors living outside the `_Chain` object entirely) that
`_link_at` cannot see — confirmed directly: using the generic functions
silently mis-contracted the bra/ket legs at the window's own edge sites
(auto-contracting them together instead of leaving them properly
dangling for a fresh bra-side mint), since `_link_at`'s neighbor lookup
returns "no such leg" for a boundary site regardless of whether the
tensor itself actually has one.

**A confirmed, non-obvious normalization finding**: `window_total_energy`
(the capped window's total energy, via `_extend_through_left` then
closing against `env_HR`) is *not* simply `e0 * (site count)` — with both
`env_HL_ket` and `env_HR_ket` left as free/dangling legs (rather than
contracted against a specific boundary vector), `<window|window>` is not 1
(confirmed numerically against an independent Lanczos-Rayleigh-quotient
computation using the identical `env_HL`/`env_HR`/`W_pL`/`W_pR`, to machine
precision) — `window_total_
energy` divides by it. Even normalized, the *absolute* total still
carries a large, `n_window`-independent baseline from whatever
`env_HL`/`env_HR` themselves represent (every macro-iteration already
absorbed before the snapshot) — so `window_energy_density(result,
n_window)`, a finite difference between window sizes `n_window` and
`n_window+1` divided by the number of physical sites actually added
(mirroring `idmrg_ground_state`'s own
`density = (energy-prev_energy)/(2*n_uc)` diagnostic), is the well-posed
sanity check this module validates against, not `window_total_energy`
directly.

**Three assumptions of exact left-canonicality, and why the window stopped
tiling `U_list`.** `window_total_energy` used to shortcut `<window|window>`
to `dim(env_HR_ket)`, which is what the chain collapses to when *every*
tiled tensor is an exact isometry. That held for `U_list` (svd()'s own `U`
factors) but not for `cell_raw`, whose second tensor
`S.V.lambda_o^{-1}` is an isometry only to ~1e-3 — and it failed
informatively, reporting an energy density *below* the true `e0`, which is
variationally impossible and is what identified the normalization rather
than the state as the culprit. `window_norm_squared` now contracts the
overlap in one pass, reproducing the shortcut to machine precision wherever
the shortcut was valid.

Two more places made the same assumption. `_close_array_chain` closed on
the left with a bare trace; `S(x=0,t=0)`, which must equal
`<Sz Sz> = 0.25` for spin-1/2, came out at -0.078. It now closes with the
transfer operator's own left fixed point at the chain's starting cell
position and calibrates against the same contraction with no operator
inserted, so a bare chain returns exactly 1. And `snapshot_correlator`/
`_padded_arrays` indexed bulk tensors by *sublattice* (`position % n_uc`);
they now index by *cell* position (`% n_cell`), while operator matrices are
still looked up by sublattice.

With all three fixed, `window_energy_density` reproduces `e0` to
1e-14..1e-9 tiling `cell_raw`, against 2.2e-4..9.5e-3 tiling `U_list` —
and the `U_list` error saturated in `n_window`, i.e. it was a genuine gauge
error rather than a window-size effect. One consequence: a window's tiling
unit is the *cell*, which is two unit cells when `n_uc=1`, so `n_window` is
rounded up to a whole number of cells and the realized size is reported
back as `IBCWindow.n_window` (`dynamical_correlator_td` takes its centre
from that, or the perturbation would sit off-centre). For `n_uc=2` nothing
changes. `build_window`'s old wraparound-dimension precondition is gone: it
existed because tiling `U_list` needed its two ends to agree, and the
cell's ends agree by construction. **A second, real bug this same finite-difference check caught**:
an earlier version compared window sizes `n_window` and `n_window+n_uc`
(not `n_window+1`) while still dividing by `n_uc` — `n_window` counts
*unit-cell copies*, so `+n_uc` copies adds `n_uc*n_uc` physical sites,
not `n_uc` — invisible for `n_uc=1` (where `n_uc*n_uc==n_uc`), but
confirmed directly to inflate every `n_uc=2` result by *exactly* a factor
of 2 for an easy, cleanly-converged test model (an exact ratio, not
convergence noise, which is what exposed it as a real off-by-`n_uc` bug
rather than an iDMRG convergence limitation).

**Real-time evolution, shifted overlaps, and `S(k,ω)`**: `window_tdvp_step`
two-site-TDVP-evolves a capped window in place (mirroring `tdvp.py`'s own
`tdvp_step`), but `tdvp.py`'s own sweep functions cannot be reused
directly — they build each bond's local effective Hamiltonian and
environments via `mpsalgebra.py`'s `_link_at`, which finds a site's
left/right Link purely by looking for a same-Index *neighbor in the chain
itself*; correct for an ordinary finite chain (whose boundary sites truly
have no further leg) but wrong for a window, whose edge sites *do* carry
one more leg each (connecting to `env_HL`/`env_HR`, tensors outside the
chain object entirely) that `_link_at` cannot see — confirmed directly,
passing a window straight through `tdvp.py`'s own sweep functions
silently mis-contracts the bra/ket legs at the window's own edge sites.
`idmrg_window.py`'s own `_window_two_site_heff`/`_window_one_site_heff`
read a site's left/right Link directly off its own tensor instead (always
present for a window, unlike an ordinary chain's boundary), and
`_all_left_environments_window`/`_all_right_environments_window` build
environments via idmrg.py's own explicit-Index `_extend_HL`/`_extend_HR`.
`apply_local_operator`/`local_expectation` add the perturbation and
readout primitives; `local_expectation` needed its own real fix, distinct
from the TDVP one: a plain observable's own boundary weighting is *not*
symmetric between the window's two ends — the *right* boundary needs the
dominant right fixed point (`_all_right_fixed_points`), exactly like
`onsite_expectation`/`two_point_correlator` already rely on; an earlier
version used a bare
trace on both ends and produced a visibly non-uniform profile for a
converged, translation-invariant ground state that should read exactly
uniform everywhere.

**A real, previously-undiscovered bug and its fix (`eshift`)**:
`window.env_HL`/`env_HR` are, like `window_total_energy`'s own docstring
documents, *not* energy-baseline-subtracted — they carry a large,
`n_window`-independent additive constant left over from every
macro-iteration `idmrg_ground_state` already absorbed before the
environment snapshot was taken, which varies run to run with iDMRG's own
unseeded convergence path even for the *same*, equally-converged physical
ground state. Since this constant multiplies the identity on the window's
own local Hilbert space at every bond (a standard property of a
canonically consistent DMRG/TDVP environment — `<θ|Heff|θ>` reproduces
the *total* energy of the whole capped window+environment system
regardless of which bond `θ` sits at), TDVP-evolving under the unshifted
Heff produced a *global*, unphysical, run-dependent phase
`exp(-i*const*t)` — confirmed directly: for the identical physical model,
`window_total_energy` differing between -33 and -96 across otherwise-
equivalent iDMRG runs left the resulting `S(x,t)` trajectory essentially
uncorrelated between runs. The fix (`window_tdvp_step`'s `eshift`
parameter, subtracted from every local effective Hamiltonian before
exponentiating) is not a novel design choice — it restores consistency
with the rest of this codebase's own established real-time correlator
convention: an ordinary finite chain's `mpscpp3::quench_tdvp` already
computes `EGS=<wf0|H|wf0>` and evolves under `Hshift=H-EGS*Id` (never raw
`H`), and `edtk/timedependent.py::evolution_DC` (the ED reference) does
the same. `dynamical_correlator_td` passes `eshift=window_total_energy
(window_B)` (measured once, on the unperturbed ground window, before any
perturbation) into every `window_tdvp_step` call. Verified two ways: (1)
a dense-matrix exact cross-check — a 2-site window's own (env_HL, window,
env_HR) system is small enough to diagonalize directly, giving a literal
`exp(-i(H-EGS)t)` to compare `window_tdvp_step`'s own result against, no
approximation involved — matches to ~1e-6 (Krylov/SVD-truncation-limited,
not a real discrepancy); (2) two independent, well-converged iDMRG solves
of the same physical model now give closely-agreeing `S(x,t)`, not just a
coincidentally-agreeing `e0` (`tests/test_idmrg_window_free_fermion.py`).
The native ITensor v3 port (below) had the identical bug. It was fixed
there first with a post-hoc `exp(+i*eshift*t)` correction (ITensorTDVP's
own `tdvp()` is a vendored black box that cannot be wrapped at the source
like pyitensor's own matvec), and then, with the 2026-08-29 gauge fix,
replaced by the same vacuum normalization described below: once the
snapshot closes the window's boundary legs with transfer-matrix fixed
points rather than a bare trace, an `eshift` measured with those legs
*traced* is no longer the right constant to divide out.

`dynamical_correlator_td`/`snapshot_correlator` reconstruct
`S(x,t)=<ψ|A_x e^{-i(H-EGS)t}B_0|ψ>` for *every* `x` from a *single* window
evolution (the paper's own headline efficiency result), by inserting
operator `A` directly into the *bra* side of the overlap at the shifted
absolute position, padding whichever side's explicit tensor list falls
short with unevolved `cell_raw` copies (`_padded_arrays`) — valid because,
outside each window's own causal cone, its tensors already equal the
converged cell exactly. This is a simplification of the paper's own Eq. 7 to
`t1=0` (the naive Eq. 3 form) rather than the full two-branch trick
(evolving a *second* window backward by `t1` too, doubling the accessible
total time for the same TDVP cost) — `window_tdvp_step` already supports
backward evolution via a negative `dt`, so only the second branch and its
overlap bookkeeping are missing, a documented follow-up. `connected=True`
(the default) subtracts the disconnected background `<A><B>`
(`onsite_expectation`) — confirmed directly to matter: the *raw*
correlator approaches `<A><B>` (not 0) at large `|x|`, and summing that
non-decaying background over `x` with `e^{-ikx}` produces a spurious,
dominant `k=0` contribution with no discernible dispersion.
`dynamical_correlator_komega` does the spatial DFT
(`S(k,t)=Σ_x e^{-ikx}S(x,t)`) then reuses `timedependent.py`'s own
`_fourier_transform_correlator` unchanged (it was explicitly factored out
for other time-domain submodes to reuse, see `tdz.py`'s own use of it),
so `delta` means the same Lorentzian-broadening thing here as in every
other dynamical-correlator submode in this codebase.

Validated: `S(x,t=0)` matches idmrg.py's own exact static
`two_point_correlator` to machine precision (1e-16) for every `x` tested
(`tests/test_idmrg_window.py`); `S(k,ω)` (with `connected=True`) shows a
clear, k-dependent peak instead of a featureless background; cross-checked
against `kpm_finite` (an independent approximation scheme — open-boundary
window + KPM, vs. this method's IBC window + TDVP) for the same physical
correlator — an exact match isn't expected, but both agree on where the
dominant spectral weight sits (`tests/test_infinite_chain.py`,
`examples/idmrg/td_dynamical_correlator/main.py`). Also cross-checked at
`t>0` against an *exact* non-interacting reference: the XX spin chain
(`H=J*Σ(Sx_iSx_{i+1}+Sy_iSy_{i+1})`) maps exactly to free fermions under
Jordan-Wigner (`tests/_free_fermion_reference.py`), whose connected Sz-Sz
correlator has a closed form computable via one N×N diagonalization —
polynomial cost, so exact for N far beyond many-body ED's exponential
reach (`tests/test_idmrg_window_free_fermion.py`). The reference itself
is validated against exact ED to ~1e-12; a *dimerized* (SSH-like, still
exactly free-fermion) variant is used for the quantitative iDMRG
comparison since it opens a gap and converges far better than the
gapless uniform chain (state_overlap ~0.98-0.999 vs ~0.87-0.96 at the
same bond dimension) — this is the check that found the `eshift` bug
above.

**Fermionic operators (2026-08-29)**: a parity-odd pair on this path
needs the Jordan-Wigner string threaded on both sides, and did not have
it. `apply_local_operator` now dresses the ket with `F` on every window
site left of the perturbation *before* the evolution (so the string is
evolved along with it), `snapshot_correlator` dresses the bra the same
way up to its own measurement site, and both truncate at the same left
reference — beyond it the two strings are identical and cancel (`F²=Id`),
which is what makes the `t=0` reduction to the finite between-the-
endpoints string of `two_point_correlator` exact. Truncating there is the
same boundary approximation the IBC window already makes (the evolution
sees `P_L H P_L`, differing from `H` only in terms crossing the window's
left edge). `snapshot_correlator` also had a latent convention bug the
same work exposed: `_close_array_chain` conjugates the bra, so applying
`A` there computes `<ψ|A†…>`, not the documented `<ψ|A…>` — invisible for
every Hermitian name (only `Sz` appeared in the tests and examples),
exactly wrong for `C`/`Cdag`/`Sp`/`Sm`. It now applies the adjoint.

**Vacuum normalization (2026-08-29, all operators)**: `eshift`
(`window_total_energy`) is a Rayleigh quotient with the window's two
dangling boundary legs *traced*, while `snapshot_correlator` measures
with them closed by the transfer-matrix fixed points (identity on the
left, `rho` on the right — the closure `local_expectation`'s own
docstring derives). The two closings see different energies, so the
evolved branch carried a spurious factor relative to the un-evolved one:
measured on a dimerized spinless-fermion chain, the *unperturbed*
window's own `<ψ|ψ(t)>` stayed at arg = -0.0014 rad at t=0.4 through the
trace closure while rotating at 3.21 rad per unit time and growing to
|·| = 1.05 through the fixed-point one. `dynamical_correlator_td` now
evolves a second, unperturbed copy of the same window alongside and
divides every `S(x,t)` by its `<ψ|ψ(t)>`, which is exact for an
eigenstate (both branches carry the same `exp(-i(E0-eshift)t)`) and costs
one extra TDVP evolution. This was not a small correction: the residual
of the dimerized-XX free-fermion check dropped from ~0.07 — which the
test had attributed to iDMRG convergence quality — to ~1e-5, and the
fermionic Green-function check from 7.1e-1 to 4.2e-2 (that one is
window-limited, not phase-limited, at the parameters used).

`Infinite_Many_Body_Chain.td_dynamical_correlator` (`infinitechain.py`)
wires this up as the public API, mirroring `kpm_finite`'s own docstring
conventions (`p_i` sublattice position, an `n_window` convergence
caveat), calling `gs_energy()` automatically like `vev`/`correlator`.

**Native ITensor v3 port (`Chain::td_dynamical_correlator_window`,
`mpscpp3/chain_session.h`)**: `itensor_version=3` is also supported now,
via a genuinely different (and simpler) implementation than the
pyitensor path above, made possible by a real ITensor v3 capability
pyitensor's own from-scratch tensor engine has no equivalent of: the
vendored ITensorTDVP library (`mpscpp3/TDVP/tdvp.h`) already ships an
overload taking explicit left/right boundary tensors —
`tdvp(psi,H,t,LH,RH,sweeps,args)`, backed by `LocalMPO(H,LH,RH,args)` —
so the C++ port needs no window-aware sweep/environment machinery of its
own at all: it tiles the converged unit cell's own per-sublattice ket
tensors and automaton MPO rows into an ordinary window MPS/MPO, then
hands `idmrg_ground_state`'s own converged `HL`/`HR` straight to that
existing overload as `LH`/`RH`. `idmrg_ground_state` was extended to
snapshot `HL`/`HR` (and the per-sublattice ket tensors) as private
`Chain` state — mirroring pyitensor's own `IDMRGResult` snapshot, but
kept off the public `idmrg_ground_state()` pybind11 binding so its
existing 3-value return stays unchanged; `td_dynamical_correlator_window`
requires `idmrg_ground_state` to have been run first on the same `Chain`
(`infinitechain.py` keeps that `Chain` instance alive across calls in
`self._session3` for exactly this reason). One genuine bookkeeping
adaptation was needed: `idmrg_extend_HL`/`HR`'s own "bra" leg convention
(minted via `sim()`, an independent Index unrelated to the "ket" leg) is
not what `LocalMPO::makeL`/`makeR`'s own `dag(prime(psi))` convention
expects (a boundary tensor's bra leg being literally `prime()` of its own
ket leg) — a one-time `replaceInds` relabeling reconciles the two.
**Gauge (fixed 2026-08-29)**: this port originally tiled `idmrg_U_`, the
raw per-micro-step iDMRG factors, in both `idmrg_build_window` and
`idmrg_window_snapshot_correlator`'s own bra — not the gauge-consistent
unit cell that `idmrg_onsite_expectation`/`idmrg_two_point_correlator`
were moved onto for exactly this reason, and that `idmrg_window.py`'s own
`_window_cell` uses on the Python side. That failed the exact
`S(x,t=0) == two_point_correlator` identity by up to 1.7e-1 on a plain
spin chain, with no fermions involved; `x=0` stayed exact and the error
grew with `|x|`, the signature of a bond-basis mismatch, on a model whose
energy density the two backends agreed on to 6.7e-11 — which is why an
energy-only comparison never noticed. Both now tile the *raw* (not
re-gauged) cell, kept as `idmrg_cell_raw_`: only there are the two outer
legs still literally `idmrg_HL_ket_`/`idmrg_HR_ket_`, so the cell attaches
to the window's environment caps with no assumption at all — the same
choice, for the same reason, as `IDMRGResult.cell_raw` on the Python side.
Because the raw cell is not exactly left-canonical, the snapshot's
closures moved with it: a bare left trace became the cell's own left
transfer fixed point, and the result is divided by the same contraction
run over the ground state, so an operator-free chain returns exactly 1
(`iw_build_cache`, `idmrg_close_array_chain`).

Two further bugs surfaced in the same pass. First, that closure change
invalidated this backend's `exp(+i*eshift*t)` phase correction: `eshift`
is measured with the window's dangling boundary legs *traced*, and the
two closings genuinely see different energies, so `S(x,t)` came out exact
at `t=0` and drifted linearly after it (5.2e-2 by `t=0.15`). v3 now
co-evolves an unperturbed vacuum window and divides by its
`<psi|psi(t)>`, measured through the identical contraction — the same
exact construction `dynamical_correlator_td` uses, which needs no energy
measurement at all, so the whole `eshift` machinery is gone. Second, a
latent index-stride bug in `idmrg_close_array_chain`'s own chain
composition indexed the bra's left axis with the bra's *right* extent:
silent for as long as every tiled tensor was square, and wrong the moment
it is not — a dimerized spinless chain at `maxm=20` converges to a
(18,20)/(20,18) cell, where it put `S(x,0)` at 0.83 against an exact 0.51
for every operator. See `tests/test_idmrg_window_fermionic.py` (the t=0
identity and the free-fermion Green function, both on v3),
`tests/test_idmrg_window_v3.py`, and
`examples/idmrg/td_dynamical_correlator_python_VS_v3/main.py`, now an
agreement check again — and evaluated at every `x`, not just at `x=0`,
which is the one point a gauge error cannot touch and the reason the old
form of that example missed this for so long.

The window's own static overlap/measurement machinery (`S(x,t)` itself,
closed by the converged unit cell's own left and right dominant
transfer-matrix fixed points and calibrated against the operator-free
contraction, mirroring
`idmrg_window.py`'s own `_close_array_chain`/`_all_right_fixed_points`)
is, by contrast, implemented with plain dense arrays extracted from
ITensor tensors rather than ITensor's own Index-based auto-contraction —
for the same reason pyitensor's own version is: the window's
TDVP-evolved ket tensors have their own, freely SVD-re-minted internal
bond Indices with no shared identity to the statically-tiled bra tensors
used for measurement, so only their *dimensions* (not Index identity)
are guaranteed to align.

Scope differences versus the `"python"` backend: `x_values` may not
extend beyond the window's own explicit range (no beyond-window padding
is implemented — increase `n_window` instead); and, independently,
`itensor_version=3`'s own `idmrg_ground_state` has a known,
pre-existing divergence bug for Hamiltonians carrying an onsite ("field")
term, unrelated to this feature but inherited by it. See
`docs/user_guide.md`'s own iDMRG section for the user-facing version.

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

`get_mode()` is a two-step decision, and the split is deliberate.
`resolve_mode()` picks the solver on the availability grounds above
alone; only then does `get_mode()` check that answer against the chain's
conserved sector, if it has one. Written the other way round — the sector
guard first, as it originally was — the guard had to *predict* which
solver would answer, and it predicted wrong in both directions: it
refused ED outright (which now targets a sector perfectly well, §4.3a)
and it refused the extension-missing and `n<3` fallbacks (which land on
ED, i.e. on the one solver that can still answer correctly). This is
§4.10's pattern exactly: a precondition placed ahead of the branch that
qualifies it.

#### 4.3a Conserved sectors on the ED backend

ED is the third implementation of `Many_Body_Chain.set_conserved_sector`,
alongside `mpscpp3`'s block-sparse QN tensors (§4.4) and `pyitensor`'s
charge grading plus penalty (§4.5), and it is the simplest: a conserved
charge is diagonal in the ED product basis, so a sector is a *set of
basis states* and confining a calculation to it is taking a submatrix.
This is the same construction the sector tests have always used as their
independent reference (`ed_sector_energy`), promoted to a supported path.

Three pieces:

- **The charge operators come from the chain, as `MultiOperator`s.**
  `Many_Body_Chain.get_sector_charge_operators()` returns
  `{name: MultiOperator}` — `{"Nf": sum(self.N)}` for spinless fermions,
  `{"Sz": 2*sum(self.Sz)}` for spin chains (the factor of two is what
  puts `Sz` in ITensor's integer $2S^z$ units, the same units every other
  backend takes), `Nf`/`Sz` for both spinful chains, `{"Nb": sum(self.N)}`
  for bosons, and nothing at all for parafermions, which is what makes
  `set_conserved_sector` refuse them. Deriving the charges from the
  chain's own operators rather than from a per-ED-class hook is what keeps
  the two spinful chains — with quite different mode layouts — from each
  needing a special case.
- **`EDchain.set_conserved_sector(qns, charges)`** assembles each charge
  operator on the *full* Hilbert space (its own sector is still off at
  that point, so no recursion), checks it is diagonal and real
  (`_diagonal_of`), and intersects the per-quantity masks into
  `self._sector_mask`. An empty intersection raises rather than producing
  a chain every later call would fail in.
- **`EDchain.sector_restrict(A)`** is applied to every *assembled*
  many-body operator on its way out — `get_operator(MultiOperator)`,
  `MO2matrix`, `MBFermion.get_hamiltonian` — and to nothing else. The
  single-site matrices in `self.operators` stay full-space on purpose:
  they are the building blocks `multioperator.MO2matrix` multiplies
  together, and restricting a factor is not the same as restricting the
  product.

`sector_restrict` *refuses* an operator that does not commute with the
charge (checked on the operator's own nonzero entries — the charge is
diagonal, so this is just "does A connect two basis states of different
charge", with a relative tolerance so that an exactly-cancelling
$S^+S^+$/$S^-S^-$ pair reads as the cancellation it is). Refusing rather
than projecting matters: $PAP$ is exact for a static expectation value but
identically *zero* for a charge-changing operator, so a dynamical
correlator of `C` or `S+` inside a fixed-charge sector would come back as
a clean, wrong zero. That also matches what the DMRG backends do with a
flux-violating operator, so the same script behaves the same way on all
three.

Two asymmetries against the DMRG backends, both deliberate:

- **`promote_to_dense()`/`promote_mps()` stay DMRG-only.** They exist
  because a sector re-solve is expensive on a DMRG backend; in ED it is
  not, so clearing the sector and re-solving is the whole answer.
- **`Spinful_Fermionic_Chain` (Jordan--Wigner) can target `Sz` only under
  ED**, because its DMRG session is $2n$ *spinless* fermionic sites that
  carry `Nf` and nothing else, while its ED object has an explicit
  up/down mode per orbital. `Many_Body_Chain._ed_only_quantum_numbers()`
  names that case explicitly rather than having
  `_apply_conserved_sector` infer it from a swallowed exception — which
  would have silently thrown away the session's eager validation of
  *every* request, including unreachable targets on ordinary chains.
  `_sector_on_session` records whether the session actually carries the
  sector, and `get_mode()` refuses DMRG when it does not, so an ED-only
  sector can never be answered by a session that is not in it.

Regression coverage: `tests/test_sector_conservation_ed.py` (sector
energies against an independent full-space-then-restrict reference, ED
against `itensor_version=3` sector by sector, both spinful chains,
bosons, the enforcement errors, and the DMRG refusal on an ED-only
sector).

Above `algebra.maxsize` (2000) that diagonalization is iterative, and
`algebra.lowest_states` does **not** hand the Hermitian case to a single
`eigsh` call. A plain `eigsh(h,k=n)` silently loses members of a degenerate
level — it stops once the Ritz pairs it holds have converged, which one copy
of a degenerate eigenvalue already satisfies — so the ED path, the very one
used to judge whether a DMRG result is right, returned 1 of 4, 3 of 8 and 4
of 13 ground-level copies on a ferromagnetic spin-1/2 chain at dim 4096.
`algebra._deflated_lowest_hermitian` instead collects levels one at a time,
each round solving for the lowest eigenvalue of `P h P + sigma|V><V|` with
`P = 1 - |V><V|` over the eigenvectors found so far, which makes any
remaining copy the extremal eigenvalue rather than an already-converged one.
No setting of `tol`/`ncv`/over-requesting substitutes for this; measured,
each fixes some (N,n) and breaks others. The non-Hermitian branch (`eigs`)
is *not* deflated — the projector needs orthogonal eigenvectors — and can
still drop copies. See `docs/ed_vs_dmrg_degenerate_multiplets.md` and
`tests/test_ed_degenerate_levels.py`.

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
the effective MPO bond dimension, the bond-ramp settings below), or the
session object itself actually
changed since the last send (`_session_ham_cache`). Repeated calculations on an
unchanged Hamiltonian (e.g. successive `get_dynamical_correlator` calls,
each of which re-verifies the ground state) then hit the session's
caches instead of re-running warm DMRG sweeps and band-edge solves. Code
paths that *want* a fresh solve of the same Hamiltonian either pass
through `restart()`/`set_hamiltonian` (which force DMRG via
`skip_dmrg_gs=False`) or clear `_session_ham_cache` explicitly (see
`groundstate.py`'s best-of-`n` loops).

**The ground-state sweep schedule is ramped, not flat**
(`Chain::make_sweeps_ramped()` in both `chain_session.h`s,
`pyitensor/chain.py::_make_sweeps_ramped()`; the public knobs are
`Many_Body_Chain.bond_ramp`/`bond_ramp_start`/`bond_ramp_fraction`/
`bond_ramp_noise_decay`, documented physics-side in `user_guide.md` §3).
`Chain::gs_energy()` builds its `Sweeps` through this variant rather than
the flat `make_sweeps()`: the first `floor(nsweeps*bond_ramp_fraction)`
sweeps interpolate the bond dimension geometrically from
`bond_ramp_start` up to `maxm`, the rest run at `maxm`, and the noise
term decays by `bond_ramp_noise_decay` per ramping sweep and is off
afterwards. Three things about it are load-bearing:

- **It is plumbed as its own session setter** (`Chain::set_bond_ramp`,
  bound in both `bindings.cc`), not as extra `set_sweep_params`
  arguments, because `set_sweep_params` has ~15 call sites across
  `vev.py`/`kpmdmrg.py`/`excited.py`/`nhdmrg.py`/`timedependent.py` that
  have nothing to say about the ground-state ramp. Only
  `groundstate.py::gs_energy_single` calls it (`hasattr`-guarded, so an
  extension compiled before the method existed still loads).
- **Only `Chain::gs_energy()` uses it.** `make_sweeps()` itself is
  unchanged in shape and still backs `excited_states()`, the band-edge
  solves, and `gs_energy_generalized()`'s per-outer-iteration schedules —
  the last of which would restart the ramp from its floor on every
  Lagrange-multiplier update.
- **A warm start is floored, not truncated.** `gs_energy()` passes
  `maxLinkDim(wf0_)` (v2: `maxM`) as the ramp's `floor_dim` whenever a
  wavefunction is already present, so re-entering it with a converged
  state (`set_wavefunction`, or `gs_energy_single`'s own `maxde`
  reconvergence retry, which also disables the ramp outright) cannot
  throw that state away. A fresh start instead builds its random MPS
  directly at the ramp's first bond dimension.

The ramp's shape is deliberately a *fraction of the schedule* rather than
"double until you reach `maxm`". The latter was implemented first and
measured: it minimizes the number of cheap sweeps, so on a 30-site
inhomogeneous Heisenberg-Hubbard chain at `nsweeps=20` it leaves only 4 of
20 sweeps below a `maxm` of 150 and buys almost nothing. Filling half the
schedule instead measures **2.15x** at `maxm=90` and 1.7–2.0x at
`maxm=60` on that model, for the same energy to ~1e-8 (BLAS pinned to one
thread on two dedicated cores, as any timing in this repo must be — see
§19). The speedup grows with `maxm`, since that is what makes the ramped
sweeps cheap relative to the full ones. Rerun it with
`examples/groundstate/bond_dimension_ramp`.

While rewriting that schedule, a real pre-existing off-by-one in the flat
`make_sweeps()` was fixed in both C++ backends: `Sweeps` is 1-based (its
arrays are sized `nsweep+1`, and ITensor's `dmrg()` loops
`for(sw=1; sw<=nsweeps; ++sw)`), but the "turn the noise off for the
second half" loop was written 0-based (`for (i=ns/2;i<ns;i++)`), which
left the *final* sweep running with the full noise term — i.e. every
returned state and energy was the output of a noisy sweep.
`pyitensor/chain.py`'s mirror of the same loop was already correct, so
this also removes a silent C++/Python divergence.

Notable, deliberate implementation details (not bugs to "fix"):

- `mpscpp3` builds every site with `ConserveQNs=false` *by default*, and
  starts DMRG from an actual `randomMPS`, not a plain product state —
  matching v2's real (unconstrained-search) behavior rather than ITensor
  v3's stricter, QN-conserving-from-the-start convention, at the cost of
  losing the QN block-sparsity speedup for `itensor_version=3`.
  `Many_Body_Chain.set_conserved_sector()` is the opt-in that turns the
  quantum numbers back on for a chain whose caller wants the search
  confined to one sector — see "Conserved-sector mode" below, which also
  records why the "DMRG gets stuck at the product state" observation
  behind the default does *not* generalize to sector mode. The pure-Python
  backend implements the same public API by a different mechanism (dense
  storage plus a charge penalty, §4.5); no other backend has quantum
  numbers at all.
- **Conserved-sector mode** (`Chain::set_conserved_sector`,
  `mpscpp3/chain_session.h`): rebuilds `sites_` through
  `SpinX(site_types,conserved)` (`get_sites.h`) with QN-carrying indices,
  keeps a permanently dense `dense_sites_` alongside it, and drops every
  piece of session state built on the old site set. Three parts are worth
  knowing about, each forced by a hard ITensor behavior confirmed
  directly:
  - *The start state is a sum of random in-sector product states*
    (`sector_mps`), not `randomMPS`: `randomMPS(SiteSet)` refuses to
    guess a sector under QNs and `randomMPS(InitState,m>1)` does not
    exist in this vendored copy (both hard `Error`s in `mps.cc`). The
    per-site assignment comes from a small dynamic program over reachable
    partial charges (`sector_state_plan`), exact rather than greedy
    because a Hubbard chain at fixed `Nf` *and* `Sz` dead-ends a greedy
    fill; arrangements are then shuffled within each site-type code,
    which is charge-neutral by construction. The old `default_mps()`
    comment about DMRG being "stuck at the product energy no matter how
    large a noise term" is about the *all-first-basis-state* start it was
    measured with — that state is the unique state of its own sector, so
    nothing can move it. An in-sector start converges to the identical
    energy as the dense path (−5.142090632836 on a 12-site Heisenberg
    chain, twelve digits, from either a Néel or a sum-of-products start).
  - *Terms are normalized and flux-checked before AutoMPO sees them*
    (`sector_terms`, `mo_terms.h`'s `expand_xy_terms`/`combine_terms`).
    A non-conserving operator does not fail cleanly in ITensor: building
    `Sx` on Sz-conserving sites aborts the process inside `ITensor::set`,
    and so does handing AutoMPO a flux-violating term *even when its
    coefficient cancels exactly* ("Index does not contain given QN
    block"). That last one is why the textbook Sx.Sx+Sy.Sy+Sz.Sz
    Heisenberg Hamiltonian needs the Sx/Sy → S± expansion plus explicit
    cancellation of identical strings here, rather than being left to
    AutoMPO. An operator's charge is inferred from its *dense* matrix
    elements against the QN labels of the local basis states
    (`op_charge`), never by building it on QN indices. Since ITensor's
    own `Error()` calls `abort()` (`util/error.h`), the sector code
    throws `std::invalid_argument`/`std::runtime_error` instead, so the
    failures reach Python as catchable exceptions.
  - *Every `terms → MPO` path in the session goes through
    `mpo_from_terms`/`ampo_from_terms`*, so the check covers `vev`,
    correlators, KPM vertices and excited states, not just
    `set_hamiltonian`. Both helpers are pass-through when no sector is
    set, which is what keeps the default path byte-identical.
  On the Python side the sector lives on the chain
  (`Many_Body_Chain.conserved_sector`), is re-applied whenever a session
  is rebuilt (`sites.py::initialize`, `__deepcopy__`/`clone`), is part of
  `groundstate.py`'s Hamiltonian-send cache key, and is checked in
  `mode.py` against the solver that will actually answer (§4.3): DMRG on
  `itensor_version=2`/`"julia_live"` *raises*, since those have no
  quantum numbers and would silently return the global ground state
  instead of the sector's, while a fallback to ED is fine because ED
  targets the sector too (§4.3a).
- **Promotion out of a sector** (`Chain::promote_to_dense`/`promote_mps`,
  `Many_Body_Chain.promote_to_dense()`/`promote_mps(wf)`): leaves sector
  mode *keeping* the state computed inside it, which
  `set_conserved_sector()` with no arguments cannot do (it drops
  everything built on the QN site set). The conversion is ITensor's own
  `removeQNs` (`mps.cc`) — a block scatter into dense storage plus an
  index-level QN strip, exact and truncation-free — followed by
  `replaceSiteInds` onto `dense_sites_`. The point is that
  `sector_terms`' conserving-operator rule, unavoidable while the sector
  is set, forbids exactly the operators a fixed-$N$ state is usually
  wanted for (`C`, `Sx`, a pairing term); afterwards they are ordinary
  again. `dense_sites_` is therefore now the chain's *own* site set, kept
  from construction and never rebuilt — clearing the sector puts it back
  into `sites_` too — so that anything built before a sector was set, and
  any two states promoted out of *different* sectors, live on the same
  indices and stay comparable (a photoemission weight
  `|<gs_{N-1}|c_i|gs_N>|^2` needs exactly that). What promotion does not
  keep is the Hamiltonian MPO and the band-edge/iDMRG/VUMPS caches, all
  built against the QN indices; the ground state and its energy survive,
  so `computed_gs` stays true and a bare `gs_energy()` afterwards returns
  the sector's energy rather than re-solving unconstrained. `promote_mps`
  exists because `self.wf0` is a *separate* C++ MPS from the session's own
  (`groundstate.gs_energy_single` copies `gs_wavefunction()` out), as is
  any MPS the caller kept a reference to.
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
- `tevol_method="TEBD"` (`itensor_version="python"` via
  `pyitensor/tebd.py`'s `TEBDEvolver`, `itensor_version=3` via
  `mpscpp3/tebd.h` -- a from-scratch C++ port, not a shared code path --
  or `itensor_version="julia_live"` via `mpsjulialive/tebd.jl`, a third
  independent from-scratch port -- see "TEBD on `julia_live`" below)
  runs the standard 2nd-order-Trotter, even/odd-bond algorithm instead of
  TDVP: every bond's evolution gate `exp(-i*tau*h_bond)` is the exact
  exponential of the *bare* local 2-site Hamiltonian, built once from a
  `bond_hamiltonians()` function and reused unchanged for every subsequent
  time step — no per-step Krylov/Lanczos work at all, unlike TDVP's
  per-bond, per-step environment-projected exponentiation. `"python"`
  exponentiates a small dense matrix directly (`scipy.linalg.expm`);
  `itensor_version=3` instead builds the gate as an `ITensor` via
  ITensor's own `BondGate` primitive (`itensor/mps/bondgate.h`, its own
  Taylor-series `exp`, order 100 -- no `expm()` port needed). Only valid
  for a strictly nearest-neighbor Hamiltonian (`bond_hamiltonians()`
  raises `NotImplementedError` in `"python"`, a catchable `ITError` in
  `3`, for any term spanning 3+ distinct sites; span is computed from each
  term's *resolved* per-site operator against the identity, not from raw
  op-name lists, since both `MultiOperator` placeholder-site bookkeeping
  and native-spinful-site Jordan-Wigner string pairs can otherwise look
  like spurious long-range support — see `tebd.py`'s `_true_span()`).
  `mpscpp3/tebd.h`'s own `bond_hamiltonians()` reimplements the identical
  algorithm as `pyitensor/autompo.py`'s `HTerm.resolve()` (same
  carried-parity Jordan-Wigner threading, same per-site factor composition
  order) directly against ITensor v3's own site operators, composing
  same-site factors via ITensor's `multSiteOps()` primitive rather than
  going through `AutoMPO`/`toMPO()`: a probe against real ITensor v3
  confirmed a single, unsummed `HTerm` does *not* compile to a
  bond-dimension-1 MPO (a bare one-site term already comes out with
  uniform bond dimension 2 on every link), so there is no cheap way to
  "slice" a term's own local operator out of its compiled MPO the way
  other bond-dimension-1-only extraction tricks in this file work.
  `Chain::quench_tebd()`/`evolve_and_measure_tebd()` mirror the TDVP
  counterparts one-for-one (down to the ground-state-energy-shift
  convention `quench_tebd()`/`quench_tdvp()` share), so
  `get_dynamical_correlator(submode="TD")` (a `quench_tebd()` call under
  `timedependent.py`'s `evolution_dmrg_DC()`) and
  `evolve_and_measure_dmrg()` both support `"TEBD"` on `itensor_version`
  `3` or `"python"` the same way they support `"TDVP"`. Individually
  validated against exact diagonalization (both agree to ~1e-4) on small
  nearest-neighbor systems, including a native-Hubbard chain, in
  `tests/test_time_evolution.py`; the two `"TEBD"` backends additionally
  agree with *each other* to machine precision (~1e-14) on a fermionic
  hopping+onsite cross-check, confirming the C++ port's Jordan-Wigner
  threading and onsite-splitting match the Python reference exactly.
  TEBD's cost advantage scales up dramatically since TDVP's per-step cost
  scales with the local Hilbert-space dimension squared (16 for a
  dimension-4 Hubbard site) while TEBD's gates are built once and reused
  every step — see
  `examples/dynamical_correlator/dynamical_correlator_tebd_VS_tdvp_hubbard_20site`
  for a 20-orbital native-Hubbard timing (TEBD ~2s/step vs TDVP
  ~200s/step there, a >100x per-step speedup). At that system size and
  `maxm`, though, TEBD and TDVP no longer agree tightly with *each
  other* (found directly while building that benchmark: essentially the
  whole bulk of the chain is already truncated at the `maxm` cap from
  the first step, and TDVP's per-bond Krylov-approximated update vs
  TEBD's per-bond exact-exponential update are then two genuinely
  different approximations whose errors compound differently across
  ~16 saturated bonds) — a real characteristic of running two different
  truncated real-time propagators at a saturated bond dimension on a
  strongly-interacting model, not a correctness bug in either (ruled out
  by the small-system ED cross-checks above); see that example's module
  docstring for the full diagnosis. `evolve_and_measure_tdvp[_gse]`/
  `evolve_and_measure_tebd` all copy their input `wf` before evolving
  (`wf.copy()`) rather than evolving it in place — both TDVP's
  `_tdvp_step_fn` and `TEBDEvolver.step()` mutate their input MPS's
  tensor list via `set_A()`, so without the copy, evolving the default
  `wf=self.wf0` (the common case, see `timedependent.py`'s
  `evolve_and_measure_dmrg()`) would silently corrupt the cached ground
  state itself.
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
  `itensor_version` 2, 3 and `"python"`; the MPS Arnoldi route (default
  `algebra/arpacktk.py`, an Implicitly Restarted Arnoldi Method ported
  from ARPACK; `algebra/arnolditk.py`'s explicit-restart Arnoldi remains
  available for comparison, see `mpsalgebra.py`) remains as the fallback
  for other backends.
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
- `mpscpp3`-specific: `Chain::gs_energy_generalized`
  (`mpscpp3/chain_session.h`, exposed as
  `Many_Body_Chain.gs_energy_generalized`) solves the generalized
  eigenproblem $H|\psi\rangle=\lambda A|\psi\rangle$ for a Hermitian
  positive-definite metric MPO $A$ via a self-consistent
  Lagrange-multiplier iteration: each outer iteration builds the MPO
  $H-\lambda A$ (`sum()`) from the current $\lambda$ estimate, sweeps it
  with ITensor v3's own `dmrg()` for a single length-1 `Sweeps` step
  (looping outer iterations by hand rather than handing the whole
  schedule to one `dmrg()` call, since $\lambda$ must be updated between
  sweeps — noise is halved off partway through exactly like `nhdmrg()`'s
  own per-sweep loop above), then updates $\lambda$ to the swept state's
  generalized Rayleigh quotient
  $\langle\psi|H|\psi\rangle/\langle\psi|A|\psi\rangle$. A line-for-line
  port of `pyitensor/dmrg.py`'s `dmrg_generalized` (§4.5) against real
  ITensor v3 calls instead of the hand-rolled two-site sweep pyitensor
  uses; `mpscpp2` has no analogous session method yet, so
  `itensor_version=2` still raises `NotImplementedError`. See
  `examples/groundstate/dmrg_generalized_benchmark`, which heads all
  three available routes (pyitensor, v3, and `algebra/arpacktk.py`'s
  pre-existing ARPACK-mode-2 `mpsiram_generalized`) up against each
  other on the same test problem — v3 is consistently the fastest of the
  three at the (small) sizes benchmarked there.
- `mpscpp3`-specific: `Chain::nhdmrg_generalized` (exposed as
  `Many_Body_Chain.gs_energy_generalized` for a non-Hermitian
  `self.hamiltonian`) is the non-Hermitian counterpart of the bullet
  above, generalizing `Chain::nhdmrg` (§4.4's own NH-DMRG bullet) the
  same way `gs_energy_generalized` generalizes plain `gs_energy()`.
  Since the metric $A$ is Hermitian, $(\lambda A)^\dagger=\bar\lambda A$,
  so the same self-consistent Lagrange-multiplier trick carries over
  with a *complex* $\lambda$ and *biorthogonal* (not plain) expectation
  values: each outer iteration shifts both $H$ and its precomputed
  adjoint $H^\dagger$ by $\lambda A$/$\bar\lambda A$ respectively, runs
  one ordinary NH-DMRG sweep (a new private `nhdmrg_one_sweep` helper,
  factored out of `Chain::nhdmrg`'s own loop body) against the shifted
  pair, then updates $\lambda$ to
  $\langle\psi_L|H|\psi_R\rangle/\langle\psi_L|A|\psi_R\rangle$. Because
  `nhdmrg_one_sweep`'s two-site solve is hand-rolled directly against
  `arnoldi_smallest_real`/manual ITensor contractions rather than a call
  into ITensor v3's own `dmrg()`, it does *not* need
  `gs_energy_generalized`'s own short-chain guard (confirmed directly, a
  2-site chain runs it without aborting) — a real asymmetry between the
  two `mpscpp3`-specific bullets in this section, not an oversight.
  Pyitensor gained the identical algorithm first (`nhdmrg.py::
  nhdmrg_generalized`); `mpscpp2` still has no analogous session method
  for either backend's non-Hermitian generalized solver.

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

One partial exception to that shared-surface rule:
`dmrg.py::dmrg_generalized` (exposed as `Chain.gs_energy_generalized`/
`Many_Body_Chain.gs_energy_generalized`) solves the generalized
eigenproblem $H|\psi\rangle=\lambda A|\psi\rangle$ for a Hermitian
positive-definite metric MPO $A$ — a self-consistent Lagrange-multiplier
iteration built directly on top of `dmrg()`'s own per-sweep machinery
(factored into a shared `_dmrg_one_sweep` helper): each outer iteration
rebuilds the MPO $H-\lambda A$ from the current $\lambda$ estimate,
sweeps it with the ordinary two-site solver, then updates $\lambda$ to
the swept state's generalized Rayleigh quotient
$\langle\psi|H|\psi\rangle/\langle\psi|A|\psi\rangle$. `mpscpp3` has
since gained a line-for-line port of this same algorithm against real
ITensor v3 calls (§4.4's own `mpscpp3`-specific bullet), so
`itensor_version` 3 supports it too now — only `mpscpp2`
(`itensor_version=2`) and `"julia_live"` still raise
`NotImplementedError`, not having an analogous session method yet.

`nhdmrg.py::nhdmrg_generalized` is the non-Hermitian counterpart,
generalizing NH-DMRG (§4.4's own bullet above, and pyitensor's own port
of it) exactly the way `dmrg_generalized` generalizes plain `dmrg()`:
`groundstate.gs_energy_generalized` checks `self.hamiltonian` for
Hermiticity and transparently dispatches a non-Hermitian one here instead
of raising. Since $A$ is Hermitian, $(\lambda A)^\dagger=\bar\lambda A$,
so the same trick carries over with a *complex* $\lambda$ and
*biorthogonal* (not plain) expectation values: each outer iteration
shifts both $H$ and its precomputed adjoint $H^\dagger$ by
$\lambda A$/$\bar\lambda A$ respectively, runs one ordinary NH-DMRG sweep
(`_nhdmrg_one_sweep`, factored out of `nhdmrg()`'s own loop body the same
way `_dmrg_one_sweep` was) against the shifted pair, then updates
$\lambda$ to $\langle\psi_L|H|\psi_R\rangle/\langle\psi_L|A|\psi_R\rangle$.
`mpscpp3` has since gained a line-for-line port of this algorithm too
(§4.4's own second `mpscpp3`-specific bullet), so `itensor_version` 3
supports a non-Hermitian `self.hamiltonian` here just like the Hermitian
case above — `mpscpp2` is still the only backend without an analogous
session method for either. See
`examples/non_hermitian/nhdmrg_generalized_benchmark`, which cross-checks
both implementations against `mpsiram_generalized` (ARPACK mode 2 needs
no adaptation at all for a non-Hermitian primary operator, only its own
$M$ positive-definite precondition was ever required) on a non-Hermitian
test problem: the same accuracy/speed pattern as the Hermitian benchmark
holds against ARPACK, but unlike the Hermitian case v3 isn't consistently
faster than pyitensor here — NH-DMRG's per-bond cost already pays for
*two* Arnoldi solves (right block and its adjoint) regardless of backend,
narrowing the compiled-vs-pure-Python gap relative to plain ground-state
DMRG's single local diagonalization per bond.

**Conserved-sector mode (`pyitensor/sector.py`).** The pure-Python
backend implements the same `set_conserved_sector`/`promote_to_dense`
surface as `mpscpp3` (§4.4), reached by the identical
`Many_Body_Chain.set_conserved_sector(Nf=6)`/`(Sz=0)` call and giving the
same answers, but by a different mechanism — and the difference is the
interesting part.

An `Index` here gains an optional *charge grading*: one integer charge per
basis state, per conserved quantity (`index.py`, from the per-site-type
`_QN` tables in `sites/base.py` and friends, in ITensor's own units — `Sz`
in integer $2S^z$). That is all it gains. Real ITensor v3's QN indices also
sort their basis states into contiguous per-charge blocks and make every
tensor built on them block-sparse; these only *label* the states they
already had, so storage stays dense, `Link` indices are never graded, and
no contraction or SVD in the engine changes at all. `SiteX(site_types,
conserved)` mints the graded set with fresh Index identities (mirroring
v3, where a QN site set is necessarily a different object), and the
session keeps its own permanently dense `dense_sites` alongside, so
`promote_to_dense`/`promote_mps` is exactly an index relabeling —
position-based, so it still works on a state handed back long after the
sector that produced it was cleared.

Three consequences of the dense choice:

- *The term-list normalization is still needed, for a different reason.*
  `sector_terms` transcribes `mo_terms.h`'s `expand_xy_terms` +
  `combine_terms` + charge check, so the textbook $S^xS^x+S^yS^y+S^zS^z$
  Heisenberg Hamiltonian is accepted (its $S^+S^+$/$S^-S^-$ strings cancel
  exactly once expanded) and a genuinely non-conserving operator is
  rejected by name. In `mpscpp3` this exists because ITensor *aborts the
  process* over a flux-violating term; here nothing would abort — the
  check exists so that the operator being diagonalized provably commutes
  with the conserved charges, which is what makes the answer the sector's
  answer. Every `terms → MPO` path in `chain.py` routes through
  `Chain._ampo`/`_mpo`, so the check covers `vev`, correlators, KPM
  vertices and time evolution, not just `set_hamiltonian`; both are
  pass-through outside sector mode.
- *A charge penalty is needed, which `mpscpp3` does not need.* Contraction
  preserves exact zeros, but a dense LAPACK SVD does not: Householder
  bidiagonalization mixes rows across charge blocks, so every truncation
  injects ~$10^{-16}$ of amplitude into neighboring sectors, and a
  variational sweep amplifies that geometrically toward whichever sector
  is lower in energy. So the *variational* solves (ground state, excited
  states, band edges — nothing else) minimize $H+\lambda\sum_k(\hat
  Q_k-q_k)^2$ instead of $H$, with the penalty MPO built in closed form as
  the standard bond-dimension-3 $(\sum_i B_i)^2$ automaton
  (`charge_penalty_mpo`) rather than by handing $O(N^2)$ two-site terms to
  `mpobuilder.py`. It is identically zero on the target sector and
  `gs_energy` reports $\langle H\rangle$ under the plain Hamiltonian
  regardless, so no reported number changes; $\lambda$ defaults to the sum
  of the Hamiltonian's coefficient magnitudes (enough to outweigh its
  spectral range, so the target sector's ground state is the global
  minimum of the penalized operator) and `Chain.set_sector_penalty` (`Many_Body_Chain.set_sector_penalty`)
  overrides it. This is not defensive coding: measured on a 12-site chain
  with attractive $V$ asked for $N_f=2$, an in-sector start with the
  penalty disabled converges to the *full* band ($\langle(\hat Q-q)^2
  \rangle = 144$), and with it enabled reproduces sector-restricted ED —
  `tests/test_sector_conservation_python.py::test_the_charge_penalty_is_load_bearing`
  locks that down. Every solve also checks its converged state's
  $\langle(\hat Q-q)^2\rangle$ and raises rather than reporting the wrong
  sector's energy.
- *The start state is a sum of random in-sector product states*
  (`sector_mps`), from the same exact dynamic program over reachable
  partial charges as `chain_session.h`'s `sector_state_plan` (greedy
  dead-ends on a Hubbard chain at fixed `Nf` *and* `Sz`). A single product
  state would be a poor start for the same reason `default_mps()` avoids
  one: this backend's DMRG has no noise term at all, so a start that is an
  exact eigenstate has nothing to move it. Routing every solve through
  `Chain._default_mps` means excited states and band edges inherit the
  in-sector start for free.

Two smaller differences from `mpscpp3`. `tevol_method="TEBD"` *works*
here, where v3 refuses it: v3's gate assembly sums bond gates before their
fluxes agree, while these gates are dense and simply inherit the
Hamiltonian's charge conservation (checked against dense TDVP/TEBD to
5e-7). And because mismatched dense site indices do not fail loudly — they
would silently not contract, producing an outer product — `Chain` checks
the site set of every wavefunction handed to `vev`/`apply_operator`/
`overlap`/`set_wavefunction` and raises instead. METTS refuses in sector
mode on both backends, for the physical reason on this one: it averages
over an ensemble sampled from every sector at once.

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
`examples/time_evolution/tdvp_VS_ED_time_evolution`): `evolve_and_measure` agreed with
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
this — they route through `gs_energy()`/the generic MPS Arnoldi method
(default `algebra/arpacktk.py`'s IRAM, via `mpsalgebra.mps_excited_states`;
`algebra/arnolditk.py`'s `mpsarnoldi` remains available for comparison),
neither of which is backend-specific.

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
the same tolerance `examples/dynamical_correlator/dynamical_correlator_VS_ED/main.py` uses
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

METTS (`metts_vev`, finite-temperature sampling, arXiv:1002.1305) needed
no new evolution primitive either — its imaginary-time half-beta
propagation is exactly `tdvp_step` again, with a purely imaginary rather
than purely real/complex effective step, the same generalization TDZ's
paragraph above already relies on. What it does need, unlike CVM/TDZ, is
genuinely new machinery: the sequential-sampling collapse of an MPS onto a
new classical product state (`mpsjulialive/metts.jl`'s
`metts_collapse_to_cps`), since nothing already backend-agnostic performs
that. Rather than deriving it as a fresh ITensor-level tensor primitive
(the route `mpscpp3/chain_session.h`'s C++ port took, via `diagHermitian`/
`setElt`), `metts.jl` is a value-level port of `pyitensor/metts.py`'s
already-validated NumPy implementation: at each site, extract the local
MPS tensor as a plain Julia array (`Array(T, s, right_link)`), rotate into
the eigenbasis of the requested single-site operator via a plain matrix
product, sample from the resulting marginal, and rebuild a bond-dimension-1
product state from the sampled eigenvectors by hand (mirroring
`build_product_state()`/`product_state()` on the other two backends).
`LinearAlgebra.eigen(Hermitian(...))` replaces `np.linalg.eigh`/
ITensor's `diagHermitian` for the per-site eigendecomposition (cached once
per `(opname, site)` pair, same as both other backends, since it never
changes across a sampling run), and a hand-rolled inverse-CDF categorical
sampler replaces `np.random.Generator.choice`/`std::discrete_distribution`
(no `Distributions.jl` dependency declared in `juliapkg.json`, and this is
the only place one would be needed). The whole `nwarmup+nsamples` Markov
chain runs inside one Julia call, same design as `kpm.jl`'s moment
recursion and `tdvp.jl`'s `quench_tdvp`. Passing the list of observable
MPOs to sample together across that one call needed no explicit
Python→Julia list conversion (unlike `sites.jl`'s `get_sites`/
`read_operator.jl`'s `toMPO`, which both need `to_julia_strvec` because
their signatures declare a strict `Vector{String}` parameter): confirmed
directly that `juliacall`'s default `PyList{Any}` wrapping already
transparently unwraps each element back to its real underlying Julia
object on indexing, for any Julia function that doesn't itself declare a
stricter parameter type. Validated against exact ED thermal averages the
same way `tests/test_metts_vev.py` already does for the other two
backends (several seeds, energy and local magnetization, generous
Markov-correlated-error tolerance) — no systematic bias, deviations
consistent with the reported (correlated, so likely optimistic) standard
error. `metts_dynamical_correlator` (the real-time finite-temperature
generalization, arXiv:2405.18484) is now ported too, in the same
`metts.jl` file: for each retained METTS sample, `apply_op`
(`mpsalgebra.jl`, the same primitive KPM/TDVP already share) applies the
`B` operator once to build `v(0)=B|phi>`, then `tdvp_step` — the exact
same real-time-evolution primitive `evolve_and_measure_tdvp`/`quench_tdvp`
already use, just with a purely real rather than purely imaginary `dt` —
advances `v(t)`/`w(t)=phi(t)` independently under `H`. No manual
norm-restoration step is needed the way `mpscpp3/chain_session.h`'s C++
port requires for the deliberately non-unit-norm `v(t)` (see that file's
own comment): this backend's `tdvp_step` already passes `normalize=false`
to ITensorMPS's `tdvp()` for every caller (unlike ITensorTDVP's hardcoded
`DoNormalize=true` the C++ port has to work around), so `v`'s norm simply
falls out of the (truncated but otherwise unitary) evolution itself.
`tdvp_cutoff`/`tdvp_maxdim` are wired through as ordinary per-call
arguments to Julia's `tdvp_step`/`apply_op`, unlike v3's hardcoded
restriction to this chain's own `cutoff`/`maxm`.

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
`examples/staticcorrelators/four_correlation_tensor_spinful_native`). Before this
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

A third implementation, `ctmode="sweep"` (`Chain::four_correlation_tensor_sweep`,
`mpscpp3/chain_session.h`/`bindings.cc`, plain non-native-spinful sites
only), replaces `ctmode="full"`'s `O(N^4)` *independent* AutoMPO builds
with a single-sweep, environment-reuse algorithm following the idea of
[ITensorCorrelators.jl](https://github.com/ITensor/ITensorCorrelators.jl).
Only the `N(N-1)(N-2)(N-3)` pairwise-distinct-index entries go through
the fast sweep (nested loops over four strictly increasing sites
`a<b<c<d`, each level extending a running left environment by exactly
one site rather than rebuilding it, with both `Cdag`/`C` tried at each
of the four operator sites); the remaining, subdominant `O(N^3)`
repeated-index entries fall back to the same per-tuple AutoMPO method
`ctmode="full"` uses, since those need same-site multi-operator products
(e.g. `Cdag_i C_i = N_i`) the sweep isn't built to produce. Reordering
the four fermionic operators from the abstract `(i,j,k,l)` order into
physical site order picks up the standard fermion-anticommutation sign
(parity of the reordering permutation) *plus* one further correction:
of the six possible `Cdag`/`C` type patterns across the four sorted
sites, the two that strictly alternate (`Cdag,C,Cdag,C` or
`C,Cdag,C,Cdag`) need no extra sign, but the other four (which have at
least one pair of *adjacent same-type* operators) need one more flip on
top of the plain permutation parity — this was found empirically (every
entry with such a pattern disagreed with `ctmode="full"` by exactly this
sign, for every site quadruple and chain length tried, before the
correction was added) rather than purely re-derived by hand from Jordan-
Wigner strings, though it traces to the same extra `-1` that
`jordanwigner.py`'s `CC()`/`CdagCdag()` carry but `CdagC()`/`CCdag()`
don't (reordering two *same-type* operators into JW-canonical form costs
one more local `F`/`Adag` or `F`/`A` anticommutation than a mixed pair
does). Validated to machine precision against `ctmode="full"` and to ED
solver tolerance across `n=4..8` (every index tuple, not just a sample)
and spot-checked at `n=16`; measurably faster than `ctmode="full"` at
every size tried, by a margin growing with `n` (~1.4x at `n=6` up to
over 2.5x at `n=16` for a nearest-neighbor Hubbard chain at fixed
`maxm=60` — absolute numbers are noisy across runs on a shared machine
(measured 2.8x-3.4x at `n=16` across separate runs), but the *trend*,
measured within any single run, consistently holds — see
`examples/staticcorrelators/four_correlation_tensor_sweep_VS_full`,
whose Part 3 loop actually runs and reproduces this trend, rather than
reporting only development-time measurements that would not be
reproducible from the checked-in tree; `n=20` was also checked
informally during development, ~4.4x, but is not part of the shipped
example since a single `n=20` `ctmode="full"` call alone takes ~5
minutes). `accelerate` (default `True`) only gates the subdominant
repeated-index fallback here — unlike `ctmode="full"`'s own
`accelerate`, which skips ~half of the *dominant* per-tuple AutoMPO
builds via conjugate-pair symmetry, this method's dominant pairwise-
distinct sweep pays its cost (the shared environment sweep, and six
cheap per-pattern `eltC()` evaluations) once regardless of how many
output entries a given leaf value gets scattered into, so there's no
equivalent saving to skip there — confirmed directly (`accelerate=True`
vs.\ `False` gives identical wall-clock on the dominant part).

`four_correlation_tensor_sweep()` originally renormalized its internal
copy of `wf` before the fast-sweep computation (`psi /=
innerC(psi,psi).real()`, mirroring `reduced_dm()`'s own convention), but
the repeated-index fallback loop calls `innerC(wf,...)` on the raw,
un-renormalized `wf` — a real bug, caught by code review, not by any of
the (unit-norm-only) validation above: for any non-unit-norm input MPS
(e.g. `wf*c` for a scalar `c`, as opposed to an already-converged DMRG
ground state, which is unit-norm to solver precision and so didn't
expose it), the two halves of one output tensor came back scaled
inconsistently with each other *and* with `ctmode="full"`/`"explicit"`.
Fixed by dropping the renormalization entirely, matching
`four_correlation_tensor()`'s own convention (no enforced normalization,
caller's responsibility) exactly — confirmed with a direct repro
(`wf*3.0`: `ctmode="full"` scales every entry by $3^2=9$ as expected;
`ctmode="sweep"` did too only on the repeated-index entries, and by
$1/9$ on the pairwise-distinct ones, before the fix). The same bug,
same fix, applies to the `pyitensor/chain.py` port below.

Not ported to `Spinful_Fermionic_Chain_Native`: that class's
`ctmode="full"` gets its Jordan-Wigner threading for free from ITensor's
own `AutoMPO` acting on flavor-resolved operator names, since two
flavors share one physical site there — a different threading problem
than this sweep's per-flat-mode gap rule solves, not attempted here.

The same algorithm was then ported to `pyitensor/chain.py::
Chain.four_correlation_tensor_sweep` for `itensor_version="python"` —
close to a mechanical transcription of the C++ version, including
`_four_pt_perm_table()` (the same sign/index-role table, generated with
`itertools.permutations` instead of hand-rolled `next_permutation`), with
one adaptation: this engine never primes `Link`-tagged indices to keep a
self-overlap's bra and ket from colliding the way `chain_session.h`'s
`dag(prime(psi.A(site),TagSet("Link")))` does (see `mpsalgebra.py`'s
`inner()`) — instead the bra side is built once per outer site `a` via
`mpsalgebra.py`'s own `_fresh_link_copy(psi)` (freshly relabeled `Link`
indices, physical indices untouched), and the running environment
carries one leg from the *ket*'s own unrelabeled links and one from this
fresh bra copy's links; the boundary `delta` "shortcut" tensors are built
from `commonIndex` on each side (ket vs.\ fresh-bra) rather than from a
single primed/unprimed pair. Matched the C++ version and ED to machine
precision immediately (0 mismatches out of every distinct- and
repeated-index tuple tried, `n=5`), confirming the ported sign table is
correct and the fresh-link-copy substitution is a faithful translation,
not just an analogous-looking one. Real but smaller speedup than the C++
version (~1.2–1.4x at `n=5..7`, vs.\ ~1.4–4.4x at `n=6..20` for
`itensor_version=3`): pure-Python loop/function-call overhead is a
bigger fraction of the per-step cost here than the underlying NumPy
linear algebra, so there is less for environment reuse to save on
relative to the total.

### The four-point tensor as batched GEMMs (`pyitensor/fourpoint.py`)

`ctmode="batched"` is a third implementation of the same tensor, for
`itensor_version="python"` only, and the default there. It exists because
both of the others are shaped wrong for the hardware. `"sweep"` and
`"fold"` evaluate one tuple's worth of MPS transfer at a time -- a chain of
`chi x chi` by `chi x d x chi` contractions -- and there are `O(n^4)` such
tuples, each contraction individually far too small to keep a BLAS call
busy. On top of that, `"sweep"` shares environments only across the
pairwise-*distinct* tuples; the `O(n^3)` repeated-index ones fall back to
independent per-tuple folds, and a smaller count times a costlier per-tuple
price still wins: instrumented at `maxm=20`, that fallback is 61–65% of the
sweep's own runtime at `n=12`, `16` and `20`. Fixing only the distinct half
would therefore cap the achievable speedup below 1.6x.

The reorganization is a trie. Written in site-sorted order, a tuple becomes
a sequence of at most four *ranks* -- rank `r` is the `r`-th smallest
distinct site the four factors occupy, carrying the product of whichever
factors sit there, with the Jordan-Wigner `F` folded in. Two tuples that
agree on their first `r` `(matrix, parity)` steps have the same partial
environment *whatever sites those steps landed on*, so all `n_modes^4`
tuples collapse onto a trie of a few dozen matrix sequences (75 shapes for a
spinless chain: `m! * S(4,m)` = 1, 14, 36, 24 for `m = 1..4` distinct sites,
times the mode-name combinations, minus every shape whose composed local
matrix vanishes -- which on spinless sites is most of the repeated-index
enumeration, `Cdag*Cdag = 0`). Each trie node holds one array of
environments, `(B, chi, chi)` with `B` up to `C(n, r)`, batched over the site
combinations that realize it, and one site of the sweep becomes, per node,
exactly two GEMMs:

    X[b,j,p,k] = sum_i  E[b,i,j] A_op[i,p,k]
    E'[b,k,l]  = sum_jp X[b,j,p,k] conj(A)[j,p,l]

each reshaped into a single `(B*chi, chi) x (chi, d*chi)` matmul. Distinct
and repeated tuples go through the same machinery, so the trie subsumes the
sweep's environment reuse and extends it to the half the sweep never
covered.

Every *convention* is inherited from `four_correlation_tensor_fold` rather
than re-derived, because that is the one implementation covering all tuples
under a single rule and already validated against the AutoMPO+`to_mpo`+
`inner` pipeline: operator matrices indexed `[ket, bra]` and applied as
`tensordot(A, mat, ([1],[0])).transpose(0,2,1)`; same-site factors composed
in increasing factor-position order; and the Jordan-Wigner string reduced to
one statement covering gap sites and operator sites alike — *a site carries
an extra `F` iff the number of factors placed up to and including it is
odd*, which is what the fold's incremental `carry != odd` test amounts to.
The cross-site reordering sign `_four_pt_site_sort_sign` depends only on the
rank assignment, so it is constant across every site combination realizing
one shape and is applied once per shape at scatter time.

The one thing not inherited is the gauge. The fold re-runs
`psi.position(mn)` per group so it can start from an identity left
environment; a single left-to-right batched pass cannot re-gauge per tuple,
so this calls `psi.position(1)` once and grows the operator-free left
environment `E0` alongside the trie. Closing is still a plain trace, because
`position(1)` leaves every site to the right right-canonical.

Two implementation notes that are load-bearing:

- **The arrays are read through `backend.xp().asarray`, not
  `np.ascontiguousarray`.** `chain.py::_mps_arrays_lpr` uses the latter,
  which on a device backend is a silent device→host copy per site — the
  exact round trip the GPU port exists to avoid — so `fourpoint.py` carries
  its own `_arrays_lpr`. That single conversion, `O(n chi^2 d)`, is the only
  host↔device traffic in the whole calculation apart from one leaf-value
  transfer per site; `tests/test_four_point_correlator.py::
  test_batched_stays_on_the_device` asserts it by the returned array's
  *type*, since a kernel that transferred every step would return identical
  numbers, only slower.
- **It is deliberately not routed through `backend.jit`.** The composites
  that module fuses have shapes fixed by the bond dimension, which
  `set_pad_bonds` can freeze; this one's leading axis is the environment
  batch, which changes at every site and every trie node by construction, so
  jitting would trade one dispatch for thousands of retraces. The batch is
  what replaces the fusion: one GEMM here already carries what would
  otherwise be tens of thousands of calls. That is also why this is a GPU
  target at all — at the ~0.35 ms eager dispatch floor Phase 0 measured,
  164430 separate `chi=20` contractions would cost ~57 s of pure dispatch at
  `n=30`, more than the entire calculation costs on one CPU core.

Measured single-threaded (`MKL_NUM_THREADS=1`), spinless chain, `maxm=20`,
full tensor, `"sweep"` / ITensor v3 (C++) / `"batched"`: `n=16`, 10.65 s /
6.69 s / 0.38 s; `n=30`, 117.3 s / 72.4 s / 7.4 s — 15.9x over the
pure-Python sweep it replaces and 9.8x over the compiled C++ backend, from
NumPy. (The `n=30` v3 number is a separate process on the same host with the
same pinning.) For `Spinful_Fermionic_Chain_Native` the comparison is
against `"fold"`: 1.29 s → 0.22 s at 5 sites, 3.02 s → 0.31 s at 6, with the
crossover at 4 sites — below that the trie build is a fixed cost the tensor
is too small to amortize. See `examples/staticcorrelators/
four_correlation_tensor_batched_VS_sweep`.

On a device the same kernel is the one calculation in this library that
wins below the chi ~ 120-160 crossover `docs/gpu_cpu_performance.md`
establishes for everything else, and for a reason worth stating plainly:
the arithmetic per array operation comes from the tuple batch, not from
bond dimension, so it does not need a large chi to clear the dispatch
floor. H200 against one Xeon core, n=30, warm: 0.97 / 0.95 / 0.98 s at
`maxm` = 20 / 40 / 80 against 1.89 / 9.44 / 22.81 s on the host (2.0x,
9.9x, 23.3x) and against 24.7 / 40.0 s for the compiled v3 backend at the
first two (25x, 42x). The device column does not move across the range --
it is entirely dispatch-bound, the GEMMs are free at these sizes -- so
every one of those ratios is a lower bound and the growth is the host
getting slower. The price is a ~20-25 s near-constant cold run, which is
XLA compiling one kernel per array shape and which `set_pad_bonds` cannot
touch here (what varies is the batch, not the bond): worth it for several
tensors in one process, or at maxm=80 where a single cold run is already
break-even at 0.92x. The ground state does not need to be on the device --
`_arrays_lpr` converts the MPS once, so `backend.set_backend("jax")`
immediately before the call is enough, which sidesteps DMRG's own
below-crossover device tax.

**Blocking, and why the block width is a budget rather than a constant.**
The trie as described holds `O(n^3)` live environments at level 3 --
`6*C(n,3)`, which is fine at n=30 (140 MB at chi=20) and 6.0 GB at n=100.
The sweep is therefore split on the *first occupied site*: tuples whose
smallest site falls in a given range are swept together. That split is exact
rather than an approximation, and for a reason that falls straight out of
the trie -- two tuples merge only when their `(matrix, parity, site)`
prefixes agree, so environments with different first sites never meet, and
partitioning on that site costs no arithmetic at all.

The block *width* is the knob, and both extremes are wrong. One block is the
6.0 GB above. One block per site holds only `6*C(n-a,2)` (0.19 GB at n=100,
chi=20) but issues n times as many array operations -- and on a device each
pays the ~0.35 ms dispatch floor, so n=100 would spend ~50 s in pure
dispatch, reintroducing exactly the cost the batching exists to remove. So
`_block_width()` solves a byte budget (`_LIVE_BYTES`, 2 GB) for the width:
one block for everything up to about n=40, 11 first-sites per block at
n=100/chi=20, one per site only when chi is large enough to force it. The
scatter is bounded the same way -- codes are accumulated incrementally as
the trie descends (`code*(n+1)+site`) and a leaf node's values are aligned
to `itertools.combinations` order by a single `argsort`, rather than through
a dense `(n+1)**m` table that would be 1.66 GB per m=4 node at n=100.

Measured at n=100 (maxm=20, chi=20, 10^8 tuples, 23.5M distinct-index leaf
values, a 1.6 GB output tensor): 481.6 s on one Xeon core, 19.9 s warm on an
H200 -- and, unlike n=30, 103.7 s *cold* as well, i.e. 4.65x even for a
single one-shot tensor, because the near-constant compile cost is now small
against an 8-minute host run. At maxm=40: 2171.8 s / 79.1 s warm / 328.6 s
cold. Agreement with the host stays at 1.8e-15 to 2.8e-15.




`get_four_correlation_tensor(ctmode=...)`'s default changed from a fixed
`ctmode="explicit"` to `ctmode=None`, resolved per-wavefunction by the
new `_four_correlation_tensor_default_ctmode()`
(`entropytk/correlationentropy.py`): `"sweep"` whenever it applies
(`itensor_version` in `(3,"python")`, non-native-spinful), else
`"full"` whenever it applies (`itensor_version` in `(2,3,"python")`, or
native-spinful under `itensor_version=3`), else `"explicit"`. An
explicit `ctmode=` argument is unaffected — it is still a hard request
that raises if unavailable, not a hint the resolver can override. ED-
backed wavefunctions (`edtk/edchain.py`, `pyfermion/mbfermion.py`) never
reach the resolver at all: they hardcode `ctmode="explicit"` themselves,
as before. Verified the resolver picks the right mode across every
backend actually reachable from the public API (`itensor_version` `2`
falling back to ED when the C++ extension is not compiled, `3`,
`"python"`, and `Spinful_Fermionic_Chain_Native`), and made it use
`getattr()`/`hasattr()` defensively rather than assume `wf.MBO` always
carries `.itensor_version`/`._session` — cheap insurance since the
function is easy to reach with an unusual `wf.MBO` from outside the two
call sites it is actually designed for.

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
  example exercising it, `examples/utilities/power_vev/main.py`, sets
  `sc.itensor_version = "julia"` directly (bypassing
  `setup_julia()`/`_reset_dmrg_state()` entirely) — the legacy,
  already-inert subprocess-based backend documented below, not
  `julia_live` — so it wouldn't validate this fix even if made.
- **`get_distribution`/`kpmdmrg.py::general_kpm`** (spectral distribution
  of an arbitrary operator, not just the Hamiltonian; has example
  coverage in four `examples/magnetization/magnetization_distribution*` scripts and
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

#### Generalized-eigenvalue and non-Hermitian DMRG on `julia_live`

`mpsjulialive/generalized.jl` and `mpsjulialive/nhdmrg.jl` close what
`ROADMAP.md` listed as the single biggest structural gap in this backend:
the generalized eigenproblem $H|\psi\rangle=\lambda A|\psi\rangle$ and
non-Hermitian DMRG, both of which previously existed only on
`itensor_version` `3`/`"python"` and raised `NotImplementedError` here.

Both follow the same "reuse the backend-agnostic algebra, write a Julia
primitive only for the piece that genuinely needs it" pattern as CVM and
TDZ above:

- `generalized.jl::get_gs_generalized` is the Hermitian solver — the same
  self-consistent Lagrange-multiplier iteration as
  `pyitensor/dmrg.py::dmrg_generalized` and
  `Chain::gs_energy_generalized`, i.e. rebuild `Heff = H - lambda*A`
  exactly (`add(...; cutoff=0.0)`), run one ordinary `dmrg()` sweep
  against it, reset `lambda` to the swept state's generalized Rayleigh
  quotient, repeat. The whole loop runs in one Julia call (same design as
  `kpm.jl`/`tdvp.jl`), so no MPO/MPS marshalling per outer sweep;
  `groundstate.gs_energy_generalized` gained a `julia_live` branch, and
  its pre-existing `self._session is None` guard was narrowed to exclude
  this backend, which legitimately has no session object at all.
  Cross-checked against `scipy.linalg.eigh`'s exact generalized
  eigensolver and against pyitensor: agreement to ~1e-15
  (`tests/test_dmrg_generalized.py`, now parametrized over `julia_live`
  too).
- `nhdmrg.jl::get_nhdmrg_pair` and
  `generalized.jl::get_gs_generalized_nhdmrg` are the non-Hermitian plain
  and generalized solvers. Unlike the other three backends, the local
  sweep here is *not* a dmrgpy port of ITensorNHDMRG.jl — it is
  ITensorNHDMRG.jl (already a declared dependency in `juliapkg.json`, and
  already used by `get_gs.jl::get_gs_nhdmrg` for plain non-Hermitian
  ground-state *energies*, which discards the left vector). The retry
  loop and two-sided eigen-residual certificate in `nhdmrg.py` are built
  purely on generic `MultiOperator`×MPS algebra, so they are shared
  unchanged; `mpsjulialive/nhdmrg.py` supplies only a per-attempt
  function, and `nhdmrg.py`'s `nhdmrg()`/`nhdmrg_generalized()` were
  refactored to call an `attempt()` closure instead of `self._session`
  directly.

Two ITensorNHDMRG-specific corrections were needed, both documented in
full in `nhdmrg.jl`'s header, and both **invisible on a
complex-*symmetric* Hamiltonian** — which is what every non-Hermitian
model previously in `tests/`/`examples/` happens to be, since their
hoppings are symmetric and every non-Hermitian piece is diagonal. Pinning
them down needed a genuinely non-symmetric model, hence
`tests/test_nh_dmrg.py::nh_asymmetric_hopping_chain` and
`examples/non_hermitian/nhdmrg_julia_asymmetric_VS_ED` (Hatano-Nelson
asymmetric hopping + staggered imaginary potential):

1. **Left-vector convention.** ITensorNHDMRG's `ProjNHMPO` builds its
   "adjoint" projector from `dag(swapprime(conj(H), 0 => 1))`; since
   these sites carry no QNs (`sites.jl` uses plain `Index(...)`), `dag()`
   is just a conjugation, the two conjugations cancel, and what it
   actually sweeps against is `swapprime(H, 0 => 1)` = $H^{T}$. So its
   returned `wfl` solves the *transpose* eigenvalue equation
   $H^{T}|wfl\rangle=\lambda|wfl\rangle$, whereas dmrgpy's convention
   throughout (`nhdmrg.py`, and what its residual certificate checks) is
   the adjoint one $H^{\dagger}|\psi_L\rangle=\bar\lambda|\psi_L\rangle$.
   The two differ by exactly a complex conjugation, so
   `nh_biorthogonal_pair` returns `conj(wfl)`, renormalized (`inner()`
   conjugates its first argument, so ITensorNHDMRG's own
   `inner(wfl,wfr)=1` does not survive the conjugation, and the correct
   rescaling factor is `1/conj(s)`, not `1/s`). Confirmed directly:
   without it the left residual sits at ~2 while the right one is ~1e-15.
2. **Unanchored left solve on a real-part-degenerate spectrum.** dmrgpy
   targets the eigenvalue of smallest *real* part; when a
   complex-conjugate pair ties for it — the generic PT-symmetric /
   Hatano-Nelson situation — nothing in ITensorNHDMRG ties its left solve
   to whichever member its right solve picked, so the two can converge to
   *different* eigenvalues (observed deterministically, on every attempt
   from every random start: `conj(wfl)` came back satisfying the equation
   for `lambda` rather than `conj(lambda)`, and its overlap with `wfr`
   was exactly 0). dmrgpy's own ports anchor against precisely this
   (`arnoldi_smallest_real`'s `Sel` comment in
   `mpscpp3/chain_session.h`), but that isn't reachable from outside the
   package. `nhdmrg_solve` breaks the tie instead: it re-solves against
   `exp(i*theta)*H`, which leaves every eigen*vector* untouched and maps
   every eigen*value* to `exp(i*theta)*lambda`, so
   `Re(exp(i*theta)*(a±ib)) = a*cos(theta) ∓ b*sin(theta)` is no longer
   degenerate; the eigenvalue is mapped back exactly by dividing by the
   same phase.

   The rotation is **not** harmless on its own, and an earlier version of
   this section wrongly claimed it was. It shifts every eigenvalue's real
   part by `-Im(lambda)*sin(theta)`, which grows with `|Im lambda|` and
   therefore with system size and coupling strength: ~5e-3 at n=4 with
   gamma=0.7, but ~0.07 at n=20 — comparable to a real low-lying gap. A
   large enough `theta` can make ITensorNHDMRG converge on a genuinely
   different eigenvalue, which maps back to a perfectly valid eigenpair
   of H that simply isn't the smallest-real-part one, and *neither*
   residual in `nhdmrg.py`'s certificate can catch that (both sit at
   ~1e-15 for any true eigenpair). Two things make it safe:

   - The untied solve's own eigenvalue is kept as an **anchor**. Only its
     left vector is suspect — its right vector and eigenvalue are a
     converged eigenpair — so a tie-break result is accepted only if it
     reproduces that eigenvalue (`same_eigenvalue`, a 1e-3 relative
     separator: two runs finding the same eigenvalue have differed by up
     to ~1e-5 in measurement, two different ones by O(0.1–1)). Otherwise
     it is rejected, `nh_biorthogonal_pair` raises, and `nhdmrg.py`
     treats it as a failed attempt. Refusing to answer is the right
     outcome: the alternative is a wrong eigenvalue with a clean bill of
     health.
   - `theta` runs over a ladder, **largest first** (`TIEBREAK_ANGLES`).
     That order is measured, not obvious: the induced split has to be
     large enough for the solver to resolve, and starting at 1e-4
     reproducibly returned the right eigenvalue with a ~3e-2 residual on
     a chain where 1e-2 gives ~1e-15. The two failure modes point in
     opposite, safely asymmetric directions — too small a `theta` gives a
     poorly-converged pair that the residual certificate catches and
     retries, too large gives a well-converged pair for the wrong
     eigenvalue that only the anchor check catches.

   The mechanism runs only when the bilinear overlap actually collapses,
   so the normal path costs one extra `inner()`.

One deliberate backend difference is documented rather than forced:
`get_gs_generalized` forwards `self.noise` (DMRG's density-matrix noise
term, default 1e-1) into its `Sweeps` exactly as the session backends do
through `set_sweep_params`, since without it this backend would silently
run a noise-free variant of an algorithm whose whole point is escaping
local minima. The ITensorNHDMRG-backed paths deliberately do **not**:
noise is defined against a Hermitian density matrix, and the "fidelity"
algorithm's biorthogonal `rho=(rho_l+rho_r)/2` isometry is not that.
Feeding a noisy `Sweeps` through it was tried and measurably broke
previously-converged runs — the tie-break's anchor check began rejecting
every angle on a chain that had converged to ~1e-15 without it.

A related robustness fix landed in the shared driver: `nhdmrg()` now
treats an attempt that *raises* the same way `nhdmrg_generalized()`
already did — as a bad random draw, redrawn — and only raises when every
attempt fails. Before this, a single unlucky/degenerate draw aborted the
whole call even when the very next draw would have converged to ~1e-15.
`mpsjulialive/generalized.py::metric_guard` translates the `.jl` layer's
own collapse guards into the `RuntimeError` the other backends raise for
the same condition, matching on the guard's message rather than blanket-
translating every `JuliaError` (so a genuine bug in the `.jl` code — an
`UndefVarError`, a dispatch failure — keeps propagating as itself instead
of being silently swallowed `ntries` times by that retry loop).

#### TDVP_GSE on `julia_live`

`self.tevol_method` used to be ignored entirely on this backend --
`timedependent.py` routed every `julia_live` call to
`mpsjulialive/timedependent.py`, which only ever ran plain two-site TDVP.
It is now honored: `"TDVP"`, `"TDVP_GSE"`, and `"TEBD"` (see below) all
work, and anything else (the legacy `"MPO"` path) raises
`NotImplementedError` instead of silently running a different integrator
than the caller asked for -- the failure mode that would make a
backend-comparison script quietly compare the wrong things.

`"TDVP_GSE"` needed no algorithm port at all, unlike v3
(`mpscpp3/TDVP/basisextension.h`) and pyitensor (`pyitensor/gse.py`, a
from-scratch reimplementation of the same scheme): ITensorMPS.jl already
ships the Yang-White global subspace expansion as
`expand(psi, H; alg="global_krylov")` (its own docstring cites
arXiv:2005.06104) and its `tdvp` accepts `nsite=1`. So
`mpsjulialive/tdvp.jl`'s `gse_expand`/`tdvp_step_onesite` and the
`quench_tdvp_gse`/`evolve_and_measure_tdvp_gse` loops around them are
only the wiring, in the same "expand for the first `gse_sweeps` steps,
then one-site step" structure `pyitensor/chain.py` uses, with the same
measurement/renormalization bookkeeping as the existing plain-TDVP loops
in that file.

One genuine difference, deliberately not papered over: v3 and pyitensor
hard-cap the *enlarged* bond dimension at `maxm`
(`bond_maxdim=self.maxm` / `"MaxDim",maxm_`), whereas ITensorMPS's
`expand` takes no such argument -- what bounds growth there is the bond
dimension of the Krylov vectors themselves, via `apply_kwargs`, so that
is what `gse_expand` caps at `maxm`. A `truncate!` after the fact would
be actively wrong: the whole point of the expansion is that its newly
added directions carry exactly zero weight in the current state, so they
are precisely what any truncation would discard first.

Validated against ED on a quench *from a genuine product state* (bond
dimension 1), which is the configuration that separates "the expansion
ran" from "the call fell through to the two-site path": one-site TDVP
conserves bond dimension exactly, so with the expansion switched off
(`tdvp_gse_sweeps=0`) the same integrator is structurally unable to
follow the quench. Measured on a 6-site chain: two-site TDVP 3.1e-7 from
ED, one-site+GSE 9.4e-9, one-site with GSE off 4.9e-1. Both directions
are asserted, in `tests/test_julia_live.py` and in
`examples/time_evolution/tdvp_gse_julia_VS_ED_time_evolution`. Note this
same product-state start is where `itensor_version=3` has a known,
isolated failure (see that example's own comment and
`examples/time_evolution/tdvp_gse_VS_ED_time_evolution`, which has to mix
in an XX+YY coupling to avoid it) -- the Julia route needs no such
workaround.

#### TEBD on `julia_live`

`mpsjulialive/tebd.jl` is a third, independent from-scratch port of the
same 2nd-order-Trotter algorithm `pyitensor/tebd.py`/`mpscpp3/tebd.h`
already implement (see the main TEBD bullet above) -- not a value-level
transcription of either, since real ITensors.jl gives it different
native primitives to build on:

- **Gate construction is entirely ITensor-native, not dense-matrix-based.**
  Each per-site factor of a Hamiltonian term becomes an `ITensor` directly
  via `op(name,s)`; composing several factors at the same site uses plain
  dense-matrix multiplication on the `(dim,dim)` arrays extracted via
  `Array(op(name,s),s',s)` (the same convention `metts.jl`'s
  `metts_site_operator_matrix` already established), but *placing* a
  local or 2-site operator onto a bond uses ITensor outer products
  (`op("Id",s) * local_tensor`, `tA * tB`) rather than `np.kron`/
  `scipy.linalg.expm`-style dense reshaping -- there is no row/column
  storage-order convention to get wrong, since ITensor tracks index
  identity rather than positional array layout. The gate itself is then
  built via ITensors.jl's own tensor exponential, `exp(-im*tau*Hbond)`
  (`exp(::ITensor)`, which auto-splits primed/unprimed indices into
  "output"/"input"), the same ITensor-native mechanism `mpscpp3/tebd.h`'s
  `BondGate` uses in C++ -- unlike `pyitensor/tebd.py`, which has no
  ITensor-level `exp()` of its own and falls back to a manual
  `scipy.linalg.expm` on the dense `(D,D)` bond matrix.
- **Jordan-Wigner threading is genuinely exercised here, not a
  pass-through.** `mpsjulialive/tebd.jl`'s `term_add!`/`resolve_term` are
  a direct port of `pyitensor/autompo.py`'s `HTerm.add()`/
  `HTerm.resolve()` -- but fed the *raw* `"C"/"Cdag"` term list
  (`mpo.py`'s `text_mpo()`/`MO2list()`, the same serialization the
  existing MPO path already uses to talk to Julia), not
  `MultiOperator.to_terms()`'s Jordan-Wigner-*predressed* `"A"/"Adag"/"F"`
  form the other two backends' TEBD consume. This is a deliberate,
  confirmed-necessary choice, not a stylistic one: real ITensors.jl's
  builtin `"Fermion"` site type
  (`ITensors/src/lib/SiteTypes/src/sitetypes/fermion.jl`) only defines
  `op("C",..)`/`op("Cdag",..)`/`op("F",..)` -- there is no `"A"/"Adag"`
  operator registered for it at all, unlike dmrgpy's own hand-rolled
  pyitensor/C++ site tables, which do define `"A"/"Adag"` as their own
  name for the bare, string-free annihilation/creation operator (exactly
  what ITensors.jl's `"C"`/`"Cdag"` already are -- the Jordan-Wigner
  string is inserted separately, by the automaton or, here, by
  `term_add!`/`resolve_term`'s own carry tracking, never baked into the
  operator itself on either side). So on this backend, unlike the other
  two (where the equivalent carry-tracking pass is mostly an inert
  no-op over already-dressed input, see the main TEBD bullet's
  `_true_span()` discussion), the Jordan-Wigner carry logic does the
  actual string threading.
- **The observable/quench-operator path is untouched.** Only the
  Hamiltonian goes through `bond_hamiltonians()`'s from-scratch term
  resolution; `A`/`B` operators (`evolution_ABA`, `quench_tebd`'s
  two-state overlap) are still built as ordinary MPOs via the ITensors.jl
  `AutoMPO` automaton and applied via `tdvp.jl`'s existing `apply_clean`,
  unchanged from the TDVP path. One consequence, confirmed directly and
  not a TEBD-specific bug: applying a single *bare*, fermion-parity-odd
  operator (e.g. `Cdag[0]` alone, as `tests/test_time_evolution.py`'s own
  `evolution_ABA`-based TEBD tests do for `itensor_version` `3`/
  `"python"`) fails on `julia_live` independent of `tevol_method`, since
  ITensorMPS.jl's own `OpSum`-to-`MPO` conversion does not yet support a
  parity-odd operator sum ("Parity-odd fermionic terms not yet supported
  by OpSum to MPO conversion") -- `tests/test_julia_live.py`'s TEBD
  fermion test measures `N[2]` (parity-even) against a hopping+staggered
  quench instead, mirroring
  `test_tebd_v3_matches_python_fermion_chain`'s setup rather than
  `test_tebd_matches_ed_fermion_chain`'s.
- `psi[i] = U` / `psi[i+1] = S*V` (`apply_bond_gate!`'s SVD-truncation
  step) relies on `ITensorMPS.jl`'s own `setindex!(::MPS,::ITensor,n)`
  to keep the MPS's orthogonality-center bookkeeping (`leftlim`/
  `rightlim`) consistent automatically -- unlike `mpscpp3/tebd.h`'s
  `apply_bond_gate`, which has to set `psi.leftLim(b)`/`psi.rightLim(b+2)`
  by hand after `psi.set()`, since ITensor C++'s `MPS::set()` primitive
  does no such bookkeeping on its own.

Validated against ED on both a spin (Neel-favoring field to Heisenberg)
and a fermion (staggered field to hopping+staggered) quench, and against
a direct rejection check for a term spanning 3+ sites (`bond_hamiltonians`
raises a plain Julia `error(...)`, surfacing to Python as a
`juliacall.JuliaError` -- the same propagation mechanism
`mpsjulialive/generalized.py` already relies on for its own retry-loop
guard) -- see `tests/test_julia_live.py` and
`examples/time_evolution/tebd_julia_VS_ED_time_evolution`.

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

#### 4.8a-sector Sector-resolved correlators (`sectordc.py`, `submode="SECTOR"`)

`submode="SECTOR"` is architecturally unlike every other submode: it does
not compute a spectrum on the chain it is called on. It solves *two*
quantum-number sectors on an internal clone and contracts states from one
against states of the other.

The mechanism it rests on is a property of `promote_mps`, not of anything
new: promotion rebases a wavefunction onto the chain's **original** dense
site indices (`dense_sites_` in `mpscpp3/chain_session.h`, fixed at
construction, `self.dense_sites` in `pyitensor/chain.py`), and
`set_conserved_sector` only ever replaces `sites_`. So two states solved
in *different* sectors of one live session promote onto the same indices
and can be contracted with each other — which is exactly what a
cross-sector matrix element `<n^{N+1}|Cdag_i|0^N>` is. This is
independently regression-tested in
`tests/test_sector_promotion.py::test_states_promoted_from_different_sectors_are_comparable`,
which predates the feature.

Pipeline (`sectordc.sector_poles`), in an order that is forced rather
than chosen:

1. clone the chain — every sector switch happens on the clone, since the
   pipeline ends *outside* sector mode and that must not be a side effect
   the caller sees. `__deepcopy__` already rebuilds a fresh session and
   re-applies the sector.
2. solve the reference sector, promote its ground state;
3. solve the target sector `reference + charge(B)`, promote its `nex`
   lowest states;
4. clear the sector, and only *now* build A and B — a charge-changing
   operator cannot be built while a sector is set (`sector_terms` rejects
   it, and ITensor's AutoMPO would abort the process over it);
5. matrix elements via `aMb` on the promoted handles, then the Lehmann
   sum in `dcex.py`'s convention.

Three consequences worth knowing:

- **Everything must happen in one session.** ITensor Index identity is
  per session, so a fresh `Chain` mints indices no earlier handle
  matches. The states are therefore cached *together with the clone that
  produced them* (`_sector_states_cache`, a small dict, dropped by
  `restart()`), which is also what makes a site sweep cheap: nothing in
  the expensive half depends on which operators are correlated, so one
  pair of sector solves serves the whole `(i,j)` matrix. Measured on a
  12-site chain, a full `S(q,w)` map costs about one site's correlator.
  The cache is a dict rather than a single slot because a spectral
  function alternates between the `N+1` and `N-1` targets as it sweeps —
  a one-entry cache is evicted by every call.
- **Charge inference happens in Python, before any session is touched**
  (`multioperatortk/charge.py`). It reuses `pyitensor/sector.py`'s
  `op_charge`/`expand_xy_terms`/`combine_terms` as a *library* — an
  ungraded `SiteX` built from the chain's site codes serves any backend,
  since charge tables are a property of the site types, not of the
  grading. The Sx/Sy expansion is load-bearing: the textbook Heisenberg
  Hamiltonian conserves Sz only as a sum of terms, so the "does H
  conserve this quantity" test would answer "cannot tell" without it.
  The same module's `charge_components` runs that split the other way,
  for the *operators being correlated* rather than for H: an `Sx` with no
  definite charge is decomposed into `(S+ + S-)/2`, the components whose
  charges cancel are paired, and `sectordc.sector_lehmann` sums one
  contribution per surviving channel -- which is why `name=(Sx,Sx)`,
  the way essentially every spin correlator in `examples/` is written,
  works rather than raising. `sector_poles` remains the single-sector
  primitive and still refuses a multi-channel pair, since its `info`
  names one sector and carries that sector's matrix elements.
- **Gating is explicit in the submode, not inherited.** `mode.py`'s
  sector guard keys on `self.conserved_sector`, and the sector here lives
  on the clone — the caller's chain typically has none, so that guard
  never fires. `sectordc._check_backend` therefore rejects
  `itensor_version` 2/`julia_live`, `ns<3` on v3 and non-Hermitian
  Hamiltonians itself, `dynamics.py` rejects `"SECTOR"` on the
  `julia_live` branch, and `edtk/dynamics.py` names it too — the last one
  because `mode.py` can route a well-formed DMRG call to ED *silently*
  (extension not compiled, or v3 below 3 sites), and a generic "no ED
  implementation" would misdescribe what happened.

`spectral_function`/`spin_spectral_function` in the same module are thin
assemblies on top: two (or three) of these correlators, plus the
chemical-potential bookkeeping for the particle/hole halves.

#### 4.8a Cost of the dynamical-correlator submodes, and what makes them fast

Both `submode="KPM"` and `submode="TD"` are dominated by one inner loop
repeated hundreds of times — the Chebyshev recursion
`a_{k+1} = 2 H a_k - a_{k-1}` for KPM, the TDVP step for TD — so their
wall time is set almost entirely by the cost of a single MPO application /
local Krylov solve. Three backend-internal choices dominate that cost;
each is invisible at the API level (every number reported is unchanged),
but each is worth knowing about when reading or extending these paths.
Reference point throughout: an `L=30` `S=1/2` Heisenberg chain, `delta=0.2`,
default `maxm=30`/`kpmmaxm=50`, one BLAS thread (see §5 and
`src/dmrgpy/blasthreads.py` on why thread pinning comes first).

**Real- vs complex-valued MPOs (`mpscpp3/mo_terms.h`).** ITensor decides
whether an MPO's tensors are real or complex from the AutoMPO
*coefficients* only; the per-site operator matrices go in unexamined. `Sy`
is purely imaginary (`spinhalf.h` sets ±0.5i), so the textbook dmrgpy
Heisenberg Hamiltonian `Sx·Sx + Sy·Sy + Sz·Sz` — all coefficients real —
still produces a *complex* MPO, even though the operator it represents is
real in the S<sup>z</sup> basis. DMRG's ground state then comes out complex
too, and every subsequent `applyMPO`/`sum`/`inner` runs in complex
arithmetic. Since `Sy = (S+ − S−)/(2i)` exactly (any spin S) and a real
Hamiltonian contains `Sy` an even number of times per term, `build_ampo()`
expands `Sy` — and `Sx = (S+ + S−)/2` alongside it, which is what lets
AutoMPO see the exact `Sx·Sx + Sy·Sy = (S+·S− + S−·S+)/2` cancellation
rather than growing the MPO bond dimension — into products of the *real*
ladder operators, with a real overall coefficient. Measured at bond
dimension 50: `applyMPO` 4.4x, `sum(MPS,MPS)` 3.0x, `inner` 3.2x. The
rewrite is gated (real coefficients, even `Sy` count per term, at most four
`Sx`/`Sy` factors, and the sites must actually define `S+`/`S-`, which
keeps it away from Z3/Z4 clock operators whose complexity is genuine) and
can be turned off process-wide with
`mpscpp3._dmrgcpp.set_realify_spin_terms(False)` — used by
`tests/test_dynamical_correlator.py` to check both paths agree. It is
`mpscpp3`-only: `mpscpp2` is untouched, and `pyitensor` stores every tensor
as `complex128` by design (`pyitensor/tensor.py`), so there is nothing for
it to win there.

**Krylov convergence in the pure-Python TDVP (`pyitensor/tdvp.py`).**
`_lanczos_expm_multiply` builds a Krylov subspace up to `niter` (its
callers pass 50, matching the compiled backend's `MaxIter`) and
exponentiates the projected tridiagonal matrix. `niter` is an upper
*bound*, not the iteration count: the recursion stops as soon as Saad's a
posteriori residual estimate `‖v‖·β_k·|e_k^T exp(coeff·T_k) e_1|` falls
below 1e-10 — the same convergence target ITensor's own `applyExp()`
(`iterativesolvers.h`) uses via its Expokit-style extended-T-matrix
correction. The first column of `exp(coeff·T_k)` is needed to assemble the
result anyway, so the test costs one extra `eigh_tridiagonal` of a k×k
matrix per iteration against an O(D³) MPS-level matvec. Typical
convergence is under ten iterations, so running the full 50 unconditionally
was the single largest cost in the pure-Python `submode="TD"` path.

**Truncation via the smaller Gram matrix (`pyitensor/svd.py`).** Every
truncation in the pure-Python engine is "reshape across an index partition,
factorize, keep the largest singular values". Those matrices are strongly
rectangular — MPS bond dimension × physical dimension on one side, MPS bond
dimension × MPO bond dimension on the other — and at most `maxdim` vectors
are ever kept. `_svd_truncated()` therefore eigendecomposes the smaller of
`M M†` / `M† M` and recovers the kept factor on the other side with one
matmul, the same trick ITensor's own `densityMatrixApplyMPOImpl` uses. It
falls back to the exact SVD whenever the smallest *kept* singular value
drops below 1e-7 of the largest, where the squared spectrum stops resolving
it — so this is a speedup, never a silent accuracy trade.

Net effect at the reference point (before → after):

| backend | KPM | TD |
| --- | --- | --- |
| `itensor_version=3` | 67.4 s → 19.6 s (3.4x) | 128.5 s → 110.2 s (1.17x) |
| `itensor_version="python"` | 57.4 s → 37.4 s (1.5x) | 1081.7 s → 156.0 s (6.9x) |

The modest C++ `TD` figure is expected rather than a missed opportunity:
`exp(-iHt)` makes the evolved state genuinely complex whatever the MPO's
type, so only the H-side of each contraction gets cheaper there.

#### 4.8a-bis Time-grid convention: measure before evolving

Every real-time evolution loop in the library — `Chain::quench*` /
`Chain::evolve_and_measure*` in `mpscpp2/` and `mpscpp3/`, their
counterparts in `pyitensor/chain.py`, `mpsjulialive/{tdvp,tebd}.jl`, and
`edtk/timedependent.py`'s `evolution_ABC`/`evolution_DC` — records its
observable at the **top** of the step loop. `correlator[k]` is therefore
the value at `t = k*dt`, matching the `ts = [0, dt, ..., (nt-1)*dt]` grid
`timedependent.py` returns alongside it.

This was previously the other way round (evolve, then measure), which
made `correlator[k]` the value at `(k+1)*dt` while still labelling it
`k*dt`. The consequences were real, not cosmetic:

- a spurious `exp(i*omega*dt)` phase on the Fourier transform, and
- more seriously, `C(0)` — the largest sample — was omitted from the
  Riemann sum in `_fourier_transform_correlator` altogether, leaving an
  `O(dt*C(0))` error in every `submode="TD"` spectrum that did **not**
  vanish as `nt` grew at fixed total time.

Confirmed on an L=10 Heisenberg chain: `correlator[0]` came back as
0.2409 / 0.2477 / 0.2494 for `dt` = 0.2 / 0.1 / 0.05 against an exact
`C(0) = <A B> = 0.25`, its imaginary part halving exactly with `dt`; and
the `submode="TD"` spectral weight moved 75% across `dt` = 0.1 / 0.05 /
0.025 (−0.105 / −0.159 / −0.184) for a `C(t)` that is itself
`dt`-independent to 5e-5. After the fix the same sweep gives −0.2047 /
−0.2061 / −0.2065 — converging rather than drifting — and
`correlator[0]` reproduces `C(0)` to 1e-13 on `itensor_version` 2, 3,
`"python"` and `mode="ED"` alike.

All backends were changed together, ED included, so every DMRG-vs-ED
cross-check in `tests/` and `examples/` still compares like with like.
`submode="TD"` results from before this change are not comparable with
results from after it.

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
and the feature's other scope limitations (single tip-coupled site, the
potential-scattering interference term's general-spin form is an
extrapolation from the paper's own worked example) -- as well as the
normalization conventions, which every term's prefactor is checked
against via the paper's absolutely scaled Figs. 3b/3d.

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

The potential-interference term (`U!=0`, part of `order=3`) is also
supported for `mode="DMRG"`, via `potentialdc.py`: its own `T=0` limit
collapses the excited-state sum to a convolution of the *same* `T=0`
dynamical structure factor against the `F0` kernel instead of `Theta0`'s
cumulative-sum weighting (`Theta0(eV-x)` is a step, integrated
cumulatively; `F0` is not, so this piece genuinely convolves rather than
cumulatively sums) -- so it reuses `get_dynamical_correlator` exactly
like the second-order term above, needing no excited-state enumeration
either, and no new DMRG-side primitive. Verified (`mode="ED"`,
`submode="ED"`) against the excited-state-sum
`conductance.third_order_potential_dIdV` to ~0.07% relative error
(`tests/test_kondo_spectrum_potentialdc.py`); a real `mode="DMRG"`,
`submode="KPM"` run through the full public API was tried directly but
found impractically expensive for a routine test -- a single KPM
dynamical-correlator call took ~40s on the development machine even with
a trivially small number of moments, confirmed directly -- so the
`Spin_Chain._get_kondo_spectrum_dmrg` wiring that combines this term
with the other two is instead checked with the three underlying
(expensive) calls monkeypatched out
(`tests/test_kondo_spectrum.py::test_get_kondo_spectrum_dmrg_combines_terms_correctly`).

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

### 4.10 Where "silently wrong" was structurally possible

The 2026-08 audit (`docs/audit_2026_08_hole_hunt.md`) found that most of
this codebase's silent-wrong-answer bugs shared one shape: **a dispatch
decision taken before the information that should inform it**. The
architecture now guards each of those points, and they are worth knowing
as a class, because new code can reintroduce any of them.

- **A short circuit ahead of a cache key.** `gs_energy()`/`get_gs()`
  returned on `computed_gs` alone, before `groundstate.gs_energy_single()`
  and its send-cache -- which correctly keys on
  `maxm`/`nsweeps`/`cutoff`/`noise`/ramp/sector -- was ever entered, so a
  convergence ramp over `sc.maxm` on one chain returned the first energy
  every time. The key now lives in `groundstate.solver_key()`, is recorded
  on the chain when a state is stored, and `groundstate.gs_is_current()`
  is what the short circuit tests. A state with no recorded key (one
  injected via `set_gs()`) is still returned unconditionally, by design.

- **A precondition tested ahead of the dispatch it should qualify.**
  `dynamics.get_dynamical_correlator()` tested Hermiticity *before*
  branching on `submode`, so a non-Hermitian Hamiltonian silently routed
  every submode to the explicit CVM resolvent. The check is now
  per-submode: KPM and CVM/CVM\_explicit have genuine non-Hermitian
  implementations, EX/maxent are backend-agnostic and pass through, and
  everything else raises. `edtk/dynamics.py` mirrors it. The
  `julia_live` branch had already been fixed this way and is what the
  others now follow.

- **An `else` branch doing two jobs.** Each `tevol_method` dispatch chain
  ended in a bare `else` that ran the legacy MPO-Taylor integrator. That
  branch is the documented fallback for a backend without TDVP
  (`itensor_version=2`), and it was also catching every misspelling.
  `timedependent.check_tevol_method()` now validates the *name* up front
  at all three dispatch sites (`evolution_dmrg_DC`,
  `evolve_and_measure_dmrg`, `tdz.py`); backend capability still falls
  back exactly as before.

- **`**kwargs` with no consumer.** `kpmdmrg.dynamical_correlator_moments`
  absorbed every unknown keyword; its sibling submodes always raised.
  Solver parameters are chain attributes, never call arguments, so
  anything unrecognized now raises `TypeError`.

- **Resolution done in the branches instead of before them.** The
  documented string form of `name=` (`"ZZ"`, `"cdc"`, ... plus `i=`/`j=`)
  was understood only by whichever submodes happened to call
  `operatornames.str2MO` themselves. It is now resolved once in
  `Many_Body_Chain.get_dynamical_correlator`, which is also the only
  object that *can* resolve it -- an `EDchain` has no `Sx`/`Sz`/`C`
  attributes of its own.

Two `mpscpp3` entries belong to the same audit and are specific to
QN-conserving (`set_conserved_sector`) mode:

- `Chain::reduced_dm` flattened the density matrix with `rho.visit()`,
  which enumerates only the **stored** blocks; with QN-carrying indices
  that is a handful of 1x1 blocks rather than the full `dim*dim`, and
  `bindings.cc` copied that short vector into an uninitialized array. It
  now reads the elements explicitly by index (`eltC(rho,s(a),prime(s)(b))`,
  the convention `op_charge()` documents), and the binding size-checks.
- `Chain::sector_terms` validated each term's **net** charge only, so a
  net-neutral term could still pile an impossible charge onto a single
  site (`Cdag_0 C_1 Cdag_0 C_1`, which the four-point correlator produces
  routinely) and abort the process inside AutoMPO. It now also composes
  each site's own factors: an identically-zero product drops the term
  (matching ED), and a charge the site cannot carry raises
  `std::invalid_argument` -- catchable, unlike ITensor's own `Error()`.

## 5. Backend performance: v3 vs the pure-Python backend

The pure-Python backend (`itensor_version="python"`) trades raw speed for
zero build dependencies. The tables below were originally measured with
`numba` explicitly opted in (`pyitensor.kernels.USE_NUMBA = True`); the
*default* path (`USE_NUMBA = False`) has since been rewritten to
precompute its own contraction plan once per bond/Lanczos-or-Krylov loop
rather than re-deriving it on every matvec call (see `pyitensor/
kernels.py`'s module docstring, "plain path" section), so it now performs
comparably to these numba-opted-in numbers for ground-state DMRG (spot-
checked directly: within ~10% of the §5.1 table below at n=8-20) without
any opt-in at all. `USE_NUMBA` is not just unnecessary now but frequently
counter-productive: re-measured directly, an independent-chain workload
(each chain starting from its own random MPS, e.g. a parameter sweep)
never reaches a stable steady state with numba on, because each chain's
own random start produces a different bond-growth trajectory and hence
different contraction shapes to compile every time -- see `kernels.py`'s
own docstring for the full, current measurements and guidance (numba is
now only ever worth opting into for a single long-running TDVP/METTS-like
session, and even there the win is modest, ~10%, not a large one).
Absolute times below will vary by machine and load, but the qualitative
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
| Hubbard ground state (3 sites) | 0.04s | 0.11s | 3.1x |

TDVP-based workloads (real-time evolution, and METTS finite-temperature
sampling in particular) are a materially different, larger gap than the
ground-state/KPM numbers above, because they call the same effective-
Hamiltonian matvec thousands of times more per run (many time steps x
many Lanczos/Krylov iterations x many METTS samples, vs. a handful of
DMRG sweeps) -- see `examples/finite_temperature/
backend_timing_metts_dynamical_correlator/main.py`, which measures this
directly (dynamical METTS, n=10: v3 ~9s vs. python ~51-59s for
nsamples=10-40, a 5-6x gap after the matvec-planning and Lanczos-
bookkeeping fixes referenced above -- it was ~13x before them).

### 5.4 Practical guidance

- Use `itensor_version="python"` for quick prototyping, CI/environments
  without a C++ toolchain, or small systems (n ≲ 20) where the 1.3–2x
  overhead doesn't matter.
- `USE_NUMBA`/`USE_JAX` (both off by default, see §5's note above) are
  no longer needed to get the numbers above -- the default path already
  precomputes its own contraction plan. Opting into `USE_NUMBA` is only
  ever worth it for a single long-running TDVP/METTS-like session (same
  matvec shapes recurring many times within one run), and even there the
  win is modest (~10%); it's actively counter-productive for independent-
  chain workloads like parameter sweeps (see `kernels.py`'s own
  docstring for the measurements behind both claims).
- For production-size chains, dynamical correlators, TDVP/METTS
  workloads, or anything performance-sensitive, compile the C++
  extension (`python install.py`) and use `itensor_version=2` or `3`.
- `examples/groundstate/backend_timing_gs_energy/main.py` reproduces the n=8..20
  rows of the §5.1 table directly against the current codebase and
  current machine (n=24..32 aren't included, to keep the example fast to
  run); `examples/finite_temperature/
  backend_timing_metts_dynamical_correlator/main.py` does the same for
  the TDVP/METTS gap discussed in §5.3 — re-run both rather than trusting
  these numbers verbatim on a different setup.


## GPU: the `pyitensor` array backend

`itensor_version="python"` has a second, orthogonal axis: *which array
library* its ITensors are made of. `pyitensor/backend.py` selects it
process-wide.

```python
from dmrgpy.pyitensor import backend
backend.set_backend("jax")     # or "numpy" (the default)
print(backend.device_info())   # e.g. "jax: cuda:0"
```

This is a device axis, not a new backend: `set_backend` changes nothing
about dispatch (`mode.py` still decides DMRG vs ED, `cppext` still decides
which compiled extension backs `itensor_version in (2,3)`), only what
`ITensor.array` is. The compiled C++ backends and the ED path are
unaffected.

The design rests on one measured constraint: **arrays must stay on the
device**. A single host->device->host round trip for one two-site tensor
costs *more* than the contraction it feeds above modest bond dimension
(22.0 ms vs 8.3 ms at chi=1024 on an H200), so any design that converts
per call loses however fast the device is. Consequently:

* there is exactly one conversion point in the engine, `ITensor.__init__`
  (via `backend.asarray`), and everything built from a converted array
  stays put;
* host round trips are confined to three O(1)/O(chi) places -- the
  Lanczos/Krylov alpha/beta coefficients (`dmrg.py`'s ground-state solver
  and `tdvp.py`'s time propagator, whose k x k projected matrices and
  convergence tests are host work by design), the singular-value vector
  `svd.py`'s truncation rule branches on, and measurement results;
* almost everything else needed no namespace at all, because the
  operations this engine performs on tensor data are ndarray *methods*
  (`.transpose`, `.reshape`, `.conj`, `@`), identical on NumPy and JAX.

Three places did need explicit work, and each is a trap worth knowing:

* **`np.transpose(arr, perm)` and friends silently transfer.** JAX arrays
  implement no `__array_function__`, so a NumPy free function falls back
  to `__array__` and copies the whole tensor to the host -- per call, with
  no error. `kernels.py`'s planned matvec used exactly this. Use methods.
* **JAX arrays are immutable**, so `mpsalgebra.py`'s direct-sum block
  write goes through `backend.setblock`, which mutates in place on NumPy
  and returns a new array on JAX. Callers must use the return value.
* **XLA compiles a kernel per (operation, shape)**, and DMRG mints a new
  shape whenever a bond dimension changes, so a cold run pays a large
  one-time compile tax (measured: 672 compilations, 18.4 s of a 29.1 s
  run). `backend.set_pad_bonds(K)` freezes every bond at K to collapse
  that shape zoo; it appends zero singular values after truncation, so
  the represented state is unchanged.

**The dispatch floor, and the two knobs against it.** Eager dispatch has
a per-call floor (~0.35 ms on an H200, against ~0.07 ms for the same
kernel under `jax.jit`), which is why the device cannot win below
chi ~ 120 whatever the hardware. `backend.jit(fn, static_argnums)`
registers a *hot composite* -- a function whose array arguments are
traced and whose shapes are static -- and returns a wrapper that runs it
under `jax.jit` when jitting is on and calls it directly otherwise, so
the NumPy path is the plain Python function and does not change. Four
composites are registered, covering everything the engine does per
Lanczos iteration:

| composite | eager dispatches | where |
|---|---|---|
| planned matvec chain | 3-4 per step, 8-12 per call | `kernels.py::_matvec_chain` |
| contraction (transpose+reshape+matmul) | 5 | `tensor.py::_contract_matmul` |
| Gram matrix + `eigh` + descending sort | 5 | `svd.py::_gram_spectrum` |
| Lanczos recurrence / block reorthogonalization | 3 / 2 per basis vector | `dmrg.py` |

`backend.set_jit()` decides whether they compile: `"auto"` (the default)
turns jitting on exactly when `set_pad_bonds` is set, `True`/`False`
force it. The tie to padding is not a convenience -- `jax.jit` traces
once per input *shape*, and DMRG mints a new shape at every bond
dimension (the same trap `kernels.py` documents for numba), so without
padding the compile count can rise instead of falling. `backend.
compilations()` reports how many kernels the registered composites have
traced, which is the number padding exists to collapse.

`make_matvec` additionally refuses its two opt-in host paths (the numba
contractor, the fused-einsum one) whenever `backend.is_device()`: both
convert per call, so either would reinstate exactly the round trip the
port exists to remove.

Dispatch is not the only device cost the engine has to manage;
*synchronization* is the other, and it binds hardest on the long, small-
operation loops that time evolution is made of. JAX dispatch is
asynchronous, so the host runs ahead enqueueing kernels while the device
works -- and every `bk.to_host()` drains that queue. Two places therefore
deliberately trade extra arithmetic for fewer stalls, both device-only so
the NumPy path is untouched:

* `tdvp.py`'s Krylov exponentiator keeps its alpha/beta on the device and
  brings a block of them home per checkpoint instead of two per
  iteration, *speculating* past its own stopping point and then rolling
  back to the exact Krylov dimension the per-iteration loop would have
  chosen -- so the returned vector is identical, and only a few matvecs
  are wasted. `tdvp.set_krylov_defer_sync()` forces or disables it.
* `mpsalgebra.BatchedBras` prepares a fixed family of bra MPS once
  (transposed, conjugated, zero-padded to a common per-bond width) so
  that many overlaps against one evolving ket come from a single batched
  sweep. `tdz.py`'s complex-time correlator uses it for its n_max+1
  phi^(n) overlaps per step; `tdz.set_batched_bras()` switches it off.
  Backends without `Chain.batched_bras` fall back to a loop over
  `overlap()`.

`svd()` follows the same rule: only the O(chi) spectrum comes home (the
truncation rule is data-dependent branching, which belongs on a host), it
comes home once rather than twice, and the diagonal S tensor is built from
the device-resident spectrum rather than from `np.diag()` of the host copy
-- the latter was an O(chi^2) transfer per factorization, i.e. per bond
per half-sweep per sweep.

Sampling methods inherit all of the above rather than needing their own
port. METTS (`pyitensor/metts.py`) is the case worth stating, because
what remained after TDVP, `dmrg.py`'s environments, `svd.py` and
`mpsalgebra.inner` were resident was a single loop: the collapse sweep
that turns each evolved sample back into a product state. It is O(n)
small operations per sample against thousands inside the evolution, but
it was host-bound in the silent way -- `np.sum`/`np.abs` applied to a
device array pull the whole (d, chi) amplitude block home, once per site,
per collapse, per sample, and nothing raises. The sweep now runs in the
active namespace and only the O(d) conditional-probability vector comes
back, because the multinomial draw it feeds is data-dependent branching
and belongs on a host; the collapsed-prefix amplitude stays resident.

One consequence is user-visible: `metts_vev`/`metts_dynamical_correlator`
accept an `njobs>1` that fans independent Markov chains over worker
processes, and those workers are started with `spawn` -- so each
re-imports `pyitensor` with the default NumPy backend instead of
inheriting the parent's device. That combination now raises rather than
running, on the same principle as `set_backend` itself: a GPU run that
quietly became a CPU run is worse than a failure.

Measured speedups, the crossover, the primitive-by-primitive table and the
benchmarking pitfalls are in `docs/gpu_cpu_performance.md`; the port's
design history is in `docs/pyitensor_gpu_port_plan.md`.

## 6. Directory reference

| Path | Contents |
|---|---|
| `src/dmrgpy/` | the Python library |
| `src/dmrgpy/manybodychain.py` | `Many_Body_Chain`, the shared base class |
| `src/dmrgpy/infinitechain.py`, `pyitensor/idmrg.py`, `mpscpp3/chain_session.h` | `Infinite_Many_Body_Chain` and infinite DMRG (iDMRG), pyitensor + ITensor v3 (energy density, static correlators and local gap on both) |
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
| `examples/backend_comparison/`, `v2_VS_v3_*` examples within each theme folder | backend-vs-backend correctness comparisons |
| `examples/groundstate/backend_timing_gs_energy/` | backend-vs-backend timing comparison |
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
