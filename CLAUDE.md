# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

DMRGPY is a Python library for quasi-one-dimensional spin chains, fermionic
systems, parafermion and bosonic models using matrix product states (MPS).
It exposes a single Python API (`src/dmrgpy/`) backed by interchangeable
solvers per calculation:

- **DMRG**, executed in-process via a pybind11 extension over ITensor
  (C++, either v2 or v3 -- see below), or via a separate Julia/ITensors.jl
  backend
- **ED** (exact diagonalization), a pure-Python fallback used for small
  systems, for cross-checking DMRG results, and automatically whenever the
  C++ extension isn't compiled

The Python layer is what this project actually develops. The C++/ITensor
backends under `src/dmrgpy/mpscpp2/ITensor/` (ITensor v2) and
`src/dmrgpy/mpscpp3/ITensor/` (ITensor v3) are vendored copies of the
ITensor library maintained upstream — treat them as black boxes, not code
to read or modify here.

## Install / build

There is no `setup.py`/`pyproject.toml`; the project is used from `src/`
via a path added to `PYTHONPATH` (or a symlink into `site-packages`).

```bash
python install.py                    # compiles ITensor v2 (mpscpp2) + its pybind11 extension
python install.py --itensor-version=3    # compiles ITensor v3 (mpscpp3) instead
python install.py --itensor-version=both # compiles both v2 and v3 backends
python install.py --gpp=g++-6        # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install_julia.py              # alternative: Julia backend instead of/alongside C++
```

`install.py` runs in two phases: `installtk/requirements.py` first checks
every build requirement -- OS, `make`, a C++ compiler, LAPACK/BLAS,
`pybind11` -- actually trial-compiling/trial-linking each one (not just
checking a version string), and only once everything checks out does
`installtk/install2.py` compile the vendored ITensor static library and
the `pybind` Makefile target for whichever backend(s) `--itensor-version`
selects (`2` = mpscpp2/ITensor v2, the default; `3` = mpscpp3/ITensor v3;
`both` builds both, one after the other).`install2.py` picks the right
`mpscppN` directory and C++ language standard per version (v2 needs only
C++14; v3 needs C++17 with the concepts TS extension, `-fconcepts`, see
`mpscpp3/ITensor/options.mk.sample`'s own `CCCOM` line) but otherwise runs
the identical compiler/BLAS/LAPACK configuration resolved in phase 1.
`pybind11` is a hard requirement (auto-installed via `pip` into the
current interpreter if missing, since the pybind11 extension is the only
DMRG backend besides ED). The compiler is auto-detected rather than
defaulting to plain `g++`, via a cascade of candidates rather than a
single guess: `installtk/requirements.py::_compiler_candidates()` yields
a conda-provided compiler first when running under a conda Python
(`installtk/cppversion.py::find_conda_compiler()` looks for one, e.g. the
`gxx_linux-64` package, next to the running interpreter -- preferred
because the extension is loaded into the same process as conda's own
numpy/scipy, which bundle their own libstdc++, and building against the
system compiler instead has reproducibly segfaulted in the past, see the
long comment in `mpscpp2/Makefile`), then falls back to the system
`g++`/`c++` on PATH if that candidate doesn't actually work (`--gpp`
overrides the whole cascade with a single hard choice). Each candidate is
verified for real, not just checked for existence: besides an actual
trial compile+link (`cppversion.trial_compile()`),
`cppversion.backend_missing()` proactively asks the driver
(`gpp -print-prog-name=cc1plus`) whether it can even locate its own
compilation backend, since `--version`/`shutil.which()` both pass for a
compiler binary that's present but incomplete. This isn't hypothetical:
confirmed directly on Aalto's Triton cluster, whose `scicomp-python-env`
module ships a conda-style Python distribution with exactly such a
broken `x86_64-conda-linux-gnu-g++` wrapper (no matching `cc1plus`
installed) -- the cascade transparently falls back to Triton's login-node
system compiler in that case and prints a note about it.
`installtk/cppversion.py::lmod_state()` detects environment-module
systems (Lmod, via `$LMOD_CMD`/`$MODULEPATH`) so failure messages can
suggest `module load ...` instead of the sudo-based apt-get/brew hints
that are meaningless (no sudo) or simply wrong on a cluster; `check()`
also prints a reminder to reload the same modules/conda environment used
at install time in any later session that imports `dmrgpy`, since the
compiled extension depends on them being present again at import time.
`install.py --doctor` runs only this requirement-checking phase and
exits, without starting the multi-minute ITensor build -- useful for
diagnosing a cluster environment or reporting a bug.
`installtk/blaslapack.py::find_working_config()`
similarly tries several LAPACK/BLAS candidates (conda's own libs, system
`-lblas -llapack`, OpenBLAS, and a Debian/Ubuntu-multiarch workaround that
asks the system compiler where it actually finds `libblas.so`/
`liblapack.so`) and picks the first one that actually links, rather than
assuming a fixed default. Finally, `install.py` adds the repo's `src` to
the Python path and to the user's shell rc file
(`installtk/addpythonpath.py`, `installtk/addsystem.py`, both idempotent
across reruns). There is no separate lint or CI config in this repo.

## Running / verifying changes

There is no unit test suite for the Python code (the only `test.py`/`*_test.cc`
files live under `pychain/` and inside the vendored ITensor tree). The
`examples/` directory (100+ self-contained scripts, one per physical model or
feature) is the de facto regression suite — to check that a change works,
run the relevant example(s) directly, e.g.:

```bash
cd examples/static_correlator_S12 && python <script>.py
```

`examples/readme_examples/` mirrors the snippets shown in `README.md`.
`examples/v2_VS_v3_*` and `examples/dynamical_correlator/
dynamical_correlator_v2_VS_v3` directly compare the ITensor v2 and v3 C++
backends against each other (same script, both `itensor_version`s, small
systems) — the fastest way to check a `mpscpp3` change didn't diverge from
`mpscpp2`'s numerics. After running examples, `python clean.py` recursively
removes generated working directories (`.mpsfolder`, `.pychainfolder`,
`.dmrgfolder`) and stray `ERROR`/`*.OUT` files from the tree.

## Architecture

### Entry point classes

Each physical model is a thin subclass of `Many_Body_Chain`
(`manybodychain.py`), which holds all shared state (bond dimensions,
sweep count, cutoffs, mode, temp folder, etc.) and implements the generic
operations (`set_hamiltonian`, `gs_energy`, `vev`, `get_gs`, correlators,
time evolution, entanglement, ...). Model modules just define the local
Hilbert space and convenience operators:

- `spinchain.py` → `Spin_Chain` (aliased directly to `Many_Body_Chain`)
- `fermionchain.py` → `Fermionic_Chain` (spinless fermions, `C`/`Cdag`/`N`/Fermi-string `F`)
- `spinfermionchain.py`, `bosonchain.py`, `parafermionchain.py`, `kondochain.py` — same pattern for other statistics/models

### Operator representation: `MultiOperator`

Hamiltonians and observables are built as `MultiOperator` objects
(`multioperator.py`) — sums of products of named single-site operators
(e.g. `Sx[i]*Sx[j]`). This is a backend-agnostic intermediate
representation: the *same* `MultiOperator` is later either

- written out as an AMPO-style operator file for the C++/Julia DMRG backend, or
- converted to a sparse matrix via `multioperator.MO2matrix` for the ED backend (`edtk/edchain.py`).

`multioperatortk/` holds supporting machinery (Jordan-Wigner strings for
fermionic operators, static/long-range operator construction, sympy-based
symbolic building).

### Backend dispatch (DMRG vs ED, C++ vs Julia)

`mode.py` decides, per call, whether a calculation runs via DMRG or ED
(`get_mode`/`run`): DMRG unless `self.mode` forces ED, or unless
`self.itensor_version` is `2` or `3` and the corresponding pybind11
extension isn't compiled (see `cppext.available(version)`), in which case
it silently falls back to ED. Most public methods on `Many_Body_Chain`
accept a `mode="DMRG"|"ED"` kwarg so results can be cross-validated
between solvers (see the bilinear-biquadratic example in `README.md`).

- **DMRG, C++ (`itensor_version=2` (the default) or `itensor_version=3`)**:
  entirely in-process, see "In-process pybind11 extension" below. Both
  versions expose the identical `Chain` session API and public
  `Many_Body_Chain` surface — `itensor_version` only picks which compiled
  extension (`mpscpp2._dmrgcpp` vs `mpscpp3._dmrgcpp`, via
  `cppext.get_backend(version)`) backs `self._session`.
  `Many_Body_Chain.setup_cpp(version=2)` switches an existing chain
  between them. There is no file-based/subprocess fallback for either —
  if the requested extension isn't compiled, `mode.py` routes to ED
  instead.
- **DMRG, Julia (`itensor_version="julia_live"`)**: a live, in-process
  Julia session (`mpsjulialive/`, via `pyjulia`/`juliasession.py`) with its
  own parallel set of modules (`mpsjulialive/groundstate.py`,
  `vev.py`, `mps.py`, `mpo.py`, ...) mirroring the top-level ones but
  talking to Julia's ITensors.jl instead of the C++ pybind11 extension.
  `Many_Body_Chain.setup_julia()` switches a chain to this mode.
  `itensor_version="julia"` (a separate, older subprocess-based Julia path
  via `juliarun.py`) is not reachable through the normal public API today
  (`gs_energy()` and friends only special-case `2` and `"julia_live"`)
  and should be treated as legacy/inert.
- **ED path**: `edtk/edchain.py` (`EDchain`) builds dense/sparse operators
  directly in Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`,
  `pyzn/` provide the per-statistics many-body operator construction,
  `edtk/one2many.py` promotes single-site operators to the full Hilbert
  space) and diagonalizes with `scipy.sparse.linalg`.

### In-process pybind11 extension (the C++ DMRG backend: v2 and v3)

`mpscpp2/bindings.cc` and `mpscpp3/bindings.cc` each compile a pybind11
extension module (`mpscpp2/_dmrgcpp*.so` / `mpscpp3/_dmrgcpp*.so`, built by
the `pybind` target in each backend's own `Makefile`) that runs DMRG
entirely in-process against its own vendored ITensor copy — no task files,
no operator files, no subprocess, no temp directory. `mpscpp3` is a
line-by-line port of `mpscpp2` to the ITensor v3 API (IQTensor/IQIndex
merged into a single ITensor/Index type with QN blocks living on the
Index; see the porting-notes comment at the top of `mpscpp3/chain_session.h`
for the full v2→v3 API mapping) exposing the *identical* `Chain`
class/Python surface, not an independent reimplementation. There is no
file-based/subprocess design for either version — that was removed
entirely from `mpscpp2`, and `mpscpp3` never had one.

- **`cppext.py`** lazily imports whichever compiled extension a chain
  needs and caches whether it loaded (`get_backend(version)`/
  `available(version)`, `version` is `2` or `3`).
- **`mpscppN/chain_session.h`**'s `Chain` class is the stateful, in-process
  session: one `Chain` instance per `Many_Body_Chain` (held as
  `self._session`, created in `sites.py::initialize` whenever
  `itensor_version in (2,3)` and the corresponding extension is
  available), with a method per DMRG task (`gs_energy`, `vev`,
  `apply_operator`, KPM/CVM dynamical correlators, `correlation_matrix`,
  `reduced_dm`, `excited_states`, time evolution, ...), built directly
  against ITensor's `AutoMPO`/`MPS`/`MPO` APIs. `check_task.h`,
  `get_sites.h` and `mo_terms.h` are the only other `mpscppN/*.h` headers
  it depends on; everything else it needs is self-contained in
  `chain_session.h`.
- **`mpscppN/mo_terms.h`** and `MultiOperator.to_terms()`
  (`multioperator.py`) convert a `MultiOperator` into a plain list of
  `(coefficient, [(opname, site), ...])` tuples that `Chain` turns directly
  into an ITensor `HTerm`/`AutoMPO`.
- **The `cpp_handle` pattern**: `mps.py::MPS` and
  `multioperatortk/staticoperator.py::StaticOperator` both carry an opaque
  `cpp_handle` attribute holding the pybind11-wrapped ITensor object living
  in C++ memory — these classes have no other representation now, so
  `cpp_handle` is always set for `itensor_version in (2,3)` chains (it's
  `None` only for objects that were never meant to touch the C++ backend,
  e.g. `mpsjulialive`'s own `MPS`/`MPO` classes are separate types
  entirely).
- Because `cpp_handle` objects have no pickle/deepcopy support, `MPS.copy()`
  / `StaticOperator.copy()` use a shallow `copy.copy()` (never
  `deepcopy()`). `Many_Body_Chain.__deepcopy__` never copies `_session`
  directly (a live C++ session has no well-defined "copy" semantics);
  instead it builds a fresh `Chain` for the clone when the original had
  one, which is what `clone()`/`bandwidth()`/`lowest_eigenvalue()` rely on.
- A few **pre-existing bugs in the original design are deliberately
  reproduced, not fixed**, in both `chain_session.h`s (see comments at the
  call sites for details): `evoloperator()`'s z³/6 term multiplies `H2`
  instead of the computed `H3`; the old file-based backend's `"noise"`
  tasks.in key was actually read back as `"moise"`, and `"tevol_fit"` was
  written but `"tevol_fit_td"` was read, making that fitting branch
  unreachable — both reproduced as hardcoded values rather than silently
  "fixed". If ITensor results look numerically "close but not exact" to a
  hand-derivation, check here before assuming a porting bug.
- **`mpscpp3`-specific: every site is built with `ConserveQNs=false`**
  (`mpscpp3/get_sites.h`), and `Chain`'s DMRG starting state is an actual
  `randomMPS(sites_,maxm_)`, not a plain product state
  (`chain_session.h`'s `default_mps()`). This isn't cosmetic: v2's own
  plain `MPS(sites_)` (no `InitState`) ends up, in practice, performing an
  *unconstrained* search of the full Hilbert space rather than one pinned
  to a single total-Sz/particle-number sector — confirmed by cross-checking
  against the compiled v2 extension with a symmetry-breaking field term,
  where v2 correctly finds a fully-polarized ground state. ITensor v3 is
  stricter about QN-conserving tensors needing a well-defined flux from
  the start, and a naive port (QN on, product-state start) gets DMRG
  permanently stuck at that starting product state's energy (a textbook
  "DMRG trapped at an exact eigenstate" failure, reproducible even with a
  large noise term) instead of raising an error. Turning QN off and
  starting from an actual random MPS reproduces v2's real behavior at the
  cost of the QN block-sparsity speedup for `itensor_version=3` — a
  genuine perf/memory tradeoff of this backend, not a bug to "fix" later.
- **`mpscpp3`-specific: `applyMPO()` doesn't restore the standard unprimed
  physical index** the way v2's `exactApplyMPO`/`fitApplyMPO` implicitly
  did — it just contracts whichever leg of the MPO happens to match the
  input MPS's own (possibly already-primed) physical index, so the *next*
  `sum()` of two independently-produced results can hit a hard
  `ITensor::operator+=: different index structure` abort. Every MPO
  application in `mpscpp3/chain_session.h` goes through a private
  `apply_mpo()` wrapper that always calls `.noPrime("Site")` on the result
  before returning, restoring that invariant globally instead of patching
  individual call sites. Likewise, `nmultMPO(A,B)` (used to build `H*H` for
  `custom_exp()`/`evoloperator()`) requires `nmultMPO(A,prime(B))` when `A`
  and `B` share the same physical indices (any two operators built on the
  same chain always do) — v2 tolerated the unprimed form directly — and the
  result's output leg then needs `.mapPrime(2,1)` to get back to the
  standard single-application convention before it can be summed against
  any other operator (see `mult_mpo()`'s comment, and ITensor v3's own
  `unittest/mpo_test.cc` "nmultMPO" section, which documents exactly this
  usage).
- **`mpscpp3`-specific: plain `inner()` throws on complex MPS/MPO** ("Cannot
  use inner(...) with complex MPS/MPO, use innerC(...) instead"), unlike
  v2's `overlap()` (which just silently returns the real part regardless).
  Since `MOTerm::coef` is always `Cplx`, any MPO/MPS built from
  Python-supplied terms is complex-typed even when every coefficient's
  imaginary part is exactly 0 — confirmed directly, the KPM dynamical
  correlator path aborted here in `same_mps()` until every "Re[<x|y>]" use
  in `chain_session.h` was changed to `innerC(...).real()`.
- **`mpscpp3`-specific: real-time MPS evolution defaults to TDVP.**
  `mpscpp3/TDVP/` (vendored from ITensor's own TDVP repo, see its
  `README.md`) provides a proper two-site TDVP integrator, wired into
  `chain_session.h` as `Chain::quench_tdvp()`/
  `Chain::evolve_and_measure_tdvp()` (via the private `tdvp_step()`
  helper) alongside the pre-existing `Chain::quench()`/
  `Chain::evolve_and_measure()`, which apply a hand-rolled 2nd-order
  Taylor expansion of `exp(-i dt H)` as an MPO (`evoloperator()`) instead.
  `Many_Body_Chain.tevol_method` (`manybodychain.py`, default `"TDVP"`)
  picks between them in `timedependent.py`'s `evolution_dmrg_DC()`/
  `evolve_and_measure_dmrg()` — `"TDVP"` only applies when
  `itensor_version==3` (mpscpp2 has no `TDVP/` and no `_tdvp` methods at
  all); any other combination silently falls back to the MPO-Taylor path,
  which remains the only option for `itensor_version=2`. Only the
  two-site TDVP algorithm is used (`NumCenter=2`); the global subspace
  expansion machinery in `TDVP/basisextension.h` (mainly useful for
  one-site TDVP or long-range Hamiltonians) is not wired in, since
  two-site TDVP already grows bond dimension via SVD like two-site DMRG.
- **Both backends can be loaded in the same process** (needed for anything
  that directly compares `itensor_version=2` vs `=3` results, e.g.
  `examples/v2_VS_v3_*`, or `examples/dynamical_correlator/
  dynamical_correlator_v2_VS_v3`) only because `mpscpp3/bindings.cc`
  registers its pybind11 types with `py::module_local()`. Without it,
  `mpscpp2` and `mpscpp3` define their own, ABI-incompatible
  `Chain`/`MPS`/`MPO` C++ types (same unqualified/unnamespaced names,
  different vendored ITensor copies), and libstdc++'s cross-DSO RTTI
  equality falls back to comparing the mangled type-name *string*, not the
  address — so pybind11 (keyed on `std::type_index`) sees the exact same
  "type" already registered the moment both extensions are imported into
  one process, and aborts with `generic_type: type "Chain" is already
  registered!`. `module_local()` keeps mpscpp3's registration in a
  module-scoped table instead of the shared global one, avoiding the
  collision regardless of import order (confirmed both ways). This is
  applied only on the mpscpp3 side — mpscpp2 is untouched and still uses
  the ordinary global registration.

### Supporting `*tk` packages

Functionality is generally split into a top-level module (the public method)
and a `<topic>tk/` package with the implementation, e.g.
`correlator.py`/`fermionchaintk/staticcorrelator.py`,
`dynamics.py`/`dynamicstk/`, `entropy.py`/`entropytk/`,
`mpsalgebra.py`/`mpsalgebratk/`, `entanglement.py`. When changing behavior
for a specific feature, look for both the thin dispatch module and its `tk`
counterpart.

### KPM / dynamical correlators

Dynamical correlators and generic operator distributions are computed with
the Kernel Polynomial Method (`kpmdmrg.py`, `Chain::kpm_dynamical_correlator`/
`Chain::general_kpm` in `mpscpp2/chain_session.h`/`mpscpp3/chain_session.h`
on the C++ DMRG side, `pyfermion/algebra/kpm.py` on the ED side) rather
than direct spectral decomposition, since exact diagonalization of the
full spectrum is infeasible for large chains.

### Julia vs C++ backend

The two DMRG backends are independent implementations, not a shared
protocol: the C++ path runs in-process through the pybind11 extension (see
above), while `itensor_version="julia_live"` (`mpsjulialive/`) drives a
live Julia/ITensors.jl session via `pyjulia` with its own mirrored set of
modules. A feature missing on one side simply isn't implemented there yet
(check the relevant `mpsjulialive/*.py` file, or lack thereof) rather than
falling back automatically — the only automatic fallback is DMRG → ED when
the C++ extension isn't compiled (`mode.py`).

### ITensor library
the ITensor folders (`mpscpp2/ITensor` = v2, `mpscpp3/ITensor` = v3) are a library that is developed elsewhere, hence you do not need to read them carefully. Read them only if there is some feature that does not work with the mpscpp2/mpscpp3 code, and you need to figure out why
