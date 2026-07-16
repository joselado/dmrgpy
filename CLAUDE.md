# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

DMRGPY is a Python library for quasi-one-dimensional spin chains, fermionic
systems, parafermion and bosonic models using matrix product states (MPS).
It exposes a single Python API (`src/dmrgpy/`) backed by interchangeable
solvers per calculation:

- **DMRG**, executed in-process via a pybind11 extension over ITensor
  (C++), or via a separate Julia/ITensors.jl backend
- **ED** (exact diagonalization), a pure-Python fallback used for small
  systems, for cross-checking DMRG results, and automatically whenever the
  C++ extension isn't compiled

The Python layer is what this project actually develops. The C++/ITensor
backend under `src/dmrgpy/mpscpp2/ITensor/` is a vendored copy of the
ITensor library maintained upstream — treat it as a black box, not code to
read or modify here.

## Install / build

There is no `setup.py`/`pyproject.toml`; the project is used from `src/`
via a path added to `PYTHONPATH` (or a symlink into `site-packages`).

```bash
python install.py                    # compiles ITensor + the pybind11 extension
python install.py --gpp=g++-6        # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install_julia.py              # alternative: Julia backend instead of/alongside C++
```

`install.py` calls `installtk/install2.py`, which compiles the vendored
ITensor static library and then `mpscpp2`'s `pybind` Makefile target
(needs `pybind11` installed; skipped with a warning if it isn't), then
adds the repo's `src` to the Python path and to the user's shell rc file
(`installtk/addpythonpath.py`, `installtk/addsystem.py`). There is no
separate lint or CI config in this repo.

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
After running examples, `python clean.py` recursively removes generated
working directories (`.mpsfolder`, `.pychainfolder`, `.dmrgfolder`) and stray
`ERROR`/`*.OUT` files from the tree.

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
`self.itensor_version==2` and the pybind11 extension isn't compiled (see
`cppext.available()`), in which case it silently falls back to ED. Most
public methods on `Many_Body_Chain` accept a `mode="DMRG"|"ED"` kwarg so
results can be cross-validated between solvers (see the
bilinear-biquadratic example in `README.md`).

- **DMRG, C++ (`itensor_version=2`, the default)**: entirely in-process,
  see "In-process pybind11 extension" below. There is no file-based/
  subprocess fallback for this path — if the extension isn't compiled,
  `mode.py` routes to ED instead.
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

### In-process pybind11 extension (the C++ DMRG backend)

`mpscpp2/bindings.cc` compiles a pybind11 extension module
(`mpscpp2/_dmrgcpp*.so`, built by the `pybind` target in `mpscpp2/Makefile`)
that runs DMRG entirely in-process against vendored ITensor — no task
files, no operator files, no subprocess, no temp directory. This is the
only C++ DMRG path; a previous file-based/subprocess design (writing
`tasks.in` and invoking a standalone `mpscpp.x` executable) has been
removed.

- **`cppext.py`** lazily imports the compiled extension and caches whether
  it loaded (`get_backend()`/`available()`).
- **`mpscpp2/chain_session.h`**'s `Chain` class is the stateful, in-process
  session: one `Chain` instance per `Many_Body_Chain` (held as
  `self._session`, created in `sites.py::initialize` whenever
  `itensor_version==2` and the extension is available), with a method per
  DMRG task (`gs_energy`, `vev`, `apply_operator`, KPM/CVM dynamical
  correlators, `correlation_matrix`, `reduced_dm`, `excited_states`, time
  evolution, ...), built directly against ITensor's
  `AutoMPO`/`MPS`/`MPO` APIs. `check_task.h`, `get_sites.h` and
  `mo_terms.h` are the only other `mpscpp2/*.h` headers it depends on;
  everything else it needs is self-contained in `chain_session.h`.
- **`mpscpp2/mo_terms.h`** and `MultiOperator.to_terms()`
  (`multioperator.py`) convert a `MultiOperator` into a plain list of
  `(coefficient, [(opname, site), ...])` tuples that `Chain` turns directly
  into an ITensor `HTerm`/`AutoMPO`.
- **The `cpp_handle` pattern**: `mps.py::MPS` and
  `multioperatortk/staticoperator.py::StaticOperator` both carry an opaque
  `cpp_handle` attribute holding the pybind11-wrapped ITensor object living
  in C++ memory — these classes have no other representation now, so
  `cpp_handle` is always set for `itensor_version==2` chains (it's `None`
  only for objects that were never meant to touch the C++ backend, e.g.
  `mpsjulialive`'s own `MPS`/`MPO` classes are separate types entirely).
- Because `cpp_handle` objects have no pickle/deepcopy support, `MPS.copy()`
  / `StaticOperator.copy()` use a shallow `copy.copy()` (never
  `deepcopy()`). `Many_Body_Chain.__deepcopy__` never copies `_session`
  directly (a live C++ session has no well-defined "copy" semantics);
  instead it builds a fresh `Chain` for the clone when the original had
  one, which is what `clone()`/`bandwidth()`/`lowest_eigenvalue()` rely on.
- A few **pre-existing bugs in the original design are deliberately
  reproduced, not fixed**, in `chain_session.h` (see comments at the call
  sites for details): `evoloperator()`'s z³/6 term multiplies `H2` instead
  of the computed `H3`; the old file-based backend's `"noise"` tasks.in key
  was actually read back as `"moise"`, and `"tevol_fit"` was written but
  `"tevol_fit_td"` was read, making that fitting branch unreachable — both
  reproduced as hardcoded values rather than silently "fixed". If ITensor
  results look numerically "close but not exact" to a hand-derivation,
  check here before assuming a porting bug.

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
`Chain::general_kpm` in `mpscpp2/chain_session.h` on the C++ DMRG side,
`pyfermion/algebra/kpm.py` on the ED side) rather than direct spectral
decomposition, since exact diagonalization of the full spectrum is
infeasible for large chains.

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
the ITensor folder is a library that is developed elsewhere, hence you do not need to read it carefully. Read it only if there is some feature that does not work with the mpscpp code, and you need to figure out why
