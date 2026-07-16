# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

DMRGPY is a Python library for quasi-one-dimensional spin chains, fermionic
systems, parafermion and bosonic models using matrix product states (MPS).
It exposes a single Python API (`src/dmrgpy/`) backed by two interchangeable
solvers per calculation:

- **DMRG**, executed by an external ITensor-based backend (C++ or Julia)
- **ED** (exact diagonalization), a pure-Python fallback used for small
  systems and for cross-checking DMRG results

The Python layer is what this project actually develops. The C++/ITensor
backend under `src/dmrgpy/mpscpp2/ITensor/` is a vendored copy of the
ITensor library maintained upstream — treat it as a black box invoked via a
compiled executable, not code to read or modify here.

## Install / build

There is no `setup.py`/`pyproject.toml`; the project is used from `src/`
via a path added to `PYTHONPATH` (or a symlink into `site-packages`).

```bash
python install.py                    # compiles ITensor + the C++ backend (mpscpp.x)
python install.py --gpp=g++-6        # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install_julia.py              # alternative: Julia backend instead of C++
```

`install.py` calls `installtk/install2.py` to compile the C++ program, then
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

### Dual backend dispatch (DMRG vs ED)

`mode.py` decides, per call, whether a calculation runs via DMRG or ED
(`get_mode`/`run`), based on `self.mode`, whether `mpscpp.x` is compiled, and
`self.itensor_version` (`2`/`"C++"` → C++ backend via `cpprun.py`;
`"julia"` → `juliarun.py`). Most public methods on `Many_Body_Chain` accept a
`mode="DMRG"|"ED"` kwarg so results can be cross-validated between solvers
(see the bilinear-biquadratic example in `README.md`).

- **DMRG path, in-process extension (default)**: see "In-process pybind11
  extension" below.
- **DMRG path, file-based (fallback)**: `taskdmrg.py` writes a `tasks.in`
  file describing the requested task (`GS`, `excited`, `correlator`,
  `dynamical_correlator`, `time_evolution`, ...) plus DMRG parameters
  (`maxm`, `nsweeps`, `cutoff`, etc.), then `cpprun.py`/`juliarun.py` invoke
  the compiled backend (`mpscpp2/mpscpp.x`, built from the headers in
  `mpscpp2/*.h` on top of vendored ITensor) or the Julia scripts in
  `mpsjulia/`. Each call runs in an isolated temp working directory
  (`/tmp/dmrgpy_tmp/<random-id>`, see `filesystem.py` and
  `id_generator`/`to_folder`/`execute` in `manybodychain.py`) so concurrent
  chains don't collide. This path still fully exists and is what runs
  whenever the in-process extension is unavailable or hasn't been ported
  for a given call.
- **ED path**: `edtk/edchain.py` (`EDchain`) builds dense/sparse operators
  directly in Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`,
  `pyzn/` provide the per-statistics many-body operator construction,
  `edtk/one2many.py` promotes single-site operators to the full Hilbert
  space) and diagonalizes with `scipy.sparse.linalg`.

### In-process pybind11 extension (default DMRG path)

`mpscpp2/bindings.cc` compiles a pybind11 extension module
(`mpscpp2/_dmrgcpp*.so`, built by the `pybind` target in `mpscpp2/Makefile`)
that embeds the same ITensor-based DMRG code as the file-based backend but
runs it in-process — no `tasks.in`, no operator files, no subprocess, no
temp directory. This is the default backend as of `use_cpp_extension=True`
(`manybodychain.py`); the file-based backend above is kept as a complete,
independent fallback, not dead code.

- **`cppext.py`** lazily imports the compiled extension and caches whether
  it loaded (`get_backend()`/`available()`). If the `.so` isn't built (or
  `itensor_version != 2`), the extension is simply absent and every call
  site below falls back to the file-based path — nothing needs a `--flag`
  to degrade gracefully.
- **`mpscpp2/chain_session.h`**'s `Chain` class is the stateful, in-process
  counterpart to the whole file-based task protocol: one `Chain` instance
  per `Many_Body_Chain` (held as `self._session`, created in
  `sites.py::initialize`), with a method per DMRG task (`gs_energy`, `vev`,
  `apply_operator`, KPM/CVM dynamical correlators, `correlation_matrix`,
  `reduced_dm`, `excited_states`, time evolution, ...). It is built
  directly against ITensor's `AutoMPO`/`MPS`/`MPO` APIs, not by shelling
  out to `mpscpp.x`.
- **`mpscpp2/mo_terms.h`** and `MultiOperator.to_terms()`
  (`multioperator.py`) are the in-memory replacement for the AMPO operator
  file format: a `MultiOperator` converts to a plain list of
  `(coefficient, [(opname, site), ...])` tuples that `Chain` turns directly
  into an ITensor `HTerm`/`AutoMPO`.
- **The `cpp_handle` pattern**: `mps.py::MPS` and
  `multioperatortk/staticoperator.py::StaticOperator` both carry an opaque
  `cpp_handle` attribute — `None` means "this object lives in the
  file-based world" (an `.mps`/`.mpo` file on disk), non-`None` means "this
  is a pybind11-wrapped ITensor object living only in C++ memory". Every
  dispatch site that can take either kind of object must check
  `self._session is not None` **and** that every wavefunction/operator
  involved has a non-`None` `cpp_handle` before taking the in-process path
  (see `mpsalgebra.py::_use_cpp_ext`) — a chain can have the extension
  available while still holding file-based objects produced by a
  not-yet-ported code path, and the two must never be mixed.
- Because `cpp_handle` objects have no pickle/deepcopy support, `MPS.copy()`
  / `StaticOperator.copy()` use a shallow `copy.copy()` (never
  `deepcopy()`), and `Many_Body_Chain.__deepcopy__` explicitly drops
  `_session` (set to `None`) on any clone.
- A few **pre-existing bugs in the original file-based backend are
  deliberately reproduced, not fixed**, in the in-process port (see
  comments at the call sites for details): `evoloperator()`'s z³/6 term
  multiplies `H2` instead of the computed `H3`; `tasks.in`'s `"noise"` key
  is read back as `"moise"`; `"tevol_fit"` is written but
  `"tevol_fit_td"` is read, making the `fitApplyMPO` time-evolution branch
  unreachable in both backends. If ITensor results look numerically
  "close but not exact" to a hand-derivation, check here before assuming a
  porting bug.
- Not every feature has an in-process implementation yet — anything built
  under `mpscpp2/` that doesn't route through `chain_session.h` (e.g. the
  `"dos"` task, which is dead/commented-out in `mpscpp.cc` itself) simply
  never gets a `cpp_handle`-carrying object, so it stays on the file-based
  path automatically.

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
the Kernel Polynomial Method (`kpmdmrg.py`, `mpscpp2/kpm*.h` on the DMRG
side, `pyfermion/algebra/kpm.py` on the ED side) rather than direct spectral
decomposition, since exact diagonalization of the full spectrum is
infeasible for large chains.

### Julia vs C++ backend

Both backends implement the same task protocol (read `tasks.in` + operator
files, write results back to the working directory). `mpsjulia/*.jl` is the
Julia mirror of `mpscpp2/*.h`; when a feature exists in one backend but not
the other, check `mode.py`/`juliarun.py`/`cpprun.py` for how the dispatch
falls back (currently ED is the fallback when `mpscpp.x` isn't compiled).

### ITensor library
the ITensor folder is a library that is developed elsewhere, hence you do not need to read it carefully. Read it only if there is some feature that does not work with the mpscpp code, and you need to figure out why
