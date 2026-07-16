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

- **DMRG path**: `taskdmrg.py` writes a `tasks.in` file describing the
  requested task (`GS`, `excited`, `correlator`, `dynamical_correlator`,
  `time_evolution`, ...) plus DMRG parameters (`maxm`, `nsweeps`, `cutoff`,
  etc.), then `cpprun.py`/`juliarun.py` invoke the compiled backend
  (`mpscpp2/mpscpp.x`, built from the headers in `mpscpp2/*.h` on top of
  vendored ITensor) or the Julia scripts in `mpsjulia/`. Each call runs in
  an isolated temp working directory (`/tmp/dmrgpy_tmp/<random-id>`, see
  `filesystem.py` and `id_generator`/`to_folder`/`execute` in
  `manybodychain.py`) so concurrent chains don't collide.
- **ED path**: `edtk/edchain.py` (`EDchain`) builds dense/sparse operators
  directly in Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`,
  `pyzn/` provide the per-statistics many-body operator construction,
  `edtk/one2many.py` promotes single-site operators to the full Hilbert
  space) and diagonalizes with `scipy.sparse.linalg`.

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
