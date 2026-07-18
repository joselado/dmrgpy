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

There is no `setup.py`/`pyproject.toml`; DMRGPY is used directly from
`src/`, added to `PYTHONPATH` (or symlinked into `site-packages`) by the
installer.

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
— only NumPy/SciPy (and optionally `numba`/`jax` for faster execution,
see §5.4).

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
| `bosonchain.py` | `Boson_Chain` | bosonic sites |
| `parafermionchain.py` | `Parafermion_Chain` | parafermion sites |

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
construction, and sympy-based symbolic building.

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
| `"julia_live"` | live in-process Julia session (`mpsjulialive/`), ITensors.jl | `pyjulia` + Julia install | none (feature-by-feature; missing methods simply aren't implemented) |

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

Notable, deliberate implementation details (not bugs to "fix"):

- `mpscpp3` builds every site with `ConserveQNs=false` and starts DMRG
  from an actual `randomMPS`, not a plain product state — matching v2's
  real (unconstrained-search) behavior rather than ITensor v3's stricter,
  QN-conserving-from-the-start convention, at the cost of losing the QN
  block-sparsity speedup for `itensor_version=3`.
- `mpscpp3`'s real-time MPS evolution defaults to a proper two-site TDVP
  integrator (vendored under `mpscpp3/TDVP/`), selectable via
  `Many_Body_Chain.tevol_method` (default `"TDVP"`); `mpscpp2` has no TDVP
  and always uses a hand-rolled 2nd-order Taylor expansion of
  `exp(-i dt H)` as an MPO instead.
- A few pre-existing bugs in the original (pre-refactor) file-based
  backend are deliberately reproduced rather than silently fixed —
  see the call-site comments in `chain_session.h` for both versions.

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
(via `pyjulia`) with its own set of modules mirroring the top-level ones
(`mpsjulialive/groundstate.py`, `vev.py`, `mps.py`, `mpo.py`, ...) that
talk to Julia's ITensors.jl instead of the C++ pybind11 extension. This
is an independent implementation, not a shared protocol with the C++
path — a feature missing on one side simply isn't implemented there yet.
(A separate, older subprocess-based Julia path, `itensor_version="julia"`
via `juliarun.py`, is not reachable through the normal public API and
should be treated as legacy/inert.)

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
- Install `numba` (and optionally `jax`) alongside the pure-Python
  backend *and* set `pyitensor.kernels.USE_NUMBA = True` (off by
  default, see §5's note above) — without both steps, `pyitensor` falls
  back to plain NumPy loops and is substantially slower than the numbers
  above (see `pyitensor/__init__.py`'s own docstring).
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
