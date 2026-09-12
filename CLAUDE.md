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

There is a `pyproject.toml`, but it is *only* for building the PyPI
distribution (see "Packaging / PyPI" below) -- it is not how this repo is
developed. For development the project is used from `src/` via a path
added to `PYTHONPATH` (or a symlink into `site-packages`), which is what
`install.py` sets up; there is no `setup.py` and no editable install.

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
selects (`2` = mpscpp2/ITensor v2; `3` = mpscpp3/ITensor v3, the default;
`both` builds both, one after the other).

Both `mpscppN/Makefile`s list `$(wildcard *.h)` among the extension's
prerequisites: without that, editing `chain_session.h` — where the actual
DMRG implementation lives — leaves `make pybind` reporting "Nothing to be
done" while the `.so` keeps the old code, which is a confusing way to test a
change that was never compiled.

`install2.py` picks the right
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

## Packaging / PyPI

`pyproject.toml` builds a **pure-Python** distribution: a `py3-none-any`
wheel (~450KB) containing the ED backend, the `itensor_version="python"`
pyitensor DMRG/TDVP backend, and the Julia bridge. It deliberately ships
no C++ at all. `[tool.setuptools.packages.find]` excludes `dmrgpy.mpscpp2`
and `dmrgpy.mpscpp3` — those are ~400MB of vendored ITensor plus a build
that needs a compiler, LAPACK/BLAS and pybind11, none of which belongs in
a pip install. Users who want the compiled backends clone the repo and
run `install.py`; without them everything still works via `"python"`/ED
(`mode.py`), so the wheel is fully usable as installed.

```bash
python -m build                  # -> dist/*.whl + dist/*.tar.gz
python -m twine check dist/*     # metadata/README validation
python -m twine upload --repository testpypi dist/*   # dry run on TestPyPI first
python -m twine upload dist/*    # the real thing
```

Things in here that are load-bearing and easy to break:

- **`numba` is a core dependency, not an extra.** `algebra/kpmextrapolate.py`
  imports it at module level and sits on the unconditional
  `manybodychain -> dynamics -> kpmdmrg -> algebra.kpm` import path, so
  without it even `import dmrgpy.spinchain` raises `ModuleNotFoundError`.
  This is unrelated to pyitensor's *own* opt-in numba/JAX matvec kernel
  (`pyitensor/kernels.py`), which is properly `try/except`-guarded and
  genuinely optional.
- **`requires-python = ">=3.10"` is set by numba**, which declares
  `>=3.10`. Claiming a lower floor doesn't make dmrgpy work on 3.9 — pip
  just backtracks onto an ancient, untested numba/numpy pair.
- **`setuptools>=77` in `[build-system]`**, not the more common 61: the
  `[project]` table uses PEP 639 licensing (bare SPDX `license` string plus
  `license-files`), which older setuptools rejects outright.
- **The `julia` extra is `juliacall`, not `julia`.** The live backend
  migrated from PyJulia to juliacall/PythonCall.jl. The only surviving
  `import julia` is `juliarun.py`, part of the legacy subprocess path that
  isn't reachable from the public API, so PyJulia is intentionally not
  pulled in.
- **`src/dmrgpy/juliapkg.json` must stay in `package-data`.** juliacall/
  juliapkg reads it at import time to resolve and precompile ITensors,
  ITensorMPS and ITensorNHDMRG; if it's missing from the wheel,
  `pip install dmrgpy[julia]` yields a Julia backend with no ITensors.jl.
  (`.gitignore`'s unanchored `julia*` pattern would also swallow it —
  hence the explicit `!src/dmrgpy/juliapkg.json` un-ignore. Don't remove
  either one.)
- **The classifiers say `OS Independent`**: the pip-installed package
  (ED and `itensor_version="python"` backends) works on Windows.
  `sites.py::initialize()` used to shell out to `mkdir -p` via
  `subprocess` and break there, but that call was removed in `1f7b7a9`
  ("Stop creating .mpsfolder; remove dead file-based-backend
  leftovers"), which also means `self.path`/`.mpsfolder` is never
  created on disk any more -- `Many_Body_Chain.clean()`, whose only job
  was `rm -rf`-ing that folder (and whose only caller, `bandwidth()`,
  just deleted nothing every time), was removed outright rather than
  ported to `shutil.rmtree`, since there was nothing left for it to do.
  The compiled C++ backends (`install.py`) remain POSIX-only and are
  not part of the wheel regardless.
- Every module in the wheel must at least byte-compile, or pip prints a
  `SyntaxError` on install; `python -m compileall` over the installed
  package is the cheap check.

Verify a release candidate by installing the built wheel into a clean
virtualenv (not from `src/`, so import bugs and missing package-data
actually surface) and cross-checking the two shipped backends against each
other — `itensor_version="python"` DMRG vs `mode="ED"` on a small chain
should agree to ~1e-15. `dist/`/`build/` are gitignored; build artifacts
are not committed. Note that a given version can never be re-uploaded to
PyPI (or TestPyPI) once published — bump the version rather than trying to
replace a release.

## Running / verifying changes

There is a `tests/` directory with a real pytest suite (`test_spin_chain.py`,
`test_fermion_chain.py`, `test_spinful_fermion_chain.py`,
`test_parafermion_chain.py`, `test_dynamical_correlator.py`,
`test_entanglement.py`, `test_excited_states.py`, `test_time_evolution.py`,
`test_long_chain.py`, `test_jordan_wigner.py`, `test_internal_consistency.py`,
`test_benchmarks.py`) — run it with `python run_tests.py` from the repo root
(a thin wrapper: runs `pytest tests` from the repo root, then `clean.py`) or
plain `pytest tests`. `tests/conftest.py` inserts this checkout's own `src/`
at the front of `sys.path`, so — unlike running an `examples/*/*/main.py`
script directly — pytest always exercises the current working tree's code,
never a stale `site-packages`/symlinked checkout (see the worktree-symlink
pitfall below). `tests/_helpers.py`'s `energy_ed_v2_v3`/`vev_ed_v2_v3` are the
standard way most tests cross-check ED against DMRG on the C++ backends (a
`versions=` kwarg lets a test opt a specific `itensor_version` out — e.g. v3
on an exactly-2-site chain, a real `mpscpp3` bug, see below); `versions=`
only covers `(2, 3)` today, not the pure-Python backend
(`itensor_version="python"`) — which is *not* to say that backend is
untested: dozens of files in `tests/` exercise it, they just parametrize
their own `itensor_version` instead of going through this helper.
`tests/reference_data.py` holds shared golden values pinned to a specific
historical commit (its own module docstring explains how to regenerate
them); several individual test files also just inline their own golden
values as local constants instead, which is equally fine for values that
are only used in one place. When a test's assertion is against exact
DMRG/ED agreement (not a golden regression value), prefer
`pytest.approx(..., abs=1e-6)`-style tolerances over exact equality, since
DMRG is iterative.

**`itensor_version="julia_live"` tests only need to run when a change
touches Julia code** (`mpsjulialive/`, any `*.jl` file, or a shared entry
point like `manybodychain.py`/`infinitechain.py` that dispatches into
that backend). These tests carry a large, mostly *fixed* cost unrelated
to what they're actually testing: `juliacall`'s JIT compiles each
distinct Julia function signature the first time it's called in a
process, and this compile cost dominates wall time regardless of how
small the test's own parameters are — confirmed directly by isolating a
`julia_live` `metts_dynamical_correlator` call down to `nsamples=1`,
`nwarmup=1`, `nt=1`: it still took ~73s, almost entirely JIT, versus ~27s
of genuine algorithmic work for the same call at its normal (nsamples=40)
parameters once the function was already warm. Reducing test parameters
further cannot buy back that fixed cost, and it's paid independently by
*every* distinct Julia function signature exercised across the whole
suite (`tests -k julia_live` alone still pays it many times over), which
is why the full suite's slowest handful of tests are almost all
`[julia_live]`-parametrized. For a change that doesn't touch Julia code
at all, skip them locally with `pytest tests -k "not julia_live"` (the
parametrize id `julia_live` appears in every such test's node id) rather
than paying this tax on every iteration; run the full suite including
`julia_live` before finishing work that *does* touch Julia code, or
before a release.

The `examples/` directory (100+ self-contained scripts, one per physical
model or feature, organized into thematic subfolders — `groundstate/`,
`staticcorrelators/`, `dynamical_correlator/`, `spin_models/`,
`fermion_models/`, `boson_models/`, `time_evolution/`, `excited_states/`,
`entanglement/`, `topological/`, `kondo/`, `non_hermitian/`,
`finite_temperature/`, `magnetization/`, `backend_comparison/`,
`utilities/`, `algebra/`, ...) is a second, complementary regression
surface — broader and more exploratory than `tests/` (most scripts just
print/plot rather than assert), but also covers more ground (every
physical model, every dynamical-correlator submode, every backend
combination including `"python"`). Some `examples/*/*/main.py` scripts
*do* use real `assert`s and are meant as regression tests in the same
spirit as `tests/` (e.g. `time_evolution/tdvp_VS_ED_time_evolution`,
`backend_comparison/backend_switch_consistency`,
`finite_temperature/thermal_purification_VS_exact`,
`staticcorrelators/static_correlator_VS_ED`,
`entanglement/entanglement_entropy_VS_ED`,
`dynamical_correlator/dynamical_correlator_VS_ED` — several of these were
added specifically to lock in bugs found and fixed in this codebase, and
are good templates for adding another one). To run one:

```bash
cd examples/readme_examples/static_correlator_S12 && python <script>.py
```

`examples/readme_examples/` mirrors the snippets shown in `README.md`.
`examples/backend_comparison/`, the `v2_VS_v3_*`-named examples living
inside each theme folder, and `examples/dynamical_correlator/
dynamical_correlator_v2_VS_v3` directly compare the ITensor v2 and v3 C++
backends against each other (same script, both `itensor_version`s, small
systems) — the fastest way to check a `mpscpp3` change didn't diverge from
`mpscpp2`'s numerics. One rule if you use one as a template: pin the sweep
schedule inside the script instead of inheriting the library defaults.
`boson_models/v2_VS_v3_boson` did the latter and its `assert` fired
non-deterministically on a clean tree (2 of 5 runs for the audit's
reviewer, 5 of 7 for its hunter) — v2 and v3 were converging to the same answer,
just not far enough to be compared at the tolerance it asserted; it now
pins `nsweeps=80`/`maxm=100` and asserts at 1e-5 (2026-09 audit #24).
After running examples, `python clean.py` recursively
removes generated working directories (`.mpsfolder`, `.pychainfolder`,
`.dmrgfolder`) and stray `ERROR`/`*.OUT` files from the tree.

**The 2026-08 audit.** `docs/audit_2026_08_hole_hunt.md` is the record of
a five-lens cross-backend audit: 21 confirmed findings, each with the
reproduction that was actually executed and an independent reviewer's
analysis, and each marked with what was done about it. Every one is fixed
except the `kpm_energy_truncate` window problem, which is guarded only in
its worst regime (`docs/known_issue_kpm_energy_truncation_window.md`).
The regressions live in `tests/test_audit_2026_08_regressions.py`. Two of
those fixes changed numbers rather than behavior -- `submode="TD"`/`"TDZ"`
were a factor of pi too large, and `itensor_version=3` real-time evolution
from a non-unit-norm start was scaled by `1/||psi||^2` -- so results from
before them are not comparable. Read `documentation.md`'s section 4.10
before adding a new dispatch: almost every finding was a dispatch decision
taken before the information that should inform it (a short circuit ahead
of a cache key, a precondition ahead of the branch it qualifies, an `else`
serving both "unsupported here" and "you typo'd it", a `**kwargs` with no
consumer).

**The 2026-09 audit.** `docs/audit_2026_09_hole_hunt.md` is the same
thing one lens-set wider: eight lenses, 36 confirmed findings, each with
its executed repro, its reviewer's attempt to refute it, and a `**Status**`
line saying what was done and which test now pins it. The regressions live
in eight files, `tests/test_audit_2026_09_<name>.py`: one per fix cluster
(there were seven) plus one for the `dispatch-leftovers` follow-up lane.
One entry is not a plain FIXED: **#20**'s compute half landed
and its memory half did not (`idmrg.py` still materializes `Es`, ~540 MB
peak at chi=64). **#3** was fixed on `itensor_version="python"` only at
the time, the C++ VUMPS half following on 2026-09-12
(`docs/known_issue_v3_vumps_variational_floor.md`, now FIXED — see the
`vx_choose_fixed_point` paragraph further down). Results from before this
audit are not comparable where it changed numbers rather than behavior,
which it did in ten places: all `itensor_version="python"` TDVP
real-time evolution (the default `tevol_method`, which carried an O(dt)
error); every `"python"` `applyMPO` consumer, at the 1e-5..1e-9 level
(`applyoperator`, `vev(npow>1)`, `gs_energy_fluctuation`, KPM/CVM, the
MPO-Taylor stepper) and ~1.4x slower for the repeated-application ones;
`submode="CVM_explicit"`, which was exactly 2x too large on every backend;
`submode="ED"`, `mode="DMRG"` `submode="CVM"` and `submode="ROOTN"`, which
now return the complex Lehmann density `i(G^R-G^A)/(2pi)`, as the
resolvent submodes on `mode="ED"` already did (`submode="TD"`/`"TDZ"` are
the routes still off it, deliberately — open item O1 in the audit record).
Unchanged wherever the Lehmann weights `M_n = <GS|A|n><n|B|GS>` are real —
`Im M_n == 0` is the discriminant, which a Hermitian pair `A = B^dagger`
implies but does not exhaust: a real Hamiltonian with real operators has
real `M_n` off-diagonally too. Where they are complex the two conventions
differ by an amount comparable to the correlator's own peak (measured
0.57 against a 0.61 peak; for `submode="ED"` the difference *is* the
discarded imaginary part, so it can never exceed the peak);
`vev(npow!=1)` and `gs_energy_fluctuation` on every ED
route, which ignored `npow`/`mode=` entirely; `exponential()` on a DMRG
backend, which computed `e^{-h}` or an unconverged 2-term Taylor series;
`SpinBoson_Chain(maxnb=...)`, which used to build a 4-level boson
whatever you asked for; `itensor_version="python"` parameter sweeps and
`bandwidth()` on a *reused* chain (finding #2), which the
pyitensor-evolution cluster recorded as wrong by up to the operator's
full spectral width (+1.25 against an exact -2.4936 on a 6-site
Heisenberg chain; `bandwidth()` 0.0 against an exact 2.366025) and which
`mpsalgebra`'s bandwidth-scaled Taylor expansions, `excited_states`'
Lagrange weight, `meanfield.py` and `thermal.py` all consume;
`itensor_version="python"` `promote_to_dense()` followed by `gs_energy()`
(finding #10), which re-solved unconstrained and returned the *global*
ground state instead of the sector's (-2.3399 against -1.9408 on a 3-site
Hubbard chain at `Nf=3`) — wrong, not merely different; and `mode="ED"`
boson occupation projectors `bc.D[i][k]` on a number-*breaking*
Hamiltonian (finding #4), where the ed-backend cluster recorded `P(n=1)`
on site 0 going from -0.0626 to 0.5405 and `sum_k P` from 0.031 to 1.000
(number-conserving Hamiltonians are unchanged to machine precision).
`src/dmrgpy/dynamics.py`'s module docstring is now
the single normative statement of the correlator convention — read it
before adding a submode, the way §4.10 is what to read before adding a
dispatch. Two items are open rather than fixed, found while re-measuring
that record and recorded in its own "Open items found while correcting
this record" section: `submode="TD"`/`"TDZ"` return the complex
one-sided transform `-(i/pi)G^A` rather than the density (they agree in
real part only, and only when `Im M_n = 0` — measured at 70% of the
correlator's peak in the imaginary part on a Hermitian pair), and
`mode="ED"` vs `mode="DMRG"` disagree pointwise under the default
`submode="KPM"` (both satisfy the sum rule exactly, and the disagreement
looks like a resolution difference rather than a convention one — but it
is recorded as unexplained, so read the record before acting on it).

**Examples should plot, not just print/assert.** What sets `examples/`
apart from `tests/` is that a human is expected to actually look at the
result, so every `examples/*/*/main.py` should end with a `matplotlib`
plot of whatever it computed (a quantity vs. time/site/parameter, a
backend-vs-backend overlay, a timing/error scaling curve, ...) whenever
the script produces a sequence of values that can meaningfully be
visualized — not just `print()`ed. This applies even to the scripts noted
above that *also* carry a real `assert` for regression purposes (e.g.
`time_evolution/tdvp_VS_ED_time_evolution`,
`time_evolution/tebd_VS_ED_time_evolution`,
`time_evolution/tdvp_gse_VS_ED_time_evolution`): the assert stays as the
pass/fail regression guard, but the script should still plot the two
trajectories being compared afterwards, so the same script is useful both
as an automated check and as a visual sanity check when run by hand. A
script whose entire output is a single scalar (e.g. one ground-state
energy, one overlap) is the one legitimate exception — there, sweep the
relevant parameter (system size, time, coupling strength, ...) instead of
computing just one point, if a sweep is cheap enough to be worth adding;
otherwise printing is fine. Treat that exception as a last resort, not a
default: an audit of the whole `examples/` tree found 105 scripts that had
quietly settled for print/assert-only even though an obvious cheap axis
was sitting right there (chain length, a coupling/field, time/frequency,
site/bond index, backend or `itensor_version` comparison, ED vs DMRG) —
usually because the script already computed an array or a comparison and
simply never plotted it. When adding a new example or touching an old
one, check first whether it already has such an axis before reaching for
the single-scalar exception.

**Worktree/symlink pitfall**: `~/.local/lib/python3.13/site-packages/dmrgpy`
is typically a symlink into *one specific* checkout's `src/dmrgpy/` (often
the primary working directory, not whichever worktree you're in). Since
every `examples/*/*/main.py` does `sys.path.append(...)` (append, not
prepend/insert), that symlinked package can resolve *before* the appended
path when running a script directly with plain `python3 main.py` from a
different worktree — silently testing the wrong checkout's code with no
error. Force the intended checkout explicitly, e.g.
`PYTHONPATH=/path/to/this/worktree/src python3 main.py`, or insert
`sys.path.insert(0, ...)` in a one-off driver script instead of relying on
the example's own `sys.path.append`. `tests/conftest.py` (see above) avoids
this problem entirely for anything run through pytest.

**A real, confirmed `mpscpp3` bug**: ITensor v3's two-site `dmrg()` (see
`chain_session.h`'s `dmrg_args()`) aborts the whole process (`SIGABRT`,
`"LocalOp is default constructed"`, an internal ITensor v3 check in
`itensor/mps/localop.h`) for any chain with fewer than 3 physical sites,
independent of physics type — confirmed for spin, fermionic, and other
chains. `mode.py`'s `get_mode()` now falls back to ED automatically for
`itensor_version==3` with `self.ns<3` (same mechanism as the existing
"extension not compiled" fallback), so this no longer crashes — but a few
older tests under `tests/` predate that fix and still manually exclude v3 for
2-site chains via `_helpers.py`'s `versions=` kwarg (see
`test_spin_chain.py`/`test_fermion_chain.py`'s dimer tests,
`test_entanglement.py`'s Bell-pair test, and `_helpers.py`'s own docstring)
— those exclusions are now redundant defense-in-depth rather than required
to avoid a crash, not something to "fix" by removing.

### Cross-backend performance benchmarking

**Before benchmarking anything, check BLAS thread counts.** DMRG is a great
many *small* dense linear-algebra calls (a two-site tensor at `maxm=30` on
spin-1/2 sites is 60x60), and at that size more BLAS threads is slower, not
faster: measured 1 thread vs 2, `svd` 1.6x and `eigh` 2.3x slower. End to end
that alone is minor (~1.13x for `"python"`, nothing for v3), but under
*oversubscription* — a shared cluster node, or several dmrgpy runs in
parallel, where every process's BLAS thinks it owns the machine — it
dominates every other effect in this file: on a 14-core host with another job
holding 10 cores, MKL's default thread count turned a 9.4s `"python"` solve
into 28.2s and a 1.6s v3 solve into 26.7s. Pin threads before importing numpy
(`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ...`) when timing anything, or
timings across runs are not comparable. `src/dmrgpy/blasthreads.py` has the
full measurements and an opt-in `limit()` context manager (needs
`threadpoolctl`; the env-var route is more effective, since some BLAS
libraries fix threading on first use). dmrgpy never changes thread counts on
its own.


`benchmarks/run_benchmarks.py` answers a different question than
`tests/`/`examples/`: not "is this backend correct" but "which backend is
fastest, and by how much, for each kind of calculation" — the thing to
check before deciding which mode is worth optimizing next. It runs on
every backend available in the current environment (`v2`, `v3`,
`"python"`, `"julia_live"`, auto-detected via `cppext.available()` and a
live Julia import, same pattern as
`examples/groundstate/dmrg_generalized_benchmark`), and sweeps two
families of calculation over a uniform S=1/2 Heisenberg chain:

- **finite**, swept over chain length: ground-state energy, static
  correlator, KPM dynamical correlator, TD (real-time + Fourier)
  dynamical correlator, first excited state, bond entanglement entropy,
  real-time evolution after a quench, and the ground state confined to
  the `Sz=0` conserved sector. Accuracy is cross-checked against ED at
  every size cheap enough for it (energies, correlator, entropy profile
  and evolution trajectory alike — the ED entropy reference is an SVD of
  the ED state vector, since `entropy.py` is `self._session`-only and has
  no ED path of its own).
- **infinite**, swept over bond dimension (`infinitechain.py`): energy
  density and `<Sz>` under both `gs_method="idmrg"` and
  `gs_method="vumps"`, checked against the Bethe-ansatz density and
  against `<Sz>=0`.

The catalogue of what is covered lives in one place, `benchmarks/calctk.py`'s
`CALCS` list; adding a capability to the benchmark is one entry there plus
a label in `reporttk.CALC_LABELS`. Calculations that only some backends
implement (the sector, the infinite chain — both `itensor_version=3` and
`"python"` only) are part of the sweep, restricted to those backends, and
show up as blank-by-construction rather than as failures in the report's
coverage table. Two of the sections exist specifically because a single
number would misrepresent them: the sector section tabulates the
sector-vs-dense ratio per size (block sparsity scales, its bookkeeping
doesn't), and the `<Sz>` sections catch a wrongly gauged iDMRG unit cell,
which an energy-only check cannot see.

It writes `benchmarks/output/report.tex`
(tables + matplotlib plots + an auto-generated "which backend is worth
optimizing" paragraph) and compiles it to `report.pdf` with `pdflatex` if
available; `benchmarks/output/` is gitignored, so nothing under it is
checked in. Run `python run_benchmarks.py --help` from `benchmarks/` for
the full option list (sizes, bond dimensions, `--calcs` selection, backend
selection, DMRG/KPM/evolution parameters, per-calculation size caps,
repeats-per-timing); `--from-json results.json` regenerates the report
from a previous run's raw data without rerunning the sweep. Note that
`--repeats>1` re-prepares each calculation from a fresh chain per repeat:
the DMRG session caches its ground-state energy (`groundstate.py`'s
send-cache), so timing a second `gs_energy()` on the same chain would
otherwise record that backend as instantaneous. Because
juliacall JIT-compiles each Julia method the first time it's called in a
process (see the `julia_live` testing note above), the script always runs
one untimed warm-up call on the Julia backend before the timed sweep, so
its reported numbers are steady-state rather than dominated by one-time
compilation cost.

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
- `spinfermionchain.py`, `bosonchain.py`, `parafermionchain.py` — same pattern for other statistics/models

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
it silently falls back to ED. `itensor_version="python"` never falls back
for that reason (it has no compiled-extension precondition at all), but
does below 2 sites: its DMRG is two-site, so a 1-site chain has no update
to make and `dmrg()` returned `None` before that fallback existed — the
same shape as, but not the same threshold as, ITensor v3's own `ns<3`
fallback. `groundstate.py::gs_energy_generalized` needs its own copy of
that guard (it has no ED fallback to be routed to), and there it *raises*,
because the symptom was worse than `None`: the outer self-consistent
iteration still returns the Rayleigh quotient of a state no sweep ever
touched, i.e. a silently wrong number (-0.3049 for an exact -0.5). Note
that guard sits *before* the non-Hermitian dispatch, unlike the
`itensor_version==3` one below it — NH-DMRG escapes v3's abort by never
calling `dmrg()`, but its own sweep is two-site too, so it does not escape
this.

**Which backend a chain gets when the caller names none** is
`cppext.default_backend()`, not `DEFAULT_ITENSOR_VERSION` directly: the
default C++ version when that extension is compiled, and `"python"`
otherwise. That fallback is what makes `pip install dmrgpy` usable — the
wheel ships no C++, so before it every default-backend chain in a pip
install silently ran ED, which cannot reach the sizes an MPS backend
exists for. Only the *implicit* choice moves: an explicit
`itensor_version=3` with no extension still falls back to ED (a caller who
named a version asked for that backend), and `mode="ED"` is untouched.
`Many_Body_Chain.__init__`/`Mixed_Spin_Fermion_Chain.__init__` take
`itensor_version=None` to mean "pick one for me". See
`tests/test_default_backend_without_cpp.py`. Most
public methods on `Many_Body_Chain` accept a `mode="DMRG"|"ED"` kwarg so
results can be cross-validated between solvers (see the
bilinear-biquadratic example in `README.md`).

- **DMRG, C++ (`itensor_version=2` or `itensor_version=3` (the default))**:
  entirely in-process, see "In-process pybind11 extension" below. Both
  versions expose the identical `Chain` session API and public
  `Many_Body_Chain` surface — `itensor_version` only picks which compiled
  extension (`mpscpp2._dmrgcpp` vs `mpscpp3._dmrgcpp`, via
  `cppext.get_backend(version)`) backs `self._session`.
  `cppext.DEFAULT_ITENSOR_VERSION` (currently `3`) is the single source of
  truth for the default, threaded through `Many_Body_Chain.__init__`,
  `setup_cpp()`, `get_backend()`/`available()`, and `install.py`'s own
  `--itensor-version` default (`installtk/__init__.py`'s matching
  constant) — change it there, not per call site.
  `Many_Body_Chain.setup_cpp(version=2)` switches an existing chain
  between them. There is no file-based/subprocess fallback for either —
  if the requested extension isn't compiled, `mode.py` routes to ED
  instead.
- **DMRG, pure Python (`itensor_version="python"`)**: `pyitensor/`
  (`Many_Body_Chain.setup_python()` switches an existing chain to it) is a
  from-scratch, pure-Python/NumPy/SciPy reimplementation of exactly the
  ITensor v3 API subset `mpscpp3/chain_session.h` uses — Index/ITensor
  tensor core, the same ten site types, an AutoMPO term compiler with its
  own from-scratch Jordan-Wigner threading, MPS/MPO algebra, two-site
  Lanczos-based DMRG (ground + excited states via the overlap-penalty
  method), two-site Krylov-based TDVP, and a `pyitensor/chain.py::Chain`
  exposing the identical method surface as the compiled backends'
  `self._session`, so every `itensor_version in (2,3)` call site elsewhere
  in this codebase treats `"python"` as a third, always-available option
  wherever such a site was updated to include it — no separate code path.
  Exists so DMRG/TDVP work with zero compiler/pybind11 dependency, at the
  cost of being substantially slower than compiled ITensor (no
  block-sparsity, no JIT) — see `pyitensor/__init__.py`'s docstring and
  the module-level docstrings throughout `pyitensor/` for the specific
  simplifications taken versus real ITensor v3 and why each one doesn't
  affect dmrgpy's own results, only internal performance/bond-dimension
  efficiency.
  **MPO construction is no longer one of them.** `pyitensor/mpobuilder.py`
  used to build one exact bond-dimension-1 MPO per term and compress the
  concatenation of all of them — right answer, right *final* bond
  dimension (5 on a nearest-neighbour Heisenberg chain), but reached via
  an intermediate of bond dimension ~3(L-1), so O(L) truncating SVDs on
  O(L)-sized matrices: **O(L^4)**. Measured on an S=1/2 Heisenberg chain
  at default `maxm=30`/`nsweeps=15`, one `gs_energy(mode="DMRG")` scaled
  as ~L^3.7 and at L=100 spent 95% of its time building the Hamiltonian
  and 5% solving it. `to_mpo` now assembles the finite-state machine over
  the terms' partial products directly (the shape ITensor's own
  `toMPO(...,{"Exact",true})` produces), keeping the same two truncating
  sweeps afterwards purely to honour `cutoff`/`maxdim`. Same operator,
  same final bond dimension, and a build that is O(L^2) rather than
  O(L^4) (the machine itself is O(L); `HTerm.resolve()` spells each of the
  O(L) terms out over all L sites, and that term is cheap but real):
  `to_mpo` alone went 0.12s →
  0.027s at L=20, 2.96s → 0.12s at L=60 and **24.6s → 0.235s (105x)** at
  L=100 (min of 5, threads pinned). Two rules in there are load-bearing
  and easy to get subtly wrong — the coefficient goes on the transition
  *into* the final state (so terms differing only by a coefficient still
  share partial states), and transitions into that state **accumulate**
  while structural ones are **assigned** (two syntactically identical
  terms trace the same path, so their coefficients must sum while their
  shared structure must not be written twice). The old construction is
  kept as `_sum_of_term_mpos`, as the reference the machine is tested
  against on a term zoo in `tests/test_mpo_automaton_builder.py`; that
  rewrite also fixed a real bug it had, silently dropping all but the last
  term on a **1-site** chain (no bonds to sweep, so the concatenation was
  returned as-is).
- **DMRG, Julia (`itensor_version="julia_live"`)**: a live, in-process
  Julia session (`mpsjulialive/`, via `juliacall`/`juliasession.py`) with its
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
- **`mpscpp3`-specific: every site is built with `ConserveQNs=false` by
  default** (`mpscpp3/get_sites.h`), and `Chain`'s DMRG starting state is an actual
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
  Note the "permanently stuck" observation is specific to the
  *all-first-basis-state* product start it was measured with: that state
  (all spins up) is the unique state of its own Sz=N/2 sector, so no noise
  term can move it and its energy is the correct answer for that sector.
  It says nothing about QN mode from a *non-trivial* in-sector start,
  which is what the opt-in sector mode below actually does.
- **`mpscpp3`-specific: opt-in conserved-sector (QN) mode.**
  `Many_Body_Chain.set_conserved_sector(Nf=6)` / `(Sz=0)` / `(Nf=8,Sz=0)`
  (off by default, no-args clears it) confines
  the whole calculation to one quantum-number sector: `Chain::
  set_conserved_sector` rebuilds `sites_` via `SpinX(site_types,conserved)`
  with QN-carrying indices and drops everything built on the old site set.
  The same public API (including `promote_to_dense`/`promote_mps`) is also
  implemented by `itensor_version="python"`, by a different mechanism --
  dense storage with a charge grading on the site indices plus a charge
  penalty on the variational solves, see `pyitensor/sector.py`'s module
  docstring and the pure-Python-backend section of `docs/documentation.md`;
  everything else below is v3-specific. `mode="ED"` implements the same
  API too, by a third mechanism (`edtk/edchain.py`: a conserved charge is
  diagonal in the ED product basis, so a sector is a set of basis states
  and every assembled operator is restricted to that submatrix; the
  charge operators come from `Many_Body_Chain.
  get_sector_charge_operators()`, and `promote_to_dense`/`promote_mps`
  are deliberately *not* part of it, since an ED re-solve is cheap). The
  backends that genuinely have no quantum numbers are DMRG on
  `itensor_version=2` and `"julia_live"`, and `mode.py` raises rather
  than letting one of those answer a sector-mode chain; falling back to
  ED, on the other hand, is now correct rather than forbidden. See
  `docs/documentation.md` §4.3a and
  `tests/test_sector_conservation_ed.py`.
  It pays off with size, and only with size: measured through the Python
  API on Heisenberg chains, sector vs dense is 0.6x (n=20, maxm=60), 2.0x
  (n=40, maxm=100), 4.2x (n=60, maxm=200) — block sparsity scales, its
  per-block bookkeeping doesn't. Don't quote a single speedup number. Three ITensor
  behaviors shape the implementation, each confirmed directly and none of
  them a clean failure:
  - `randomMPS(SiteSet)` refuses to guess a sector under QNs and
    `randomMPS(InitState,m>1)` does not exist here (`mps.cc`), so the
    start state is a normalized `sum()` of random in-sector product
    states (`sector_mps`), with the per-site assignment coming from an
    exact DP over reachable partial charges (`sector_state_plan` — greedy
    dead-ends on a Hubbard chain at fixed `Nf` *and* `Sz`).
  - Building a non-QN-definite operator (`Sx` on Sz-conserving sites)
    aborts the process inside `ITensor::set`, and AutoMPO aborts on a
    flux-violating term *even when its coefficient cancels exactly*
    ("Index does not contain given QN block") — which is why the textbook
    `Sx*Sx+Sy*Sy+Sz*Sz` Heisenberg Hamiltonian needs `mo_terms.h`'s
    `expand_xy_terms`+`combine_terms` normalization here rather than
    being left to AutoMPO, and why `op_charge()` infers an operator's
    charge from its *dense* matrix elements instead of building it on QN
    indices. Every `terms → MPO` path routes through
    `mpo_from_terms`/`ampo_from_terms` so the check covers `vev`/
    correlators/KPM vertices too; both are pass-through when no sector is
    set, keeping the default path byte-identical.
  - ITensor's own `Error()` calls `abort()` (`util/error.h`), so nothing
    it raises is catchable from Python — the sector code throws
    `std::invalid_argument`/`std::runtime_error` instead.
  Python side: the sector lives on the chain, is re-applied wherever a
  session is rebuilt (`sites.py::initialize`, `__deepcopy__`), is part of
  `groundstate.py`'s send-cache key, and makes `mode.py` raise rather than
  fall back to ED (a fallback would silently answer with the *global*
  ground state). `tests/test_sector_conservation.py` cross-checks every
  sector against ED restricted to that sector — note the reference must
  *restrict* (the charge operator is diagonal in the ED basis, so a sector
  is a submatrix), not filter full-spectrum eigenvectors by ⟨N⟩: sector-
  degenerate eigenvalues come back as arbitrary superpositions.
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
  `evolve_and_measure_dmrg()` — every one of those dispatches tests
  `itensor_version in (3,"python")`, so `"TDVP"` applies on the C++ side
  only at `itensor_version==3` (mpscpp2 has no `TDVP/` and no `_tdvp`
  methods at all); any other combination silently falls back to the
  MPO-Taylor path, which remains the only option for
  `itensor_version=2`. Note `Chain::tdvp_step` takes `DoNormalize` from
  the step (`dt.imag()==0.0`), not hardcoded: forcing unit norm after a
  COMPLEX-time step deletes exactly the decay `tdz.py`'s `submode="TDZ"`
  contour is built on, which is what made v3's TDZ return too much
  spectral weight — the audit measured 36% too much, 0.340676 against an
  exact 0.25 (2026-09 audit finding #7; that "before" figure is the audit
  record, not re-measurable now that the extension has been rebuilt,
  while the post-fix 0.249473 re-measures). Only the
  two-site TDVP algorithm is used (`NumCenter=2`); the global subspace
  expansion machinery in `TDVP/basisextension.h` (mainly useful for
  one-site TDVP or long-range Hamiltonians) is not wired in, since
  two-site TDVP already grows bond dimension via SVD like two-site DMRG.
- **`mpscpp3`-specific: `Chain::idmrg_ground_state` is a full port of
  `pyitensor/idmrg.py`, static observables included.** It was for a long
  time only *most* of one, and the three missing pieces are worth knowing
  about because they are exactly the parts an energy-only test cannot
  see: (1) McCulloch's wavefunction prediction (`theta = lambda_k . B_k .
  lambda_{k-1}^{-1} . A_k . lambda_k`) as each micro-step's Krylov start,
  now `idmrg_wavefunction_prediction` (the older heuristic it replaced —
  reusing a raw flat vector whenever its size matches — survives as the
  fallback when the previous step's shapes don't line up yet); (2)
  extraction of the converged unit cell from a single micro-step's theta,
  now `idmrg_theta_cell` + `ic_canonicalize_cell`, kept on the Chain as
  `idmrg_cell_` — tiling the per-micro-step `idmrg_U_` chain instead
  (whose two ends live in bond bases minted by *different* micro-steps)
  leaves the energy correct while putting `<Sz>` at -0.13 against an exact
  0 on the XX chain; (3) the per-site energy baseline subtracted from both
  growing environments each micro-step (`idmrg_subtract_energy_baseline`,
  the reference `idmrg.h`'s `HL += -energy*IL`), with the density adding
  the shift back. On top of (2), `bindings.cc` now exposes
  `idmrg_onsite_expectation`/`idmrg_two_point_correlator`/
  `idmrg_local_excitation_gap`, so `gs_method="idmrg"` supports
  `vev`/`correlator`/`local_excitation_gap` on `itensor_version=3`;
  `td_dynamical_correlator_window`'s own connected-background subtraction
  reads the same cell-based expectation value (it used to tile `idmrg_U_`,
  i.e. it had the gauge bug (2) exists to fix).
  `tests/test_infinite_chain.py::test_itensor_version3_matches_python_backend`
  now compares `vev` and a `correlator` sweep as well as `e0`, so the
  divergence is guarded rather than merely documented, and
  `tests/test_idmrg_correlator_v3.py` covers the new surface directly.
  Still `itensor_version="python"`-only: `local_excitation_gap(window>0)`
  (`local_excitation_gap_windowed` is an explicit prototype, not stable
  API) and the iMPS algebra over `IDMRGResult`
  (`imps_overlap`/`apply_mpo`/`imps_sum`), which `infinitechain.py` never
  exposed for `gs_method="idmrg"` and which the VUMPS path already covers
  on this backend.
- **`mpscpp3`-specific: the iDMRG static observables are dense arrays, not
  ITensors** — the `ic_*` helpers in `chain_session.h`, the same design
  (and for the same reasons) as the VUMPS `vx_*` ones next to them. Two
  deliberate departures from a literal transcription of `idmrg.py`, both
  pure performance and neither changing a returned number: the
  measurement-independent part (transfer tensors + both families of
  fixed points) is built once per ground-state run and cached
  (`ic_build_cache`) rather than recomputed inside every `vev`/
  `correlator` call the way `idmrg.py` does; and a two-point correlator
  applies the operator string to the right fixed point one site at a time
  (O(chi^4) per site) instead of composing transfer matrices (O(chi^6)
  per site). Measured on a critical Heisenberg chain at `maxm=24`: a
  single `r=3` correlator went from ~4.1s to ~0.01s. Above
  `ic_dense_eig_max_` (=64, matching `idmrg.py`'s own `_DENSE_EIG_MAX`)
  the transfer-matrix eigenproblems switch from a dense `zgeev` on the
  whole `chi^2 x chi^2` matrix to a matrix-free restarted Arnoldi
  (`ic_arnoldi_dominant`), the C++ counterpart of `idmrg.py`'s own ARPACK
  branch. Net effect versus `itensor_version="python"` on that same
  chain: the *growth loop* is no faster (v3 runs ~1.5-2x slower — it is
  ITensor-object-bound, see the benchmark note in
  `docs/documentation.md`), while the *observables* are 5-500x faster.
- **Both backends can be loaded in the same process** (needed for anything
  that directly compares `itensor_version=2` vs `=3` results, e.g. the
  `v2_VS_v3_*`-named examples in each theme folder, or `examples/
  dynamical_correlator/dynamical_correlator_v2_VS_v3`) only because
  `mpscpp3/bindings.cc`
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

`submode="SECTOR"` (`sectordc.py`) is the exception to that shape, and
the only submode that does not compute a spectrum on the chain it is
called on: it solves the reference sector and the *one* sector
`B|gs>` lands in (N+1 for `Cdag`, Sz+2 for `S+`) on an internal clone,
promotes both to dense indices and contracts them, summing an explicit
Lehmann representation over converged per-sector eigenstates. It rests
entirely on `promote_mps` rebasing onto the chain's *original*
`dense_sites_`, which `set_conserved_sector` never rebuilds -- already
regression-tested by
`tests/test_sector_promotion.py::test_states_promoted_from_different_sectors_are_comparable`.
`get_spectral_function`/`get_spin_spectral_function` assemble the
particle/hole and Sz channels on top, and a charge-indefinite operator
(`Sx`) is split into `(S+ + S-)/2` and summed channel by channel
(`multioperatortk/charge.py`) rather than rejected. Two things to know before touching
it: the gating is explicit inside the submode (`mode.py`'s sector guard
keys on the *caller's* `conserved_sector`, which is empty here, so it
never fires), and the per-sector solves are cached with the clone that
produced them (`_sector_states_cache`) because ITensor Index identity is
per session -- that cache is also what makes a site sweep, and so
`S(q,w)`/`A(k,w)`, cost about one site's correlator. See
`docs/sector_resolved_spectral_function_plan.md` for the design and
`tests/test_sector_spectral_function.py` for what is pinned.

### Julia vs C++ backend

The two DMRG backends are independent implementations, not a shared
protocol: the C++ path runs in-process through the pybind11 extension (see
above), while `itensor_version="julia_live"` (`mpsjulialive/`) drives a
live Julia/ITensors.jl session via `juliacall` with its own mirrored set of
modules. A feature missing on one side simply isn't implemented there yet
(check the relevant `mpsjulialive/*.py` file, or lack thereof) rather than
falling back automatically — the only automatic fallback is DMRG → ED when
the C++ extension isn't compiled (`mode.py`).

### ITensor library
the ITensor folders (`mpscpp2/ITensor` = v2, `mpscpp3/ITensor` = v3) are a library that is developed elsewhere, hence you do not need to read them carefully. Read them only if there is some feature that does not work with the mpscpp2/mpscpp3 code, and you need to figure out why

**Upstream reference for the infinite-chain code**: ITensor's own
infinite-MPS package, <https://github.com/ITensor/ITensorInfiniteMPS.jl>,
implements iDMRG, VUMPS, the infinite MPO/`InfiniteMPS` types and the
transfer-matrix fixed-point machinery that `infinitechain.py`,
`pyitensor/idmrg.py`, `pyitensor/vumps.py` and the `ic_*`/`vx_*` helpers in
`mpscpp3/chain_session.h` reimplement here. It is *not* vendored in this
repo and dmrgpy does not depend on it — consult it as the reference
implementation when an infinite-chain algorithm (unit-cell gauging,
McCulloch wavefunction prediction, environment energy subtraction,
fixed-point conventions for `<Sz>`/correlators) looks wrong or ambiguous,
the same way `mpscppN/ITensor` is consulted for the finite algorithms.

One thing has already been ported from it: `pyitensor/vumps.py`'s
`_subspace_expand` (after upstream's `subspace_expansion.jl`) is how the
`"python"` VUMPS D-ramp warm-starts each rung -- it picks the new bond
directions from the SVD of `H` projected into `AL`/`AR`'s null spaces,
instead of the noise-padding `_grow_initial_state` still does on
`itensor_version=3`. The expansion cannot change the state or the energy
(both enlarged tensors stay exactly isometric, `C` is embedded with zeros
in the new block), only how fast the next rung converges: a wash on most
models, and 4.6x faster plus 5x more accurate (and run-to-run
reproducible, which the noise start is not) on the hardest one tested, a
gapless `D=8` Heisenberg chain. `pyitensor/vumps_ms.py` (the sequential
solver) now has the per-bond version too -- a loop over every bond of the
cell, which is the shape upstream's own `subspace_expansion(psi, H)` has.
Two things generalize rather than repeat: bond `n`'s two null spaces come
from *different* tensors (`AL[n]` and `AR[n+1]`), so `AR[n]`'s new
directions come from bond `n-1`; and the cell's last bond straddles the
cell boundary, where the right environment/automaton/tensor are site 0's
of the next cell -- by periodicity the same `GR[0]`/`W_list[0]`/`AR[0]`
this cell already holds. Every bond grows by the same amount so the state
stays uniform-`D`. Measured the same way (`nrestarts=1`, median of 5): a
wash on TFIM and on the easy rows, 4.3x faster on gapless Heisenberg
`n_uc=2 D=8` (1701 -> 909 iterations, 21.4s -> 5.0s), and on `n_uc=1 D=8`
(where neither start converges within `maxiter`) 1.4x faster with the
accuracy spread ~5x tighter at both ends. `itensor_version=3`'s own `vms_ground_state` has not picked up
subspace expansion -- it noise-pads (`Chain::vms_grow_init`, the per-site
`vumps_grow_init`). It DID have neither: its ramp's `reuse` test compared
the previous rung's tensor size against the *new* `D`, so it could never
hold across a rung and every rung started from pure noise. That is not a
speed nicety on a model whose exact state is smaller than the requested
`D` -- see the redundant-bond-dimension note further down. The residual
criterion, the AC/C phase alignment and the variational safety net are now
ported too.

The sequential solver is also now the route for **couplings reaching
further than one unit cell**, which used to be rejected outright
(`infinitechain._canonicalize_hamiltonian` refused any term touching both
the previous and the next cell). A Hamiltonian with reach>1 has no grouped
reach-1 representation on that cell at all, so `vumps_ground_state` sends
it to `vumps_ms` regardless of `n_uc` -- `idmrg_excitations.is_reach_one`
read off the grouped automaton is the test, `_check_reach_one` is now just
the raising wrapper around it, kept for the excitation ansatz. That is
what makes a long-range model affordable: range-R costs the sequential
path R extra automaton channels (linear), where the only previous route --
rewrite the chain on an `n_uc >= R` cell so every coupling is reach-1 --
costs the grouped path a `d**R` supersite. `get_operator(name, i,
group=c)` now takes an integer cell offset (`c = -1, 0, 1` being `L`, `C`,
`R`) to write one, and `_canonicalize_hamiltonian` canonicalizes a term's
position by translating it whole onto the cell its leftmost site lives in
(well defined at any reach, identical to the old three-cell rules
everywhere those applied). Two things worth knowing before touching this:
`gs_method="idmrg"` needed no change at all -- its growth loop consumes
whatever `_build_periodic_mpo` builds, and `_active_channels_at` has
always carried one pending channel per site of a term's reach, confirmed
directly against exact answers, so do NOT "fix" that by adding a guard;
and `_window_hamiltonian` (`kpm_finite`'s open-boundary tiling) now tests
the fit per term rather than dropping one fixed copy, because a
longer-ranged term runs off the end of the window sooner.
`itensor_version=3` has the same dispatch (`Chain::vumps_ground_state`'s
`use_multisite`, on `vumps_is_reach_one`) and agrees with `"python"` to
2e-13. `vev`/`correlator` follow along on it too, at any `n_uc` and
any reach. They briefly did not, and the reason is worth keeping: the
C++ picks the SOLVER on `n_uc>2 || reach>1` while `infinitechain.py`
picked the READER on `n_uc>2` alone, so a reach>1 chain on a short cell
ran the sequential solver and was then handed to the grouped reader,
which had no snapshot -- even though `Chain::vms_onsite_expectation`/
`vms_two_point_correlator` could answer it and had existed since
`b81b6a9`. `Chain::vumps_onsite_expectation`/`vumps_two_point_correlator`
now fall through to the `vms_*` pair on `have_vms_snapshot_`, i.e. the
snapshot that exists decides the reader, and the Python side has no
dispatch of its own left to drift. `tests/test_infinite_long_range.py`,
`examples/idmrg/long_range_infinite_chain`.

**Redundant bond dimension used to be fatal on `itensor_version=3`.** A
gapped model's exact state often needs fewer directions than the requested
`maxm` (a field-polarized chain needs exactly one), the extra ones then
carry no weight, and the transfer matrix picks up a decoupled unimodular
block -- a genuinely DEGENERATE dominant eigenvalue, which looks identical
to the one thing `vx_check_perron_nondegenerate` exists to reject (a "cat
state": two branches with matched *nonzero* weight). So `gs_energy()`
raised "every attempt at D=... failed" for any such model at `maxm>1` on
the sequential solver, and the grouped one survived only by never landing
exactly on the degeneracy (its second eigenvalue was measured at 0.99996,
just outside the guard's 1e-9). Both environment builders now fall back to
the fixed points the state itself names, `C C^dag`/`C^dag C`
(`Chain::vx_bond_fixed_points`) -- exact under redundancy, and for a real
cat state the same branch mixture an unchecked eigensolver would give,
where the eigensolver may instead return an arbitrary single branch with
its own wrong energy. Note what was tried first and is the WRONG shape: a
threshold on `C`'s weight spectrum to tell redundancy from a cat state.
The guard trips mid-convergence, where the redundant direction is still on
its way down, so the ratio to catch is a moving number (2.7e-9, 1.2e-4 and
1.6e-2 on three cells of the same model) with no defensible cutoff. The
pure-Python reference has since gone one step further, and the ORDERING is
the point: both its environment builders
(`vumps_ms._cell_fixed_points` on the sequential path,
`vumps._transfer_fixed_points` on the grouped one) now *prefer* the bond
candidate -- `C C^dag` for the AL transfer's right fixed point,
`conj(C^dag C)` for the AR transfer's left one in this codebase's
`X[ket,bra]` ordering -- whenever it reproduces itself under the transfer
map to a residual <= 1e-6, which in mixed canonical gauge is an exact
algebraic identity and therefore a yes/no test rather than a tuned
threshold, so at convergence no eigensolve runs at all; the guarded
eigensolver runs only when the residual says the gauge relation does not
hold, and the two are cross-checked against each other by residual. That
is bond-candidate-FIRST, and **`itensor_version=3` got there on
2026-09-12** -- `Chain::vx_choose_fixed_point`, one shared selection that
both `vumps_build_environments` and `vms_environments` go through, where
each used to hold its own `catch (ITError const&)` reaching
`vx_bond_fixed_points` as a failure fallback. The difference between the
two orderings was measurable and was not in the C++'s favour: on a
field-polarized reach-1 one-site cell at `maxm=4`, `gs_energy()` on
`itensor_version="python"` raised "every attempt at D=4 failed" in 7 of
20 runs before its own change and 0 of 20 after, while
`itensor_version=3` returned energies *below* the exact variational
minimum at a rate that moved between runs and builds (13 and 17 of 80
measured at D=6 on the grouped cell, worst 2.3e-8; worst 2.8e-3 on the
sequential one) where `"python"` did so in 0 of 40. It is now 0 of 330
past 1e-12 across six cells. NUMBERS CHANGE, on the same terms the
`"python"` half's own entry sets: any `itensor_version=3` VUMPS energy
from before 2026-09-12 on a model whose exact state is smaller than the
requested `maxm` is not comparable, on either solver -- those were wrong,
not merely less converged. Energies from runs that were already fine move
only at ~1e-15. Quote the threshold with any rate, as the
known-issue file itself instructs -- see
`docs/known_issue_v3_vumps_variational_floor.md`, now FIXED, and
`tests/test_vumps_redundant_bond_dimension.py`, whose narrow `xfail` came
out and which gained a repeated-run guard (a 6-in-30 defect is invisible
~80% of the time to a one-run test).

**Porting that ordering is what caught a latent defect in the primitive
itself**, and the lesson generalizes: `vx_bond_fixed_points`' LEFT
candidate was `C^dag C` where the `X[ket,bra]` ordering needs
`conj(C^dag C)` -- the transpose. Nothing had caught it because nothing
had ever *measured* that candidate; adding the residual read it
immediately (4e-15 for the correct orientation against 0.38-0.53 for the
transpose, at convergence on the polarized cell at D=6). It survived
because the only two models this function had been exercised on -- a
field-polarized chain, whose converged `C` is real diagonal, and AKLT,
whose `C` is real -- both have a real symmetric `C^dag C`, where the two
orientations coincide. A fallback that is only reached on failure is a
fallback nothing validates.

**Reading that led to a bigger, unrelated find, and it is the one to know
about**: `pyitensor/dmrg.py`'s `_lanczos_ground_state` stops when the
lowest Ritz *value* stops moving. For a Hermitian operator the value
error is quadratic in the eigenvector error, so a value tolerance of
1e-12 caps the eigen*vector* at ~1e-6 -- fine for finite DMRG, which
reports energies, and fatal for VUMPS, whose convergence criterion is a
norm difference between `AC` and `C`, two independently-solved
eigenvectors. VUMPS's gauge mismatch therefore floored at ~1e-6 and
`tol=1e-10` was unreachable at any number of outer iterations, while the
*energy* was already correct -- which is exactly why nothing in `tests/`
caught it. `_lanczos_ground_state` now has an opt-in `residual_tol`
(Ritz residual `beta*|s_k|`, free from the tridiagonal problem) that
`vumps.py`/`vumps_ms.py` pass as `tol/10`; every other caller is
byte-identical. Alongside it, `AC`/`C` are now sign-aligned against the
vectors they were warm-started from -- solved independently, their
relative sign was arbitrary, and a flip made the mismatch read ~4 instead
of ~1e-6 in 11 of every 100 iterations. Measured through the public
driver: `D=8` TFIM went from `converged` in **0/3** runs to **3/3**, and
32.6s -> 6.7s (g=1.5) / 38.4s -> 4.6s (critical g=1.0), energies
unchanged to every printed digit. One consequence looked alarming and,
**measured, is not**: finite DMRG's own local solves carry the same ~1e-6
per-solve eigenvector cap, so anything consuming a DMRG *wavefunction*
rather than its energy was flagged here as worth a look. It has had one
(2026-09-12) and the answer is that the cap does not reach any consumer.
On an 8-site Heisenberg chain at FULL bond dimension -- where the MPS can
represent the exact state, so the local solver is the only error left --
energies, `<Sz_i>`, two-site correlators, the bond entropy, `<wf|wf>` and
the first three excited states all agree with ED to 1e-13..1e-15, and
forcing `residual_tol=1e-13` into finite DMRG's own two
`_lanczos_ground_state` call sites (`dmrg.py:421`/`:726`) moves every one
of them by <= 5e-15. Below full bond dimension the error is
truncation-dominated and the criterion makes no systematic difference
either way: swept over `maxm` in 4..32 against ED, the loose/tight
difference in the excited-state errors is non-monotonic noise (tight is
slightly *worse* at maxm=8 and 12), not an improvement. The mechanism is
the sweep: each local solve is warm-started from the previous sweep's own
tensor and the sweep is iterated to convergence, so a per-solve vector
error contracts instead of accumulating -- which is exactly what VUMPS
lacks, its criterion being a difference between two *independently*
solved eigenvectors with no outer contraction. So do not "fix" finite
DMRG's criterion on this reasoning; it costs Krylov iterations and buys
nothing measurable. The C++ `itensor_version=3` port carried the eigenvalue criterion for
a while after this, so `Chain::vumps_ground_state` floored at ~1e-6 and
reported `converged=False` at `D>=8` -- a documented divergence between the
two VUMPS backends. `Chain::vx_lanczos_ground_state` now has the same
`residual_tol` (both VUMPS drivers pass `tol/10`; every other caller is
byte-identical) and `vx_align_phase` the same alignment, removing the whole
phase rather than just the sign, since a vector from the dense
`vx_dense_eig_max_` branch comes out of `zheev` with no fixed phase at all.
Re-measured through the public driver at `D=8`, `nrestarts=6`: `converged`
in **3/3** runs at both `g=1.5` and critical `g=1.0`, from 0/3.
`tests/test_lanczos_residual_criterion.py` covers both backends. See also
`docs/documentation.md`'s D-ramp section,
`tests/test_vumps_subspace_expansion.py`, and the measurement tables in
`vumps.py`'s own section comments.

That fix also **overturned a previously-recorded diagnosis**, which is
worth knowing before trusting a similar one: the tangent-space excitation
ansatz's accuracy loss at redundant bond dimension (excitation error
1.8e-15 at `maxm=1` up to 3.5e-3 at `maxm=8` on an exactly-`D=1` state,
and a ~33% flake in the dispersion tests) had been investigated at length
and attributed to conditioning alone -- `_random_raw_tensor`'s docstring
in `vumps.py` concluded "the ground-state solver is not the thing to
change". The amplification was real, but what it amplified was the ~1e-6
eigenvector error; the ground-state solver's *stopping criterion* was
exactly the thing to change. Re-measured, the error is now 1e-13-1e-12
across `maxm` = 1..16 and no longer monotonic. The regression that pinned
the degradation now pins its absence
(`test_redundant_bond_dimension_no_longer_costs_excitation_accuracy`),
as its own predecessor's docstring had asked for; both docstrings keep the
superseded reasoning rather than deleting it.

## Documentation

`docs/documentation.{md,tex}` covers the architecture (backends, dispatch,
directory layout); `docs/user_guide.{md,tex}` covers the physics-facing
API (what calculation each method performs, with the formula behind it).
Whenever a new feature is implemented — a new physical model, a new
`Many_Body_Chain` method, a new dynamical-correlator submode, a new
post-processing tool, etc. — update `docs/user_guide.md` **and**
`docs/user_guide.tex` accordingly (and `docs/documentation.{md,tex}` too,
if the change affects architecture/backend dispatch rather than just
physics). Keep the two formats of each document in sync; verify the
`.tex` file still compiles cleanly with `pdflatex` (no errors, ideally no
overfull-hbox warnings) after editing it.

## Git workflow

Working directly on `master` is fine at all times, unless the user says
otherwise for a given task or the change is a really major one with real
potential to conflict with other concurrent development (in which case
use a feature branch/PR instead).

When a feature branch is warranted, new work (a new branch, a new PR)
should start from `master`, not from another feature/staging branch,
unless explicitly told otherwise. This repo previously reused a
long-lived branch (`mpo-algebra-tdvp-gse`) as a staging area for several
stacked PRs; once `master` and the staging branch diverged (each gained
commits the other lacked), a PR built on the staging branch picked up
unrelated content and produced conflicting/misleading diffs against
`master`. Branching from `master` (`git fetch origin master && git
checkout -b <branch> origin/master`) avoids this.
