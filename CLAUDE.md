# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with
code in this repository. It holds what every session needs: how to build, how
to test, the architecture, and the rules the audits taught. Three parts of
the codebase have their own notes under `.claude/rules/`, loaded only when a
file they cover is read or edited: `cpp-backend.md` (the pybind11 extension,
`mpscpp2/` and `mpscpp3/`), `infinite-chain.md` (iDMRG and VUMPS on both
backends) and `kondo.md` (`kondospectrumtk/`). They load when the Read,
Write or Edit tool touches a matching file, not when a shell command prints
it, so a session that reads with `cat` or `sed` opens the matching one by
hand before touching those files; a pointer at the place each used to live
says where it went. The history of every fix, with its repro and its reviewer, is in
`docs/audit_*.md`, and the ledger of every result that changed is
`docs/numerical_changes.md`.

## What this is

DMRGPY is a Python library for quasi-one-dimensional spin chains, fermionic
systems, parafermion and bosonic models using matrix product states (MPS).
It exposes a single Python API (`src/dmrgpy/`) backed by interchangeable
solvers per calculation:

- **DMRG**, executed in-process via a pybind11 extension over ITensor
  (C++, either v2 or v3, see below), via a pure-Python reimplementation
  (`itensor_version="python"`), or via a live Julia/ITensors.jl backend
- **ED** (exact diagonalization), a pure-Python fallback used for small
  systems, for cross-checking DMRG results, and automatically whenever the
  C++ extension is not compiled

The Python layer is what this project actually develops. The C++/ITensor
backends under `src/dmrgpy/mpscpp2/ITensor/` (ITensor v2) and
`src/dmrgpy/mpscpp3/ITensor/` (ITensor v3) are vendored copies of the
ITensor library maintained upstream: treat them as black boxes, not code to
read or modify here. Read them only if some feature does not work with the
`mpscpp2`/`mpscpp3` code and you need to figure out why.

## Install / build

There is a `pyproject.toml`, but it is only for building the PyPI
distribution (see "Packaging / PyPI" below); it is not how this repo is
developed. For development the project is used from `src/` via a path added
to `PYTHONPATH` (or a symlink into `site-packages`), which is what
`install.py` sets up; there is no `setup.py` and no editable install.

```bash
python install.py                    # compiles ITensor v2 (mpscpp2) + its pybind11 extension
python install.py --itensor-version=3    # compiles ITensor v3 (mpscpp3) instead
python install.py --itensor-version=both # compiles both v2 and v3 backends
python install.py --gpp=g++-6        # use a specific compiler (needs g++ >= 6, LAPACK/BLAS)
python install.py --doctor           # only the requirement checks, no multi-minute ITensor build
python install_julia.py              # alternative: Julia backend instead of/alongside C++
```

`install.py` runs in two phases: `installtk/requirements.py` first checks
every build requirement (OS, `make`, a C++ compiler, LAPACK/BLAS,
`pybind11`), actually trial-compiling and trial-linking each one rather
than checking a version string, and only once everything checks out does
`installtk/install2.py` compile the vendored ITensor static library and the
`pybind` Makefile target for whichever backend(s) `--itensor-version`
selects (`2` = mpscpp2/ITensor v2; `3` = mpscpp3/ITensor v3, the default;
`both` builds both, one after the other). `install2.py` picks the right
`mpscppN` directory and C++ language standard per version (v2 needs only
C++14; v3 needs C++17 with the concepts TS extension, `-fconcepts`, see
`mpscpp3/ITensor/options.mk.sample`'s own `CCCOM` line) but otherwise runs
the identical compiler/BLAS/LAPACK configuration resolved in phase 1.
`pybind11` is a hard requirement (auto-installed via `pip` into the current
interpreter if missing, since the pybind11 extension is the only DMRG
backend besides ED and `"python"`).

The compiler is auto-detected via a cascade of candidates rather than a
single guess: `installtk/requirements.py::_compiler_candidates()` yields a
conda-provided compiler first when running under a conda Python
(`installtk/cppversion.py::find_conda_compiler()` looks for one, e.g. the
`gxx_linux-64` package, next to the running interpreter; preferred because
the extension is loaded into the same process as conda's own numpy/scipy,
which bundle their own libstdc++, and building against the system compiler
instead has reproducibly segfaulted in the past, see the long comment in
`mpscpp2/Makefile`), then falls back to the system `g++`/`c++` on PATH if
that candidate does not actually work (`--gpp` overrides the whole cascade
with a single hard choice). Each candidate is verified for real: besides an
actual trial compile and link (`cppversion.trial_compile()`),
`cppversion.backend_missing()` asks the driver (`gpp -print-prog-name=cc1plus`)
whether it can even locate its own compilation backend, since `--version`
and `shutil.which()` both pass for a compiler binary that is present but
incomplete. This is not hypothetical: Aalto's Triton cluster ships a
conda-style Python distribution with exactly such a broken
`x86_64-conda-linux-gnu-g++` wrapper (no matching `cc1plus`), and the
cascade falls back to the system compiler there and prints a note.
`installtk/cppversion.py::lmod_state()` detects environment-module systems
(Lmod, via `$LMOD_CMD`/`$MODULEPATH`) so failure messages can suggest
`module load ...` instead of sudo-based apt-get/brew hints; `check()` also
prints a reminder to reload the same modules or conda environment in any
later session that imports `dmrgpy`, since the compiled extension depends
on them being present again at import time.
`installtk/blaslapack.py::find_working_config()` similarly tries several
LAPACK/BLAS candidates (conda's own libs, system `-lblas -llapack`,
OpenBLAS, and a Debian/Ubuntu-multiarch workaround that asks the system
compiler where it actually finds `libblas.so`/`liblapack.so`) and picks the
first one that actually links. Finally, `install.py` adds the repo's `src`
to the Python path and to the user's shell rc file
(`installtk/addpythonpath.py`, `installtk/addsystem.py`, both idempotent
across reruns). There is no separate lint or CI config in this repo.

Both `mpscppN/Makefile`s list `$(wildcard *.h)` among the extension's
prerequisites: without that, editing `chain_session.h`, where the actual
DMRG implementation lives, leaves `make pybind` reporting "Nothing to be
done" while the `.so` keeps the old code. A C++ change is not tested until
the `.so` has been rebuilt.

## Packaging / PyPI

`pyproject.toml` builds a **pure-Python** distribution: a `py3-none-any`
wheel (~450KB) containing the ED backend, the `itensor_version="python"`
pyitensor DMRG/TDVP backend, and the Julia bridge. It deliberately ships no
C++ at all. `[tool.setuptools.packages.find]` excludes `dmrgpy.mpscpp2` and
`dmrgpy.mpscpp3`, which are ~400MB of vendored ITensor plus a build that
needs a compiler, LAPACK/BLAS and pybind11, none of which belongs in a pip
install. Users who want the compiled backends clone the repo and run
`install.py`; without them everything still works via `"python"`/ED
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
  This is unrelated to pyitensor's own opt-in numba/JAX matvec kernel
  (`pyitensor/kernels.py`), which is properly `try/except`-guarded and
  genuinely optional.
- **`requires-python = ">=3.10"` is set by numba**, which declares
  `>=3.10`. Claiming a lower floor does not make dmrgpy work on 3.9; pip
  just backtracks onto an ancient, untested numba/numpy pair.
- **`setuptools>=77` in `[build-system]`**, not the more common 61: the
  `[project]` table uses PEP 639 licensing (bare SPDX `license` string plus
  `license-files`), which older setuptools rejects outright.
- **The `julia` extra is `juliacall`, not `julia`.** The live backend
  migrated from PyJulia to juliacall/PythonCall.jl. The only surviving
  `import julia` is `juliarun.py`, part of the legacy subprocess path that
  is not reachable from the public API, so PyJulia is intentionally not
  pulled in.
- **`src/dmrgpy/juliapkg.json` must stay in `package-data`.** juliacall/
  juliapkg reads it at import time to resolve and precompile ITensors,
  ITensorMPS and ITensorNHDMRG; if it is missing from the wheel,
  `pip install dmrgpy[julia]` yields a Julia backend with no ITensors.jl.
  `.gitignore`'s unanchored `julia*` pattern would also swallow it, hence
  the explicit `!src/dmrgpy/juliapkg.json` un-ignore. Do not remove either
  one.
- **The classifiers say `OS Independent`**: the pip-installed package (ED
  and `itensor_version="python"` backends) works on Windows. Nothing in the
  Python layer shells out any more (`sites.py::initialize()` used to run
  `mkdir -p` via `subprocess`, removed in `1f7b7a9` along with the
  `.mpsfolder` it created and `Many_Body_Chain.clean()`, whose only job was
  removing that folder). The compiled C++ backends (`install.py`) remain
  POSIX-only and are not part of the wheel regardless.
- Every module in the wheel must at least byte-compile, or pip prints a
  `SyntaxError` on install; `python -m compileall` over the installed
  package is the cheap check.

Verify a release candidate by installing the built wheel into a clean
virtualenv (not from `src/`, so import bugs and missing package-data
actually surface) and cross-checking the two shipped backends against each
other: `itensor_version="python"` DMRG against `mode="ED"` on a small chain
should agree to ~1e-15. `dist/`/`build/` are gitignored; build artifacts are
not committed. A given version can never be re-uploaded to PyPI (or
TestPyPI) once published, so bump the version rather than trying to replace
a release.

## Running / verifying changes

`tests/` is a real pytest suite. Run it with `python run_tests.py` from the repo root (a thin wrapper: `pytest tests`
from the repo root, then `clean.py`) or plain `pytest tests`.
`tests/conftest.py` inserts this checkout's own `src/` at the front of
`sys.path`, so pytest always exercises the current working tree's code,
never a stale `site-packages`/symlinked checkout (see the worktree pitfall
below). `tests/_helpers.py`'s `energy_ed_v2_v3`/`vev_ed_v2_v3` are the
standard way most tests cross-check ED against DMRG on the C++ backends (a
`versions=` kwarg lets a test opt a specific `itensor_version` out; it
covers `(2, 3)` only, and the pure-Python backend is tested by the dozens
of files that parametrize their own `itensor_version` instead).
`tests/reference_data.py` holds shared golden values pinned to a specific
historical commit (its own module docstring explains how to regenerate
them); several test files inline their own golden values as local
constants instead, which is equally fine for values used in one place. When
a test's assertion is against exact DMRG/ED agreement (not a golden
regression value), prefer `pytest.approx(..., abs=1e-6)`-style tolerances
over exact equality, since DMRG is iterative. The audits' regressions live
in `tests/test_audit_<date>_<cluster>.py`, one file per fix cluster; each
`**Status**` line in the matching `docs/audit_*.md` record names the test
that pins it.

**`itensor_version="julia_live"` tests only need to run when a change
touches Julia code** (`mpsjulialive/`, any `*.jl` file, or a shared entry
point like `manybodychain.py`/`infinitechain.py` that dispatches into that
backend). These tests carry a large, mostly fixed cost unrelated to what
they are testing: `juliacall`'s JIT compiles each distinct Julia function
signature the first time it is called in a process, and this compile cost
dominates wall time regardless of how small the test's own parameters are.
Confirmed directly by isolating a `julia_live` `metts_dynamical_correlator`
call down to `nsamples=1`, `nwarmup=1`, `nt=1`: it still took ~73s, almost
entirely JIT, against ~27s of genuine algorithmic work for the same call at
its normal (nsamples=40) parameters once the function was already warm.
Reducing test parameters further cannot buy back that fixed cost, and it is
paid independently by every distinct Julia function signature exercised
across the whole suite, which is why the full suite's slowest handful of
tests are almost all `[julia_live]`-parametrized. For a change that does not
touch Julia code at all, skip them locally with
`pytest tests -k "not julia_live"` (the parametrize id `julia_live` appears in every such test's
node id); run the full suite including `julia_live` before finishing work
that does touch Julia code, or before a release.

The `examples/` directory (100+ self-contained scripts, one per physical
model or feature, organized into thematic subfolders) is a second,
complementary regression surface: broader and more exploratory than
`tests/` (most scripts just print and plot rather than assert), but
covering more ground (every physical model, every dynamical-correlator
submode, every backend combination including `"python"`). Some
`examples/*/*/main.py` scripts do use real `assert`s and are meant as
regression tests in the same spirit as `tests/` (e.g.
`time_evolution/tdvp_VS_ED_time_evolution`,
`backend_comparison/backend_switch_consistency`,
`finite_temperature/thermal_purification_VS_exact`,
`staticcorrelators/static_correlator_VS_ED`,
`entanglement/entanglement_entropy_VS_ED`,
`dynamical_correlator/dynamical_correlator_VS_ED`). To run one:

```bash
cd examples/readme_examples/static_correlator_S12 && python <script>.py
```

`examples/readme_examples/` mirrors the snippets shown in `README.md`.
`examples/backend_comparison/` and the `v2_VS_v3_*`-named examples inside
each theme folder directly compare the ITensor v2 and v3 C++ backends
against each other (same script, both `itensor_version`s, small systems),
the fastest way to check a `mpscpp3` change did not diverge from
`mpscpp2`'s numerics. One rule if you use one as a template: pin the sweep
schedule inside the script instead of inheriting the library defaults
(`boson_models/v2_VS_v3_boson` inherited them and its `assert` fired
non-deterministically; it now pins `nsweeps=80`/`maxm=100`, 2026-09 audit
#24). After running examples, `python clean.py` recursively removes
generated working directories (`.mpsfolder`, `.pychainfolder`,
`.dmrgfolder`) and stray `ERROR`/`*.OUT` files from the tree.

**Examples should plot, not just print/assert.** What sets `examples/`
apart from `tests/` is that a human is expected to look at the result, so
every `examples/*/*/main.py` ends with a `matplotlib` plot of whatever it
computed whenever the script produces a sequence of values; the scripts
that also carry a regression `assert` keep it and plot afterwards. A script
whose entire output is a single scalar is the one exception, and a last
resort: sweep a cheap axis (chain length, a coupling, time, bond dimension,
backend against backend) instead. An audit found 105 scripts that had
settled for print-only output with such an axis already computed.
**Examples are short and simple** (the author's rule since 2026-10-23): 40
lines, 50 at most, written in plain English the way a condensed matter
physicist who is not a programmer would write them, about the physics and the
one capability shown, nothing else. The
`new-example` skill has the full rules and the two approved models. The
`new-example` skill (`.claude/skills/new-example/SKILL.md`) has the full
rules, the path line, and how to run a script against this worktree.

**Worktree/symlink pitfall**: `site-packages/dmrgpy` is typically a symlink
into one specific checkout's `src/dmrgpy/` (usually the primary working
directory, not whichever worktree you are in;
`python3 -c 'import dmrgpy; print(dmrgpy.__file__)'` shows which). Since every
`examples/*/*/main.py` does `sys.path.append(...)` (append, not insert),
that symlinked package resolves before the appended path when running a
script directly with plain `python3 main.py` from a different worktree,
silently testing the wrong checkout's code with no error. Force the intended
checkout explicitly with
`PYTHONPATH=/path/to/this/worktree/src python3 main.py`, or
`sys.path.insert(0, ...)` in a one-off driver script.
`tests/conftest.py` avoids this problem entirely for anything run through
pytest.

**A real, confirmed `mpscpp3` bug**: ITensor v3's two-site `dmrg()` (see
`chain_session.h`'s `dmrg_args()`) aborts the whole process (`SIGABRT`,
`"LocalOp is default constructed"`, an internal ITensor v3 check in
`itensor/mps/localop.h`) for any chain with fewer than 3 physical sites,
independent of physics type. `mode.py`'s `get_mode()` falls back to ED
automatically for `itensor_version==3` with `self.ns<3` (same mechanism as
the "extension not compiled" fallback), so this no longer crashes; a few
older tests under `tests/` predate that fix and still manually exclude v3
for 2-site chains via `_helpers.py`'s `versions=` kwarg (the dimer tests in
`test_spin_chain.py`/`test_fermion_chain.py`, `test_entanglement.py`'s
Bell-pair test). Those exclusions are redundant defense-in-depth now, not
something to "fix" by removing.

### Cross-backend performance benchmarking

**Before benchmarking anything, check BLAS thread counts.** DMRG is a great
many small dense linear-algebra calls (a two-site tensor at `maxm=30` on
spin-1/2 sites is 60x60), and at that size more BLAS threads is slower, not
faster: measured 1 thread against 2, `svd` 1.6x and `eigh` 2.3x slower. End
to end that alone is minor (~1.13x for `"python"`, nothing for v3), but
under oversubscription (a shared cluster node, or several dmrgpy runs in
parallel, where every process's BLAS thinks it owns the machine) it
dominates every other effect in this file: on a 14-core host with another
job holding 10 cores, MKL's default thread count turned a 9.4s `"python"`
solve into 28.2s and a 1.6s v3 solve into 26.7s. Pin threads before
importing numpy (`MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 ...`) when
timing anything, or timings across runs are not comparable.
`src/dmrgpy/blasthreads.py` has the full measurements and an opt-in
`limit()` context manager (needs `threadpoolctl`; the env-var route is more
effective, since some BLAS libraries fix threading on first use). dmrgpy
never changes thread counts on its own. The same cap applies to agents: at
most three Python processes at once across all subagents of a hunt or
workflow (`run3.sh` in the 2026-09-24b record's "Shared helpers").

`benchmarks/run_benchmarks.py` answers a different question than
`tests/`/`examples/`: not "is this backend correct" but "which backend is
fastest, and by how much, for each kind of calculation". It runs on every
backend available in the current environment (`v2`, `v3`, `"python"`,
`"julia_live"`, auto-detected via `cppext.available()` and a live Julia
import) and sweeps two families of calculation over a uniform S=1/2
Heisenberg chain: **finite**, swept over chain length (ground-state energy,
static correlator, KPM and TD dynamical correlators, first excited state,
bond entanglement entropy, real-time evolution after a quench, and the
ground state confined to the `Sz=0` sector), cross-checked against ED at
every size cheap enough for it (the ED entropy reference is an SVD of the
ED state vector, since `entropy.py` is `self._session`-only and has no ED
path); and **infinite**, swept over bond dimension (`infinitechain.py`):
energy density and `<Sz>` under both `gs_method="idmrg"` and
`gs_method="vumps"`, checked against the Bethe-ansatz density and against
`<Sz>=0`. The catalogue of what is covered is `benchmarks/calctk.py`'s
`CALCS` list; adding a capability is one entry there plus a label in
`reporttk.CALC_LABELS`. Calculations only some backends implement (the
sector, the infinite chain, both `itensor_version=3` and `"python"` only)
are restricted to those backends and show up as blank-by-construction
rather than as failures in the coverage table. The sector section tabulates
the sector-against-dense ratio per size (block sparsity scales, its
bookkeeping does not), and the `<Sz>` sections catch a wrongly gauged iDMRG
unit cell, which an energy-only check cannot see.

It writes `benchmarks/output/report.tex` (tables, matplotlib plots and an
auto-generated "which backend is worth optimizing" paragraph) and compiles
it to `report.pdf` with `pdflatex` if available; `benchmarks/output/` is
gitignored. Run `python run_benchmarks.py --help` from `benchmarks/` for
the full option list; `--from-json results.json` regenerates the report
from a previous run's raw data. `--repeats>1` re-prepares each calculation
from a fresh chain per repeat, because the DMRG session caches its
ground-state energy (`groundstate.py`'s send-cache) and timing a second
`gs_energy()` on the same chain would record that backend as instantaneous.
The script always runs one untimed warm-up call on the Julia backend before
the timed sweep, for the JIT reason above.

## The audits, and what to know before touching the code

Eight hole hunts, a fix pass and one review ran between 2026-08 and
2026-10-23, each hunt a multi-lens search for numbers that are plausible and
wrong, with every finding reproduced, handed to a reviewer briefed to refute
it, and fixed:
`docs/audit_2026_08_hole_hunt.md` (21 findings), `docs/audit_2026_09_hole_hunt.md`
(36), `docs/audit_2026_09_24_hole_hunt.md` (16), `docs/audit_2026_09_24b_hole_hunt.md`
(18), `docs/audit_2026_09_24c_hole_hunt.md` (18), `docs/audit_2026_09_25_open_items.md`
(ten open items, a fix pass rather than a hunt) and
`docs/audit_2026_09_25b_hole_hunt.md` (29, all fixed on 2026-09-26 together
with its leads), followed by the 2026-09-26 review of `"python"` real-time
evolution (regressions in `tests/test_pyitensor_time_evolution_review.py`,
measurements in `docs/numerical_changes.md`, leads ranked in
`docs/time_evolution_improvement_plan.md`), and by
`docs/audit_2026_10_hole_hunt.md` (7, over the TASK correlator and the
Kondo spectrum on fermion chains, all fixed on 2026-10-23), and by
`docs/audit_2026_10_23_hole_hunt.md` (four bugs and one performance item,
over the least-covered areas rather than a commit; the bugs fixed the same
day). Four things to read before
changing the kind of code they cover:

- **`docs/numerical_changes.md`** before comparing any number against
  a result from before 2026-09-26. It is the ledger of every fix that moved
  numbers rather than behaviour, by call and by size, and most of the fixes
  did: the TD/TDZ pi factor and the KPM calibration, every injected state
  (normalized now), every constructor keyword (applied now), every
  finite-temperature result, every `gs_energy_fluctuation`, every `CVM`
  spectrum, and more. The records' `**Status**` lines carry the same
  entries finding by finding.
- **`docs/development.md` section 4.10** before adding a dispatch. Almost
  every finding of the first hunt was a dispatch decision taken before the
  information that should inform it: a short circuit ahead of a cache key,
  a precondition ahead of the branch it qualifies, an `else` serving both
  "unsupported here" and "you typo'd it", a `**kwargs` with no consumer.
- **`src/dmrgpy/dynamics.py`'s module docstring** before adding a
  correlator submode. It is the single normative statement of the
  convention: every submode returns the complex Lehmann density
  `i(G^R-G^A)/(2pi)`, unchanged wherever the Lehmann weights
  `M_n = <GS|A|n><n|B|GS>` are real (`Im M_n == 0` is the discriminant,
  which a Hermitian pair implies but does not exhaust).
- **The `hole-hunt` skill** (`.claude/skills/hole-hunt/SKILL.md`) before
  running or extending an audit; the record shape and the evidence standard
  are fixed there.

The durable rules the hunts left, each with the test or file that pins it.
Do not undo one without reading its finding first.

- **A state becomes the chain's normalized.** `groundstate.unit_copy` is
  applied in `mark_injected`, in the `wf0=` copies and in `set_gs`'s ED
  branch, and a zero state raises; an explicit `wf=` to
  `evolve_and_measure`/`evolution_ABA` is still evolved as given
  (`test_audit_2026_09_25b_session.py`).
- **Every correlator starts from `ground_state_on_session`**, with
  `groundstate.mark_injected` marking a state the caller supplied, so that
  a DMRG correlator after `set_gs()`, `set_initial_wf` or
  `set_initial_wf_guess` measures that state and not the session's own
  solved one (on `julia_live` too, since 2026-09-25); the mark's third field,
  `supplied`, read through `groundstate.state_supplied()`, is what SECTOR
  uses to refuse a set state that is not its reference sector's ground
  state (`test_audit_2026_09_24b_groundstate.py`,
  `test_audit_2026_09_24c_session.py`).
- **`groundstate.stored_answer_holds` is the one condition** under which
  `gs_energy`/`get_gs` return a stored answer; a keyword the stored answer
  would ignore (`maxde=`, `reconverge=True`, `wf0=`) raises or re-solves.
  After `gs_energy_generalized`, `e0` is the state's own energy `<wg|H|wg>`
  (biorthogonal on the NH route) and lambda is kept as `lam_generalized`
  (`test_audit_2026_09_25b_session.py`).
- **The session's lower band edge, `minimum_energy()`, never sweeps a state
  the session holds without its energy**: it solves from a fresh start and
  puts the state back (`test_audit_2026_09_24c_session.py`).
- **`vev.energy_variance` is how an energy variance is computed**,
  `||(H-<H>)|psi>||^2` with the `maxm` cap lifted. `<H^2>-<H>^2` through
  `vev(npow=2)` is a small difference of two numbers of order E0^2 whose
  error is set by the truncation of H|psi>, and H^2 as an MPO is O(L^3) on
  `"python"` (`test_audit_2026_09_24c_fluctuation.py`).
- **Every threshold is relative to a scale of the problem it acts on.**
  Where a rule exists it is: multiply by the power of two that brings a
  largest coefficient below 1 into [1,2), and leave everything
  byte-identical at 1 or more. It covers v2/v3's `to_mpo_unit`/
  `solver_hamiltonian` (`hscale_up_`, read from the merged AutoMPO),
  `pyitensor`'s `to_mpo`, the NH-DMRG entry in `nhdmrg.py` (once, for all
  three session backends, so do NOT also scale inside `Chain::nhdmrg`) and
  ED's deflated eigensolver. The others: `multioperator.clean_threshold` is
  1e-12 relative to the largest coefficient (in `canonical_dict`, to the
  larger of that and the summed signature magnitudes); CVM stops relative
  to `||b||` and its early exits read the CG functional; NH-DMRG's Ritz tie
  window is `1e-6*(max Re - min Re) + 100*eps*max|ev|` with `<=`,
  identically in `pyitensor/nhdmrg.py` and both `arnoldi_select_kbest`; the
  VUMPS Lanczos reference is capped by the Hamiltonian's unit; the KPM
  same-vector shortcut and ROOTN's Lanczos breakdown compare against the
  vectors' own norms (`test_audit_2026_09_25_scale.py`,
  `test_audit_2026_09_25b_{cvm,nh,cpp,ednum,mpobuilder}.py`). Not covered,
  and held by strict xfails: a bond whose strongest term is far below the
  largest coefficient still loses channels on every backend.
- **A chain constructor keyword is a setting assigned after `initialize()`**,
  exactly as an assignment afterwards, `mode` included (`sites.check_settings`,
  the explicit allowlist `sites.SETTINGS`); anything else, a Hamiltonian
  accumulator such as `hubbard=` included, raises
  (`test_audit_2026_09_25_construction.py`,
  `test_audit_2026_09_25b_construction.py`).
- **`mode.resolve_mode` lets only a chain pinned to `"ED"` override a call's
  `mode=`**; a `"DMRG"` pin is the same as None, and both names are
  validated before any fallback can return (`test_audit_2026_09_25b_session.py`).
- **KPM broadening is calibrated.** Every KPM route goes through
  `algebra/kpm.py::polynomials_for_broadening`, which solves the Jackson
  kernel's resolution relation so that `delta` is the broadening it is
  everywhere else (FWHM `2*delta` at the band centre); the DMRG routes
  reconstruct from exactly those n moments, the ED route uses the DMRG
  rescaling so the two share the same `x` at the same physical energy (on
  the bandwidth-centred window only: under `kpm_energy_truncate`,
  `mode="ED"` raises), `kpm_n_scale` defaults to 1 and must be a positive
  integer, and
  every backend raises once a moment exceeds 1.5 times the exact bound
  `||vi|| ||vj||`. Two residuals are documented rather than removed: the
  width tightens as `sqrt(1-x^2)` away from the band centre (0.70 of the
  requested width for an intensive excitation of a long chain), and the
  Jackson line peaks about 1.5x higher than a Lorentzian of equal FWHM.
  `get_distribution()`'s own KPM path is deliberately not on the
  calibration (`test_audit_2026_09_24_kpm.py`, `test_audit_2026_09_24b_kpm.py`).
- **`submode="TD"`/`"TDZ"` return the two-sided transform** written in
  terms of forward runs only, `(F[(A,B)] + conj(F[(B^dagger,A^dagger)]))/2`
  (`timedependent.lehmann_density_from_one_sided`), collapsing to `Re F`
  exactly when `canonical.is_dagger_pair` proves `A = B^dagger`; it refuses
  rather than guesses for a name with no known adjoint. Taking `np.real`
  instead measured worse than leaving the transform whole. The Fourier sum
  is evaluated at each requested frequency, not interpolated off an FFT
  grid; `sxt_to_skomega`, the infinite-chain `S(k,omega)` reduction, and
  the two short-window direct callers stay on the FFT stage, and
  `sxt_to_skomega` on the raw transform, since it never sees the pair
  (`test_audit_2026_09_24_realtime.py`).
- **Real-time evolution on `"python"`** (`tdvp._lanczos_expm_multiply`)
  stops on `|coeff|*beta_k*|c_k| < 1e-10` relative to the vector, with
  sub-stepping when `niter` runs out; the host loop and the device path
  share one test, `_krylov_verdict`, and return their basis rather than a
  vector, so the sub-stepping above them cannot drift between the two. Keep
  it that way if either changes. Every `evolve_and_measure_*` loop on every
  backend restores the input norm `norm0` after each step, and
  `gse.global_subspace_expand` normalizes each Krylov vector
  (`test_pyitensor_time_evolution_review.py`). `TDVP_GSE` stops growing the
  bond dimension after `tdvp_gse_sweeps` steps, so a growing-entanglement
  quench saturates; recorded, not changed.
- **`mode="ED"` real-time routes are `expm_multiply`**, start from the
  cached ground state when called without `wf=`, and
  `evolve_and_measure`/`evolution_ABA` raise on an unknown keyword on DMRG,
  take `h=` on ED and default to `nt=1000` on both
  (`test_audit_2026_09_24c_realtime.py`).
- **`julia_live` starts every DMRG solve from a random MPS** of link
  dimension min(maxm,10), with the chain's noise on the first half of the
  sweeps, so its solves vary run to run by solver noise and a warm 16-site
  solve is about 40% slower than before; the KPM window there trusts `e0`
  as the lower band edge only when `groundstate.mark_lower_edge` says a
  plain solve produced it (`test_audit_2026_09_25b_julia.py`).
- **`Thermal_Spin_Chain.anneal()` is order-8 Taylor steps** of
  `exp(-dtau(H-E_mid))` at `dtau*W <= 2`, about `2*beta*W` factor
  applications (`test_audit_2026_09_25b_thermal.py`).
- **TASK's error bar is two numbers, not one.** The two-site variance of
  a pole does not see the part of `B|GS>` outside the one-site tangent
  space, which is dropped before the recursion starts: zero for a sum of
  one-site operators and for a fermion operator on a state of definite
  parity, up to all of the weight for `C_i` on a parity mixture (a
  degenerate Majorana pair) at truncated `maxm`. `taskdmrg.get_poles`
  returns it as `missing_weight` and warns above 1e-6 of the weight
  (`test_audit_2026_10_task.py`).
- **A collapse basis METTS cannot sample with is refused.** `metts_vev`
  raises on a non-Hermitian basis operator, and on an all-diagonal basis
  (`("N",)`, `("Sz",)`) when H conserves `sum_i O_i` for one of its names,
  tested on a random witness, since every sample would then stay in the
  first sample's charge sector with an error bar of 0
  (`test_audit_2026_10_23.py`).
- **The DMRG `vev` has no `**kwargs`**: an unknown keyword raises, and
  `T>0` raises naming `mode="ED"` and `metts_vev()`
  (`test_audit_2026_10_23.py`).
- **Out of scope by construction** for every hunt: the vendored ITensor and
  TDVP folders, the three legacy quirks `.claude/rules/cpp-backend.md` lists
  as deliberately reproduced (`evoloperator()`'s z^3/6 term on `H2`, the
  `"moise"` key, the unreachable `"tevol_fit_td"` branch), the open
  `docs/known_issue_*.md` items, and what `ROADMAP.md` marks as absent by
  design.

Left open and recorded as such, in the records' "Left open" and "New leads"
sections: the `kpm_energy_truncate` window problem outside its worst regime
(`docs/known_issue_kpm_energy_truncation_window.md`); the sliver just below
`kpm_scale=1/2` where the elastic line is distorted while the moments are
still under the bound; the truncated CVM regime (warned, about 2x the cost,
still unconverged); `Thermal_Spin_Chain`'s cost and its variational band
width; the bond-local truncation above; the MPO-set Hamiltonian and
small-unit real-time evolution; v3's Krylov exponentiator's absolute error
goal (`applyExp`) and v3's run-to-run `emax` noise, which is the whole of
the ED-against-v3 KPM residual; `kpm_finite`'s `window_chain_kwargs` check;
`edtk/edchain.py`'s silent `None` on an unknown operator; the reason for
`julia_live`'s slowdown and its KPM guard raising `juliacall.JuliaError`
rather than `RuntimeError`; `ShiftInv`'s `applyinverse(delta=)`; the lower-level route's pair `name=` next to
`i=`/`j=`; and the short `nt` defaults of `td_dynamical_correlator`.

## Architecture

### Entry point classes

Each physical model is a thin subclass of `Many_Body_Chain`
(`manybodychain.py`), which holds all shared state (bond dimensions, sweep
count, cutoffs, mode, etc.) and implements the generic operations
(`set_hamiltonian`, `gs_energy`, `vev`, `get_gs`, correlators, time
evolution, entanglement, ...). Model modules just define the local Hilbert
space and convenience operators:

- `spinchain.py` -> `Spin_Chain` (aliased directly to `Many_Body_Chain`)
- `fermionchain.py` -> `Fermionic_Chain` (spinless fermions, `C`/`Cdag`/`N`/Fermi-string `F`)
- `spinfermionchain.py`, `bosonchain.py`, `parafermionchain.py`: same pattern for other statistics/models

### Operator representation: `MultiOperator`

Hamiltonians and observables are built as `MultiOperator` objects
(`multioperator.py`), sums of products of named single-site operators (e.g.
`Sx[i]*Sx[j]`). This is a backend-agnostic intermediate representation: the
same `MultiOperator` is later either converted to ITensor terms for the
C++/Julia/`"python"` DMRG backends (`MultiOperator.to_terms()`), or
converted to a sparse matrix via `multioperator.MO2matrix` for the ED
backend (`edtk/edchain.py`). `multioperatortk/` holds supporting machinery
(Jordan-Wigner strings for fermionic operators, static/long-range operator
construction, the canonical form below).

**Canonical form and the Hermiticity proof.**
`multioperatortk/canonical.py` rewrites every term into one spelling
(factors sorted by site, stable so same-site factors keep their order, the
coefficient carrying the sign of every fermionic exchange the sort makes,
identity factors dropped) so that terms which are the same operator have
the same signature and collect in a dict.
`MultiOperator.simplify()`/`is_zero()`/`is_hermitian()`/`is_antihermitian()`
all read off it. It replaced a sympy round trip (`sympymultioperator.py`,
deleted with it, taking sympy out of the dependencies) that collected only
terms already spelled identically: 0.834 s on an n=20 spin-1/2 chain to
report `False` for an ordinary Heisenberg Hamiltonian, now proven Hermitian
in 0.5 ms. The rewrite is exact, so the proof is **one-sided**: an empty
canonical form of `H-H^dagger` means Hermitian, full stop, while a nonempty
one means not proven, since a same-site identity (`Sx Sx = 1/4`) or a name
alias (`Sp = Sx + i Sy`, `Sz = (Nup-Ndn)/2`) needs the site type, which
lives on the chain and not on the operator. Two factors sharing a site are
reordered only when both are diagonal (`_DIAGONAL`), which is what makes the
Hubbard `Nup[i]*Ndn[i]` provable while leaving `Sx[i]*Sz[i]` alone; do NOT
widen that to a general same-site sort, those two are both even and do not
commute. Nor is a term reordered when it names both a pre-transform
`C`-type name and a bare `A`-type one (`_LOCAL_LADDER`): with
`C_j = F_0...F_{j-1} A_j`, `A_i` anticommutes with `C_j` for j>i and
commutes for j<i, so their exchange sign depends on the site order and is
not a parity (2026-09-24 audit finding 1). Do NOT fix that by giving the
`A`-type names no parity: they are also the boson ladder operators, and the
boson and Jordan-Wigner-transformed proofs need them even. So
`Many_Body_Chain.is_hermitian()` takes the proof when it lands and falls
back to its random-witness probe otherwise; do not turn that fallback into
a rejection. Only the names in `canonical.py`'s `_PARITY` table are
reordered; anything else (parafermionic `Sig`/`Tau`, which reorder with a
Z_n phase, or a caller's own name) is left spelled as written, and the
Hermiticity proof refuses outright when any name is off the table, since
`get_dagger()` leaves an unrecognized name untouched (and carries a phase
for a name whose adjoint is a multiple of itself, `_dagger_phase`, so far
only `ISy` -> `-ISy`) and such an operator would otherwise cancel against
its own "dagger". See `docs/development.md` 4.2a,
`tests/test_multioperator_canonical.py` and
`examples/algebra/operator_canonical_form`.

### Backend dispatch (DMRG vs ED, C++ vs Julia)

`mode.py` decides, per call, whether a calculation runs via DMRG or ED
(`get_mode`/`resolve_mode`): DMRG unless `self.mode` pins the chain to ED,
or unless `self.itensor_version` is `2` or `3` and the corresponding
pybind11 extension is not compiled (see `cppext.available(version)`), in
which case it falls back to ED with a printed note. Two more automatic
fallbacks exist, both structural: ITensor v3's two-site `dmrg()` aborts
below 3 sites (see "A real, confirmed `mpscpp3` bug" above), and
`"python"`'s two-site DMRG has no update to make on 1 site, where `dmrg()`
returned `None` before the fallback existed, so `itensor_version==3` with
`ns<3` and `"python"` with `ns<2` route to ED. `itensor_version="python"`
never falls back for a missing extension (it has none).
`groundstate.py::gs_energy_generalized` needs its own copy of the short-chain
guard (it has no ED fallback to be routed to), and there it raises, because
the symptom was worse than `None`: the outer self-consistent iteration
returned the Rayleigh quotient of a state no sweep ever touched, a silently
wrong number (-0.3049 for an exact -0.5); that guard sits before the
non-Hermitian dispatch, since NH-DMRG escapes v3's abort by never calling
`dmrg()` but its own sweep is two-site too. A conserved sector
(`set_conserved_sector`) makes `get_mode` raise on DMRG with
`itensor_version=2` or `"julia_live"`, the two backends with no quantum
numbers, rather than answer with the global ground state; falling back to
ED is correct there, since ED restricts to the sector itself.

**Which backend a chain gets when the caller names none** is
`cppext.default_backend()`, not `DEFAULT_ITENSOR_VERSION` directly: the
default C++ version when that extension is compiled, and `"python"`
otherwise. That fallback is what makes `pip install dmrgpy` usable; before
it every default-backend chain in a pip install silently ran ED, which
cannot reach the sizes an MPS backend exists for. Only the implicit choice
moves: an explicit `itensor_version=3` with no extension still falls back
to ED, and `mode="ED"` is untouched.
`Many_Body_Chain.__init__`/`Mixed_Spin_Fermion_Chain.__init__` take
`itensor_version=None` to mean "pick one for me". See
`tests/test_default_backend_without_cpp.py`. Most public methods on
`Many_Body_Chain` accept a `mode="DMRG"|"ED"` kwarg so results can be
cross-validated between solvers.

- **DMRG, C++ (`itensor_version=2` or `3`, the default)**: entirely
  in-process through a pybind11 extension over a vendored ITensor, one
  `Chain` session per `Many_Body_Chain` held as `self._session`. Both
  versions expose the identical `Chain` session API and public
  `Many_Body_Chain` surface; `itensor_version` only picks which compiled
  extension (`mpscpp2._dmrgcpp` or `mpscpp3._dmrgcpp`, via
  `cppext.get_backend(version)`) backs the session.
  `cppext.DEFAULT_ITENSOR_VERSION` (currently `3`) is the single source of
  truth for the default, threaded through `Many_Body_Chain.__init__`,
  `setup_cpp()`, `get_backend()`/`available()`, and `install.py`'s own
  `--itensor-version` default (`installtk/__init__.py`'s matching
  constant); change it there, not per call site.
  `Many_Body_Chain.setup_cpp(version=2)` switches an existing chain between
  them. Everything else about this backend, the `cpp_handle` pattern, the
  unit-scale rule, the three deliberately reproduced legacy quirks, the
  `mpscpp3`-specific behaviour (QNs off by default and a random-MPS start,
  the opt-in conserved-sector mode, `apply_mpo()`'s `.noPrime("Site")`,
  `innerC`, TDVP as the default real-time method, `py::module_local()`), is
  in `.claude/rules/cpp-backend.md`, formerly this file's "In-process
  pybind11 extension" section.
- **DMRG, pure Python (`itensor_version="python"`)**: `pyitensor/`
  (`Many_Body_Chain.setup_python()` switches an existing chain to it) is a
  from-scratch, pure-Python/NumPy/SciPy reimplementation of exactly the
  ITensor v3 API subset `mpscpp3/chain_session.h` uses: Index/ITensor tensor
  core, the same eleven site types, an AutoMPO term compiler with its own
  from-scratch Jordan-Wigner threading, MPS/MPO algebra, two-site
  Lanczos-based DMRG (ground and excited states via the overlap-penalty
  method), two-site Krylov-based TDVP, and a `pyitensor/chain.py::Chain`
  exposing the identical method surface as the compiled backends'
  `self._session`, so every `itensor_version in (2,3)` call site elsewhere
  treats `"python"` as a third, always-available option wherever such a
  site was updated to include it. Exists so DMRG/TDVP work with zero
  compiler/pybind11 dependency, at the cost of being substantially slower
  than compiled ITensor (no block-sparsity, no JIT); see
  `pyitensor/__init__.py`'s docstring and the module docstrings throughout
  `pyitensor/` for the specific simplifications taken and why each one does
  not affect dmrgpy's own results. MPO construction is no longer one of
  them: `pyitensor/mpobuilder.py`'s `to_mpo` assembles the finite-state
  machine over the terms' partial products directly (the shape ITensor's
  own `toMPO(...,{"Exact",true})` produces), O(L^2) instead of the old
  O(L^4) sum-and-compress, keeping two sweeps afterwards purely to honour
  `cutoff`/`maxdim`, the first exact and only the return one truncating (a
  truncating first sweep compressed any operator whose coefficients were
  all below about 2e-7 into a different one), and an operator whose largest
  coefficient is below 1 is built at unit scale and scaled back, like
  v2/v3's `to_mpo_unit()`. Two rules in there are load-bearing and easy to
  get subtly wrong: the coefficient goes on the transition into the final
  state (so terms differing only by a coefficient still share partial
  states), and transitions into that state accumulate while structural
  ones are assigned (two syntactically identical terms trace the same path,
  so their coefficients must sum while their shared structure must not be
  written twice). The old construction is kept as `_sum_of_term_mpos`, the
  reference the machine is tested against in
  `tests/test_mpo_automaton_builder.py`.
- **DMRG, Julia (`itensor_version="julia_live"`)**: a live, in-process
  Julia session (`mpsjulialive/`, via `juliacall`/`juliasession.py`) with
  its own parallel set of modules (`mpsjulialive/groundstate.py`, `vev.py`,
  `mps.py`, `mpo.py`, ...) mirroring the top-level ones but talking to
  Julia's ITensors.jl instead of the C++ pybind11 extension.
  `Many_Body_Chain.setup_julia()` switches a chain to this mode.
  `itensor_version="julia"` (a separate, older subprocess-based Julia path
  via `juliarun.py`) is not reachable through the normal public API today
  and should be treated as legacy/inert.
- **ED path**: `edtk/edchain.py` (`EDchain`) builds dense/sparse operators
  directly in Python/NumPy/SciPy (`pyfermion/`, `pyspin/`, `pyboson/`,
  `pyzn/` provide the per-statistics many-body operator construction,
  `edtk/one2many.py` promotes single-site operators to the full Hilbert
  space) and diagonalizes with `scipy.sparse.linalg`.

**Conserved sectors.** `Many_Body_Chain.set_conserved_sector(Nf=6)`/
`(Sz=0)`/`(Nf=8,Sz=0)` (off by default, no-args clears it) confines the
whole calculation to one quantum-number sector, by three mechanisms behind
one API: QN-carrying ITensor indices on `itensor_version=3`, dense storage
with a charge grading plus a charge penalty on `"python"`
(`pyitensor/sector.py`), and a restriction of every assembled operator to
the basis states carrying the charge on `mode="ED"` (`edtk/edchain.py`,
where `promote_to_dense`/`promote_mps` are deliberately not implemented
since an ED re-solve is cheap). The sector lives on the chain, is re-applied
wherever a session is rebuilt (`sites.py::initialize`, `__deepcopy__`) and
is part of `groundstate.py`'s send-cache key. It pays off with size, and
only with size (0.6x at n=20, 4.2x at n=60 on Heisenberg chains); do not
quote a single speedup number. `tests/test_sector_conservation.py` and
`tests/test_sector_conservation_ed.py` cross-check every sector against ED
restricted to that sector; the reference must restrict, not filter
full-spectrum eigenvectors by `<N>`, since sector-degenerate eigenvalues
come back as arbitrary superpositions. See `docs/development.md` 4.3a and
the C++ rule file for the three ITensor behaviours that shape the v3
implementation.

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
on the C++ DMRG side, `algebra/kpm.py` on the ED side) rather than direct
spectral decomposition, since exact diagonalization of the full spectrum is
infeasible for large chains. The convention and the calibration are in the
rules above.

`submode="SECTOR"` (`sectordc.py`) is the exception to that shape, and the
only submode that does not compute a spectrum on the chain it is called on:
it solves the reference sector and the one sector `B|gs>` lands in (N+1 for
`Cdag`, Sz+2 for `S+`) on an internal clone, promotes both to dense indices
and contracts them, summing an explicit Lehmann representation over
converged per-sector eigenstates. It rests entirely on `promote_mps`
rebasing onto the chain's original `dense_sites_`, which
`set_conserved_sector` never rebuilds
(`tests/test_sector_promotion.py::test_states_promoted_from_different_sectors_are_comparable`).
`get_spectral_function`/`get_spin_spectral_function` assemble the
particle/hole and Sz channels on top, and a charge-indefinite operator
(`Sx`) is split into `(S+ + S-)/2` and summed channel by channel
(`multioperatortk/charge.py`) rather than rejected. Two things to know
before touching it: the gating is explicit inside the submode (`mode.py`'s
sector guard keys on the caller's `conserved_sector`, which is empty here,
so it never fires), and the per-sector solves are cached with the clone
that produced them (`_sector_states_cache`) because ITensor Index identity
is per session; that cache is also what makes a site sweep, and so
`S(q,w)`/`A(k,w)`, cost about one site's correlator. See
`docs/sector_resolved_spectral_function_plan.md` for the design and
`tests/test_sector_spectral_function.py` for what is pinned.

### STM/Kondo tunneling spectra (`kondospectrumtk/`)

`get_kondo_spectrum` (on `Spin_Chain` and both spinful fermion chains,
`kondospectrumtk/kondomethod.py`) is the third-order perturbation theory of
Ternes (NJP 17 063016, arXiv:1505.04430); `mode="ED"` is the reference
implementation at any T and `mode="DMRG"` reproduces it at T=0 via the
dynamical correlator and a two-time TDVP construction. Every prefactor was
pinned against the paper's Fig. 7, except that the Kondo term's exchange
diagram deliberately departs from the paper's eq. 25, where a brute-force
T-matrix puts the log elsewhere; do NOT "restore" eq. 25 to get back on the
figure. The notes, formerly this file's section of the same name, are in
`.claude/rules/kondo.md`.

### Julia vs C++ backend

The two DMRG backends are independent implementations, not a shared
protocol: the C++ path runs in-process through the pybind11 extension,
while `itensor_version="julia_live"` (`mpsjulialive/`) drives a live
Julia/ITensors.jl session via `juliacall` with its own mirrored set of
modules. A feature missing on one side simply is not implemented there yet
(check the relevant `mpsjulialive/*.py` file, or lack thereof) rather than
falling back automatically; the only automatic fallbacks are the DMRG to ED
ones `mode.py` makes (missing extension, too few sites), none of which
substitutes one DMRG backend for another.

### Infinite chains (iDMRG and VUMPS)

`infinitechain.py` is the public API, `pyitensor/idmrg.py`, `vumps.py` and
`vumps_ms.py` the `"python"` implementation, and the `idmrg_*`/`ic_*`/
`vx_*`/`vms_*` helpers in `mpscpp3/chain_session.h` the `itensor_version=3`
one; the upstream reference is ITensor's own
<https://github.com/ITensor/ITensorInfiniteMPS.jl>, not vendored here. The
rules, with the reasoning and the measurements in
`.claude/rules/infinite-chain.md` (formerly this file's notes under the
"ITensor library" header):

- Both environment builders on both backends try the bond candidate first
  (`C C^dag` and `conj(C^dag C)` in this codebase's `X[ket,bra]` ordering,
  `vumps._transfer_fixed_points`, `vumps_ms._cell_fixed_points`,
  `Chain::vx_choose_fixed_point`), a yes/no residual test, and run the
  guarded eigensolver only when it fails; a threshold on `C`'s weight
  spectrum is the WRONG shape, since the ratio to catch is a moving number
  mid-convergence. The LEFT candidate is the transpose, `conj(C^dag C)`,
  and only a measured residual ever caught that.
- VUMPS's Lanczos stops on the Ritz residual (`residual_tol`, passed as
  `tol/10` by `vumps.py`/`vumps_ms.py` and by `Chain::vumps_ground_state`),
  since an eigenvalue criterion caps the eigenvector at ~1e-6 and the
  VUMPS criterion is a difference of two independently solved eigenvectors;
  `AC`/`C` are phase-aligned against their warm starts. Do not "fix" finite
  DMRG's criterion on this reasoning: measured, its sweep contracts the
  per-solve error and the tight criterion moves nothing above 5e-15.
- Any Hamiltonian with reach>1 goes to the sequential solver (`vumps_ms`,
  `Chain::vumps_ground_state`'s `use_multisite`), and the snapshot that
  exists decides the reader; `gs_method="idmrg"`'s `_active_channels_at`
  carries one pending channel per site of a term's reach by design, so do
  NOT add a guard there.
- `idmrg.py`'s transfer chain is `_TransferChain`, deliberately not
  list-like, so a leftover `Es[p]` raises instead of materializing a chi^4
  array; v3's static observables come from the extracted unit cell
  (`idmrg_cell_`), never from tiling the per-micro-step `idmrg_U_`.
- Any `itensor_version=3` VUMPS energy from before 2026-09-12 on a model
  whose exact state is smaller than the requested `maxm` is not comparable
  (`docs/known_issue_v3_vumps_variational_floor.md`, FIXED;
  `tests/test_vumps_redundant_bond_dimension.py`, with a repeated-run guard).

## Documentation

`docs/development.md` is the development document: the architecture
(backends, dispatch, directory layout); the user guide covers the physics-facing API (what
calculation each method performs, with the formula behind it), kept brief
and physics-first, with no history in it; `docs/numerical_changes.md` is
the ledger of what raises and what changed numbers.

**The user guide is a website.** Its only source is `docs/guide/*.md`, one
page per topic, built with MkDocs Material from `mkdocs.yml` (repo root,
`docs_dir: docs/guide`, so nothing else under `docs/` is published) and
deployed to <https://joselado.github.io/dmrgpy/> by
`.github/workflows/guide.yml` on every push to `master` that touches it.
The PDF is generated from the same pages (`docs/build_guide_pdf.py`, pandoc
and xelatex, also run by the workflow); there is no handwritten `.tex` any
more, and `docs/user_guide.md` is only a pointer. Math is `$...$`/`$$...$$`
(MathJax through `pymdownx.arithmatex`); a `|` inside math in a table row
is written `\vert`/`\lvert`/`\rvert`, or the table splits. Pages link to
each other by page name (`[Ground state](ground-state.md)`), never by a
section number. A new page needs an entry in `nav`, or it is neither in
the menu nor in the PDF. `mkdocs build` is strict and fails on a broken
link or anchor; run it (`pip install mkdocs==1.6.1 mkdocs-material==9.7.7`)
before committing a guide change. Every figure in the guide is computed:
`docs/figures/<id>/compute.py` writes `data.npz` (committed), `plot.py`
reads only that and writes the light and dark PNGs into `docs/guide/img/`;
`docs/figures/README.md` says how to recompute, replot and run on a remote
host, and `docs/figures/plan.md` is the physics each figure checks.
Under each figure a collapsed block shows its `compute.py` and `plot.py`,
copied in at build time by `pymdownx.snippets` with `base_path:
docs/figures`, a path relative to the current directory, so `mkdocs build`
runs from the repo root.

Whenever a new feature is implemented (a new physical model, a new
`Many_Body_Chain` method, a new dynamical-correlator submode, a new
post-processing tool) update the matching page under `docs/guide/`, and
`docs/development.md` too if the change affects architecture or
backend dispatch rather than just physics; the
`docs-sync` skill (`.claude/skills/docs-sync/SKILL.md`) is the procedure.

**The user guide's style is fixed.** Its readers are condensed matter
physicists who want to compute things, not to learn how the code works. Each
capability gets a short physics explanation (the quantity, its formula, what
it is used for, what measures it), then the code, then the few keywords that
control the physics or the accuracy as a one-line-per-item list. Short plain
sentences, one idea each. No audit history, no "used to", no old numbers, no
timings, no implementation internals (caches, penalties, block sparsity,
stopping rules), and a backend caveat only when a number would otherwise be
wrong, in one sentence. Every snippet in it has been run. The
"Ground state" page (`docs/guide/ground-state.md`, formerly section 3,
"Ground-state properties") is the reference: a new or edited page matches
it in length per topic and in tone. History and changed numbers go to
`docs/numerical_changes.md`, never into the guide.

## Git workflow

Working directly on `master` is fine at all times, unless the user says
otherwise for a given task or the change is a really major one with real
potential to conflict with other concurrent development (in which case use
a feature branch/PR instead).

When a feature branch is warranted, new work (a new branch, a new PR) should
start from `master`, not from another feature/staging branch, unless
explicitly told otherwise. This repo previously reused a long-lived branch
(`mpo-algebra-tdvp-gse`) as a staging area for several stacked PRs; once
`master` and the staging branch diverged (each gained commits the other
lacked), a PR built on the staging branch picked up unrelated content and
produced conflicting/misleading diffs against `master`. Branching from
`master` (`git fetch origin master && git checkout -b <branch>
origin/master`) avoids this.
