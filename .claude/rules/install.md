---
paths:
  - "install.py"
  - "install_julia.py"
  - "installtk/**"
  - "src/dmrgpy/buildinfo.py"
  - "src/dmrgpy/doctor.py"
  - "src/dmrgpy/cppext.py"
  - "tests/test_installer.py"
  - "tests/test_backend_status.py"
---

# Install and build (`install.py`, `installtk/`)

These are the notes CLAUDE.md used to carry under "Install / build"; they
load when `install.py`, a file under `installtk/`, `cppext.py`,
`buildinfo.py` or `doctor.py` is read or edited. The commands themselves
are still in CLAUDE.md. The first part is the requirement check, the
second what the 2026-11-22 rework added around it.

## The requirement check

`install.py` checks before it builds: `installtk/requirements.py` first checks
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
first one that actually links. Finally, `install.py` links the repo's
`src/dmrgpy` into the interpreter's site-packages
(`installtk/addpythonpath.py`) and, only with `--shellrc`, adds `DMRGROOT`
to the user's shell rc file (`installtk/addsystem.py`); both are idempotent
across reruns. There is no lint config, and the only CI is the build of
the user guide.

Both `mpscppN/Makefile`s list the headers of the backend among the
extension's prerequisites: without that, editing `chain_session.h`, where
the actual DMRG implementation lives, leaves `make pybind` reporting
"Nothing to be done" while the `.so` keeps the old code. A C++ change is
not tested until the `.so` has been rebuilt.

## What the 2026-11-22 rework added

The check above was kept as it was. What changed is what happens after
it, and what is left behind for a later session. The measurements are
from this workstation (GCC 13.3, system LAPACK/BLAS, conda Python 3.14).

- **The flags are in one module**, `installtk/buildflags.py`, read by the
  check and by the build. On x86-64 with GCC the compiler command is the
  one every earlier build used, `g++ -m64 -std=c++17 -fconcepts -fPIC` for
  v3 and `-std=c++14` for v2, and `tests/test_installer.py` pins the
  string. Do NOT add `-march=native` or change `-O2`: the repository has
  bit-identical tests and a ledger, and a build made on a login node runs
  on other nodes. The order of the compiler and BLAS candidates is pinned
  by the same test file for the same reason.
- **The check compiles with the real flags** of every backend asked for,
  on a small program that needs its language standard
  (`buildflags.TRIAL_SOURCE`). Before, it compiled with `-std=c++14`
  whatever the backend, so a compiler that cannot build v3 passed.
- **Clang and non-x86 are written, not tested.** Clang gets ITensor's own
  flags from `options.mk.sample` (`-std=c++17 -Wno-gcc-compat`, no
  `-fconcepts`), is not asked for a `cc1plus` (only GCC has one, so the
  old check rejected every Clang), other architectures get no `-m64`, and
  `cppversion.parse_version` reads the version after the word "version"
  on a Clang line (the old parser read the build number in parentheses
  and called an Apple Clang too old). No Clang and no ARM machine was at
  hand, so say so if a user reports a failure there.
- **ITensor is rebuilt only when it has to be.** `install2.itensor_key`
  hashes the git tree of the vendored folder with the compiler command and
  the BLAS configuration, and the key of the last build is in
  `build_info.json`. Without git, or with local changes in the folder, the
  key is None and ITensor is always rebuilt; `--force` rebuilds too. The
  extension is cleaned only when its configuration changed, otherwise make
  decides. A second `python install.py` takes 4 s instead of 183 s (v3; a first v2
  build took 167 s).
- **`--jobs`** runs the ITensor build in parallel, 4 at most by default
  (`buildflags.default_jobs`): 76 s with four jobs, and each compiler
  process stays below 1.7 GB. The extension is one translation unit,
  100 to 160 s for v3 and 39 s for v2, and no flag makes it faster. Keep the default small: the build
  often runs on a shared login node.
- **One compiler choice.** `install.py` writes `mpscppN/build_config.mk`
  (`DMRGPY_PYTHON`, `DMRGPY_PYBIND_CCCOM`) and both Makefiles include it,
  so `make pybind` by hand uses the compiler of the install. The
  `CONDA_GXX` detection of the Makefiles is still there for a folder with
  no such file; do NOT remove it, it is the fix of the segfault described
  in `mpscpp2/Makefile`. `--extension-only` is `make pybind` with that
  file, for the step after a C++ edit. It also runs in a checkout built
  before 2026-11-22, which has no such file, with a note; it needs only
  the ITensor library, `options.mk` and `this_dir.mk`. The key of the
  ITensor library and the test for cleaning the extension both hold the
  compiler version, since a system compiler upgraded in place keeps its
  path.
- **The build record**, `mpscppN/build_info.json`, is written by
  `src/dmrgpy/buildinfo.py`, which has no import from dmrgpy because the
  Makefile runs it as a script (`--stamp`) after every link and
  `install2.py` loads it from its file (importing the package there would
  import the extension being built). Its `sources_hash` covers
  `bindings.cc`, `*.h`, `extra/*.h` and `TDVP/*.h`, the same files the
  `pybind` target lists as prerequisites. Staleness is decided by that
  hash, never by file times, which a git checkout changes.
- **`cppext.py` says why.** `_record_failure` sorts a failed import into
  `not_compiled`, `other_python` (an extension for another interpreter is
  in the folder) and `load_failed`, and warns once for the last two. It
  sets `_backends[version] = None` before it writes the warning: `explain`
  asks `get_backend`, and without that line the first real test recursed
  without end (the unit tests had preset `_backends`, so only the real
  extension showed it; `test_a_failed_first_import_does_not_import_again`
  pins it now). `dmrgpy/__init__.py` keeps its own warning for the
  not-compiled case only. Do NOT warn for `not_compiled` in `cppext`:
  that is every pip install.
- **The self-test** is `python -m dmrgpy.doctor --selftest`, run by
  `install2.selftest` in a fresh process with this checkout's `src` first
  on `PYTHONPATH`, so that it cannot test another copy of dmrgpy the
  interpreter knows. `install.py` returns 1 when it fails and leaves the
  build in place.
- **Testing the installer.** Never run a plain `python install.py` in the
  primary checkout to test a change: it rebuilds the working `.so` and
  links site-packages to that checkout. Use a
  worktree and `--build-only`. When a test moves the extension file by
  hand, quote the path: `ls` is aliased to mark executables with `*` here,
  and an unquoted `$(ls ...)` renamed the file.
- **The author's picks of 2026-11-22.** The symlink into site-packages
  stays; a `.pth` file and `pip install -e .` were offered and not taken
  (the second was never tested against the `pyproject.toml` that leaves
  the C++ folders out). The `DMRGROOT` line in the shell rc file is no
  longer written by default: `--shellrc` asks for it, the one example
  that read the variable (`topological/transverse_2d_ising_model`) adds
  its path like the others, and lines already in a user's rc file stay.
  Do NOT make the installer edit a shell rc file by default again.
  `install_julia.py` is still a separate script; folding it into
  `install.py` is the author's decision.
