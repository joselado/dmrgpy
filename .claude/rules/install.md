---
paths:
  - "install.py"
  - "install_julia.py"
  - "installtk/**"
  - "addsystem.py"
---

# Install and build (`install.py`, `installtk/`)

These are the notes CLAUDE.md used to carry under "Install / build"; they
load when `install.py` or a file under `installtk/` is read or edited. The
commands themselves are still in CLAUDE.md.

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
across reruns). There is no lint config, and the only CI is the build of
the user guide.

Both `mpscppN/Makefile`s list `$(wildcard *.h)` among the extension's
prerequisites: without that, editing `chain_session.h`, where the actual
DMRG implementation lives, leaves `make pybind` reporting "Nothing to be
done" while the `.so` keeps the old code. A C++ change is not tested until
the `.so` has been rebuilt.
