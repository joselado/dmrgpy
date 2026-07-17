"""Preflight requirement checks for install.py.

Everything the compilation step (installtk/install2.py) needs is resolved
and *verified* here first -- compiler, BLAS/LAPACK, pybind11, make -- so a
broken environment is reported clearly before a several-minutes-long
ITensor build is attempted, instead of failing at the last link step.
"""
import os
import platform
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Optional

from . import cppversion
from . import blaslapack


@dataclass
class BuildConfig:
    gpp: str              # C++ compiler to use for ITensor AND the pybind extension
    platform: str         # ITensor PLATFORM= value (lapack/openblas)
    libflags: str         # BLAS_LAPACK_LIBFLAGS
    includeflags: str     # BLAS_LAPACK_INCLUDEFLAGS
    python_exe: str       # interpreter to build the extension against
    is_conda: bool


def _fail(msg):
    print("\n### Requirement check failed ###")
    print(msg)
    sys.exit(1)


def _check_os():
    system = platform.system()
    if system not in ("Linux", "Darwin"):
        _fail("Unsupported OS: "+system+". This project only builds on "
              "Linux and Mac (see README's Windows section for using a "
              "Linux virtual machine instead).")
    return system


def _check_make():
    if not shutil.which("make"):
        _fail("'make' was not found on PATH.\n"
              "Install it, e.g.:\n"
              "  Ubuntu/Debian: sudo apt-get install build-essential\n"
              "  Mac:           xcode-select --install")


def _compiler_candidates(gpp_arg, is_conda, python_exe):
    """Yield (label, path) compiler candidates to try, in priority order.
    An explicit --gpp is a hard override, not a suggestion to cascade past,
    so it's the only candidate yielded when given.

    Otherwise a conda-provided compiler (if any) is tried first -- to keep
    a single, consistent libstdc++/BLAS ABI with conda's own numpy/scipy,
    see find_conda_compiler()'s docstring -- but the system compiler on
    PATH is *also* offered as a fallback candidate, since a conda-adjacent
    compiler can itself be broken: e.g. Aalto Triton's `scicomp-python-env`
    module ships a conda-style Python distribution whose bundled
    `x86_64-conda-linux-gnu-g++` wrapper can't find its own `cc1plus`
    (confirmed directly), while plain system `g++` on that same login node
    works fine without any module."""
    if gpp_arg is not None:
        yield "explicit --gpp", gpp_arg
        return
    conda_gpp = cppversion.find_conda_compiler(python_exe) if is_conda \
            else None
    if conda_gpp is not None:
        yield "conda-provided compiler", conda_gpp
    for name in ("g++", "c++"):
        path = shutil.which(name)
        if path:
            yield "system "+name, path
            break


def _try_compiler(gpp, check_gpp):
    """Check a single compiler candidate without exiting on failure, so
    the caller can cascade to the next one. Returns (ok, reason)."""
    if not shutil.which(gpp) and not os.path.isfile(gpp):
        return False, "not found"
    if not check_gpp:
        return True, None
    version = cppversion.compiler_version(gpp)
    if version is None or not cppversion.correct_version(gpp):
        return False, ("could not be used or is too old (found version "
                +str(version)+", need >= 6)")
    if cppversion.backend_missing(gpp):
        return False, ("the compiler driver can't locate its own cc1plus "
                "backend -- a partial/broken toolchain, not a version or "
                "PATH problem (common for module-provided Python "
                "distributions that bundle a compiler wrapper without the "
                "matching backend package)")
    ok, out = cppversion.trial_compile(gpp)
    if not ok:
        return False, "failed to compile a trivial test program:\n"+out
    return True, None


def _environment_hints(is_conda, python_exe):
    lmod_present, loaded = cppversion.lmod_state()
    hints = []
    if lmod_present:
        hint = ("This looks like an HPC cluster using environment modules "
                "(Lmod). No usable C++ compiler may be on PATH until one "
                "is loaded, e.g. on Aalto's Triton:\n"
                "  module avail gcc\n"
                "  module load <stack> gcc/<version>   "
                "# e.g. triton/2024.1-gcc gcc/12.3.0\n"
                "then rerun this script in the same shell (or pass "
                "--gpp=$(which g++) once one is loaded).")
        if loaded:
            hint += "\nCurrently loaded modules: "+", ".join(loaded)
        hints.append(hint)
    if is_conda:
        conda_prefix = os.environ.get("CONDA_PREFIX")
        writable = bool(conda_prefix) and os.access(conda_prefix, os.W_OK)
        if writable:
            hints.append("Your conda environment can install its own "
                    "compiler:\n"
                    "  conda install -c conda-forge gxx_linux-64   "
                    "# Linux\n"
                    "  conda install -c conda-forge clangxx_osx-64 "
                    "# Mac")
        else:
            hints.append("Note: "+python_exe+" comes from a read-only or "
                    "centrally-managed environment (e.g. an HPC module), "
                    "so 'conda install' won't work here -- use --gpp to "
                    "point at a working compiler instead (a module-loaded "
                    "one, or the system compiler).")
    if not lmod_present and not is_conda:
        hints.append("  Ubuntu/Debian: sudo apt-get install g++\n"
                "  Mac:           xcode-select --install")
    return "\n".join(hints)


def _resolve_compiler(gpp_arg, check_gpp, is_conda, python_exe):
    candidates = list(_compiler_candidates(gpp_arg, is_conda, python_exe))
    if not candidates:
        _fail("No C++ compiler was found on PATH.\n"
              +_environment_hints(is_conda, python_exe)
              +"\nOr rerun with 'python install.py --gpp=<compiler>'.")

    tried = []
    for label, gpp in candidates:
        ok, reason = _try_compiler(gpp, check_gpp)
        if ok:
            if tried:
                print("Note: "+tried[-1][0]+" ("+tried[-1][1]+") did not "
                      "work, falling back to "+label+" ("+gpp+").")
            return gpp
        tried.append((label, gpp, reason))

    msg = "No working C++ compiler found. Tried:\n"
    for label, gpp, reason in tried:
        msg += "  - "+label+" ("+gpp+"): "+reason+"\n"
    msg += "\n"+_environment_hints(is_conda, python_exe)
    msg += "\nOr rerun with 'python install.py --gpp=<compiler>' once " \
            "you have a working one."
    _fail(msg)


def _resolve_blas(gpp, openblas, openblas_libdir, openblas_includedir):
    cfg = blaslapack.find_working_config(gpp, openblas=openblas,
            openblas_libdir=openblas_libdir,
            openblas_includedir=openblas_includedir)
    if cfg is None:
        msg = ("Could not find a working LAPACK/BLAS installation "
              "(tried system LAPACK/BLAS, conda's own libs, OpenBLAS in "
              "common locations, and the compiler's own multiarch lib "
              "dir).\nInstall LAPACK/BLAS, e.g.:\n"
              "  Ubuntu/Debian: sudo apt-get install liblapack-dev "
              "libblas-dev\n"
              "  Mac:           brew install openblas   # then rerun with "
              "--openblas --openblas_libdir=/usr/local/opt/openblas/lib "
              "--openblas_includedir=/usr/local/opt/openblas/include\n"
              "  conda:         conda install -c conda-forge liblapack "
              "libblas\n"
              "or pass --openblas_libdir/--openblas_includedir explicitly "
              "if it's installed in a non-standard location.")
        lmod_present, _ = cppversion.lmod_state()
        if lmod_present:
            msg += ("\nOn an HPC cluster, LAPACK/BLAS is usually provided "
                  "by a module too, e.g.:\n"
                  "  module avail openblas\n"
                  "  module load <stack> openblas/<version>\n"
                  "then pass its lib/include dirs explicitly with "
                  "--openblas --openblas_libdir=... "
                  "--openblas_includedir=... (module-provided libraries "
                  "usually aren't on the default linker search path).")
        _fail(msg)
    print("Using "+cfg["description"])
    return cfg


def _check_pybind11(python_exe):
    if not cppversion.has_pybind11():
        print("pybind11 not found, installing it into the current "
              "environment ("+python_exe+")...")
        proc = subprocess.run([python_exe, "-m", "pip", "install",
                "pybind11"])
        if proc.returncode != 0 or not cppversion.has_pybind11():
            _fail("Could not install pybind11 automatically.\n"
                  "Install it manually with:\n"
                  "  "+python_exe+" -m pip install pybind11\n"
                  "or:\n"
                  "  conda install -c conda-forge pybind11")

    proc = subprocess.run([python_exe, "-m", "pybind11", "--includes"],
            capture_output=True, text=True)
    if proc.returncode != 0:
        _fail("'"+python_exe+" -m pybind11 --includes' failed:\n"
              +proc.stdout+proc.stderr)

    ext_suffix = subprocess.run([python_exe, "-c",
            "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"],
            capture_output=True, text=True)
    if ext_suffix.returncode != 0 or not ext_suffix.stdout.strip():
        _fail("Could not determine the Python extension suffix for '"
              +python_exe+"'.")

    config_tool = python_exe+"-config"
    if not os.path.isfile(config_tool):
        _fail("'"+config_tool+"' was not found -- needed by the pybind11 "
              "extension Makefile to compute its own include/link flags.\n"
              "This usually means the interpreter's *-dev/*-devel package "
              "isn't installed, e.g.:\n"
              "  Ubuntu/Debian: sudo apt-get install python3-dev")


def check(args):
    """Run all preflight checks and return a validated BuildConfig, or
    exit(1) with an actionable message if the environment can't build."""
    print("### Checking build requirements ###")
    _check_os()
    _check_make()

    python_exe = sys.executable
    conda_prefix = os.environ.get("CONDA_PREFIX")
    is_conda = bool(conda_prefix) or os.path.isdir(
            os.path.join(os.path.dirname(os.path.dirname(
            os.path.realpath(python_exe))), "conda-meta"))
    print("Python interpreter: "+python_exe+(" (conda)" if is_conda
          else ""))

    system_python3 = shutil.which("python3")
    if system_python3 and os.path.realpath(system_python3) != \
            os.path.realpath(python_exe):
        print("Note: 'python3' on PATH ("+system_python3+") differs from "
              "the interpreter running this script ("+python_exe+"). The "
              "extension will still be built against "+python_exe+".")

    gpp = _resolve_compiler(args.gpp, args.check_gpp, is_conda, python_exe)
    print("C++ compiler: "+gpp+" (version "
          +str(cppversion.compiler_version(gpp))+")")

    blas_cfg = _resolve_blas(gpp, args.openblas, args.openblas_libdir,
            args.openblas_includedir)

    _check_pybind11(python_exe)
    print("pybind11: found")

    lmod_present, loaded = cppversion.lmod_state()
    if lmod_present and loaded:
        print("Loaded modules: "+", ".join(loaded))
        print("Remember to load these same modules whenever you use "
              "dmrgpy later (e.g. in job scripts) -- the compiled "
              "extension depends on them being present at import time too.")
    if is_conda:
        print("Conda/mamba environment: "+(conda_prefix or python_exe))
        print("Remember to activate this same environment whenever you "
              "use dmrgpy later.")

    print("### All requirements satisfied, proceeding to compilation ###\n")

    return BuildConfig(gpp=gpp, platform=blas_cfg["platform"],
            libflags=blas_cfg["libflags"],
            includeflags=blas_cfg["includeflags"],
            python_exe=python_exe, is_conda=is_conda)
