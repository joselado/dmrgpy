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


def _resolve_compiler(gpp_arg, check_gpp, is_conda, python_exe):
    if gpp_arg is not None:
        gpp = gpp_arg
    else:
        conda_gpp = cppversion.find_conda_compiler(python_exe) if is_conda \
                else None
        if conda_gpp is not None:
            gpp = conda_gpp
        elif is_conda:
            print("Detected a conda Python ("+python_exe+") but no "
                  "conda-provided C++ compiler was found next to it.")
            print("This risks a crash once the compiled extension is "
                  "loaded alongside conda's own numpy/scipy (they bundle "
                  "their own libstdc++). Recommended fix:")
            print("  conda install -c conda-forge gxx_linux-64   "
                  "# Linux")
            print("  conda install -c conda-forge clangxx_osx-64 "
                  "# Mac")
            print("Falling back to the system compiler for now.\n")
            gpp = "g++"
        else:
            gpp = "g++"

    if not shutil.which(gpp) and not os.path.isfile(gpp):
        _fail("C++ compiler '"+gpp+"' was not found.\n"
              "Install a GNU C++ compiler (version 6 or higher), e.g.:\n"
              "  Ubuntu/Debian: sudo apt-get install g++\n"
              "  Mac:           xcode-select --install\n"
              "  conda:         conda install -c conda-forge gxx_linux-64\n"
              "then rerun with 'python install.py --gpp=<compiler>' if "
              "needed.")

    if check_gpp:
        version = cppversion.compiler_version(gpp)
        if version is None or not cppversion.correct_version(gpp):
            _fail("C++ compiler '"+gpp+"' could not be used or is too old "
                  "(found version "+str(version)+", need >= 6).\n"
                  "Install a newer GNU C++ compiler, e.g.:\n"
                  "  Ubuntu/Debian: sudo add-apt-repository "
                  "ppa:ubuntu-toolchain-r/test -y && sudo apt-get update -y "
                  "&& sudo apt-get install g++-11 -y\n"
                  "  conda:         conda install -c conda-forge "
                  "gxx_linux-64\n"
                  "then rerun with 'python install.py --gpp=<compiler>'.")
        ok, out = cppversion.trial_compile(gpp)
        if not ok:
            _fail("C++ compiler '"+gpp+"' failed to compile a trivial "
                  "test program:\n"+out)
    return gpp


def _resolve_blas(gpp, openblas, openblas_libdir, openblas_includedir):
    cfg = blaslapack.find_working_config(gpp, openblas=openblas,
            openblas_libdir=openblas_libdir,
            openblas_includedir=openblas_includedir)
    if cfg is None:
        _fail("Could not find a working LAPACK/BLAS installation "
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

    print("### All requirements satisfied, proceeding to compilation ###\n")

    return BuildConfig(gpp=gpp, platform=blas_cfg["platform"],
            libflags=blas_cfg["libflags"],
            includeflags=blas_cfg["includeflags"],
            python_exe=python_exe, is_conda=is_conda)
