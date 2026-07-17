#!/usr/bin/python3
"""Compile ITensor and the in-process pybind11 DMRG extension, then wire
dmrgpy into the current Python interpreter's path.

Runs in two phases: first every requirement (C++ compiler, LAPACK/BLAS,
pybind11, make) is checked -- and the compiler/BLAS choice is validated
with an actual trial compile+link, auto-detecting a conda-provided
compiler when run from an Anaconda/conda Python -- and only once all of
that succeeds does compilation begin.
"""
import argparse

parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--gpp", default=None,
        help="C++ compiler to use. Default: auto-detect (a conda-provided "
             "compiler when run from a conda Python, else the system g++).")
parser.add_argument("--check-gpp", "--check_gpp",
        dest="check_gpp", default=True,
        action=argparse.BooleanOptionalAction,
        help="Verify the C++ compiler's version and that it can actually "
             "compile (default: on).")
parser.add_argument("--openblas", default=False,
        action=argparse.BooleanOptionalAction,
        help="Use OpenBLAS instead of auto-detecting LAPACK/BLAS.")
parser.add_argument("--openblas_libdir", default=None,
        help="Path to the directory containing libopenblas (only used "
             "together with --openblas), for clusters/modules where "
             "OpenBLAS is not on the default linker search path.")
parser.add_argument("--openblas_includedir", default=None,
        help="Path to the directory containing OpenBLAS headers (only "
             "used together with --openblas).")
parser.add_argument("--itensor-version", "--itensor_version",
        dest="itensor_version", default="2", choices=["2", "3", "both"],
        help="Which C++ DMRG backend(s) to compile: 2 = ITensor v2 "
             "(mpscpp2, itensor_version=2 in Python), 3 = ITensor v3 "
             "(mpscpp3, itensor_version=3), both = compile both, "
             "one after the other (default: 2, matching the historical "
             "single-backend behavior).")
args = parser.parse_args()

versions = [2, 3] if args.itensor_version == "both" else [int(args.itensor_version)]

from installtk import requirements
config = requirements.check(args) # phase 1: check everything first

from installtk import install2
for version in versions:
    install2.compile(config, version=version) # phase 2: only reached once phase 1 passed

from installtk import addpythonpath
addpythonpath.addpath() # add the library to the Python path

from installtk import addsystem
addsystem.addrc() # add DMRGROOT to the shell rc file
