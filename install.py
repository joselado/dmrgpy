#!/usr/bin/python3
"""Compile ITensor and the in-process pybind11 DMRG extension, check that
it works, and make dmrgpy importable from the current Python interpreter.

It runs in four steps. First every requirement (C++ compiler, LAPACK/BLAS,
pybind11, make) is checked, with a trial compile and link that uses the
flags of the build itself; a conda-provided compiler is preferred when the
interpreter is a conda one. Only then does the compilation start. After
it, the extension is imported in a fresh process and a small chain is
solved against exact diagonalization. Last, dmrgpy is linked into the
interpreter's site-packages and DMRGROOT is added to the shell rc file
(--no-link and --no-shellrc skip these two, --build-only skips both).

    python install.py                    # build the default backend (v3)
    python install.py --doctor           # only the requirement checks
    python install.py --extension-only   # after a C++ edit: rebuild only
                                         # the extension, not ITensor
    python -m dmrgpy.doctor              # later: which backends work here
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time

from installtk import DEFAULT_ITENSOR_VERSION
from installtk import buildflags


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--itensor-version", "--itensor_version",
            dest="itensor_version", default=str(DEFAULT_ITENSOR_VERSION),
            choices=["2", "3", "both"],
            help="Which C++ DMRG backend(s) to compile: 2 = ITensor v2 "
                 "(mpscpp2, itensor_version=2 in Python), 3 = ITensor v3 "
                 "(mpscpp3, itensor_version=3), both = compile both, "
                 "one after the other (default: %(default)s, matching the "
                 "default backend in Python).")
    parser.add_argument("--doctor", "--check-only", dest="doctor",
            action="store_true", default=False,
            help="Only run the requirement checks (compiler, LAPACK/BLAS, "
                 "pybind11) and report what would be used, without "
                 "compiling anything. Useful for diagnosing a broken "
                 "environment (e.g. an HPC cluster) before the ITensor "
                 "build.")

    build = parser.add_argument_group("what to build")
    build.add_argument("--extension-only", "--extension_only",
            dest="extension_only", action="store_true", default=False,
            help="Rebuild only the pybind11 extension, with the compiler "
                 "of the last full install and the ITensor library that is "
                 "already there. The step after editing chain_session.h or "
                 "another C++ header of the backend.")
    build.add_argument("-j", "--jobs", type=int, default=None,
            help="Parallel compile jobs for ITensor (default: the cores "
                 "available, at most 4; use 1 on a busy login node).")
    build.add_argument("--force", action="store_true", default=False,
            help="Rebuild ITensor even when it is already built with the "
                 "same sources and configuration.")
    build.add_argument("--verbose", action="store_true", default=False,
            help="Print the whole output of the build instead of one line "
                 "that follows it (it is always written to "
                 "install_v<version>.log).")

    comp = parser.add_argument_group("compiler and LAPACK/BLAS")
    comp.add_argument("--gpp", default=None,
            help="C++ compiler to use. Default: auto-detect (a "
                 "conda-provided compiler when run from a conda Python, "
                 "else the system g++).")
    comp.add_argument("--check-gpp", "--check_gpp",
            dest="check_gpp", default=True,
            action=argparse.BooleanOptionalAction,
            help="Verify the C++ compiler's version and that it can "
                 "actually compile (default: on).")
    comp.add_argument("--openblas", default=False,
            action=argparse.BooleanOptionalAction,
            help="Use OpenBLAS instead of auto-detecting LAPACK/BLAS.")
    comp.add_argument("--openblas-libdir", "--openblas_libdir",
            dest="openblas_libdir", default=None,
            help="Path to the directory containing libopenblas (only used "
                 "together with --openblas), for clusters/modules where "
                 "OpenBLAS is not on the default linker search path.")
    comp.add_argument("--openblas-includedir", "--openblas_includedir",
            dest="openblas_includedir", default=None,
            help="Path to the directory containing OpenBLAS headers (only "
                 "used together with --openblas).")

    after = parser.add_argument_group("after the build")
    after.add_argument("--selftest", default=True,
            action=argparse.BooleanOptionalAction,
            help="Import the new extension in a fresh process and solve a "
                 "small chain against exact diagonalization (default: on).")
    after.add_argument("--link", default=True,
            action=argparse.BooleanOptionalAction,
            help="Link dmrgpy into the site-packages of this interpreter "
                 "(default: on).")
    after.add_argument("--shellrc", default=True,
            action=argparse.BooleanOptionalAction,
            help="Add DMRGROOT to the shell rc file (default: on).")
    after.add_argument("--build-only", "--build_only", dest="build_only",
            action="store_true", default=False,
            help="Build and test only: the same as --no-link --no-shellrc. "
                 "Nothing outside this folder is changed.")
    return parser


def bare_import_path(python_exe):
    """Where `import dmrgpy` resolves for `python_exe` from a neutral
    folder and with no PYTHONPATH, or None when it does not import."""
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    code = ("import importlib.util, os; s = importlib.util.find_spec("
            "'dmrgpy'); print(os.path.realpath(os.path.dirname(s.origin)) "
            "if s and s.origin else '')")
    with tempfile.TemporaryDirectory(prefix="dmrgpy_install_") as tmp:
        try:
            proc = subprocess.run([python_exe, "-c", code], cwd=tmp, env=env,
                    capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError):
            return None
    return proc.stdout.strip() or None


def main(argv=None):
    args = build_parser().parse_args(argv)
    try: # in a job script the output is a file: show each step as it starts
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass
    if args.build_only:
        args.link = args.shellrc = False
    if args.jobs is not None and args.jobs < 1:
        sys.exit("--jobs must be at least 1")
    jobs = args.jobs or buildflags.default_jobs()
    versions = ([2, 3] if args.itensor_version == "both"
            else [int(args.itensor_version)])
    t0 = time.time()

    from installtk import install2
    if args.extension_only:
        python_exe = sys.executable
        for version in versions:
            recorded = install2.recorded_python(version)
            if recorded and os.path.realpath(recorded) != \
                    os.path.realpath(sys.executable):
                print("Note: the extension of ITensor v"+str(version)
                      +" is built for "+recorded+", the interpreter of the "
                      "last full install, not for "+sys.executable+".")
                python_exe = recorded
            install2.compile_extension_only(version, force=args.force,
                    verbose=args.verbose)
    else:
        from installtk import requirements
        config = requirements.check(args, versions) # check everything first
        if args.doctor:
            print("### --doctor: requirements OK, stopping before "
                  "compilation ###")
            return 0
        python_exe = config.python_exe
        for version in versions: # only reached once the check passed
            install2.compile(config, version=version, jobs=jobs,
                    force=args.force, verbose=args.verbose)

    passed = None
    if args.selftest:
        passed = install2.selftest(python_exe, versions)

    if not args.extension_only:
        if args.link:
            from installtk import addpythonpath
            addpythonpath.addpath() # add the library to the Python path
            src = os.path.realpath(os.path.join(install2.ROOT, "src",
                    "dmrgpy"))
            found = bare_import_path(python_exe)
            if found is not None and found != src:
                print("Note: with this interpreter a bare 'import dmrgpy' "
                      "finds "+found+", not this checkout ("+src+"). "
                      "Remove that copy, or set PYTHONPATH="
                      +os.path.dirname(src)+".")
        else:
            print("Not linking dmrgpy into site-packages (--no-link). To "
                  "use this build: PYTHONPATH="
                  +os.path.join(install2.ROOT, "src"))
        if args.shellrc:
            from installtk import addsystem
            addsystem.addrc() # add DMRGROOT to the shell rc file

    print("")
    what = "ITensor v"+" and v".join(str(v) for v in versions)
    if passed is False:
        print("### "+what+" compiled in %.0f s, but the self-test FAILED "
              "(see above). python -m dmrgpy.doctor shows the state of the "
              "backends. ###"%(time.time()-t0))
        return 1
    print("### "+what+" ready in %.0f s%s ###"%(time.time()-t0,
          ", self-test passed" if passed else ""))
    print("Later, 'python -m dmrgpy.doctor' says which backends work in a "
          "session and how each was built.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
