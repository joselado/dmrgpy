import os
import re
import subprocess
import sys


def _run(cmd, cwd):
    """Run a build command, streaming its output, and stop with a clear
    error (instead of silently continuing, as bare os.system() calls used
    to) if it fails."""
    print("$ "+" ".join(cmd)+"   (in "+cwd+")")
    proc = subprocess.run(cmd, cwd=cwd)
    if proc.returncode != 0:
        print("\n### Build step failed: "+" ".join(cmd)+" ###")
        sys.exit(proc.returncode)


def compile(config):
    """Compile ITensor and the in-process pybind11 extension using an
    already-validated installtk.requirements.BuildConfig -- every choice
    made here (compiler, BLAS/LAPACK flags) was already trial-compiled and
    trial-linked during the preflight check, so nothing is re-guessed."""
    root = os.path.dirname(os.path.realpath(__file__))+"/.." # main path
    itensor_dir = root+"/src/dmrgpy/mpscpp2/ITensor"
    mpscpp2_dir = root+"/src/dmrgpy/mpscpp2"

    writemk(itensor_dir, config)
    # CCCOM must carry the same "-m64 -std=c++14 -fPIC" suffix writemk()
    # put in options.mk -- without -fPIC in particular, the final pybind
    # extension link fails with a TLS relocation error, since a shared
    # object needs position-independent code throughout.
    cccom = config.gpp+" -m64 -std=c++14 -fPIC"

    print("### Compiling ITensor (this is the slow step, several minutes) ###")
    _run(["make", "clean"], itensor_dir)
    _run(["make"], itensor_dir) # CCCOM already correct via options.mk

    print("### Compiling the in-process pybind11 extension ###")
    _run(["make", "clean"], mpscpp2_dir)
    _run(["make", "pybind",
          "PYBIND_CCCOM="+cccom,
          "PYTHON="+config.python_exe],
         mpscpp2_dir)


def writemk(itensor_dir, config):
    """Write options.mk from options.save, with the resolved compiler and
    BLAS/LAPACK configuration substituted in."""
    out = open(itensor_dir+"/options.save").read()
    out = re.sub(r"(?m)^CCCOM=.*$",
            "CCCOM="+config.gpp+" -m64 -std=c++14 -fPIC", out, count=1)
    out = re.sub(r"(?m)^PLATFORM=.*$", "PLATFORM="+config.platform, out,
            count=1)
    replacement = "BLAS_LAPACK_LIBFLAGS="+config.libflags
    if config.includeflags:
        replacement += "\nBLAS_LAPACK_INCLUDEFLAGS="+config.includeflags
    out = re.sub(r"(?m)^BLAS_LAPACK_LIBFLAGS=.*$", replacement, out,
            count=1)
    open(itensor_dir+"/options.mk", "w").write(out)
