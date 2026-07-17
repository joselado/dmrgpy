import os
import re
import subprocess
import sys

from installtk import DEFAULT_ITENSOR_VERSION


# Per-backend-version knobs: which mpscpp<version> directory to build, and
# the language-standard suffix CCCOM needs. ITensor v2 only needs C++14; v3
# requires C++17 with the concepts TS extension (-fconcepts), see
# mpscpp3/ITensor/options.mk.sample's own CCCOM line.
_STD_FLAGS = {
    2: "-std=c++14",
    3: "-std=c++17 -fconcepts",
}


def _run(cmd, cwd):
    """Run a build command, streaming its output, and stop with a clear
    error (instead of silently continuing, as bare os.system() calls used
    to) if it fails.

    Passes an explicit PWD in the child's environment matching `cwd`:
    subprocess's own cwd= argument changes the child's actual working
    directory (chdir()), but does NOT update the inherited PWD environment
    variable, and ITensor's own Makefile writes this_dir.mk (the include
    every other ITensor Makefile relies on for its own source/include path)
    from GNU Make's $(PWD) -- which, absent a fresh env var, falls back to
    whatever PWD this Python process itself inherited (e.g. the caller's
    shell cwd at the time install.py was invoked), not `cwd`. Confirmed
    directly: running this from an unrelated directory produced a
    this_dir.mk pointing at that directory instead of ITensor/, breaking
    every subsequent #include.
    """
    print("$ "+" ".join(cmd)+"   (in "+cwd+")")
    env = dict(os.environ, PWD=cwd)
    proc = subprocess.run(cmd, cwd=cwd, env=env)
    if proc.returncode != 0:
        print("\n### Build step failed: "+" ".join(cmd)+" ###")
        sys.exit(proc.returncode)


def compile(config, version=DEFAULT_ITENSOR_VERSION):
    """Compile ITensor and the in-process pybind11 extension for the given
    C++ DMRG backend version (2 = ITensor v2/mpscpp2, 3 = ITensor v3/
    mpscpp3), using an already-validated installtk.requirements.BuildConfig
    -- every choice made here (compiler, BLAS/LAPACK flags) was already
    trial-compiled and trial-linked during the preflight check, so nothing
    is re-guessed."""
    if version not in _STD_FLAGS:
        raise ValueError("Unknown ITensor backend version: "+str(version))
    std_flag = _STD_FLAGS[version]

    root = os.path.dirname(os.path.realpath(__file__))+"/.." # main path
    mpscpp_dir = root+"/src/dmrgpy/mpscpp"+str(version)
    itensor_dir = mpscpp_dir+"/ITensor"

    writemk(itensor_dir, config, std_flag)
    # CCCOM must carry the same "-m64 <std_flag> -fPIC" suffix writemk()
    # put in options.mk -- without -fPIC in particular, the final pybind
    # extension link fails with a TLS relocation error, since a shared
    # object needs position-independent code throughout.
    cccom = config.gpp+" -m64 "+std_flag+" -fPIC"

    print("### Compiling ITensor v"+str(version)+
          " (this is the slow step, several minutes) ###")
    _run(["make", "clean"], itensor_dir)
    _run(["make"], itensor_dir) # CCCOM already correct via options.mk

    print("### Compiling the in-process pybind11 extension (v"+
          str(version)+") ###")
    _run(["make", "clean"], mpscpp_dir)
    _run(["make", "pybind",
          "PYBIND_CCCOM="+cccom,
          "PYTHON="+config.python_exe],
         mpscpp_dir)


def writemk(itensor_dir, config, std_flag):
    """Write options.mk from options.save, with the resolved compiler,
    language standard and BLAS/LAPACK configuration substituted in."""
    out = open(itensor_dir+"/options.save").read()
    out = re.sub(r"(?m)^CCCOM=.*$",
            "CCCOM="+config.gpp+" -m64 "+std_flag+" -fPIC", out, count=1)
    out = re.sub(r"(?m)^PLATFORM=.*$", "PLATFORM="+config.platform, out,
            count=1)
    replacement = "BLAS_LAPACK_LIBFLAGS="+config.libflags
    if config.includeflags:
        replacement += "\nBLAS_LAPACK_INCLUDEFLAGS="+config.includeflags
    out = re.sub(r"(?m)^BLAS_LAPACK_LIBFLAGS=.*$", replacement, out,
            count=1)
    open(itensor_dir+"/options.mk", "w").write(out)
