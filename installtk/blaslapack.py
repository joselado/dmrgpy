import os
import shutil
import subprocess

from . import cppversion

# a minimal program that forces the linker to resolve a BLAS symbol,
# without needing to actually call it with valid arguments (the trial
# binary is never executed, only linked)
_TRIAL_SOURCE = """
extern "C" void dgemm_();
int main(){ void (*p)() = dgemm_; (void)p; return 0; }
"""

_OPENBLAS_DIRS = [
    "/usr/local/opt/openblas",   # Homebrew, Intel Mac
    "/opt/homebrew/opt/openblas", # Homebrew, Apple Silicon
    "/usr/local",
    "/usr",
]


def _openblas_flags(libdir=None, includedir=None):
    libflags = "-lpthread"
    if libdir:
        libflags += " -L" + libdir
    libflags += " -lopenblas"
    includeflags = "-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE"
    if includedir:
        includeflags = "-I" + includedir + " " + includeflags
    return dict(platform="openblas", libflags=libflags,
            includeflags=includeflags, description="OpenBLAS"
            + (" at "+libdir if libdir else ""))


def _lapack_flags(libdir=None):
    libflags = "-lpthread"
    if libdir:
        libflags += " -L" + libdir
    libflags += " -lblas -llapack"
    return dict(platform="lapack", libflags=libflags, includeflags="",
            description="system LAPACK/BLAS"
            + (" at "+libdir if libdir else ""))


def _system_libdir_for(libname, sys_gpp="g++"):
    """Ask the system compiler where it would actually find `libname`
    (e.g. libblas.so) -- this is the directory the linker resolves
    -lblas/-llapack against, which is the only reliable way to find it
    on Debian/Ubuntu multiarch systems where it lives under
    /usr/lib/x86_64-linux-gnu rather than directly under /usr/lib."""
    if not shutil.which(sys_gpp):
        return None
    try:
        out = subprocess.run([sys_gpp, "-print-file-name="+libname],
                capture_output=True, text=True, timeout=15).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    if not out or out == libname:
        return None # compiler didn't find it either
    d = os.path.dirname(out)
    return d or None


def candidates(gpp, openblas=False, openblas_libdir=None,
        openblas_includedir=None):
    """Yield candidate BLAS/LAPACK configurations to try, in order,
    for the given (already resolved) compiler `gpp`."""
    if openblas:
        # explicit user request: only try this, no cascading past it
        yield _openblas_flags(openblas_libdir, openblas_includedir)
        return

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        libdir = os.path.join(conda_prefix, "lib")
        yield _lapack_flags(libdir)

    yield _lapack_flags() # plain system default, e.g. /usr/lib

    for d in _OPENBLAS_DIRS:
        yield _openblas_flags(os.path.join(d, "lib"),
                os.path.join(d, "include"))

    # Debian/Ubuntu multiarch workaround: gpp may be a conda-provided
    # compiler whose linker doesn't share the system compiler's default
    # library search path, so plain "-lblas -llapack" above can fail even
    # though the libraries are present. Ask the *system* g++ where it
    # would find them and pass that directory explicitly.
    sysdir_blas = _system_libdir_for("libblas.so")
    sysdir_lapack = _system_libdir_for("liblapack.so")
    if sysdir_blas or sysdir_lapack:
        libdirs = {d for d in (sysdir_blas, sysdir_lapack) if d}
        libflags = "-lpthread " + " ".join("-L"+d for d in sorted(libdirs)) \
                + " -lblas -llapack"
        yield dict(platform="lapack", libflags=libflags, includeflags="",
                description="system LAPACK/BLAS (multiarch lib dir "
                + ", ".join(sorted(libdirs))+")")


def find_working_config(gpp, openblas=False, openblas_libdir=None,
        openblas_includedir=None):
    """Try each candidate config in order, actually compiling+linking a
    trial program against it with `gpp`. Returns the first config dict
    that links successfully, or None if none of them do."""
    for cfg in candidates(gpp, openblas=openblas,
            openblas_libdir=openblas_libdir,
            openblas_includedir=openblas_includedir):
        link_flags = cfg["libflags"].split()
        extra_flags = cfg["includeflags"].split() if cfg["includeflags"] \
                else []
        ok, _ = cppversion.trial_compile(gpp, extra_flags=extra_flags,
                link_flags=link_flags, source=_TRIAL_SOURCE)
        if ok:
            return cfg
    return None
