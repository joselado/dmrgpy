import glob
import os
import subprocess
import sys
import tempfile


def compiler_version(gpp="g++"):
    """Return the (major, minor, ...) version tuple of a C++ compiler,
    or None if it cannot be run at all."""
    try:
        out = subprocess.run([gpp, "--version"], capture_output=True,
                text=True, timeout=15).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    firstline = out.split("\n")[0].split()
    if not firstline: return None
    token = firstline[-1] # last token of the first line, e.g. "13.2.0"
    try:
        return tuple(int(n) for n in token.split(".") if n.isdigit())
    except ValueError:
        return None


def correct_version(gpp="g++", minimum=6):
    v = compiler_version(gpp=gpp)
    if v is None: return False
    return v[0] >= minimum


def has_pybind11():
    """Check whether pybind11 is importable, needed for the in-process
    extension build (see mpscpp2/bindings.cc)"""
    try:
        import pybind11
        return True
    except ImportError:
        return False


def find_conda_compiler(python_exe=None):
    """If `python_exe` (default: the running interpreter) lives inside a
    conda environment, return the path to a conda-provided C++ compiler
    sitting next to it (the `gxx_linux-64`/`clangxx_osx-64` packages install
    e.g. `x86_64-conda-linux-gnu-c++`). Returns None if not running under
    conda, or if no such compiler package is installed.

    Using this compiler (rather than the system one) for the *whole* build
    keeps a single, consistent libstdc++/BLAS ABI throughout: the resulting
    extension is loaded into the same Python process as conda's own
    numpy/scipy, which bundle their own libstdc++ -- building against the
    system compiler instead has reproducibly segfaulted in the past (see
    CLAUDE.md and the long comment in mpscpp2/Makefile).
    """
    python_exe = python_exe or sys.executable
    conda_prefix = os.environ.get("CONDA_PREFIX")
    bindir = os.path.dirname(os.path.realpath(python_exe))
    if not conda_prefix:
        # fall back to detecting conda via a conda-meta directory next to
        # the interpreter's prefix, in case CONDA_PREFIX isn't set (e.g.
        # invoked with an absolute interpreter path outside an activated env)
        prefix = os.path.dirname(bindir)
        if os.path.isdir(os.path.join(prefix, "conda-meta")):
            conda_prefix = prefix
        else:
            return None
    candidates = sorted(glob.glob(os.path.join(conda_prefix, "bin",
            "*-conda*-g++")))
    candidates += sorted(glob.glob(os.path.join(conda_prefix, "bin",
            "*-conda*-clang++")))
    for c in candidates:
        if os.access(c, os.X_OK):
            return c
    return None


def trial_compile(gpp, extra_flags=(), link_flags=(), source=None):
    """Try to compile (and link) a small C++ program with `gpp`, using
    `extra_flags` for compilation and `link_flags` for linking. Returns
    (success, output) where `output` is combined stdout/stderr, useful for
    diagnostics on failure. Runs entirely inside a throwaway temporary
    directory (never a fixed /tmp path, so parallel installs don't race).
    """
    if source is None:
        source = "int main(){return 0;}\n"
    with tempfile.TemporaryDirectory(prefix="dmrgpy_install_") as tmp:
        srcfile = os.path.join(tmp, "trial.cc")
        binfile = os.path.join(tmp, "trial.out")
        with open(srcfile, "w") as f:
            f.write(source)
        cmd = [gpp, "-std=c++14", "-fPIC"] + list(extra_flags) + \
                [srcfile, "-o", binfile] + list(link_flags)
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                    timeout=60)
        except (OSError, subprocess.SubprocessError) as e:
            return False, str(e)
        ok = proc.returncode == 0 and os.path.exists(binfile)
        return ok, proc.stdout + proc.stderr


if __name__ == "__main__":
    print("version:", compiler_version())
    print("correct version:", correct_version())
    print("conda compiler:", find_conda_compiler())
