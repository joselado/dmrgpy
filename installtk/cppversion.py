import glob
import os
import re
import subprocess
import sys
import tempfile


def parse_version(first_line):
    """The version tuple in the first line of a compiler's --version
    output, or None. GCC prints it after an optional distribution note in
    parentheses ("g++ (Ubuntu 13.2.0-23ubuntu4) 13.2.0"); Clang prints it
    after the word "version" and puts a build number in parentheses after
    it ("Apple clang version 15.0.0 (clang-1500.3.9.4)"), which must not
    be read as the version."""
    m = re.search(r"version\s+(\d+(?:\.\d+)*)", first_line)
    if m is None:
        outside = re.sub(r"\([^)]*\)", " ", first_line)
        m = (re.search(r"(\d+(?:\.\d+)+)", outside)
                or re.search(r"(\d+(?:\.\d+)+)", first_line))
    if m is None:
        return None
    return tuple(int(n) for n in m.group(1).split("."))


def _version_output(gpp):
    try:
        return subprocess.run([gpp, "--version"], capture_output=True,
                text=True, timeout=15).stdout
    except (OSError, subprocess.SubprocessError):
        return None


def compiler_version(gpp="g++"):
    """Return the (major, minor, ...) version tuple of a C++ compiler,
    or None if it cannot be run at all."""
    out = _version_output(gpp)
    if not out:
        return None
    return parse_version(out.split("\n")[0])


def compiler_family(gpp="g++"):
    """"clang" or "gcc", from the compiler's own --version output (on a
    Mac, g++ is Clang). Anything that does not say clang counts as gcc,
    which leaves the flags as they have always been."""
    out = _version_output(gpp) or ""
    return "clang" if "clang" in out.lower() else "gcc"


def correct_version(gpp="g++", minimum=6):
    """Whether a GCC is recent enough. Clang numbers its versions
    differently, so for Clang only "it runs" is asked here, and the trial
    compile with the real flags of the backend decides."""
    v = compiler_version(gpp=gpp)
    if v is None: return False
    if compiler_family(gpp) == "clang": return True
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


def backend_missing(gpp, prog="cc1plus"):
    """True if `gpp`'s own driver can't locate its compilation backend
    (`cc1plus`, the actual C++ compiler -- `gpp` itself is just a driver
    that execs it). Uses `gpp -print-prog-name=`, the same lookup gcc/clang
    use internally, so this is a real check, not a guess -- and it's much
    cheaper and more specific than waiting for a full trial_compile() to
    fail with a raw "posix_spawnp: No such file or directory".

    This exists because `--version` and `shutil.which()` both pass for a
    compiler that's just a driver front-end with no matching backend
    package installed -- confirmed on Aalto's Triton cluster, where the
    `scicomp-python-env` module's bundled conda-style
    `x86_64-conda-linux-gnu-g++` wrapper resolves and runs `--version`
    fine but has no `cc1plus` behind it.
    """
    try:
        out = subprocess.run([gpp, "-print-prog-name="+prog],
                capture_output=True, text=True, timeout=15).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return True
    if not out or out == prog:
        # gcc/clang echo the bare program name back unchanged when they
        # can't resolve it to a path at all
        return True
    return not os.path.isfile(out)


def lmod_state():
    """Return (lmod_present, loaded_module_names) using the environment
    variables Lmod itself sets ($LMOD_CMD/$MODULEPATH, $LOADEDMODULES).
    Many HPC clusters (e.g. Aalto's Triton) provide no compiler at all on
    PATH until a module is loaded -- detecting this lets error messages
    point at 'module load ...' instead of the apt-get/brew hints that are
    meaningless (no sudo) or simply wrong on a cluster."""
    present = bool(os.environ.get("LMOD_CMD")) or \
            bool(os.environ.get("MODULEPATH"))
    loaded = [m for m in os.environ.get("LOADEDMODULES", "").split(":") if m]
    return present, loaded


def trial_compile(gpp, extra_flags=(), link_flags=(), source=None,
        std_flags=("-std=c++14",)):
    """Try to compile (and link) a small C++ program with `gpp`, using
    `std_flags` for the language standard, `extra_flags` for compilation
    and `link_flags` for linking. Returns
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
        cmd = [gpp] + list(std_flags) + ["-fPIC"] + list(extra_flags) + \
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
