"""The build step of install.py: ITensor, then the pybind11 extension.

Every choice made here (compiler, flags, BLAS/LAPACK) comes from an
already validated installtk.requirements.BuildConfig, so nothing is
guessed again. What the build used is written to two files in the backend
folder: build_config.mk, which the Makefile reads so that a later
"make pybind" by hand uses the same compiler and interpreter, and
build_info.json, the record that dmrgpy reads (src/dmrgpy/buildinfo.py).
"""
import collections
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import time

from installtk import DEFAULT_ITENSOR_VERSION
from installtk import buildflags

ROOT = os.path.realpath(os.path.join(os.path.dirname(
        os.path.realpath(__file__)), ".."))
TAIL_LINES = 40 # lines of a failed step shown on the screen


def backend_dir(version, root=ROOT):
    return os.path.join(root, "src", "dmrgpy", "mpscpp"+str(version))


def log_path(version, root=ROOT):
    return os.path.join(root, "install_v"+str(version)+".log")


def _buildinfo(root=ROOT):
    """src/dmrgpy/buildinfo.py as a module, loaded from its file: importing
    the dmrgpy package here would import the very extension being built."""
    path = os.path.join(root, "src", "dmrgpy", "buildinfo.py")
    spec = importlib.util.spec_from_file_location("_dmrgpy_buildinfo", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fail(msg):
    print("\n### Build failed ###")
    print(msg)
    sys.exit(1)


def _run(cmd, cwd, log=None, verbose=False):
    """Run a build command and stop with a clear error if it fails.

    The output goes to the log file. On the screen it is one line that
    follows the build (on a terminal) or nothing (in a job script), and
    with `verbose` all of it; when the step fails, its last lines are
    printed together with the path of the log.

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
    shown = " ".join(cmd)
    print("$ "+shown+"   (in "+os.path.relpath(cwd, ROOT)+")")
    env = dict(os.environ, PWD=cwd)
    t0 = time.time()
    follow = (not verbose) and sys.stdout.isatty()
    width = shutil.get_terminal_size((80, 20)).columns - 1
    tail = collections.deque(maxlen=TAIL_LINES)
    if log is not None:
        log.write("\n$ "+shown+"   (in "+cwd+")\n")
        log.flush()
    proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, errors="replace")
    for line in proc.stdout:
        tail.append(line)
        if log is not None:
            log.write(line)
        if verbose:
            sys.stdout.write(line)
        elif follow and line.strip():
            sys.stdout.write("\r"+line.strip()[:width].ljust(width))
            sys.stdout.flush()
    proc.wait()
    if follow:
        sys.stdout.write("\r"+" "*width+"\r")
    if log is not None:
        log.flush()
    if proc.returncode != 0:
        if not verbose:
            sys.stdout.write("".join(tail))
        msg = "Step failed: "+shown
        if log is not None:
            msg += "\nThe full output is in "+log.name
        _fail(msg)
    print("  done in %.0f s"%(time.time()-t0))


def _git(args, cwd):
    """Output of a git command, or None when git is missing or fails."""
    try:
        proc = subprocess.run(["git"]+args, cwd=cwd, capture_output=True,
                text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return proc.stdout.strip() if proc.returncode == 0 else None


def itensor_key(itensor_dir, cccom, config):
    """A hash of everything the ITensor static library depends on: the
    vendored sources (the git tree of the folder) and the configuration
    written into options.mk. None when the sources cannot be pinned down
    (no git, or the folder has local changes), in which case the library
    is always rebuilt."""
    tree = _git(["rev-parse", "HEAD:./"], itensor_dir)
    dirty = _git(["status", "--porcelain", "--", "."], itensor_dir)
    if not tree or dirty is None or dirty != "":
        return None
    payload = json.dumps([tree, cccom, config.platform, config.libflags,
            config.includeflags, _version_string(config)])
    return hashlib.sha256(payload.encode()).hexdigest()


def _version_string(config):
    """The compiler version as text. It is part of what a build depends
    on: a system compiler upgraded in place keeps its path."""
    return ".".join(str(n) for n in
            (getattr(config, "compiler_version", None) or ())) or None


def options_mk(template, cccom, config):
    """The text of ITensor's options.mk: the template options.save with the
    resolved compiler command and BLAS/LAPACK configuration substituted."""
    out = re.sub(r"(?m)^CCCOM=.*$", lambda m: "CCCOM="+cccom, template,
            count=1)
    out = re.sub(r"(?m)^PLATFORM=.*$",
            lambda m: "PLATFORM="+config.platform, out, count=1)
    replacement = "BLAS_LAPACK_LIBFLAGS="+config.libflags
    if config.includeflags:
        replacement += "\nBLAS_LAPACK_INCLUDEFLAGS="+config.includeflags
    out = re.sub(r"(?m)^BLAS_LAPACK_LIBFLAGS=.*$", lambda m: replacement,
            out, count=1)
    return out


def _write_if_changed(path, text):
    try:
        if open(path).read() == text:
            return False
    except OSError:
        pass
    with open(path, "w") as f:
        f.write(text)
    return True


def writemk(itensor_dir, config, cccom):
    """Write options.mk from options.save. Left alone when it already says
    the same, so that its time stamp does not move."""
    template = open(os.path.join(itensor_dir, "options.save")).read()
    return _write_if_changed(os.path.join(itensor_dir, "options.mk"),
            options_mk(template, cccom, config))


def build_config_mk(cccom, python_exe):
    """The text of build_config.mk, which the backend's Makefile includes."""
    return ("# Written by install.py. The Makefile of this folder reads it, "
            "so that\n# \"make pybind\" uses the interpreter and the "
            "compiler install.py validated.\n# Delete it to go back to the "
            "Makefile's own detection.\n"
            "DMRGPY_PYTHON="+python_exe+"\n"
            "DMRGPY_PYBIND_CCCOM="+cccom+"\n")


def _record(version, config, cccom, key):
    """The fields of build_info.json that install.py knows."""
    return dict(
        itensor_version=version,
        compiler=config.gpp,
        compiler_version=_version_string(config),
        compiler_family=config.family,
        cccom=cccom,
        blas=config.blas_description or config.platform,
        platform=config.platform,
        libflags=config.libflags,
        includeflags=config.includeflags,
        python=config.python_exe,
        python_version=sys.version.split()[0],
        modules=list(config.modules),
        conda_prefix=config.conda_prefix,
        commit=_git(["rev-parse", "HEAD"], ROOT),
        itensor_key=key,
    )


def compile(config, version=DEFAULT_ITENSOR_VERSION, jobs=1, force=False,
        verbose=False):
    """Compile ITensor and the in-process pybind11 extension for the given
    C++ DMRG backend version (2 = ITensor v2/mpscpp2, 3 = ITensor v3/
    mpscpp3).

    The ITensor library is the slow part and depends only on the vendored
    sources and on the configuration, so it is rebuilt only when one of
    them changed since the build recorded in build_info.json, or with
    `force`. The extension is cleaned and rebuilt when the compiler, the
    interpreter or the BLAS/LAPACK flags changed, or with `force`;
    otherwise make decides from its sources, as after a C++ edit."""
    if version not in buildflags.VERSIONS:
        raise ValueError("Unknown ITensor backend version: "+str(version))
    mpscpp_dir = backend_dir(version)
    itensor_dir = os.path.join(mpscpp_dir, "ITensor")
    buildinfo = _buildinfo()
    cccom = buildflags.cccom(config.gpp, version, config.family)

    key = itensor_key(itensor_dir, cccom, config)
    previous = buildinfo.read(mpscpp_dir) or {}
    library = os.path.join(itensor_dir, "lib", "libitensor.a")
    up_to_date = (not force and key is not None
            and previous.get("itensor_key") == key
            and os.path.isfile(library)
            and os.path.isfile(os.path.join(itensor_dir, "this_dir.mk")))

    with open(log_path(version), "w") as log:
        log.write("install.py, ITensor v%s, %s\nCCCOM=%s\nBLAS/LAPACK: %s "
                "%s\nPython: %s\n"%(version, time.strftime("%Y-%m-%d %H:%M:%S"),
                cccom, config.platform, config.libflags, config.python_exe))
        writemk(itensor_dir, config, cccom)
        if up_to_date:
            print("### ITensor v"+str(version)+" is already built with this "
                  "configuration, not rebuilding it (--force rebuilds) ###")
        else:
            print("### Compiling ITensor v"+str(version)+" with "+str(jobs)
                  +" job"+("s" if jobs != 1 else "")
                  +" (this is the slow step) ###")
            _run(["make", "clean"], itensor_dir, log, verbose)
            # ITensor/Makefile's own "itensor" target runs a bare,
            # targetless "cd itensor && $(MAKE)" -- which defaults to
            # itensor/Makefile's own first target, "install: build debug",
            # compiling a second, unused libitensor-g.a debug library
            # alongside the real one (ITENSOR_LIBNAMES=itensor, never
            # itensor-g, so nothing here ever links against it). Run
            # "configure" (writes this_dir.mk/config.h, same as the
            # "itensor" target's prerequisite) ourselves, then invoke
            # itensor/Makefile directly with an explicit "build" goal to
            # skip that duplicate compile. CCCOM already correct via
            # options.mk.
            _run(["make", "configure"], itensor_dir, log, verbose)
            _run(["make", "-j"+str(jobs), "build"],
                 os.path.join(itensor_dir, "itensor"), log, verbose)

        print("### Compiling the in-process pybind11 extension (v"
              +str(version)+") ###")
        _write_if_changed(os.path.join(mpscpp_dir, "build_config.mk"),
                build_config_mk(cccom, config.python_exe))
        same_config = (not force
                and previous.get("cccom") == cccom
                and previous.get("compiler_version")
                        == _version_string(config)
                and previous.get("python") == config.python_exe
                and previous.get("libflags") == config.libflags
                and previous.get("includeflags") == config.includeflags)
        if not same_config:
            _run(["make", "clean"], mpscpp_dir, log, verbose)
        _run(["make", "pybind"], mpscpp_dir, log, verbose)

    buildinfo.stamp(mpscpp_dir, **_record(version, config, cccom, key))
    print("Build log: "+os.path.relpath(log_path(version), os.getcwd()))


def compile_extension_only(version=DEFAULT_ITENSOR_VERSION, force=False,
        verbose=False):
    """Rebuild only the pybind11 extension against the ITensor library
    that is already there, with the compiler and the interpreter of the
    last full install (build_config.mk). This is the step after a C++
    edit of chain_session.h or of another header of the backend.

    A checkout built before install.py wrote build_config.mk has none; the
    Makefile's own detection then decides, as "make pybind" always did."""
    mpscpp_dir = backend_dir(version)
    itensor_dir = os.path.join(mpscpp_dir, "ITensor")
    missing = [name for name, path in (
            ("the ITensor library", os.path.join(itensor_dir, "lib",
                    "libitensor.a")),
            ("options.mk", os.path.join(itensor_dir, "options.mk")),
            ("this_dir.mk", os.path.join(itensor_dir, "this_dir.mk")))
            if not os.path.isfile(path)]
    if missing:
        _fail("--extension-only needs a full install of ITensor v"
              +str(version)+" first ("+", ".join(missing)+" not found). "
              "Run python install.py --itensor-version="+str(version)+".")
    if not os.path.isfile(os.path.join(mpscpp_dir, "build_config.mk")):
        print("Note: this backend was built before install.py recorded its "
              "compiler (no build_config.mk), so the Makefile picks the "
              "compiler itself. A full python install.py records it.")
    with open(log_path(version), "w") as log:
        log.write("install.py --extension-only, ITensor v%s, %s\n"
                %(version, time.strftime("%Y-%m-%d %H:%M:%S")))
        print("### Compiling the in-process pybind11 extension (v"
              +str(version)+") ###")
        if force:
            _run(["make", "clean"], mpscpp_dir, log, verbose)
        _run(["make", "pybind"], mpscpp_dir, log, verbose)
    _buildinfo().stamp(mpscpp_dir, commit=_git(["rev-parse", "HEAD"], ROOT))
    print("Build log: "+os.path.relpath(log_path(version), os.getcwd()))


def recorded_python(version):
    """The interpreter build_config.mk of this backend names, or None."""
    try:
        text = open(os.path.join(backend_dir(version),
                "build_config.mk")).read()
    except OSError:
        return None
    m = re.search(r"(?m)^DMRGPY_PYTHON=(.*)$", text)
    return m.group(1).strip() if m else None


def selftest(python_exe, versions):
    """Run python -m dmrgpy.doctor --selftest for `versions` in a fresh
    process of `python_exe`, with this checkout's src first on the path
    (never another copy of dmrgpy the interpreter may know). Returns True
    when it passed."""
    src = os.path.join(ROOT, "src")
    env = dict(os.environ)
    env["PYTHONPATH"] = src+(os.pathsep+env["PYTHONPATH"]
            if env.get("PYTHONPATH") else "")
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    print("### Self-test: importing the extension in a fresh process and "
          "solving a small chain ###")
    proc = subprocess.run([python_exe, "-W", "ignore", "-m",
            "dmrgpy.doctor", "--selftest"]+[str(v) for v in versions],
            env=env, cwd=ROOT)
    return proc.returncode == 0
