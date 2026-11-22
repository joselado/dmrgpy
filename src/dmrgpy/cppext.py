"""
Single import point for the optional in-process pybind11 extensions
(mpscpp2/bindings.cc built against ITensor v2, mpscpp3/bindings.cc built
against ITensor v3), which replace the file-based Python<->C++ protocol
(see the migration plan) one task at a time. Both expose the same Chain
session API (see chain_session.h in each backend) -- callers elsewhere
never need to know which ITensor version actually backs self._session.

Also the import point for the third, non-C++ backend: itensor_version=
"python" selects pyitensor.chain.Chain, a pure-Python reimplementation of
exactly the ITensor v3 API subset mpscpp3/chain_session.h uses (see
pyitensor/__init__.py). It exposes the identical Chain(site_types)
constructor and method surface as the compiled extensions -- callers
never need to know self._session is a plain Python object rather than a
pybind11 handle either (see mps.py's MPS/multioperatortk/staticoperator.py's
StaticOperator, which only ever treat cpp_handle as opaque). Unlike the
C++ backends, "python" has no compiled-extension precondition at all (no
compiler, no pybind11, nothing to "make pybind" for), so available("python")
is always True.

Callers should use get_backend(version)/available(version) rather than
importing mpscpp2._dmrgcpp/mpscpp3._dmrgcpp/pyitensor.chain directly, since
a given C++ extension may not be built (pybind11 wasn't available at
install time, or "make pybind" was never run for that backend) -- in that
case every caller should fall back to ED, not fail outright (see mode.py).
"""

# Single source of truth for dmrgpy's default C++ DMRG backend version.
# Every other default (setup_cpp, get_backend, available) is derived from
# this one constant -- change it here only. Note it is the default *C++*
# version, not unconditionally the default backend: a chain built without
# an explicit itensor_version goes through default_backend() below, which
# picks the pure-Python one when this extension isn't compiled.
DEFAULT_ITENSOR_VERSION = 3

_backends = {} # version -> compiled _dmrgcpp module (or pyitensor.chain), or None if unavailable


def default_backend():
    """The backend a chain gets when its caller names no itensor_version:
    the default C++ version when that extension is actually compiled, and
    the pure-Python one (pyitensor) otherwise.

    The fallback is what makes a `pip install dmrgpy` usable. The wheel
    deliberately ships no C++ at all (see CLAUDE.md's "Packaging / PyPI"),
    so before this existed *every* default-backend chain in a pip install
    silently ran exact diagonalization -- mode.py's extension-not-compiled
    fallback -- which cannot reach the sizes an MPS backend exists to
    serve. pyitensor implements the same Chain API with no compiler,
    pybind11 or BLAS requirement of any kind, so it is always available
    and is the right answer there.

    Only the *implicit* choice moves. An explicit `itensor_version=3` on a
    machine with no extension still falls back to ED (mode.py), because a
    caller who named a version asked for that backend specifically; and
    `mode="ED"` stays exactly what it always was, a deliberate choice and
    the cross-check every test in tests/ is built on.
    """
    if available(DEFAULT_ITENSOR_VERSION):
        return DEFAULT_ITENSOR_VERSION
    return "python"


def get_backend(version=DEFAULT_ITENSOR_VERSION):
    """Return the module implementing the given DMRG backend version
    (2 = ITensor v2, 3 = ITensor v3, "python" = pyitensor's pure-Python
    backend), or None if it isn't available (never happens for "python")"""
    if version not in _backends:
        _backends[version] = _load(version)
    return _backends[version]


def _load(version):
    # numpy first: the compiled extension links the libstdc++ of the
    # compiler it was built with, and the first libstdc++ a process loads
    # is the one every later library gets. Under environment modules (on
    # Triton, scicomp-python-env with gcc) that compiler's libstdc++ is
    # older than the one of the Python environment, and loading the
    # extension first made every later scipy.optimize import fail
    # (CXXABI_1.3.15 not found), so `import dmrgpy` as the first import of
    # a script failed there. numpy brings the environment's own libstdc++,
    # which the extension also accepts.
    import numpy # noqa: F401
    try:
        if version==2:
            from .mpscpp2 import _dmrgcpp
        elif version==3:
            from .mpscpp3 import _dmrgcpp
        elif version=="python":
            from .pyitensor import chain as _dmrgcpp
        else:
            raise ValueError("Unknown DMRG backend version: "+str(version))
    except ImportError as e:
        _record_failure(version, e)
        return None
    _states[version] = dict(state="ok", error=None)
    if version in (2,3): _warn_if_stale(version)
    return _dmrgcpp


# ---------------------------------------------------------------------
# Why a backend is or is not available. An ImportError alone does not say
# whether the extension was never compiled (every pip install: nothing to
# report) or was compiled and cannot be loaded here (another interpreter,
# an environment module that is not loaded, a libstdc++ that does not
# match), and the second case used to be as silent as the first: a chain
# on the default backend then ran on "python" with no message at all.
# ---------------------------------------------------------------------

_states = {} # version -> dict(state=..., error=...), filled by _load

NOT_COMPILED = "not_compiled"   # no _dmrgcpp*.so in the backend folder
OTHER_PYTHON = "other_python"   # compiled, but for another interpreter
LOAD_FAILED = "load_failed"     # compiled for this one, and it does not load


def extension_dir(version):
    """Folder of the compiled backend `version` (2 or 3)."""
    import os
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
            "mpscpp"+str(version))


def extension_files(version):
    """(file for this interpreter or None, every _dmrgcpp*.so present)."""
    import glob, os, sysconfig
    d = extension_dir(version)
    present = sorted(glob.glob(os.path.join(d, "_dmrgcpp*.so")))
    suffix = sysconfig.get_config_var("EXT_SUFFIX") or ".so"
    mine = os.path.join(d, "_dmrgcpp"+suffix)
    return (mine if mine in present else None), present


def _record_failure(version, error):
    """Classify a failed import of backend `version` and warn, once per
    process, when an extension is there and cannot be used."""
    if version not in (2,3):
        _states[version] = dict(state=LOAD_FAILED, error=str(error))
        return
    mine, present = extension_files(version)
    if mine is not None: state = LOAD_FAILED
    elif present: state = OTHER_PYTHON
    else: state = NOT_COMPILED
    _states[version] = dict(state=state, error=str(error))
    # explain() below asks get_backend(); the answer has to be there
    # already, or it would start this very import again and recurse
    _backends[version] = None
    if state != NOT_COMPILED:
        import warnings
        msg = explain(version)
        if version == DEFAULT_ITENSOR_VERSION:
            msg += (" Until then a chain with no itensor_version runs on "
                    "the slower pure-Python backend, and one that asks for "
                    "itensor_version=%s runs on exact diagonalization."
                    %version)
        warnings.warn(msg, RuntimeWarning, stacklevel=4)


def _warn_if_stale(version):
    """Warn when the extension was built from other sources than the ones
    in the backend folder now: a C++ change is not tested until the
    extension has been rebuilt."""
    try:
        from . import buildinfo
        stale = buildinfo.is_stale(extension_dir(version))
    except Exception:
        return # a broken record must never stop the backend from loading
    _states[version]["stale"] = stale
    if stale:
        import warnings
        warnings.warn("The compiled ITensor v%s extension is older than its "
                "C++ sources in %s. Results come from the old code until it "
                "is rebuilt: python install.py --extension-only"
                %(version, extension_dir(version)),
                RuntimeWarning, stacklevel=4)


def status(version=DEFAULT_ITENSOR_VERSION):
    """What is known about backend `version` in this process, as a dict:
    `available`, `state` ("ok", "not_compiled", "other_python" or
    "load_failed"), `error` (the import error, if any), `file` (the loaded
    extension), `stale` (True when the extension is older than its
    sources, None when unknown) and `build` (the record install.py left,
    or None)."""
    backend = get_backend(version)
    st = dict(_states.get(version) or dict(state=NOT_COMPILED, error=None))
    out = dict(version=version, available=backend is not None,
            state=st["state"], error=st.get("error"),
            file=getattr(backend, "__file__", None),
            stale=st.get("stale"), build=None)
    if version in (2,3):
        from . import buildinfo
        out["build"] = buildinfo.read(extension_dir(version))
    return out


def explain(version=DEFAULT_ITENSOR_VERSION):
    """One or two sentences saying why backend `version` cannot be used
    here and what to do about it (or that it can)."""
    if version not in (2,3):
        if get_backend(version) is not None:
            return "The %r backend is available."%(version,)
        return "The %r backend could not be imported: %s"%(version,
                (_states.get(version) or {}).get("error"))
    name = "ITensor v%s extension"%version
    if get_backend(version) is not None:
        return "The "+name+" is available."
    st = _states.get(version) or {}
    state = st.get("state")
    if state == OTHER_PYTHON:
        import os, sys
        _, present = extension_files(version)
        return ("The "+name+" is compiled for another Python ("
                +", ".join(os.path.basename(f) for f in present)
                +"), not for the one running now ("+sys.executable
                +"). Run python install.py of the repository with this "
                "interpreter, or use the one it was built with.")
    if state == LOAD_FAILED:
        msg = ("The "+name+" is compiled but could not be loaded: "
                +str(st.get("error"))+".")
        from . import buildinfo
        build = buildinfo.read(extension_dir(version)) or {}
        import os, sys
        loaded = os.environ.get("LOADEDMODULES", "").split(":")
        missing = [m for m in build.get("modules") or [] if m not in loaded]
        if missing:
            msg += (" It was built with these environment modules, which "
                    "are not loaded now: "+", ".join(missing)+".")
        conda = build.get("conda_prefix")
        here = (os.environ.get("CONDA_PREFIX"), sys.prefix)
        if conda and not any(p and os.path.realpath(p)
                == os.path.realpath(conda) for p in here):
            msg += (" It was built in the conda environment "+conda
                    +", which is not the active one.")
        return msg+" python -m dmrgpy.doctor shows the details."
    return ("The "+name+" is not compiled. It is built by python "
            "install.py in a clone of the dmrgpy repository (a pip install "
            "has no C++ part); itensor_version=\"python\" runs DMRG "
            "without a compiler.")


def available(version=DEFAULT_ITENSOR_VERSION):
    """Whether the in-process extension for this backend version can be used"""
    return get_backend(version) is not None
