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
    try:
        if version==2:
            from .mpscpp2 import _dmrgcpp
        elif version==3:
            from .mpscpp3 import _dmrgcpp
        elif version=="python":
            from .pyitensor import chain as _dmrgcpp
        else:
            raise ValueError("Unknown DMRG backend version: "+str(version))
        return _dmrgcpp
    except ImportError:
        return None


def available(version=DEFAULT_ITENSOR_VERSION):
    """Whether the in-process extension for this backend version can be used"""
    return get_backend(version) is not None
