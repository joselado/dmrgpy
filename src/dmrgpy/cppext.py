"""
Single import point for the optional in-process pybind11 extensions
(mpscpp2/bindings.cc built against ITensor v2, mpscpp3/bindings.cc built
against ITensor v3), which replace the file-based Python<->C++ protocol
(see the migration plan) one task at a time. Both expose the same Chain
session API (see chain_session.h in each backend) -- callers elsewhere
never need to know which ITensor version actually backs self._session.

Callers should use get_backend(version)/available(version) rather than
importing mpscpp2._dmrgcpp/mpscpp3._dmrgcpp directly, since a given
extension may not be built (pybind11 wasn't available at install time, or
"make pybind" was never run for that backend) -- in that case every caller
should fall back to ED, not fail outright (see mode.py).
"""

_backends = {} # version -> compiled _dmrgcpp module, or None if unavailable


def get_backend(version=2):
    """Return the compiled _dmrgcpp module for the given C++ DMRG backend
    version (2 = ITensor v2, 3 = ITensor v3), or None if it isn't available"""
    if version not in _backends:
        _backends[version] = _load(version)
    return _backends[version]


def _load(version):
    try:
        if version==2:
            from .mpscpp2 import _dmrgcpp
        elif version==3:
            from .mpscpp3 import _dmrgcpp
        else:
            raise ValueError("Unknown C++ DMRG backend version: "+str(version))
        return _dmrgcpp
    except ImportError:
        return None


def available(version=2):
    """Whether the in-process extension for this backend version can be used"""
    return get_backend(version) is not None
