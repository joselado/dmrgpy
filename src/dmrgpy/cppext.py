"""
Single import point for the optional in-process pybind11 extension
(mpscpp2/bindings.cc), which is being introduced to replace the file-based
Python<->C++ protocol (see the migration plan) one task at a time.

Callers should use get_backend()/available() rather than importing
mpscpp2._dmrgcpp directly, since the extension may not be built (pybind11
wasn't available at install time, or "make pybind" was never run) -- in
that case every caller should fall back to the existing subprocess-based
backend, not fail outright.
"""

_backend = None
_tried = False


def get_backend():
    """Return the compiled _dmrgcpp module, or None if it isn't available"""
    global _backend,_tried
    if not _tried:
        _tried = True
        try:
            from .mpscpp2 import _dmrgcpp
            _backend = _dmrgcpp
        except ImportError:
            _backend = None
    return _backend


def available():
    """Whether the in-process extension can be used"""
    return get_backend() is not None
