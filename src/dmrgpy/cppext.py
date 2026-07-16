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

import os


def _preload_system_libstdcxx():
    """
    The compiled extension is built with the system g++ and needs a
    libstdc++ new enough to provide its GLIBCXX_* symbols. On systems where
    Python's own environment bundles an older libstdc++ (e.g. conda), the
    *first* C++ shared library loaded into the process determines which
    libstdc++ copy every C++ extension in that process ends up sharing --
    having two different libstdc++ builds loaded at once crashes (a hard
    segfault, not a clean missing-symbol error): this was discovered
    empirically when scipy.linalg/scipy.sparse (both already imported by
    dmrgpy itself) loaded their own libstdc++ before this extension got a
    chance to. Preloading the newer system libstdc++ globally, as early as
    possible -- this function is called at the top of dmrgpy/__init__.py,
    before any submodule that touches scipy -- makes sure only one copy
    is ever loaded. Harmless no-op if none of these paths exist (e.g. on a
    non-Debian/Ubuntu system) or if the extension is never used.
    """
    candidates = [
        "/usr/lib/x86_64-linux-gnu/libstdc++.so.6",
        "/usr/lib64/libstdc++.so.6",
        "/usr/lib/libstdc++.so.6",
    ]
    import ctypes
    for path in candidates:
        if os.path.isfile(path):
            try:
                ctypes.CDLL(path,mode=ctypes.RTLD_GLOBAL)
                return
            except OSError:
                continue


_preload_system_libstdcxx()


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
