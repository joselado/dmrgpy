import warnings as _warnings

from . import cppext as _cppext

if not _cppext.available(_cppext.DEFAULT_ITENSOR_VERSION):
    # Two very different audiences hit this notice, so it has to serve both.
    # From a git checkout, the fix is to compile the extension with
    # install.py. From a PyPI install there *is* no install.py -- the wheel
    # deliberately ships no C++ at all (see pyproject.toml) -- and nothing
    # needs fixing at all: chains built without an explicit itensor_version
    # now run on the pure-Python DMRG backend (cppext.default_backend()),
    # which is a real MPS solver, not the ED fallback this used to warn
    # about. Mentioning only install.py would send pip users chasing a file
    # they don't have; warning about ED would now be simply untrue.
    _warnings.warn(
        "ITensor v%s (dmrgpy's default C++ DMRG backend) is not compiled, so "
        "chains default to itensor_version=\"python\", the pure-Python DMRG "
        "backend. That is a real MPS solver and needs no compiler, but it is "
        "substantially slower than compiled ITensor; to build the C++ backend, "
        "run `python install.py --itensor-version=%s` from a clone of the "
        "dmrgpy repository. Passing itensor_version=%s explicitly on this "
        "machine falls back to exact diagonalization instead, which does not "
        "scale past small systems."
        % (_cppext.DEFAULT_ITENSOR_VERSION, _cppext.DEFAULT_ITENSOR_VERSION,
           _cppext.DEFAULT_ITENSOR_VERSION),
        stacklevel=2)
