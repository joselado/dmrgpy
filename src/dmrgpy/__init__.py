import warnings as _warnings

from . import cppext as _cppext

if not _cppext.available(_cppext.DEFAULT_ITENSOR_VERSION):
    # Two very different audiences hit this warning, so it has to serve both.
    # From a git checkout, the fix is to compile the extension with
    # install.py. From a PyPI install there *is* no install.py -- the wheel
    # deliberately ships no C++ at all (see pyproject.toml) -- and the real
    # answer is itensor_version="python", the pure-Python DMRG backend, which
    # always works and is far faster than the ED fallback on larger chains.
    # Mentioning only install.py would send pip users chasing a file they
    # don't have.
    _warnings.warn(
        "ITensor v%s (dmrgpy's default C++ DMRG backend) is not compiled, so "
        "chains fall back to exact diagonalization, which does not scale past "
        "small systems. For real DMRG without compiling anything, use the "
        "pure-Python backend: Spin_Chain(..., itensor_version=\"python\") or "
        "chain.setup_python(). To build the (faster) C++ backend instead, run "
        "`python install.py --itensor-version=%s` from a clone of the dmrgpy "
        "repository."
        % (_cppext.DEFAULT_ITENSOR_VERSION, _cppext.DEFAULT_ITENSOR_VERSION),
        stacklevel=2)
