import warnings as _warnings

from . import cppext as _cppext

if not _cppext.available(_cppext.DEFAULT_ITENSOR_VERSION):
    _warnings.warn(
        "ITensor v%s (dmrgpy's default C++ DMRG backend) is not compiled. "
        "Chains will fall back to another compiled backend (itensor_version="
        "2 or 3) or exact diagonalization. Run `python install.py` (or "
        "`python install.py --itensor-version=%s`) to compile it."
        % (_cppext.DEFAULT_ITENSOR_VERSION, _cppext.DEFAULT_ITENSOR_VERSION),
        stacklevel=2)
