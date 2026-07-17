"""pyitensor: a pure-Python implementation of the (small) subset of the
ITensor v3 API that dmrgpy's mpscpp3 backend actually uses, so DMRG/TDVP
can run without a compiled C++ extension. See mpscpp3/chain_session.h,
get_sites.h, mo_terms.h and TDVP/tdvp.h for the API surface this mirrors --
this package is not, and is not meant to be, a general ITensor port.

Phase 1 (this file's current contents): the Index/ITensor tensor core.
Later phases add SiteSet types, AutoMPO, MPS/MPO, DMRG, TDVP, and finally a
Chain facade with the same method surface as mpscpp3/bindings.cc's pybind11
class, so it can be registered as a new itensor_version in cppext.py with
no changes needed anywhere else in dmrgpy.
"""

from .index import Index, sim
from .tensor import ITensor, prime, noPrime, mapPrime, swapPrime, dag, commonIndex, delta
from .svd import svd, Spectrum
from .mpscontainer import MPS, MPO
from .mpsalgebra import inner, innerC, traceC, sum, applyMPO, nmultMPO, randomMPS
from .sweeps import Sweeps
from .autompo import HTerm, AutoMPO
from .mpobuilder import to_mpo

__all__ = [
    "Index", "sim",
    "ITensor", "prime", "noPrime", "mapPrime", "swapPrime", "dag", "commonIndex", "delta",
    "svd", "Spectrum",
    "MPS", "MPO",
    "inner", "innerC", "traceC", "sum", "applyMPO", "nmultMPO", "randomMPS",
    "Sweeps",
    "HTerm", "AutoMPO", "to_mpo",
]
