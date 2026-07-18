"""pyitensor: a pure-Python implementation of the (small) subset of the
ITensor v3 API that dmrgpy's mpscpp3 backend actually uses, so DMRG/TDVP
can run without a compiled C++ extension. See mpscpp3/chain_session.h,
get_sites.h, mo_terms.h and TDVP/tdvp.h for the API surface this mirrors --
this package is not, and is not meant to be, a general ITensor port.

Index/ITensor tensor core, site types, AutoMPO, MPS/MPO algebra, DMRG
(ground + excited states), TDVP, and a Chain facade (chain.Chain) with the
same method surface as mpscpp3/bindings.cc's pybind11 class -- registered
as itensor_version="python" in cppext.py, with no other changes needed
anywhere else in dmrgpy (mps.py's MPS/staticoperator.py's StaticOperator
only ever treat the session's return values as an opaque cpp_handle).
kernels.py is an optional, JAX-accelerated fast path for the single
hottest loop (the DMRG/TDVP effective-Hamiltonian matvec) -- never a hard
dependency; everything here works with only NumPy/SciPy installed.
"""

from .index import Index, sim
from .tensor import ITensor, prime, noPrime, mapPrime, swapPrime, dag, commonIndex, delta
from .svd import svd, Spectrum
from .mpscontainer import MPS, MPO
from .mpsalgebra import inner, innerC, traceC, sum, applyMPO, nmultMPO, randomMPS
from .sweeps import Sweeps
from .autompo import HTerm, AutoMPO
from .mpobuilder import to_mpo
from .chain import Chain
from . import kernels

__all__ = [
    "Index", "sim",
    "ITensor", "prime", "noPrime", "mapPrime", "swapPrime", "dag", "commonIndex", "delta",
    "svd", "Spectrum",
    "MPS", "MPO",
    "inner", "innerC", "traceC", "sum", "applyMPO", "nmultMPO", "randomMPS",
    "Sweeps",
    "HTerm", "AutoMPO", "to_mpo",
    "Chain",
    "kernels",
]
