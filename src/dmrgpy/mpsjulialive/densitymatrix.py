import numpy as np


def reduced_dm(self,wf,i):
    """Single-site reduced density matrix on the Julia backend, via
    mpsjulialive/densitymatrix.jl's reduced_dm (a native Julia port of
    pyitensor.chain.Chain.reduced_dm, itself a port of
    mpscpp3/chain_session.h's Chain::reduced_dm). i is 0-indexed, matching
    densitymatrix.py::reduced_dm's own convention."""
    from .juliasession import Main as Mainjl
    rho = Mainjl.reduced_dm(wf.jlmps,i+1)
    return np.array(rho)
