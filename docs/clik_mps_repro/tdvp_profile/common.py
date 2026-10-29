import sys, time
sys.path.insert(0, "/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import bosonchain
from dmrgpy.lindblad import Lindblad_Chain
from dmrgpy.pyitensor import mpsalgebra, tdvp

def lindblad_case(L=5, N=4, maxm=64):
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
    h = sum(-(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1)) + 0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L))
    bc.set_hamiltonian(h); bc.maxm = maxm
    jumps = [(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
    lc = Lindblad_Chain(bc, jumps)
    S = lc.doubled._session; G = S._mpo(lc.L.to_terms())
    psi = mpsalgebra.randomMPS(S.sites, 10, rng=np.random.default_rng(0))
    return S, G, psi

def heis_case(n=20, maxm=64):
    from dmrgpy import spinchain
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python")
    h = 0
    for i in range(n-1):
        h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
    sc.set_hamiltonian(h); sc.maxm = maxm
    sc.setup_python() if hasattr(sc, "setup_python") else None
    return sc
