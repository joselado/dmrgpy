import sys; sys.path.insert(0, "/u/40/ladovj1/data/Documents/programs/dmrgpy/src")
import numpy as np
from dmrgpy import spinchain, bosonchain
from dmrgpy.lindblad import Lindblad_Chain

def xxz(n=4, maxm=64):
    sc = spinchain.Spin_Chain([2]*n, itensor_version="python"); sc.maxm = maxm
    H = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + 0.5*sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    sc.set_hamiltonian(H)
    sp = lambda i: sc.Sx[i] + 1j*sc.Sy[i]; sm = lambda i: sc.Sx[i] - 1j*sc.Sy[i]
    jumps = [0.6*sp(0), 0.6*sm(n-1)] + [0.4*sc.Sz[i] for i in range(n)]
    return sc, Lindblad_Chain(sc, jumps)

def bh(L=3, N=2, maxm=64):
    bc = bosonchain.Bosonic_Chain(L, maxnb=[N+1]*L, itensor_version="python")
    h = sum(-(bc.Adag[j+1]*bc.A[j] + bc.Adag[j]*bc.A[j+1]) for j in range(L-1))
    h = h + 0.5*sum(bc.N[j]*(bc.N[j]-1) for j in range(L))
    bc.set_hamiltonian(h); bc.maxm = maxm
    jumps = [(bc.Adag[j+1]+bc.Adag[j])*(bc.A[j+1]-bc.A[j]) for j in range(L-1)]
    return bc, Lindblad_Chain(bc, jumps)

def bh_start(L, N):
    psi = [1]*N + [0]*(L-N); phi = [0]*(L-N) + [1]*N
    return psi, phi

def errors(lam, lamE, k=6):
    """distance of each of the first k ED eigenvalues to the nearest Krylov one"""
    return np.array([np.min(np.abs(lam - x)) for x in lamE[:k]])
