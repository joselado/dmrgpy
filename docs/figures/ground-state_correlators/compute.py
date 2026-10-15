# <Sz_i Sz_{i+r}> in the open spin-1/2 and spin-1 Heisenberg chains, from the
# site i0 = n/2-16 outwards, so both sites stay 16 or more sites from the ends.
import numpy as np
from dmrgpy import spinchain

n = 64
rs = np.arange(1, 33)
i0 = n//2 - 16


def correlator(sites, maxm):
    sc = spinchain.Spin_Chain(sites)
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                           + sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)))
    sc.maxm, sc.nsweeps, sc.cutoff = maxm, 20, 1e-12
    sc.set_conserved_sector(Sz=0)
    return np.array([sc.vev(sc.Sz[i0]*sc.Sz[i0+r]).real for r in rs])


czz_half = correlator(["S=1/2"]*n, 120)
# A spin 1/2 at each end binds the free edge spins of the spin-1 chain into
# a singlet, so the ground state is unique and the bulk is unchanged.
czz_one = correlator(["S=1/2"] + ["S=1"]*(n-2) + ["S=1/2"], 100)

np.savez("data.npz", r=rs, czz_half=czz_half, czz_one=czz_one)
