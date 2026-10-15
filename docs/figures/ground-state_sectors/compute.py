# Ground-state energies E0(M) of the open spin-1/2 XX and Heisenberg chains in
# every sector of total Sz = M; in a field h the magnetization steps from M to
# M+1 at h_M = E0(M+1) - E0(M).
import numpy as np
from dmrgpy import spinchain

n = 32
Ms = np.arange(0, n//2 + 1)


def sector_energy(heisenberg, M):
    sc = spinchain.Spin_Chain(["S=1/2"]*n)
    h = sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] for i in range(n-1))
    if heisenberg:
        h = h + sum(sc.Sz[i]*sc.Sz[i+1] for i in range(n-1))
    sc.set_hamiltonian(h)
    sc.maxm, sc.nsweeps, sc.cutoff = 60, 12, 1e-12
    sc.set_conserved_sector(Sz=2*int(M))  # Sz in units of 1/2
    return sc.gs_energy()


e_xx = np.array([sector_energy(False, M) for M in Ms])
e_heis = np.array([sector_energy(True, M) for M in Ms])

# XX chain as free fermions with hopping J/2: E0(M) is the sum of the lowest
# n/2+M levels J cos(pi k/(n+1))
eps = np.sort(np.cos(np.pi*np.arange(1, n+1)/(n+1)))
e_xx_ff = np.array([np.sum(eps[:n//2 + M]) for M in Ms])

np.savez("data.npz", n=n, M=Ms, e_xx=e_xx, e_heis=e_heis, e_xx_ff=e_xx_ff)
