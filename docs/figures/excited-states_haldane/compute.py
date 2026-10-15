# Edge states and bulk gap of the open spin-1 Heisenberg chain, from the
# ground states of the sectors Sz_tot = 0, 1 and 2 (the two spin-1/2 edge
# states form a singlet and a triplet, the gap is a magnon on the triplet).
import numpy as np
from dmrgpy import spinchain

ns = [16, 24, 32, 48, 64, 96, 128]
n_profile = 48  # length of the <Sz_i> profile


def sector(n, sz):
    sc = spinchain.Spin_Chain(["S=1"]*n)
    sc.set_hamiltonian(sum(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1]
                           + sc.Sz[i]*sc.Sz[i+1] for i in range(n-1)))
    sc.maxm, sc.nsweeps, sc.cutoff = 150, 30, 1e-12
    sc.set_conserved_sector(Sz=2*sz)  # Sz is given in units of 1/2
    return sc, sc.gs_energy()


e0, e1, e2 = [], [], []
for n in ns:
    e0.append(sector(n, 0)[1])
    sc1, e = sector(n, 1)
    e1.append(e)
    e2.append(sector(n, 2)[1])
    if n == n_profile:
        profile = [sc1.vev(sc1.Sz[i]).real for i in range(n)]

np.savez("data.npz", ns=ns, e0=e0, e1=e1, e2=e2, profile=profile)
