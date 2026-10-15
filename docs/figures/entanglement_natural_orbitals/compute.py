# Natural-orbital occupations n_alpha and correlation entropy S of the
# half-filled open Hubbard chain as U grows, with ED of a 4-site chain as the
# check of S(U).
import numpy as np
from dmrgpy import fermionchain

n, n_small = 12, 4
Us = [0, 1, 2, 4, 8, 16]
U_fine = np.linspace(0, 16, 33)


def hubbard(nsite, U):
    fc = fermionchain.Spinful_Fermionic_Chain(nsite)
    h = -sum(cd[i]*c[i+1] + cd[i+1]*c[i] for i in range(nsite-1)
             for c, cd in ((fc.Cup, fc.Cdagup), (fc.Cdn, fc.Cdagdn)))
    if U != 0:
        h = h + U*sum(fc.Nup[i]*fc.Ndn[i] for i in range(nsite))
    fc.set_hamiltonian(h)
    fc.set_conserved_sector(Nf=nsite)  # half filling
    fc.maxm, fc.nsweeps = 100, 24
    return fc


occ, s, s_small_dmrg = [], [], []
for U in Us:
    fc = hubbard(n, U)
    occ.append(np.sort(fc.get_correlation_eigenvalues().real)[::-1])
    s.append(fc.get_correlation_entropy())
    s_small_dmrg.append(hubbard(n_small, U).get_correlation_entropy())
s_fine_ed = [hubbard(n_small, U).get_correlation_entropy(mode="ED")
             for U in U_fine]

np.savez("data.npz", n=n, n_small=n_small, U=Us, occ=occ, s=s,
         s_small_dmrg=s_small_dmrg, U_fine=U_fine, s_fine_ed=s_fine_ed)
