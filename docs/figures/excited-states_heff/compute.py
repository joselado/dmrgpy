# Exchange J of the half-filled Hubbard dimer, fitted by get_heff on the
# four lowest states, with and without projecting them on the singly
# occupied states.
import numpy as np
from dmrgpy import effectivehamiltonian, fermionchain

t = 1.0
lam = 20.0  # lam*(N-2)^2 keeps the four lowest states at half filling
Us = np.geomspace(1, 20, 25)

jfit, jfit_noproj = [], []
for U in Us:
    fc = fermionchain.Spinful_Fermionic_Chain(2)
    h = t*fc.Cdagup[0]*fc.Cup[1] + t*fc.Cdagdn[0]*fc.Cdn[1]
    h = h + h.get_dagger()
    h = h + sum(U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) for i in range(2))
    ntot = fc.Ntot[0] + fc.Ntot[1]
    fc.set_hamiltonian(h + lam*(ntot-2)*(ntot-2))
    ops = [fc.Sx[0]*fc.Sx[1], fc.Sy[0]*fc.Sy[1], fc.Sz[0]*fc.Sz[1]]
    proj = effectivehamiltonian.single_occupancy_projector(fc)
    jfit.append(np.real(fc.get_heff(operators=ops, n=4, mode="ED",
                                    projector=proj)))
    jfit_noproj.append(np.real(fc.get_heff(operators=ops, n=4, mode="ED")))

np.savez("data.npz", t=t, U=Us, jfit=jfit, jfit_noproj=jfit_noproj)
