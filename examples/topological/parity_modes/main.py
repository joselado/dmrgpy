# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fermionic parity of the many-body ground state, which is +1 or -1 for an
# even or an odd number of electrons. We take a spinful chain with
# hopping, pairing, Hubbard interaction and a magnetic impurity.
# The parity changes when the chemical potential mu is changed.
import numpy as np
from dmrgpy import fermionchain

L = 4 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(L)
fc.maxm, fc.nsweeps = 20, 10 # bond dimension and sweeps
U, Bz = 0.5, 1.0 # Hubbard and impurity field

def parity(mu):
    h = 0
    for i in range(L-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(L): # chemical potential and pairing
        h = h + mu*(fc.Cdagup[i]*fc.Cup[i] + fc.Cdagdn[i]*fc.Cdn[i])
        h = h + mu*fc.Cup[i]*fc.Cdn[i]
    for i in range(L-1): # Hubbard
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + Bz*fc.Sz[0] # magnetic impurity
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    wf = fc.get_gs(mode="DMRG")
    return wf.dot(sum(fc.N)*wf).real, wf.get_fermionic_parity(fpmode="full")

mus = np.linspace(-1., 0., 20) # chemical potentials
data = np.array([parity(mu) for mu in mus])

import matplotlib.pyplot as plt
plt.plot(mus, data[:, 0], "o-", label="electrons")
plt.plot(mus, data[:, 1], "s-", label="parity")
plt.xlabel("chemical potential")
plt.legend()
plt.show()
