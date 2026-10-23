# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Fermion parity of the ground state of a superconducting chain with a
# magnetic impurity. The parity is even or odd depending on whether the
# number of electrons is even or odd. Changing the chemical potential
# makes the ground state switch parity, a signature of a quantum phase
# transition that can be detected with Yu-Shiba-Rusinov states in STM.
import numpy as np
from dmrgpy import fermionchain

L = 4 # number of sites

def number_and_parity(mu):
    fc = fermionchain.Spinful_Fermionic_Chain(L)
    h = 0
    for i in range(L-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(L): # s-wave superconductivity
        h = h + 0.2*fc.Cup[i]*fc.Cdn[i]
    h = h + h.get_dagger()
    for i in range(L): # chemical potential
        h = h + mu*(fc.Cdagup[i]*fc.Cup[i] + fc.Cdagdn[i]*fc.Cdn[i])
    for i in range(L-1): # Hubbard interaction
        h = h + 0.5*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + 2.0*fc.Sz[0] # magnetic impurity
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    wf = fc.get_gs(mode="DMRG")
    return wf.dot(sum(fc.N)*wf).real, wf.get_fermionic_parity()

mus = np.linspace(-1., 0., 20)
res = np.array([number_and_parity(mu) for mu in mus])

import matplotlib.pyplot as plt
plt.plot(mus, res[:,0], "o", label="number of electrons")
plt.plot(mus, res[:,1], "x", label="parity")
plt.xlabel("chemical potential")
plt.legend()
plt.show()
