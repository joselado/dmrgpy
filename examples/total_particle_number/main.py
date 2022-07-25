# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
L = 4 # number of sites

def get_N(mu):
    fc = fermionchain.Fermionic_Chain(L) # create the chain
    h = 0 # initialize
    for i in range(L-1): h = h + fc.Cdag[i]*fc.C[i+1] # first neighbor hopping
    h = h + h.get_dagger() # add Hermitian conjugate
    V = 1.0 # strength of many-body repulsion
    for i in range(L-1): h = h + V*fc.N[i]*fc.N[i+1] # first neighbor repulsion
    for i in range(L): h = h + mu*fc.N[i] # chemical potential
    fc.set_hamiltonian(h) # set Hamiltonian
    wf = fc.get_gs() # compute ground state
    Nop = 0 
    for i in range(L): Nop = Nop + fc.N[i] # total number operator
    nf = wf.dot(Nop*wf).real # expectation value of the number operator
    print("Particle number for mu =",mu," N = ",nf)
    return nf

mus = np.linspace(-5.,5.,30) # chemical potentials
ns = [get_N(mu) for mu in mus] # particle number operator for different mu

import matplotlib.pyplot as plt

plt.plot(mus,ns,marker="o")
plt.xlabel("mu")
plt.ylabel("# of electrons")
plt.show()










