# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

def get_den_par(mu=0.,U=0.5,Bz=1.0,delta=0.5):
    """Return parity and number of electrons"""
    L = 4 # number of fermionic sites
    fc = fermionchain.Spinful_Fermionic_Chain(L) # create the chain
    h = 0
    for i in range(L-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1]
        h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(L): # hopping
        h = h + mu*fc.Cdagup[i]*fc.Cup[i]
        h = h + mu*fc.Cdagdn[i]*fc.Cdn[i]
    for i in range(L): # s-wave superconductivity
        h = h + mu*fc.Cup[i]*fc.Cdn[i]
    for i in range(L-1): # Hubbard
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + Bz*fc.Sz[0] # add magnetic impurity
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    wf = fc.get_gs(mode="DMRG")
    from dmrgpy.fermionicparity import explicit_parity
    p = explicit_parity(wf) # parity of the state
    n = wf.dot(sum(fc.N)*wf) # total number of electrons 
    return n,p


mus = np.linspace(-1.,0.,20)
nps = np.array([get_den_par(mu=mu) for mu in mus])

ns = nps[:,0] # densities
ps = nps[:,1] # parities

import matplotlib.pyplot as plt

plt.subplot(1,2,1)
plt.scatter(mus,ns,c="blue")
plt.xlabel("Chemical potential")
plt.ylabel("Particle number")

plt.subplot(1,2,2)
plt.scatter(mus,ps,c="red")
plt.xlabel("Chemical potential")
plt.ylabel("Particle parity")

plt.tight_layout()

plt.show()
