# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

##################
# this example shows how to compute the parity of a many-body fermionic
# state, meaning whether if it made of an odd or even number of electrons
##################


def get_par(mu=0.,U=0.5,Bz=1.0,delta=0.5,fpmode="full",L=4):
    """Return parity and number of electrons"""
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
    wf = fc.get_gs(mode="DMRG") # compute ground state
    return  wf.get_fermionic_parity(fpmode=fpmode) # parity of the state

p = get_par(L=20,fpmode="iterative")
print("Parity",p)

