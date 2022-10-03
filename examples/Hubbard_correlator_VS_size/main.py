# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np

def get(n):
    from dmrgpy import fermionchain
    U = 2.0 # value of U
    fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdagup[i]*fc.Cup[i+1]
        h = h + fc.Cdagdn[i]*fc.Cdn[i+1]
    for i in range(n): # Hubbard
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    h = h + h.get_dagger()
    ##############################
    # Setup the Many Body Hamiltonian
    fc.maxm = 40
    fc.nsweeps = 40
    fc.set_hamiltonian(h) # set the hoppings
    O = fc.Cdagup[n//2]*fc.Cup[n//2+1] + fc.Cdagdn[n//2]*fc.Cdn[n//2+1]
    print("Doing",n)
    return np.abs(fc.vev(O,mode="DMRG"))



Ls = range(4,10,1)
cs = [get(L) for L in Ls]

import matplotlib.pyplot as plt
plt.plot(Ls,cs)
plt.ylim([0.,1.0])
plt.xlabel("L")
plt.ylabel("Correlator")
plt.show()



