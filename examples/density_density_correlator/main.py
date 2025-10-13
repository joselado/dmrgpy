# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain


def get_correlator(U=0.):
    L = 6 # number of sites
    fc = fermionchain.Spinful_Fermionic_Chain(L) # create the fermion chain object
    H = 0 # initialize Hamiltonian
    for i in range(L-1): # nearest neighbor hopping
        H = H + fc.Cdagup[i]*fc.Cup[i+1] # add hopping
        H = H + fc.Cdagdn[i]*fc.Cdn[i+1] # add hopping
    H = H + H.get_dagger() # Hermitian conjugate

    for i in range(L): # onsite interaction
        H = H + U*(fc.Nup[i]-0.5)*(fc.Ndn[i]-0.5) # add interaction

    fc.set_hamiltonian(H) # set the Hamiltonian
    wf0 = fc.get_gs(mode="ED") # get ground state
    return [wf0.dot((fc.Cdagup[0]*fc.Cup[i])*wf0).real for i in range(L)]

plt.figure(figsize=(8,4))

Us = [0,2,10] # interactions
for U in Us:
    cij = get_correlator(U=U)
    plt.plot(range(len(cij)),cij,label="$U =$"+str(U)+"t",marker="o",markersize=10)

plt.legend(ncol=len(Us),fontsize=15)
plt.xlabel("Distance") ; plt.ylabel("$\\langle c^{\\dagger}_{0\\uparrow} c_{N\\uparrow} \\rangle $") #; plt.ylim([0.,2.0]) ; plt.xlim([2,max(Ls)])
plt.show()
