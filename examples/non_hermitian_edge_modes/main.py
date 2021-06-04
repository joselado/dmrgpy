# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain,mpsalgebra
import numpy as np

# this is a minimal example of an interacting non-Hermitian system
# with edge zero modes

def get(V=0.0,c=1.0):
    """Return non-Hermitian eigenvalues
    V in the interaction
    c is the coupling between first and last sites"""
    n = 2 # number of tetramers
    ns = 4*n # total number of sites
    mh = np.zeros((ns,ns),dtype=np.complex) # TB matrix
    for i in range(ns-1):
        mh[i,i+1] = 1.0
        mh[i+1,i] = 1.0
    
    mh[ns-1,0] = c
    mh[0,ns-1] = c
    
    eta = 1.0
    for i in range(n):
        mh[4*i,4*i] = 1j*eta
        mh[4*i+1,4*i+1] = -1j*eta
        mh[4*i+2,4*i+2] = -1j*eta
        mh[4*i+3,4*i+3] = 1j*eta
    
    
    
    fc = fermionchain.Fermionic_Chain(ns) # create the fermion chain
    h = 0 # initialize Hamiltonian
    # put all the hoppings
    for i in range(ns):
        for j in range(ns):
            h = h + mh[i,j]*fc.Cdag[i]*fc.C[j]
    for i in range(ns-1): h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    
    
    fc.set_hamiltonian(h)
    nex = 10 # lowest states
    # the following methods target the states with minimum Re(E)
    # use this if you wanted to use MPS
#    es,wfs = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,verbose=1,n=nex,maxit=7)
    # this method is just for ED
    es,wfs = fc.get_excited_states(mode="ED",n=nex) # get excited states
    return es - es[0] # return energies

# open file
f = open("EIGS.OUT","w")

# compute for different interactions
for V in np.linspace(0.0,1.0,20):
    esed = get(c=0.0,V=V)
    print("Computing ",V)
    for e in esed:
        f.write(str(e.real)+"  ")
        f.write(str(e.imag)+"  ")
        f.write(str(V)+"\n")
    f.flush()
f.close()

