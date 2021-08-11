# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain

n = 1 # number of tetramers
ns = 4*n # number of unit cells

# create the tight binding matrix
mh = np.zeros((ns,ns),dtype=np.complex) # TB matrix
for i in range(ns-1):
    mh[i,i+1] = 1.0
    mh[i+1,i] = 1.0

eta = 1.0 # imaginary part
for i in range(n): 
    mh[4*i,4*i] = 1j*eta
    mh[4*i+1,4*i+1] = -1j*eta
    mh[4*i+2,4*i+2] = -1j*eta
    mh[4*i+3,4*i+3] = 1j*eta

# create the many-body object
fc = fermionchain.Fermionic_Chain(ns) # create the fermion chain

# add the single particle term
h = 0 # initialize Hamiltonian
for i in range(ns):
    for j in range(ns):
        h = h + mh[i,j]*fc.Cdag[i]*fc.C[j]

# now add the density-density interactions
V = 1.0
for i in range(ns-1): h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)


# Compute the states with lowest Re(E) with exact diagonalization
# Note that this will only work for small systems
fc.set_hamiltonian(h)
esed,wfsed = fc.get_excited_states(mode="ED",
                         n=3 # number of excited states
                         )
# Compute the states with lowest Re(E) with MPS non-hermitian arnoldi
from dmrgpy import mpsalgebra
fc.maxm = 20 # increase it if you need more accuracy
es,wfs = mpsalgebra.lowest_energy_non_hermitian_arnoldi(fc,h,
         verbose=1, # print info on how it is progressing
         n=3, # number of excited states
         maxit=7 # this sets a cutoff in the number of itnerrations
         )


# now define the function to compute the correlators
def proj(wf1,wf2,i):
    """Compute |<0|c^\dagger_n|1>|^2 + |<0|c_n|1>|^2"""
    return np.abs(wf1.dot(fc.Cdag[i]*wf2))**2 + np.abs(wf1.dot(fc.C[i]*wf2))**2


def get_proj(wfs):
    """Compute the projection with the lowest energy state"""
    ps = [proj(wfs[0],wfs[1],i) for i in range(ns)]
    return ps

# compute the correlators with both methods
pss = get_proj(wfs) # get the projection for the MPS wavefunctions
pssed = get_proj(wfsed) # get the projection for the ED wavefunctions

import matplotlib.pyplot as plt
plt.scatter(range(len(pss)),pss,c="blue",s=80,label="MPS")
plt.scatter(range(len(pss)),pssed,c="red",s=40,label="ED")
plt.legend()
plt.xlabel("Site")
plt.ylabel("Correlator")
plt.show()






