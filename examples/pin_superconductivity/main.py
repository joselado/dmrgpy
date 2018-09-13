import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain

n = 9 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain

####### Input matrices #######

# Array with the hoppings and with hubbard couplings
# These are the matrices that you have to modify
hopping = np.zeros((n,n))
hubbard = np.zeros((n,n))
pairing = np.zeros((n,n))
for i in range(n-1):  hopping[i,i+1] = 1. ; hopping[i+1,i] = 1.
for i in range(n-3):  hopping[i,i+3] = 1./3. ; hopping[i+3,i] = 1./3.
for i in range(n): U = -10. ; hubbard[i,i] = U/2. ; hopping[i,i] = -U
for i in range(n): pairing[i,i] = 3.

# The implemented Hamiltonian is
# H = \sum_ij hopping[i,j] c^dagger_i c_j + hubbard[i,j] n_i n_j
# with n_i = c^\dagger_{i,up} c_{i,up} + c^\dagger_{i,dn} c_{i,dn}

# the previous matrices are for a half filled Hubbard chain

##############################


# Setup the Many Body Hamiltonian
fc.set_hoppings(lambda i,j: hopping[i,j]) # set the hoppings
fc.set_hubbard(lambda i,j: hubbard[i,j]) # set the hubbard constants
fc.set_pairing(lambda i,j: pairing[i,j]) # set the hubbard constants
#fc.set_fields(lambda i: [0.,0.,0.2]) # set the hubbard constants

#fc.nsweeps = 7


# Compute the dynamical correlator defined by
# <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>

delta = 0.1 # energy resolution (approximate)
fc.maxm = 30 # maximum bond dimension
fc.nsweeps = 7

wf0 = fc.get_gs() # ground state wavefunction

fc.set_pairing(lambda i,j: 0.*pairing[i,j]) # set to zero

#fc.get_gs(wf0=wf0) # use an initial guess

e = fc.gs_energy(wf0=wf0)
print("Energy",e)

# The result will be written in a file called DYNAMICAL_CORRELATOR.OUT
ds = fc.get_delta()
print(fc.get_delta())



import matplotlib.pyplot as plt

plt.plot(range(len(ds)),ds)
plt.show()

