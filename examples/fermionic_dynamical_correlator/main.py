import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain

n = 10 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain

####### Input matrices #######

# Array with the hoppings and with hubbard couplings
# These are the matrices that you have to modify
hopping = np.zeros((n,n))
hubbard = np.zeros((n,n))
for i in range(n-1):  hopping[i,i+1] = 1. ; hopping[i+1,i] = 1.
for i in range(n): U = 6.0 ; hubbard[i,i] = U/2. ; hopping[i,i] = -U

# The implemented Hamiltonian is
# H = \sum_ij hopping[i,j] c^dagger_i c_j + hubbard[i,j] n_i n_j
# with n_i = c^\dagger_{i,up} c_{i,up} + c^\dagger_{i,dn} c_{i,dn}

# the previous matrices are for a half filled Hubbard chain

##############################


# Setup the Many Body Hamiltonian
fc.set_hoppings(lambda i,j: hopping[i,j]) # set the hoppings
fc.set_hubbard(lambda i,j: hubbard[i,j]) # set the hubbard constants



# Compute the dynamical correlator defined by
# <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>

i = 0 # first index of the dynamical correlator
j = 0 # second index of the dynamical correlator
delta = 0.1 # energy resolution (approximate)
fc.kpmmaxm = 20 # maximum bond dimension in KPM

# The result will be written in a file called DYNAMICAL_CORRELATOR.OUT


# compute the dynamical correlator using KPM DMRG
(x,y) = fc.get_dynamical_correlator(i=i,j=j,delta=delta,name="cdc")

#import matplotlib.pyplot as plt
# plot the result
#plt.plot(x,y.real,marker="o")
#plt.show()






