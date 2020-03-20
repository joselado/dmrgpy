# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import fermionchain
n = 30 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Chain(n) # create the chain
####### Input matrices #######
# Array with the hoppings and with hubbard couplings
# These are the matrices that you have to modify
hopping = np.zeros((n,n))
hubbard = np.zeros((n,n))
for i in range(n-1):  hopping[i,i+1] = 1. ; hopping[i+1,i] = 1.
for i in range(n-3):  hopping[i,i+3] = 1./3. ; hopping[i+3,i] = 1./3.
for i in range(n): U = -0.5 ; hubbard[i,i] = U/2. ; hopping[i,i] = -U
# The implemented Hamiltonian is
# H = \sum_ij hopping[i,j] c^dagger_i c_j + hubbard[i,j] n_i n_j
# with n_i = c^\dagger_{i,up} c_{i,up} + c^\dagger_{i,dn} c_{i,dn}
# the previous matrices are for a half filled Hubbard chain
##############################
# Setup the Many Body Hamiltonian
fc.set_hoppings(lambda i,j: hopping[i,j]) # set the hoppings
fc.set_hubbard(lambda i,j: hubbard[i,j]) # set the hubbard constants
#fc.set_fields(lambda i: [0.,0.,0.2]) # set the hubbard constants
#fc.nsweeps = 7
# Compute the dynamical correlator defined by
# <0|c_i^dagger \delta(H-E_0-\omega) c_j |0>
i = 0 # first index of the dynamical correlator
j = 0 # second index of the dynamical correlator
delta = 0.1 # energy resolution (approximate)
fc.maxm = 60 # maximum bond dimension
fc.nsweeps = 4
# The result will be written in a file called DYNAMICAL_CORRELATOR.OUT
print(np.mean(fc.get_delta()))
exit()
import matplotlib.pyplot as plt
for i in range(1):
# compute the dynamical correlator using KPM DMRG
  (x,y) = fc.get_dynamical_correlator(i=i,j=j,delta=delta,name="deltadelta")
#(x1,y1) = fc.get_dynamical_correlator(i=i,j=j,delta=delta,name="cdcup")
#(x2,y2) = fc.get_dynamical_correlator(i=i,j=j,delta=delta,name="cdcdn")
# plot the result
  plt.plot(x,y.real,label=str(i))
#plt.plot(x1,y1.real,label="Up")
#plt.plot(x2,y2.real,label="Dn")
plt.legend()
plt.show()


