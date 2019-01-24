# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import fermionchain
n = 4 # number of spinful fermionic sites
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
####### Input matrices #######
# Array with the hoppings and with hubbard couplings
# These are the matrices that you have to modify
hopping = np.zeros((n,n))
hubbard = np.zeros((n,n))
for i in range(n-1):  hopping[i,i+1] = 1. ; hopping[i+1,i] = 1.
for i in range(n): 
#    if i<n//2: U = 8.0 
#    else: U = 8.0
#    hubbard[i,i] = U/2. 
    hopping[i,i] = 5.0 #-1.0
# The implemented Hamiltonian is
# H = \sum_ij hopping[i,j] c^dagger_i c_j + hubbard[i,j] n_i n_j
# with n_i = c^\dagger_{i,up} c_{i,up} + c^\dagger_{i,dn} c_{i,dn}
# the previous matrices are for a half filled Hubbard chain
##############################
# Setup the Many Body Hamiltonian
fc.set_hoppings(lambda i,j: hopping[i,j]) # set the hoppings
#print(fc.get_excited())
#exit()
(xf,yf) = fc.get_gr_free()
(x,y) = fc.get_gr()
#print(x,y)
import matplotlib.pyplot as plt
#x = x - 2.
plt.plot(x,y.real,c="blue")
plt.scatter(xf,yf.real,c="blue")
plt.plot(x,y.imag,c="red")
plt.scatter(xf,yf.imag,c="red")
plt.show()


