# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain

# random matrix for the hoppings
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H
m *= 5.
# random matrix for hubbard
mu = np.matrix(np.random.random((n,n)))
mu = mu + mu.H
mu *= 0.1
#mu *= 0. ; mu[0,1] = 1.0 ; mu[1,0] = 1.0 ; m = m*0. ; m[0,0] =-2.0;m[1,1]=-2.0
# create hoppings and hubbard
fc.set_hoppings(lambda i,j: m[i,j]) # hoppings
fc.set_hubbard(lambda i,j: mu[i,j]) # density-density interaction

e = fc.gs_energy(mode="DMRG") # energy with DMRG
e1 = fc.gs_energy(mode="ED") # energy with ED
print("Ground state energy DMRG",e)
print("Ground state energy ED",e1)
exit()


import matplotlib.pyplot as plt

# Now compute the ground state correlator
pairs = [(0,i) for i in range(n-1)]
cs = fc.get_correlator(pairs=pairs)
x = np.array(range(len(cs))) # distances
plt.plot(x,cs,c="blue",label="DMRG")
plt.xlabel("N")
plt.xlabel("Correlator")
plt.legend()
plt.show()


