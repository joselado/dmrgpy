# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 5
fc = fermionchain.Fermionic_Hamiltonian(n,spinful=False) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # Make it Hermitian


def ft(i,j):
    return m[i,j]

def fu(i,j):
    if i==n//2 and j==n//2+1: return 1.0
    if j==n//2 and i==n//2+1: return 1.0
    else: return 0.0

# Initialize the Hamiltonian
fc.set_hoppings(ft) # hoppings
fc.set_hubbard(fu) # hoppings
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


### Compute the dyamical correlator ###

name = "density" # name of the correlator
x0,y0 = fc.get_dynamical_correlator(i=n//2,j=n//2+1,mode="ED",name=name,
        delta=1e-2)
x1,y1 = fc.get_dynamical_correlator(i=n//2,j=n//2+1,mode="DMRG",name=name,
        delta=1e-3)



### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x0,y0.real,label="ED")
plt.plot(x1,y1.real,label="DMRG")
plt.ylabel("density-density correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()

