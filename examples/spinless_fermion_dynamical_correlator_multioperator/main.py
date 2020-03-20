# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # Make it Hermitian


def ft(i,j):
    if abs(i-j)==1: return 1.0
    return 0.0
    return m[i,j]

def fu(i,j):
    if abs(i-j)==1: return 1.0
    else: return 0.0

h = 0
for i in range(n):
  for j in range(n):
    h = h + ft(i,j)*fc.Cdag[i]*fc.C[j]
   # h = h + fu(i,j)*fc.N[i]*fc.N[j]
#    h = h + fu(i,j)*fc.C[i]*fc.Cdag[i]*fc.C[j]*fc.Cdag[j]

# Initialize the Hamiltonian
fc.set_hamiltonian(h) # set the Hamiltonian

e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


### Compute the dyamical correlator ###
i,j = 1,1
mi = fc.Cdag[i]
mj = fc.C[j]

x0,y0 = fc.get_dynamical_correlator(mode="ED",name=(mi,mj),delta=1e-2)

x1,y1 = fc.get_dynamical_correlator(mode="DMRG",name=(mi,mj))
#x0,y0 = fc.get_dynamical_correlator(i=i,j=j,mode="DMRG",name=name,delta=1e-2)



### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x0,y0.real,label="ED")
plt.plot(x1,y1.real,label="DMRG")
plt.ylabel("density-density correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()

