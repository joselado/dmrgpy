import numpy as np
from dmrgpy import fermionchain

nf = 6 # number of different spinless fermionic orbitals
fc = fermionchain.Fermionic_Hamiltonian(nf) # create the object
C = fc.C # annihilation
Cdag = fc.Cdag # creation

H = 0 # initialize Hamiltonian 

# random Hamiltonian
for i in range(nf):
  for j in range(nf):
    H = H + Cdag[i]*C[j]*np.random.random() # random hopping

for i in range(nf-1):
    H = H + Cdag[i]*C[i]*Cdag[i+1]*C[i+1]*np.random.random() 

H = H + H.get_dagger() # make it Hermitian

fc.set_hamiltonian(H) # set the Hamiltonian
print("Energy with DMRG",fc.gs_energy(mode="DMRG")) # energy with DMRG
print("Energy with ED",fc.gs_energy(mode="ED")) # energy with exact diag.

nc = 6 # number of random correlators
def r(): return np.random.randint(0,nf-1)
ops = [fc.Cdag[r()]*fc.C[r()] for i in range(nc)] # random pairs of operators

# now compute the static correlators with ED and DMRG
c1 = np.array([fc.vev(o,mode="DMRG") for o in ops])
c2 = np.array([fc.vev(o,mode="ED") for o in ops])

import matplotlib.pyplot as plt
plt.scatter(range(nc),c1.real,label="DMRG",c="blue")
plt.plot(range(nc),c2.real,label="ED",c="red")
plt.legend()

plt.show()


