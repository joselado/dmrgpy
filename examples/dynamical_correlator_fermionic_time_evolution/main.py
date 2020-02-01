# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 3
fc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # Make it Hermitian

h = 0

for i in range(n-1):
      h = h + np.random.random()*fc.Cdag[i]*fc.C[i+1]
      h = h  + (np.random.random()-.5)*fc.N[i]*fc.N[i+1]

h = h + h.get_dagger()

fc.set_hamiltonian(h) # Hamiltonian
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


### Compute the dyamical correlator ###

i,j = 1,1
name = (fc.N[i],fc.N[j])

es = np.linspace(-0.5,6.0,100) # energies of the correlator
delta = 1e-1 # smearing of the correlator
import time
t1 = time.time()
x1,y1 = fc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name,
        es=es,delta=delta)
t2 = time.time()
x2,y2 = fc.get_dynamical_correlator(mode="DMRG",submode="TD",name=name,
        es=es,delta=delta)
t3 = time.time()

print("Time in KPM",t2-t1)
print("Time in TD",t3-t2)


### Plot the result ###

import matplotlib.pyplot as plt

#plt.plot(x0,y0.real,label="ED",marker="o")
plt.plot(x1,y1.real,label="DMRG-KPM",marker="o")
plt.plot(x2,y2.real,label="DMRG-TD",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()

