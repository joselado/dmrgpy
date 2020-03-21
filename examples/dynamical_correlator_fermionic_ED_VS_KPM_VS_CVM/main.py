# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 4
fc = fermionchain.Fermionic_Chain(n) # create the chain
h = 0

for i in range(n):
  for j in range(n):
      t = np.random.random() + 1j*np.random.random()
      h = h + fc.Cdag[i]*fc.C[j]*t
      h = h + fc.N[i]*fc.N[j]*np.random.random()


h = h + h.get_dagger()

fc.set_hamiltonian(h)

e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


### Compute the dyamical correlator ###

i,j = 1,1
name = (fc.C[i],fc.Cdag[j])

es = np.linspace(-0.5,6.0,100) # energies of the correlator
delta = 3e-2 # smearing of the correlator
x0,y0 = fc.get_dynamical_correlator(mode="ED",name=name,
        es=es,delta=delta)
x1,y1 = fc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name,
        es=es,delta=delta)
x2,y2 = fc.get_dynamical_correlator(mode="DMRG",submode="CVM",name=name,
        es=es,delta=delta)



### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x0,y0.real,label="ED",marker="o")
plt.plot(x1,y1.real,label="DMRG-KPM",marker="o")
plt.plot(x2,y2.real,label="DMRG-CVM",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()

