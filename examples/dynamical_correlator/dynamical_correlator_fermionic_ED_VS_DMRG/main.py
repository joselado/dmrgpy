# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

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
#      h = h + -1*fc.N[i]*fc.N[j]*np.random.random()


h = h + h.get_dagger()

fc.set_hamiltonian(h)

e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)

print("Excited",fc.get_excited(mode="ED",n=10))

### Compute the dyamical correlator ###

def randop():
    i = np.random.randint(n)
    j = np.random.randint(n)
    print(i)
    ci = np.exp(2*1j*np.random.random()*np.pi)
    cj = np.exp(2*1j*np.random.random()*np.pi)
    return ci*fc.C[i] #+ cj*fc.Cdag[j]

name = (randop().get_dagger(),randop())

es = np.linspace(-0.5,6.0,300) # energies of the correlator
delta = 1e-1 # smearing of the correlator
x0,y0 = fc.get_dynamical_correlator(mode="ED",submode="KPM",name=name,
        es=es,delta=delta)
x1,y1 = fc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name,
        es=es,delta=delta)
x2,y2 = fc.get_dynamical_correlator(mode="ED",submode="INV",name=name,
        es=es,delta=delta)



### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x0,y0.real,label="ED-KPM",marker="o")
plt.plot(x2,y2.real,label="ED-INV",marker="o")
plt.plot(x1,y1.real,label="DMRG-KPM",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()










