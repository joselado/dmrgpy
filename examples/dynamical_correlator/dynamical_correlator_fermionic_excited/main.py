# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6
fc = fermionchain.Fermionic_Chain(n) # create the chain
m = np.matrix(np.random.random((n,n)) + 1j*np.random.random((n,n)))
m = m + m.H # Make it Hermitian
h = 0

for i in range(n):
    for j in range(n):
        h = h + fc.Cdag[i]*fc.C[j]

for i in range(n-1):
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)


# Initialize the Hamiltonian
fc.set_hamiltonian(h) # Hamiltonian
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)


### Compute the dyamical correlator ###

i,j = 1,2
name = (fc.Cdag[i],fc.C[j])

es = np.linspace(-0.5,6.0,100) # energies of the correlator
eskpm = np.linspace(min(es),max(es),len(es)*8) # finer grid for KPM


delta = 3e-2 # smearing of the correlator
x2,y2 = fc.get_dynamical_correlator(mode="DMRG",submode="EX",name=name,
        purify = False,
        nex=10, # number of excited states
        es=eskpm,delta=delta)
x0,y0 = fc.get_dynamical_correlator(mode="ED",name=name,submode="KPM",
        es=eskpm,delta=delta)
x1,y1 = fc.get_dynamical_correlator(mode="DMRG",submode="KPM",name=name,
        es=eskpm,delta=delta)



### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x0,y0.real,label="ED",marker="o")
plt.plot(x1,y1.real,label="DMRG-KPM",marker="o")
plt.plot(x2,y2.real,label="DMRG-Excited",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()










