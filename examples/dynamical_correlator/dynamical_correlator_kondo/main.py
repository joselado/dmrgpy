# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
nb = 10
n  = nb*2 + 1 + 1 # total number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n) # create the chain
ic = nb # central site
h = 0
# Kondo Hamiltonian
U = 4.0 # Hubbard
h = U*(fc.Nup[ic]-0.5)*(fc.Ndn[ic]-0.5) # Kondo site
# now lets add the bath
eps = np.linspace(0.,1.0,nb) # energies for the bath
lamb = 0.6
eps = np.array([0.4**(nb-i-1) for i in range(nb)])
eps = np.sort(np.concatenate([-eps,[0.],eps])) # sorted energies
inds = [i for i in range(nb)] + [nb+1] + [nb+i+2 for i in range(nb)]
print(eps)
print(inds)
for i in range(len(eps)): # loop over sites in the bath
    t = 1./nb#*eps[i]
    e,ii = eps[i],inds[i]
    h = h + t*fc.Cdagup[ii]*fc.Cup[ic]
    h = h + t*fc.Cdagdn[ii]*fc.Cdn[ic]
    h = h + e*(fc.Nup[ii] + fc.Ndn[ii]) # energy spacing

h = h + h.get_dagger()

fc.set_hamiltonian(h)
fc.maxm = 20
fc.kpmmaxm = 20

#e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
name = (fc.Cdagup[ic],fc.Cup[ic])
print("Density",np.round(fc.get_density(),2))

es = np.linspace(-0.5,10.0,300) # energies of the correlator
delta = 1e-1 # smearing of the correlator
x1,y1 = fc.get_dynamical_correlator(name=name,
        es=es,delta=delta)

np.savetxt("SPECTRAL_FUNCTION.OUT",np.array([x1,y1]).T)


### Plot the result ###

import matplotlib.pyplot as plt

plt.plot(x1,y1.real,label="DMRG-KPM",marker="o")
plt.ylabel("Dynamical correlator")
plt.xlabel("Frequency")
plt.legend()
plt.show()










