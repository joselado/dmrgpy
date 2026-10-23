# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S_i(w) = sum_n |<n|Sz_i|GS>|^2 delta(w - E_n + E_0)
# of an S=1/2 Heisenberg chain with an S=1 impurity in the middle.
# The impurity changes the spectrum on the sites close to it.
# This is the signal that spin-polarized STM would measure at each site.
import numpy as np
from dmrgpy import spinchain

n = 7 # number of sites, odd so the impurity is in the middle
spins = ["S=1/2" for i in range(n)]
spins[n//2] = "S=1" # the impurity
sc = spinchain.Spin_Chain(spins)
h = 0
for i in range(n-1):
    h = h + sc.SS(i, i+1) # S_i . S_{i+1}
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
zs = []
for i in range(n): # loop over sites
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", i=i, j=i, name="ZZ", es=es, delta=0.1)
    zs.append(s.real)
zs = np.array(zs).T # energy along the first axis

import matplotlib.pyplot as plt
plt.contourf(range(n), es, zs, 100, cmap="hot")
plt.xlabel("site i")
plt.ylabel("energy w")
plt.colorbar(label="S_i(w)")
plt.show()
