# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Map of the local spin spectrum S_i(w) = sum_n |<n|Sz_i|GS>|^2
# delta(w - E_n + E_0) along a Heisenberg chain with one S=1 impurity in
# the middle. The impurity changes the spin excitations close to it.
# Spin-polarized STM measures S_i(w) site by site.
import numpy as np
from dmrgpy import spinchain

spins = ["S=1/2" for i in range(4)]
spins = spins + ["S=1"] + spins # S=1 impurity in the middle
n = len(spins)
sc = spinchain.Spin_Chain(spins)
# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
zs = []
for i in range(n): # loop over sites
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", name=(sc.Sz[i], sc.Sz[i]),
                                         es=es, delta=0.05)
    zs.append(s.real)
zs = np.array(zs).transpose()

import matplotlib.pyplot as plt
plt.contourf(range(n), e, zs/np.max(zs), 100, cmap="magma", vmax=0.8)
plt.xlabel("site")
plt.ylabel("energy w")
plt.show()
