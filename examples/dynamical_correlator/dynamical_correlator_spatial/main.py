# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spatially resolved spin spectrum of the S=1/2 Heisenberg chain,
# S_i(w) = sum_n |<n|Sz_i|GS>|^2 delta(w - E_n + E_0).
# It tells how much each site contributes to each spin excitation,
# as measured by spin-polarized STM on a chain of magnetic atoms.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
delta = 0.05 # broadening of the peaks
zs = [] # spectrum of each site
for i in range(n):
    name = (sc.Sz[i], sc.Sz[i])
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=delta)
    zs.append(s.real)

import matplotlib.pyplot as plt
plt.contourf(range(n), e, np.array(zs).T, 100, cmap="magma")
plt.xlabel("site i")
plt.ylabel("energy w")
plt.colorbar(label="S_i(w)")
plt.show()
