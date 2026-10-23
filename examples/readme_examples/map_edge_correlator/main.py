# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_i|GS>|^2 delta(w - E_n + E_0)
# at every site i of an open Heisenberg chain. The map shows how the
# spin excitations change from the edge to the bulk, which spin-polarized
# STM would measure site by site.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 4.0, 200) # energies
zs = []
for i in range(n): # loop over sites
    name = (sc.Sz[i], sc.Sz[i])
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=0.05)
    zs.append(s.real)
zs = np.array(zs).transpose()

import matplotlib.pyplot as plt
plt.contourf(range(n), e, zs, 100, cmap="magma")
plt.colorbar(label="S(w)")
plt.xlabel("site i")
plt.ylabel("energy w [J]")
plt.show()
