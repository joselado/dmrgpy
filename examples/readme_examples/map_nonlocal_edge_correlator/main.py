# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin correlator <Sz_0 Sz_j>(w) between the first site and site j of a
# Heisenberg chain, as a function of the energy w. The map shows how a
# spin excitation created at the edge spreads along the chain.
# It is measured with inelastic neutron scattering or spin-flip STM.
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
for j in range(n): # correlator between site 0 and site j
    (e, s) = sc.get_dynamical_correlator(mode="DMRG", name=(sc.Sz[0], sc.Sz[j]), es=es, delta=0.05)
    zs.append(s.real)

import matplotlib.pyplot as plt
vmax = np.max(np.abs(zs))
plt.contourf(range(n), e, np.array(zs).T, 100, cmap="bwr", vmin=-vmax, vmax=vmax)
plt.xlabel("site j")
plt.ylabel("energy w")
plt.colorbar()
plt.show()
