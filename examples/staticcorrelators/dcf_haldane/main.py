# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator of Sz_i in the S=1 Heisenberg chain (Haldane chain),
# S_i(w) = sum_n |<n|Sz_i|GS>|^2 delta(w - E_n + E_0).
# It is computed site by site with the kernel polynomial method. The map
# of site against energy is what inelastic neutron scattering measures.
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain([3 for i in range(n)]) # 3 means S=1

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 30 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
sc.kpmmaxm = 10 # bond dimension used in the KPM

data = []
for i in range(n): # loop over sites
    (xs, ys) = sc.get_dynamical_correlator(mode="DMRG", name=(sc.Sz[i], sc.Sz[i]))
    data.append(np.abs(ys))

import matplotlib.pyplot as plt
plt.imshow(data, aspect="auto", origin="lower", extent=[xs[0], xs[-1], 0, n-1])
plt.colorbar(label="|S_i(w)|")
plt.xlabel("energy w")
plt.ylabel("site i")
plt.show()
