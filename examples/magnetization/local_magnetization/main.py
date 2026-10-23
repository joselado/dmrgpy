# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local magnetization <Sz_i> of the S=1/2 Heisenberg chain with a magnetic
# field acting on the first site, H = sum_i S_i . S_{i+1} + 0.1 Sz_0.
# The field polarizes the edge and the response decays into the chain,
# which a spin-polarized STM would measure site by site.
# We compute it with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
h = h + 0.1*sc.Sz[0] # field on the first site
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

mx1, my1, mz1 = sc.get_magnetization(mode="DMRG")
mx2, my2, mz2 = sc.get_magnetization(mode="ED")
print("Max difference DMRG-ED =", np.max(np.abs(np.array(mz1)-mz2)))

import matplotlib.pyplot as plt
plt.plot(range(n), mz2, "o", markersize=10, label="ED")
plt.plot(range(n), mz1, "x", markersize=10, label="DMRG")
plt.xlabel("site i")
plt.ylabel("<Sz_i>")
plt.legend()
plt.show()
