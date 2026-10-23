# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin-spin correlation <Sz_0 Sz_j> in the ground state of the S=1/2
# Heisenberg chain. The sign alternates (antiferromagnetic order) and the
# amplitude decays slowly with distance. We compute it with DMRG and
# check it against exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 12 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

# correlation between the first site and site j
cdmrg = [sc.vev(sc.Sz[0]*sc.Sz[j], mode="DMRG").real for j in range(n)]
ced = [sc.vev(sc.Sz[0]*sc.Sz[j], mode="ED").real for j in range(n)]
print("Max difference DMRG-ED =", np.max(np.abs(np.array(cdmrg)-ced)))

import matplotlib.pyplot as plt
plt.plot(range(n), ced, "o", markersize=10, label="ED")
plt.plot(range(n), cdmrg, "x", markersize=10, label="DMRG")
plt.axhline(0, color="gray", linewidth=0.5)
plt.xlabel("site j")
plt.ylabel("<Sz_0 Sz_j>")
plt.legend()
plt.show()
