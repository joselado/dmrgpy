# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Energy per site of the infinite S=1/2 Heisenberg chain with infinite
# DMRG (iDMRG), for a unit cell of one and of two sites. The exact result
# from the Bethe ansatz is 1/4 - log(2). We compare the pure Python
# backend with the compiled ITensor v3 backend.
import numpy as np
from dmrgpy import infinitechain

exact = 0.25 - np.log(2) # Bethe ansatz energy per site
ncs = [1, 2] # sites in the unit cell
e = {"python": [], 3: []}
for version in e:
    for nc in ncs:
        ic = infinitechain.Infinite_Spin_Chain(["1/2"]*nc, itensor_version=version)
        ic.gs_method = "idmrg"
        h = 0 # bonds inside the cell, and from the last site to the next cell
        for i in range(nc):
            if i < nc-1: (A, B, C) = (ic.SxC[i+1], ic.SyC[i+1], ic.SzC[i+1])
            else: (A, B, C) = (ic.SxR[0], ic.SyR[0], ic.SzR[0])
            h = h + ic.SxC[i]*A + ic.SyC[i]*B + ic.SzC[i]*C
        ic.set_hamiltonian(h)
        ic.maxm = 30 # bond dimension
        e[version].append(ic.gs_energy())
print("Max difference python-v3 =", np.max(np.abs(np.array(e["python"]) - e[3])))

import matplotlib.pyplot as plt
plt.plot(ncs, e["python"], "o", markersize=12, label="python")
plt.plot(ncs, e[3], "x", markersize=12, label="v3")
plt.axhline(exact, color="gray", linestyle="--", label="exact")
plt.xticks(ncs)
plt.xlabel("sites in the unit cell")
plt.ylabel("energy per site")
plt.legend()
plt.show()
