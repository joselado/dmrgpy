# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Magnetization Mz = sum_i <Sz_i> of the S=1/2 Heisenberg chain in a
# magnetic field B, H = sum_i S_i . S_{i+1} + B sum_i Sz_i.
# The magnetization grows in steps, and it saturates once the Zeeman
# energy wins over the exchange. We use DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 10 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
sc.maxm = 30 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

h0 = 0 # Heisenberg Hamiltonian
Mz = 0 # total magnetization
for i in range(n-1):
    h0 = h0 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    Mz = Mz + sc.Sz[i]

bs = np.linspace(0., 4.0, 20) # magnetic fields
m1, m2 = [], []
for b in bs:
    sc.set_hamiltonian(h0 + b*Mz)
    m1.append(sc.vev(Mz, mode="ED").real)
    m2.append(sc.vev(Mz, mode="DMRG").real)
print("Max difference DMRG-ED =", np.max(np.abs(np.array(m1)-m2)))

import matplotlib.pyplot as plt
plt.plot(bs, m1, linewidth=4, label="ED")
plt.plot(bs, m2, "--", linewidth=2, label="DMRG")
plt.xlabel("magnetic field B")
plt.ylabel("magnetization Mz")
plt.legend()
plt.show()
