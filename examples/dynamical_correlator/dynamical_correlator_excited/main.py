# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0)
# of a small Heisenberg chain with random exchange. We compute it with
# DMRG using excited states (submode "EX"), with DMRG using the kernel
# polynomial method (submode "KPM"), and with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

np.random.seed(1)
n = 5 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1): # random first-neighbor exchange
    h = h + np.random.random()*(sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1])
sc.set_hamiltonian(h)
sc.maxm = 20 # bond dimension
sc.kpmmaxm = 20 # bond dimension of the KPM vectors
sc.nsweeps = 30 # number of DMRG sweeps

es = sc.get_excited(n=14, mode="DMRG")
print("Mean difference of excited energies DMRG-ED =",
      np.mean(np.abs(es - sc.get_excited(n=14, mode="ED"))))
name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
(x1, y1) = sc.get_dynamical_correlator(mode="ED", submode="ED", name=name)
(x2, y2) = sc.get_dynamical_correlator(mode="DMRG", submode="EX", name=name)
(x3, y3) = sc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name)

import matplotlib.pyplot as plt
plt.plot(x1, np.abs(y1), linewidth=4, label="ED")
plt.plot(x2, np.abs(y2), "--", linewidth=2, label="DMRG excited states")
plt.plot(x3, np.abs(y3), ":", linewidth=2, label="DMRG KPM")
plt.xlim([-0.5, 4.5])
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
