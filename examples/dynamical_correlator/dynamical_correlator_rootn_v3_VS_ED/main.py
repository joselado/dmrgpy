# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum S(w) of the S=1/2 Heisenberg chain, computed with
# the root-N Krylov correction vector method (submode="ROOTN") in DMRG,
# and compared with exact diagonalization (ED). With a Krylov space of
# 25 vectors the root-N result must agree with ED.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 60 # bond dimension
sc.cutoff = 1e-10 # truncation error
sc.nsweeps = 10 # number of DMRG sweeps

name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
es = np.linspace(0.5, 3.0, 20) # energies
delta = 0.2 # broadening of the peaks
(x1, y1) = sc.get_dynamical_correlator(mode="ED", submode="ED",
        name=name, es=es, delta=delta)
(x2, y2) = sc.get_dynamical_correlator(mode="DMRG", submode="ROOTN",
        name=name, es=es, delta=delta, N=6, nkry=25)
err = np.max(np.abs(y2.real - y1.real))
print("Max difference =", err)
assert err < 1e-8

import matplotlib.pyplot as plt
plt.plot(x1, y1.real, linewidth=4, label="ED")
plt.plot(x2, y2.real, "--", linewidth=2, label="DMRG root-N")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
