# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Density-density spectrum of a small interacting fermion chain,
# S(w) = sum_n |<n|N_1|GS>|^2 delta(w - E_n + E_0).
# Peaks are charge excitations, seen for example in electron
# energy loss spectroscopy. We compute it with DMRG, with KPM and with
# time evolution (TD), and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
np.random.seed(1)
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1): # random hopping and interaction
    h = h + np.random.random()*fc.Cdag[i]*fc.C[i+1]
    h = h + (np.random.random()-.5)*fc.N[i]*fc.N[i+1]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 6.0, 100) # energies
delta = 0.1 # broadening of the peaks
name = (fc.N[1], fc.N[1])
(x1, y1) = fc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name, es=es, delta=delta)
(x2, y2) = fc.get_dynamical_correlator(mode="DMRG", submode="TD", name=name, es=es, delta=delta)
(x3, y3) = fc.get_dynamical_correlator(mode="ED", name=name, es=es, delta=delta)
print("Max difference KPM-ED =", np.max(np.abs(y1-y3)))

import matplotlib.pyplot as plt
plt.plot(x3, y3.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG-KPM")
plt.plot(x2, y2.real, "--", linewidth=2, label="DMRG-TD")
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
