# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Density-density spectrum of a small interacting spinless fermion chain,
# S(w) = sum_n |<n|N_0|GS>|^2 delta(w - E_n + E_0).
# It is computed with exact diagonalization (matrix inversion), and with
# DMRG using the kernel polynomial method (KPM) and correction vectors (CVM).
# The KPM peaks have the Jackson shape, so they are higher than Lorentzians.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
fc = fermionchain.Fermionic_Chain(n)

# hopping and nearest-neighbor interaction
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h = h + fc.N[i]*fc.N[i+1]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 5.0, 100) # energies
delta = 0.1 # broadening of the peaks
name = (fc.N[0], fc.N[0]) # the two operators of the correlator
(x0, y0) = fc.get_dynamical_correlator(mode="ED", submode="INV", name=name, es=es, delta=delta)
(x1, y1) = fc.get_dynamical_correlator(mode="DMRG", submode="KPM", name=name, es=es, delta=delta)
(x2, y2) = fc.get_dynamical_correlator(mode="DMRG", submode="CVM", name=name, es=es, delta=delta)

import matplotlib.pyplot as plt
plt.plot(x0, y0.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG-KPM")
plt.plot(x2, y2.real, ":", linewidth=2, label="DMRG-CVM")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
