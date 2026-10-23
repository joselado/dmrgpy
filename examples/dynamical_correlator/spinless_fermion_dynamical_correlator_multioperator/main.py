# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spectral function A(w) = sum_n |<n|Cdag_i|GS>|^2 delta(w - E_n + E_0)
# of a chain of spinless fermions, H = sum_i Cdag_i C_{i+1} + h.c.
# It is the local density of states that STM would measure.
# We compute it with DMRG and compare it with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)

# first-neighbor hopping
h = 0
for i in range(n-1):
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

name = (fc.Cdag[1], fc.C[1]) # the two operators of the correlator
(x0, y0) = fc.get_dynamical_correlator(mode="ED", name=name, delta=1e-2)
(x1, y1) = fc.get_dynamical_correlator(mode="DMRG", name=name, delta=1e-2)

import matplotlib.pyplot as plt
plt.plot(x0, y0.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("A(w), site 1")
plt.legend()
plt.show()
