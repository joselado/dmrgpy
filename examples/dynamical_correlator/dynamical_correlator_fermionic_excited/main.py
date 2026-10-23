# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spectral function <Cdag_1 C_2> of a small interacting spinless fermion
# chain, G(w) = sum_n <GS|Cdag_1|n><n|C_2|GS> delta(w - E_n + E_0).
# We compute it with three methods: exact diagonalization (ED), DMRG with
# the kernel polynomial method (KPM), and DMRG with explicit excited states.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)

# hopping between all pairs of sites, and nearest neighbor interaction
h = 0
for i in range(n):
    for j in range(n):
        h = h + fc.Cdag[i]*fc.C[j]
for i in range(n-1):
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

name = (fc.Cdag[1], fc.C[2]) # the two operators of the correlator
es = np.linspace(-0.5, 6.0, 800) # energies
delta = 0.03 # broadening of the peaks
(x0, y0) = fc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED", submode="KPM")
(x1, y1) = fc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG", submode="KPM")
(x2, y2) = fc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG", submode="EX", nex=10)

import matplotlib.pyplot as plt
plt.plot(x0, y0.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG-KPM")
plt.plot(x2, y2.real, ":", linewidth=2, label="DMRG-excited")
plt.xlabel("energy w")
plt.ylabel("G(w)")
plt.legend()
plt.show()
