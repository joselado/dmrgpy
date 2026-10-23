# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin spectrum of the S=1/2 Heisenberg chain,
# S(w) = sum_n |<n|Sx_0|GS>|^2 delta(w - E_n + E_0).
# The peaks are the spin excitations seen in inelastic neutron scattering.
# We compute it with the correction vector method (CVM) and compare it
# with the kernel polynomial method (KPM).
import numpy as np
from dmrgpy import spinchain

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 10 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps
sc.kpmmaxm = 10 # bond dimension of KPM
sc.cvm_tol = 1e-3 # tolerance of CVM

es = np.linspace(-0.5, 4.5, 100) # energies
name = (sc.Sx[0], sc.Sx[0]) # the two operators of the correlator
(e1, s1) = sc.get_dynamical_correlator(submode="KPM", name=name, es=es, delta=0.1)
(e2, s2) = sc.get_dynamical_correlator(submode="CVM", name=name, es=es, delta=0.1)

import matplotlib.pyplot as plt
plt.plot(e1, np.abs(s1), linewidth=4, label="KPM")
plt.plot(e2, np.abs(s2), "--", linewidth=2, label="CVM")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
