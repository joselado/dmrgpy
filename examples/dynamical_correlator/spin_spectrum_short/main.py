# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Local spin excitation spectrum of the S=1/2 Heisenberg chain,
# S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0).
# Each peak is a spin excitation, at the energy that inelastic
# neutron scattering or spin-polarized STM would measure.
# We compute it with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import spinchain

n = 8 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 3.0, 300) # energies
delta = 0.05 # broadening of the peaks
name = (sc.Sz[0], sc.Sz[0]) # the two operators of the correlator
(e1, s1) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG")
(e2, s2) = sc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED")

import matplotlib.pyplot as plt
plt.plot(e2, np.real(s2), linewidth=4, label="ED")
plt.plot(e1, np.real(s1), "--", linewidth=2, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
