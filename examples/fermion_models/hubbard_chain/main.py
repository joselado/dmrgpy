# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin spectrum of the Hubbard chain with n spinful sites,
# H = sum t c^+_i c_i+1 + h.c. + U (n_up - 1/2)(n_dn - 1/2).
# S(w) = sum_n |<n|Sz_0|GS>|^2 delta(w - E_n + E_0) is the spin
# excitation spectrum, measured by neutron scattering or spin-polarized STM.
# We compute it with DMRG and compare it with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger()
for i in range(n): # Hubbard interaction
    h = h + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 40 # number of DMRG sweeps
fc.kpmmaxm = 20 # bond dimension in KPM

name = (fc.Sz[0], fc.Sz[0])
es = np.linspace(-0.5, 4.0, 300) # energies
delta = 0.1 # broadening of the peaks
x1, y1 = fc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="DMRG")
x2, y2 = fc.get_dynamical_correlator(name=name, es=es, delta=delta, mode="ED", submode="INV")

import matplotlib.pyplot as plt
plt.plot(x2, y2.real, linewidth=4, label="ED")
plt.plot(x1, y1.real, "--", linewidth=2, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("S(w), first site")
plt.legend()
plt.show()
