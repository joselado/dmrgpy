# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator of a Z3 parafermion chain with random couplings,
# S(w) = sum_n <GS|A|n><n|B|GS> delta(w - E_n + E_0), with A and B
# parafermion operators. We compute it with the kernel polynomial method
# (KPM) in DMRG and in exact diagonalization (ED), and with the matrix
# inversion method (INV) in ED. The peaks are at the same energies.
import numpy as np
from dmrgpy import parafermionchain

np.random.seed(1) # same random couplings every run
n = 3 # number of sites
pc = parafermionchain.Parafermionic_Chain(n, Z=3)
h = 0
for i in range(n):
    for j in range(n): # random couplings between all pairs
        h = h + pc.N[i]*pc.N[j]*np.random.random()
        h = h + pc.Sig[i]*pc.Sig[j]*np.random.random()
        h = h + 1j*pc.Sig[i]*pc.Sigd[j]*np.random.random()
        h = h + pc.Tau[i]*pc.Tau[j]*np.random.random()
        h = h + 1j*pc.Tau[i]*pc.Taud[j]*np.random.random()
h = h + h.get_dagger() # make it Hermitian
pc.set_hamiltonian(h)
pc.maxm = 60 # bond dimension
pc.kpmmaxm = 60 # bond dimension of the KPM vectors

name = [pc.Chi[0], pc.Psi[1]] # the two operators of the correlator
delta = 0.1 # broadening of the peaks
e1, d1 = pc.get_dynamical_correlator(name=name, mode="ED", submode="INV", delta=delta)
e2, d2 = pc.get_dynamical_correlator(name=name, mode="ED", submode="KPM", delta=delta)
e3, d3 = pc.get_dynamical_correlator(name=name, mode="DMRG", submode="KPM", delta=delta)

import matplotlib.pyplot as plt
plt.plot(e1, d1.real, linewidth=4, label="ED-INV")
plt.plot(e2, d2.real, linewidth=2, label="ED-KPM")
plt.plot(e3, d3.real, "--", linewidth=2, label="DMRG-KPM")
plt.xlabel("energy w")
plt.ylabel("S(w)")
plt.legend()
plt.show()
