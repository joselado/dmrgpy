# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state of an interacting non-Hermitian chain of spinless fermions,
# H = sum_ij M_ij c^dag_i c_j + V sum_i (n_i - 1/2)(n_{i+1} - 1/2),
# where M has first-neighbor hopping and an imaginary onsite potential
# i*gamma*cos(2 pi sqrt(2) i). We compare ED with the Julia DMRG backend.
# Needs Julia and juliacall installed.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites

def get_energy(gamma, julia):
    if julia: fc = fermionchain.Fermionic_Chain(n, itensor_version="julia_live")
    else: fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1): # hopping
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    for i in range(n): # imaginary onsite potential
        h = h + 1j*gamma*np.cos(2*np.pi*np.sqrt(2)*i)*fc.Cdag[i]*fc.C[i]
    for i in range(n-1): # interaction
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
    fc.set_hamiltonian(h)
    return fc.gs_energy(mode="DMRG" if julia else "ED")

gammas = np.linspace(0.0, 0.6, 5)
e_ed = np.array([get_energy(g, False) for g in gammas])
e_jul = np.array([get_energy(g, True) for g in gammas])
print("Max difference =", np.max(np.abs(e_ed - e_jul)))

import matplotlib.pyplot as plt
plt.plot(gammas, e_ed.real, linewidth=4, label="ED")
plt.plot(gammas, e_jul.real, "--", linewidth=2, label="DMRG Julia")
plt.xlabel("imaginary potential gamma")
plt.ylabel("ground state energy (real part)")
plt.legend()
plt.show()
