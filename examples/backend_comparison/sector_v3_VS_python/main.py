# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Addition spectrum E0(N) of a spinless t-V chain, H = sum_i c^+_i c_{i+1}
# + h.c. + V n_i n_{i+1}. With set_conserved_sector(Nf=N) each DMRG
# calculation is confined to N particles. We compare the backend
# itensor_version=3 and the pure-Python one with each other.
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of sites
V = 2.0 # interaction

def addition_spectrum(version):
    fc = fermionchain.Fermionic_Chain(n)
    if version == "python": fc.setup_python()
    else: fc.setup_cpp(version=version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
        h = h + V*fc.N[i]*fc.N[i+1]
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    es = []
    for nf in range(n+1):
        fc.set_conserved_sector(Nf=nf) # fix the particle number
        es.append(fc.gs_energy())
    return np.array(es)

e3 = addition_spectrum(3)
e2 = addition_spectrum("python")
print("Max difference =", np.max(np.abs(e3-e2)))
assert np.max(np.abs(e3-e2)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(range(n+1), e3, linewidth=4, label="itensor_version=3")
plt.plot(range(n+1), e2, "--", linewidth=2, label="python")
plt.xlabel("particle number N")
plt.ylabel("E0(N)")
plt.legend()
plt.show()
