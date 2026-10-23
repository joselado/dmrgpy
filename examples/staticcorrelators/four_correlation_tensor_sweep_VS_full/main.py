# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Four-point correlator <Cdag_i C_j Cdag_k C_l> of the ground state of
# interacting spinless fermions. It is the building block of the
# correlation entropy of a correlated state. We compute it with the
# fast single sweep (ctmode="sweep") and with one operator per entry
# (ctmode="full"), which must give the same tensor.
import numpy as np
from dmrgpy import fermionchain

n = 5 # number of sites
fc = fermionchain.Fermionic_Chain(n, itensor_version=3)
h = 0
for i in range(n-1): # hopping and first-neighbor interaction
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] + 0.5*fc.N[i]*fc.N[i+1]
for i in range(n): # onsite energy
    h = h + 0.1*fc.N[i]
fc.set_hamiltonian(h)
fc.maxm = 60 # bond dimension
fc.nsweeps = 8 # number of DMRG sweeps
wf = fc.get_gs(mode="DMRG")

c1 = wf.get_four_correlation_tensor(ctmode="sweep").flatten()
c2 = wf.get_four_correlation_tensor(ctmode="full").flatten()
print("Max difference sweep-full =", np.max(np.abs(c1-c2)))
assert np.max(np.abs(c1-c2)) < 1e-8

import matplotlib.pyplot as plt
plt.plot(c2.real, linewidth=4, label="full")
plt.plot(c1.real, "--", linewidth=2, label="sweep")
plt.xlabel("index (i,j,k,l)")
plt.ylabel("Re <Cdag_i C_j Cdag_k C_l>")
plt.legend()
plt.show()
