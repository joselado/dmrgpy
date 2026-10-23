# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlation matrix <c^dag_i c_j> of the ground state of a spinless
# fermion chain with complex hopping and interaction U. It gives the
# density and the coherence between sites. We compute it with the
# "fast" (Python) and the "full" (C++) algorithms, which must agree.
import numpy as np
from dmrgpy import fermionchain

n = 12 # number of sites
U = -0.95 # nearest-neighbor interaction
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n-1):
    h = h + 1j*fc.Cdag[i]*fc.C[i+1] # complex hopping
    h = h + U*(fc.N[i]-.5)*(fc.N[i+1]-.5) # interaction
h = h + h.get_dagger()
fc.maxm = 20 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps
fc.set_hamiltonian(h)
wf = fc.get_gs(mode="DMRG")

m1 = wf.get_correlation_matrix(dmmode="fast") # Python algorithm
m2 = wf.get_correlation_matrix(dmmode="full") # C++ algorithm
print("Max difference =", np.max(np.abs(m1-m2)))

import matplotlib.pyplot as plt
plt.plot(np.abs(m2[0]), linewidth=4, label="full (C++)")
plt.plot(np.abs(m1[0]), "--", linewidth=2, label="fast (Python)")
plt.xlabel("site j")
plt.ylabel("|<c^dag_0 c_j>|")
plt.legend()
plt.show()
