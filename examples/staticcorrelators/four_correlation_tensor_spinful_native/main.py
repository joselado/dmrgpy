# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Four-point correlator <Cdag_i C_j Cdag_k C_l> in the ground state of a
# spinful Hubbard chain, built with native Electron sites. It is the
# object needed to compute the correlation entropy. We compute it with
# DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 3 # number of sites
U = 1.3 # Hubbard interaction
fc = fermionchain.Spinful_Fermionic_Chain_Native(n)
h = 0
for i in range(n-1): # complex hopping
    h = h + (0.9+0.2j)*fc.Cdagup[i]*fc.Cup[i+1]
    h = h + (0.9+0.2j)*fc.Cdagdn[i]*fc.Cdn[i+1]
for i in range(n): # Hubbard term
    h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
h = h + h.get_dagger()
h = h + 0.3*(fc.Nup[0]-fc.Ndn[0]) # lifts the spin degeneracy
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 14 # number of DMRG sweeps
fc.setup_cpp(3) # native sites need ITensor v3

ct_dmrg = fc.get_gs(mode="DMRG").get_four_correlation_tensor(ctmode="full")
ct_ed = fc.get_gs(mode="ED").get_four_correlation_tensor()
diff = np.max(np.abs(ct_dmrg - ct_ed))
print("Max difference DMRG-ED =", diff)
assert diff < 1e-4

import matplotlib.pyplot as plt
plt.plot(ct_ed.flatten().real, linewidth=4, label="ED")
plt.plot(ct_dmrg.flatten().real, "--", linewidth=2, label="DMRG")
plt.xlabel("index (i,j,k,l)")
plt.ylabel("<Cdag_i C_j Cdag_k C_l>")
plt.legend()
plt.show()
