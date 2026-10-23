# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Four-body correlation tensor of a ground state, M_ijkl = <Cdag_i C_j Cdag_k C_l>,
# for interacting spinless fermions with complex hopping. It is the input
# needed for two-particle quantities, for example the energy of any
# interaction written in the fermion operators. We compute it from the
# DMRG state and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
fc = fermionchain.Fermionic_Chain(n)

# random complex hopping plus random nearest-neighbor interaction
h = 0
for i in range(n-1):
    h = h + np.random.random()*(1.+1j)*fc.Cdag[i]*fc.C[i+1]
    h = h + np.random.random()*(fc.N[i]-.5)*(fc.N[i+1]-.5)
h = h + h.get_dagger() # make it Hermitian
fc.set_hamiltonian(h)
fc.maxm = 20 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

m1 = fc.get_gs(mode="DMRG").get_four_correlation_tensor(ctmode="full")
m2 = fc.get_gs(mode="ED").get_four_correlation_tensor()
print("Max difference DMRG-ED =", np.max(np.abs(m1-m2)))

# slice of the tensor, |<Cdag_0 C_j Cdag_0 C_1>| as a function of j
import matplotlib.pyplot as plt
plt.plot(range(n), np.abs(m2[0,:,0,1]), "o-", linewidth=4, label="ED")
plt.plot(range(n), np.abs(m1[0,:,0,1]), "x--", linewidth=2, label="DMRG")
plt.xlabel("site j")
plt.ylabel("|<Cdag_0 C_j Cdag_0 C_1>|")
plt.legend()
plt.show()
