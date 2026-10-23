# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../../src')

# Correlation entropy density of the ground state of interacting
# spinless fermions, on each site. It measures how far the state is from
# a single Slater determinant, so it is zero without interactions.
from dmrgpy import fermionchain

def get(mode):
    n = 6 # number of sites
    U = 1.0 # nearest-neighbor interaction
    fc = fermionchain.Fermionic_Chain(n)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + U*(fc.N[i]-.5)*(fc.N[i+1]-.5)
    h = h + 0.5*fc.N[0] # a local perturbation
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    wf = fc.get_gs(mode=mode)
    return wf.get_correlation_entropy_density()

Sed = get("ED")
Sdmrg = get("DMRG")
print("Max difference DMRG-ED =", max(abs(Sed - Sdmrg)))

import matplotlib.pyplot as plt
plt.plot(Sed, linewidth=4, label="ED")
plt.plot(Sdmrg, "--", linewidth=2, label="DMRG")
plt.xlabel("site")
plt.ylabel("correlation entropy density")
plt.legend()
plt.show()
