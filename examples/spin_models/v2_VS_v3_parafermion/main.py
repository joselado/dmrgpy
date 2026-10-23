# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground-state energy of a random Z3 parafermion chain, computed with the
# ITensor v2 and v3 backends of DMRG and with exact diagonalization (ED).
# The three results must agree, for every chain length.
import numpy as np
from dmrgpy import parafermionchain

def get_energy(n, version, mode="DMRG"):
    pc = parafermionchain.Parafermionic_Chain(n)
    pc.setup_cpp(version) # ITensor version
    np.random.seed(3) # same random Hamiltonian for every version
    h = 0
    for i in range(n):
        for j in range(n):
            h = h + pc.N[i]*pc.N[j]*np.random.random()
            h = h + pc.Sig[i]*pc.Sigd[j]*np.random.random()
            h = h + pc.Tau[i]*pc.Taud[j]*np.random.random()
    h = h + h.get_dagger() # make it Hermitian
    pc.set_hamiltonian(h)
    pc.maxm = 40 # bond dimension
    pc.nsweeps = 10 # number of DMRG sweeps
    return pc.gs_energy(mode=mode)

ns = [3, 4, 5, 6] # number of sites
e2 = [get_energy(n, 2) for n in ns]
e3 = [get_energy(n, 3) for n in ns]
eed = [get_energy(n, 2, mode="ED") for n in ns]
print("Max difference v2-v3 =", np.max(np.abs(np.array(e2)-e3)))
assert np.max(np.abs(np.array(e2)-eed)) < 1e-2

import matplotlib.pyplot as plt
plt.plot(ns, eed, linewidth=4, label="ED")
plt.plot(ns, e2, "--", linewidth=2, label="DMRG v2")
plt.plot(ns, e3, ":", linewidth=2, label="DMRG v3")
plt.xlabel("number of sites")
plt.ylabel("ground-state energy")
plt.legend()
plt.show()
