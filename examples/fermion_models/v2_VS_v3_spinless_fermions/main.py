# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spinless fermion hopping chain, H = sum_i Cdag_i C_{i+1} + h.c.
# We compute the correlation <Cdag_0 C_j> in the ground state with the
# two DMRG backends, ITensor v2 and v3, and compare them.
# This correlation is the one-body density matrix of the chain.
import numpy as np
from dmrgpy import fermionchain

n = 8 # number of sites

def correlation(version):
    fc = fermionchain.Fermionic_Chain(n, itensor_version=version)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1]
    h = h + h.get_dagger()
    fc.set_hamiltonian(h)
    fc.maxm = 20 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    wf = fc.get_gs()
    return fc.gs_energy(), wf.get_correlation_matrix()[0, :]

e2, c2 = correlation(2)
e3, c3 = correlation(3)
print("Energy v2 and v3 =", e2, e3)
print("Max difference v2-v3 =", np.max(np.abs(c2-c3)))

import matplotlib.pyplot as plt
plt.plot(range(n), c2.real, linewidth=4, label="ITensor v2")
plt.plot(range(n), c3.real, "--", linewidth=2, label="ITensor v3")
plt.xlabel("site j")
plt.ylabel("<Cdag_0 C_j>")
plt.legend()
plt.show()
