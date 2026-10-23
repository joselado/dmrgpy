# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Ground state energy of a Hubbard chain with two interleaved spinless
# fermion sites (up and down) per physical site, as a function of the
# interaction U. H = -t sum c^dag_i c_j + h.c. + U sum (n_up-1/2)(n_dn-1/2).
# The energy sets the charge gap and the double occupancy. We compare the
# ITensor v2 and v3 backends.
import numpy as np
from dmrgpy import fermionchain

ns = 3 # number of physical sites

def get_energy(version, U):
    fc = fermionchain.Fermionic_Chain(2*ns, itensor_version=version)
    C = fc.C ; Cdag = fc.Cdag ; N = fc.N
    h = 0
    for i in range(ns-1):
        for j in range(2): # spin up and down
            h = h + Cdag[2*i+j]*C[2*(i+1)+j]
    h = h + h.get_dagger() # hopping
    for i in range(ns):
        h = h + U*(N[2*i]-0.5)*(N[2*i+1]-0.5) # Hubbard interaction
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    return fc.gs_energy()

Us = [0.0, 1.0, 2.0, 4.0, 8.0]
e2 = np.array([get_energy(2, U) for U in Us])
e3 = np.array([get_energy(3, U) for U in Us])
print("Max difference v2-v3 =", np.max(np.abs(e2-e3)))

import matplotlib.pyplot as plt
plt.plot(Us, e2, linewidth=4, label="ITensor v2")
plt.plot(Us, e3, "--", linewidth=2, label="ITensor v3")
plt.xlabel("Hubbard U")
plt.ylabel("ground state energy")
plt.legend()
plt.show()
