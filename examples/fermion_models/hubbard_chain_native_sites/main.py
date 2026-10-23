# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Hubbard chain with two kinds of sites: native spinful sites (one site of
# dimension 4 per orbital) and two spinless sites per orbital. The physics is
# the same, so the ground state energy and the spectral function
# A(w) = -Im G(w)/pi of the electron operator Cup_0 have to agree.
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
U = 2.0 # Hubbard interaction

def build(chain_class):
    fc = chain_class(n)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
    h = h + h.get_dagger() # make it Hermitian
    for i in range(n):
        h = h + U*(fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
    fc.set_hamiltonian(h)
    fc.maxm = 60 # bond dimension
    fc.nsweeps = 14 # number of DMRG sweeps
    fc.setup_cpp(3)
    return fc

es = np.linspace(-4, 4, 100) # energies
fc1 = build(fermionchain.Spinful_Fermionic_Chain_Native)
fc2 = build(fermionchain.Spinful_Fermionic_Chain)
print("Energies", fc1.gs_energy(), fc2.gs_energy())
(x, y1) = fc1.get_dynamical_correlator(name=(fc1.Cup[0], fc1.Cdagup[0]), es=es, delta=0.3)
(x, y2) = fc2.get_dynamical_correlator(name=(fc2.Cup[0], fc2.Cdagup[0]), es=es, delta=0.3)
print("Max difference =", np.max(np.abs(y1 - y2)))
assert np.max(np.abs(y1 - y2)) < 5e-3

import matplotlib.pyplot as plt
plt.plot(x, y2.real, linewidth=4, label="two spinless sites")
plt.plot(x, y1.real, "--", linewidth=2, label="native spinful site")
plt.xlabel("energy w")
plt.ylabel("A(w), first site")
plt.legend()
plt.show()
