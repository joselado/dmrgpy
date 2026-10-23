# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator <N_0(z) N_0> of a non-Hermitian fermion chain
# (staggered imaginary potential), computed with the non-Hermitian kernel
# polynomial method (NH-KPM) on the biorthogonal ground state.
# We compare exact diagonalization (ED) with DMRG (ITensor v3).
import numpy as np
from dmrgpy import fermionchain

n = 4 # number of sites
def chain():
    fc = fermionchain.Fermionic_Chain(n)
    fc.setup_cpp(version=3)
    h = 0
    for i in range(n-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
        h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    for i in range(n): h = h + 1j*(-1)**i*fc.Cdag[i]*fc.C[i] # imaginary potential
    h = h + 0.03*fc.N[0] # small field, removes a degeneracy
    fc.set_hamiltonian(h)
    return fc

es = np.linspace(0., 4., 12) # frequencies
kw = dict(E_max=10, n=80, delta=0.3, es=es, submode="KPM")
fc = chain()
e1, y1 = fc.get_dynamical_correlator(mode="ED", name=[fc.N[0], fc.N[0]], **kw)
fc = chain()
e2, y2 = fc.get_dynamical_correlator(mode="DMRG", name=[fc.N[0], fc.N[0]], **kw)
diff = np.max(np.abs(y1-y2))
print("Max difference =", diff)
assert diff < 1e-6

import matplotlib.pyplot as plt
plt.plot(e1, y1.real, linewidth=4, label="ED")
plt.plot(e2, y2.real, "--", linewidth=2, label="DMRG")
plt.xlabel("frequency")
plt.ylabel("Re <N_0(z) N_0>")
plt.legend()
plt.show()
