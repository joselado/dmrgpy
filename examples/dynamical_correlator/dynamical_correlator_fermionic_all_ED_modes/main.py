# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Dynamical correlator of a small interacting fermionic chain with random
# complex hopping, S(w) = sum_n <GS|A^+|n><n|A|GS> delta(w - E_n + E_0)
# with A = c^+_i c_j. Each peak is an excitation that a spectroscopy sees.
# With mode="ED" there are several submodes, and all must give the same
# spectrum: CVM, INV and ED (exact Lehmann sum).
import numpy as np
from dmrgpy import fermionchain

np.random.seed(1)
n = 4 # number of sites
fc = fermionchain.Fermionic_Chain(n)
h = 0
for i in range(n):
    for j in range(n):
        t = np.random.random() + 1j*np.random.random()
        h = h + t*fc.Cdag[i]*fc.C[j] # random hopping
        h = h + 0.2*np.random.random()*fc.N[i]*fc.N[j] # random interaction
h = h + h.get_dagger()
fc.set_hamiltonian(h)

op = fc.Cdag[0]*fc.C[1] # operator A
name = (op.get_dagger(), op)
es = np.linspace(-0.5, 6.0, 500) # energies
import matplotlib.pyplot as plt
for submode, style in zip(["ED", "CVM", "INV"], ["-", "--", ":"]):
    x, y = fc.get_dynamical_correlator(mode="ED", name=name, submode=submode,
                                       es=es, delta=2e-2)
    plt.plot(x, y.real, style, linewidth=3, label=submode)
plt.xlabel("energy w")
plt.ylabel("dynamical correlator")
plt.legend()
plt.show()
