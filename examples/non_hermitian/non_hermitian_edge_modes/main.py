# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Interacting non-Hermitian chain made of tetramers with onsite gain and loss,
# +i, -i, -i, +i, hopping 1 and interaction V. The complex many-body
# energies E-E0 show the edge zero modes of the chain. We follow the
# lowest 10 states as a function of V, using exact diagonalization.
import numpy as np
from dmrgpy import fermionchain

ns = 8 # number of sites (two tetramers)
fc = fermionchain.Fermionic_Chain(ns)

def spectrum(V):
    h = 0
    for i in range(ns-1):
        h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i] # hopping
        h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5) # interaction
    for i in range(ns):
        gain = 1j*(-1)**(i//2 + i) # +i,-i,-i,+i pattern
        h = h + gain*fc.Cdag[i]*fc.C[i]
    fc.set_hamiltonian(h)
    es, wfs = fc.get_excited_states(mode="ED", n=10)
    return es - es[0]

import matplotlib.pyplot as plt
for V in np.linspace(0.0, 1.0, 20):
    es = spectrum(V)
    plt.scatter(es*0+V, es.real, c="black", s=8)
    plt.scatter(es*0+V, es.imag, c="red", s=8)
plt.xlabel("interaction V")
plt.ylabel("E-E0 (black: real part, red: imaginary part)")
plt.show()
