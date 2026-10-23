# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Many-body spectrum of an interacting chain with non-reciprocal hopping,
# t_{i,i+1} = 1 + gamma and t_{i+1,i} = 1 - gamma, plus a nearest-neighbor
# interaction V. The non-Hermitian hopping gives a complex spectrum, and the
# energies change strongly with gamma (the non-Hermitian skin effect).
# We show the real and the imaginary part of the lowest energies (exact diagonalization).
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
V = -0.3 # interaction
fc = fermionchain.Fermionic_Chain(n)

def spectrum(gamma):
    h = 0
    for i in range(n-1):
        h = h + V*(fc.N[i]-0.5)*(fc.N[i+1]-0.5)
        h = h + (1.0+gamma)*fc.Cdag[i]*fc.C[i+1]
        h = h + (1.0-gamma)*fc.Cdag[i+1]*fc.C[i]
    fc.set_hamiltonian(h)
    es = fc.get_excited(mode="ED", n=10)
    return es - es[0].real

import matplotlib.pyplot as plt
for g in np.linspace(0.0, 1.5, 20):
    es = spectrum(g)
    plt.scatter(es*0+g, es.real, c="black", s=8)
    plt.scatter(es*0+g, es.imag, c="red", s=8)
plt.xlabel("non-reciprocity gamma")
plt.ylabel("E (black: real part, red: imaginary part)")
plt.show()
