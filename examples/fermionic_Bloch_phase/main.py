# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain

def get_energy(phi):
    n = 6
    fc = fermionchain.Fermionic_Chain(n) # create the chain
    h = 0
    for i in range(n-1):
        h = h + fc.C[i]*fc.Cdag[i+1]
    h = h + np.exp(1j*phi*np.pi*2)*fc.C[n-1]*fc.Cdag[0]
    h = h + h.get_dagger() # hermitian
    fc.set_hamiltonian(h)
    e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
    return e0


phis = np.linspace(0.,1.,30) # phases
es = [get_energy(phi) for phi in phis]

plt.plot(phis,es,marker="o")
plt.xlabel("$\phi/(2\pi)$")
plt.ylabel("Energy")
plt.show()






