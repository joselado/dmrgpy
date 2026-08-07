# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 2
spins = ["S=1" for i in range(n)] # spin 1/2 heisenberg chain

def get_spectrum(c0):
    sc = spinchain.Spin_Chain(spins) # create the spin chain
    h = 0
    h = h +sc.Sx[0]*sc.Sx[1]
    h = h +sc.Sy[0]*sc.Sy[1]
    h = h +sc.Sz[0]*sc.Sz[1]
    h = -h
    S111 = sc.Sx[0] + sc.Sy[0] + sc.Sz[0]
    L111 = sc.Sx[1] + sc.Sy[1] + sc.Sz[1]
    h = h + c0*S111*S111
    sc.set_hamiltonian(h)
    return sc.get_excited(n=9,mode="ED")

c0 = -0.1
es = get_spectrum(c0)
print(es)

# sweep the single-ion anisotropy c0 and plot how the 9 lowest energy
# levels of the S=1 dimer evolve/split as a function of it
c0s = np.linspace(-0.5,0.5,7)
spectra = np.array([get_spectrum(c) for c in c0s]) # shape (len(c0s), 9)

for level in range(spectra.shape[1]):
    plt.plot(c0s,spectra[:,level],marker="o")
plt.xlabel("c0")
plt.ylabel("Energy levels")
plt.show()
