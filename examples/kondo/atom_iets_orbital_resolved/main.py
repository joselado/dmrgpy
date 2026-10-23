# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin-flip IETS spectrum of a d^2 atom (V III) with a crystal field,
# resolved by the d orbital that the STM tip couples to (`iorb`).
# The t2g orbitals sit below the eg ones, so each orbital gives a different
# peak weight. Peaks at the spin-flip energies are what STM measures in dI/dV.
import numpy as np
from dmrgpy import atom

orbs = ["dxy", "dyz", "dz2", "dxz", "dx2y2"] # d orbitals
tij = np.diag([0.0, 0.05, 0.6, 0.05, 0.6]) # crystal field in eV
fc = atom.generate_atom(orbs=orbs, tij=tij, Ne=2, soc=0.03, U=4., J=0.5)

es = np.linspace(0., 0.1, 301) # bias window in eV
import matplotlib.pyplot as plt
for iorb in range(len(orbs)):
    x, y = atom.get_spinflip(fc, es=es, iorb=iorb, delta=1e-3, dex=1e-3)
    plt.plot(x*1e3, y, label=orbs[iorb])
plt.xlabel("bias (meV)")
plt.ylabel("spin-flip signal")
plt.legend()
plt.show()
