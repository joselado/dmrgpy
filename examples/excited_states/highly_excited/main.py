# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Excited states at a chosen energy in the middle of the spectrum,
# found with the shift-invert method, which returns the eigenstate
# of H whose energy is closest to a target energy e.
# We use a short S=1/2 chain with random anisotropic couplings
# and plot the energy we recover against the target.
import numpy as np
from dmrgpy import spinchain, mpsalgebra

n = 4 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])
h = 0
for i in range(n-1):
    h = h + np.random.random()*sc.Sx[i]*sc.Sx[i+1]
    h = h + np.random.random()*sc.Sy[i]*sc.Sy[i+1]
    h = h + np.random.random()*sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)

targets = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6] # target energies
found = []
for e in targets:
    (ei, wf) = mpsalgebra.mpsarnoldi(sc, h, e=e, mode="ShiftInv")
    found.append(np.array(ei).flatten()[0])
    print("Target", e, "recovered energy", found[-1])

import matplotlib.pyplot as plt
plt.plot(targets, found, "o-", label="recovered energy")
plt.plot(targets, targets, "--", color="gray", label="target energy")
plt.xlabel("target energy")
plt.ylabel("recovered energy")
plt.legend()
plt.show()
