# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Low energy manifold of a Heisenberg chain. We write the Hamiltonian H
# and the total spin S^2 as matrices in the basis of the lowest excited
# states. H is diagonal by construction, and S^2 tells which states
# belong to the same total spin multiplet (S=0, 1, 2, ...).
import numpy as np
from dmrgpy import spinchain
from dmrgpy.effectivehamiltonian import get_representation

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

# total spin operator S^2
s2 = 0
for s in [sc.Sx, sc.Sy, sc.Sz]:
    stot = sum(s[i] for i in range(n))
    s2 = s2 + stot*stot

es, ws = sc.get_excited_states(n=4) # lowest states
hm = get_representation(ws, h) # H in that basis
sm = get_representation(ws, s2) # S^2 in that basis
print("Energies", np.round(np.linalg.eigvalsh(hm), 3))
print("Eigenvalues of S^2", np.round(np.linalg.eigvalsh(sm), 3))

import matplotlib.pyplot as plt
plt.imshow(sm.real)
plt.colorbar(label="<i|S^2|j>")
plt.xlabel("state j")
plt.ylabel("state i")
plt.show()
