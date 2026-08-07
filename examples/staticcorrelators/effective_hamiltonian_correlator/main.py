# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h0 = 0
for i in range(n-1): # nearest-neighbor Heisenberg exchange
    h0 = h0 + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h0)
e = sc.gs_energy() # compute the ground state energy
print("Ground state energy",e)
print(sc.get_excited(n=4))

# Represent the Hamiltonian and the total spin S^2 in the low-energy
# excited-state manifold (the "effective Hamiltonian" correlator matrices):
# h is the Hamiltonian itself in that basis (diagonal, by construction),
# b is S^2, whose eigenvectors group the low-energy states into total-spin
# multiplets.
from dmrgpy.effectivehamiltonian import get_representation
ne = 4
es,ws = sc.get_excited_states(n=ne)
sxtot = sum(sc.Sx[i] for i in range(n))
sytot = sum(sc.Sy[i] for i in range(n))
sztot = sum(sc.Sz[i] for i in range(n))
s2 = sxtot*sxtot + sytot*sytot + sztot*sztot
h = get_representation(ws,sc.hamiltonian)
b = get_representation(ws,s2)

import scipy.linalg as lg
print(np.round(lg.eigvalsh(b),3))
es_b,vs = lg.eigh(-b)
vs = vs.transpose()
print(vs[0].real)
print(np.round(lg.eigvalsh(h),3))

# heatmaps of the Hamiltonian and S^2 representations in the low-energy
# excited-state manifold
import matplotlib.pyplot as plt
fig,axes = plt.subplots(1,2,figsize=(9,4))
im0 = axes[0].imshow(h.real,aspect="auto")
axes[0].set_title("Hamiltonian representation")
axes[0].set_xlabel("Excited state j")
axes[0].set_ylabel("Excited state i")
fig.colorbar(im0,ax=axes[0])
im1 = axes[1].imshow(b.real,aspect="auto")
axes[1].set_title("S^2 representation")
axes[1].set_xlabel("Excited state j")
axes[1].set_ylabel("Excited state i")
fig.colorbar(im1,ax=axes[1])
plt.tight_layout()
plt.show()
