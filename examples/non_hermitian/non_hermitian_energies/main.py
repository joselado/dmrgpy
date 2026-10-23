# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Lowest energies of a non-Hermitian Heisenberg chain, where a staggered
# imaginary field i*0.6*(-1)^i Sz_i (gain and loss) makes the energies
# complex. We compute them with DMRG and compare with exact
# diagonalization (ED), and plot them in the complex plane.
import numpy as np
from dmrgpy import spinchain

n = 6 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian plus the imaginary staggered field
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
for i in range(n):
    h = h + (-1)**i*0.6j*sc.Sz[i]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

e1 = np.array(sc.get_excited(mode="ED", n=4))
e2 = np.array(sc.get_excited(mode="DMRG", n=4))
print("Energies ED  :", np.round(e1, 3))
print("Energies DMRG:", np.round(e2, 3))

import matplotlib.pyplot as plt
plt.scatter(e1.real, e1.imag, s=200, marker="o", label="ED")
plt.scatter(e2.real, e2.imag, s=80, marker="x", label="DMRG")
plt.xlabel("Re(E)")
plt.ylabel("Im(E)")
plt.legend()
plt.show()
