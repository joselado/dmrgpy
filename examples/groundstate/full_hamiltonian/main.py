# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Full matrix of the Hamiltonian of a Heisenberg dimer,
# H = Sx_0 Sx_1 + Sy_0 Sy_1 + Sz_0 Sz_1, in the basis of the four spin
# states. It is useful to check small systems and to diagonalize by hand.
import numpy as np
from dmrgpy import spinchain

n = 2 # number of sites
sc = spinchain.Spin_Chain(["S=1/2" for i in range(n)])

# Heisenberg Hamiltonian
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]

hm = sc.get_full_matrix(h) # operator as a sparse matrix
dense = np.array(hm.todense()) # operator as a dense matrix
print(dense)
print("Eigenvalues =", np.linalg.eigvalsh(dense))

import matplotlib.pyplot as plt
plt.imshow(np.abs(dense))
plt.colorbar(label="|H_ij|")
plt.xlabel("column")
plt.ylabel("row")
plt.show()
