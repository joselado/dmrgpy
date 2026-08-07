# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import spinchain
n = 2
spins = ["S=1/2" for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0 # initialize Hamiltonian
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]

hm = sc.get_full_matrix(h) # get the operator as an sparse matrix
print("Operator as sparse matrix")
print(hm)
print("Operator as dense numpy matrix")
dense = hm.todense()
print(dense)

# Visualize the full Hamiltonian matrix -- sparsity pattern and the
# magnitude of every entry -- instead of only printing the raw arrays.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.2))
ax1.spy(hm, markersize=10)
ax1.set_title("Sparsity pattern")
im = ax2.imshow(np.abs(np.asarray(dense)), cmap="viridis")
ax2.set_title(r"$|H_{ij}|$")
fig.colorbar(im, ax=ax2, fraction=0.046)
fig.suptitle("Heisenberg dimer Hamiltonian matrix (n=%d)" % n)
fig.tight_layout()
plt.show()










