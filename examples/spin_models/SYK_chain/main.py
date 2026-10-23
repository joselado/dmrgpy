# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Chain of Sachdev-Ye-Kitaev (SYK) dots with random quartic couplings,
# H = sum J_ijkl c_i c_j c^+_k c^+_l + h.c., inside each dot of n fermionic
# sites and between neighboring dots. There is no hopping, so the
# interaction is the only energy scale. We compute the density of
# fermions with DMRG and compare with exact diagonalization (ED).
import numpy as np
from dmrgpy import fermionchain

np.random.seed(1)
n = 3 # sites per dot
L = 3 # number of dots
fc = fermionchain.Fermionic_Chain(n*L)
C, Cdag = fc.C, fc.Cdag

h = 0
for l in range(L):
    for i in range(n):
        for j in range(i+1, n):
            for k in range(n):
                for m in range(k+1, n):
                    A = C[i+l*n]*C[j+l*n]
                    h = h + (np.random.random()-0.5)*A*Cdag[k+l*n]*Cdag[m+l*n]
                    if l < L-1: # coupling to the next dot
                        h = h + (np.random.random()-0.5)*A*Cdag[k+(l+1)*n]*Cdag[m+(l+1)*n]
h = h + h.get_dagger()
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

d1 = fc.get_density(mode="DMRG")
d2 = fc.get_density(mode="ED")
print("Max difference =", np.max(np.abs(np.array(d1)-d2)))

import matplotlib.pyplot as plt
plt.plot(range(n*L), d2, "o", markersize=10, label="ED")
plt.plot(range(n*L), d1, "x", markersize=10, label="DMRG")
plt.xlabel("site i")
plt.ylabel("density <n_i>")
plt.legend()
plt.show()
