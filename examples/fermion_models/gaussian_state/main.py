# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The Fermi sea of a tight-binding chain written as an MPS (Fishman-White
# circuit). A free-fermion state is fixed by its correlation matrix
# <Cdag_i C_j>, which we know exactly from the single-particle orbitals,
# so we can check the MPS at a size far beyond exact diagonalization.
# Such a state is a good start for DMRG of an interacting model.
import numpy as np
from dmrgpy import fermionchain

n = 40 # number of sites
nf = 20 # number of fermions, half filling
fc = fermionchain.Fermionic_Chain(n)
fc.maxm = 100 # bond dimension of the MPS

# hopping matrix, H = sum_ij h_ij Cdag_i C_j
h = np.zeros((n,n))
for i in range(n-1):
    h[i,i+1] = h[i+1,i] = -1.0
fc.set_hamiltonian(sum(h[i,j]*fc.Cdag[i]*fc.C[j] for i in range(n) for j in range(n) if h[i,j]!=0))

wf = fc.get_gaussian_state(h=h, nf=nf, eigval_cutoff=1e-12) # the Fermi sea as an MPS

# exact correlation matrix from the lowest nf orbitals
e,u = np.linalg.eigh(h)
cm = u[:,:nf].conj() @ u[:,:nf].T

i0 = n//2 # a site in the bulk
cmps = [fc.vev(fc.Cdag[i0]*fc.C[j], wf=wf).real for j in range(n)]
cex = [cm[i0,j].real for j in range(n)]
print("Max difference MPS-exact =", np.max(np.abs(np.array(cmps)-cex)))
print("Energy MPS =", fc.vev(fc.hamiltonian, wf=wf).real, " exact =", np.sum(e[:nf]))
assert np.max(np.abs(np.array(cmps)-cex)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(range(n), cex, linewidth=4, label="exact")
plt.plot(range(n), cmps, "--", linewidth=2, label="MPS")
plt.xlabel("site j")
plt.ylabel("<Cdag_i C_j>, i = "+str(i0))
plt.legend()
plt.show()
