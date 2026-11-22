# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# The ground state of a superconducting chain without interactions
# (Kitaev chain) written as an MPS. With pairing the state has anomalous
# correlations <C_i C_j>, which we know exactly from the Bogoliubov-de
# Gennes (BdG) matrix, so we can check the MPS on a long chain.
import numpy as np
from dmrgpy import fermionchain, gaussianmps

n = 40 # number of sites
fc = fermionchain.Fermionic_Chain(n)
fc.maxm = 60 # bond dimension of the MPS

# H = sum_ij h_ij Cdag_i C_j + 1/2 sum_ij (D_ij Cdag_i Cdag_j + h.c.)
mu = 2.5 # chemical potential, trivial superconductor for |mu| > 2
h = -mu*np.eye(n)
D = np.zeros((n,n))
for i in range(n-1):
    h[i,i+1] = h[i+1,i] = -1.0 # hopping
    D[i,i+1], D[i+1,i] = 1.0, -1.0 # p-wave pairing, D_ij = -D_ji

wf = fc.get_gaussian_state(h=h, delta=D, eigval_cutoff=1e-12)

# exact correlations <gamma_a^dag gamma_b>, gamma = (C_0, Cdag_0, C_1, ...)
G = gaussianmps.bdg_correlation_matrix(fc, h, D)

i0 = n//2 # a site in the bulk
pmps = [fc.vev(fc.C[i0]*fc.C[j], wf=wf).real for j in range(n)]
pex = [G[2*i0+1,2*j].real for j in range(n)] # <C_i0 C_j>
print("Max difference MPS-exact =", np.max(np.abs(np.array(pmps)-pex)))
assert np.max(np.abs(np.array(pmps)-pex)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(range(n), pex, linewidth=4, label="exact")
plt.plot(range(n), pmps, "--", linewidth=2, label="MPS")
plt.xlabel("site j")
plt.ylabel("<C_i C_j>, i = "+str(i0))
plt.legend()
plt.show()
