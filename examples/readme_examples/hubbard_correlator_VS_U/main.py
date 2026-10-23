# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spin correlator <Sz_0 Sz_i> of the half-filled Hubbard chain,
# H = sum_i (c_i^dagger c_{i+1} + h.c.) + U sum_i (n_up - 1/2)(n_dn - 1/2).
# Increasing U builds up antiferromagnetic correlations, which neutron
# scattering or spin-polarized STM can detect. We use DMRG.
import numpy as np
from dmrgpy import fermionchain

n = 12 # number of sites
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger()
hU = 0
for i in range(n): # Hubbard term
    hU = hU + (fc.Nup[i]-.5)*(fc.Ndn[i]-.5)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

import matplotlib.pyplot as plt
for U in [0.0, 2.0, 4.0]:
    fc.set_hamiltonian(h + U*hU)
    zz = [fc.vev(fc.Sz[0]*fc.Sz[i]).real for i in range(n)]
    plt.plot(range(n), zz, "o-", label="U=%s" % U)
plt.xlabel("site i")
plt.ylabel("<Sz_0 Sz_i>")
plt.legend()
plt.show()
