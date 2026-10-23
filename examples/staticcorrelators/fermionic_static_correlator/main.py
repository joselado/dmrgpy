# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Correlator <Cdag_0 C_j> of a spinful tight-binding chain with chemical
# potential mu, summed over spin. It measures how a fermion can hop from
# site 0 to site j, and it decays with the distance in the ground state.
import numpy as np
from dmrgpy import fermionchain

n = 10 # number of sites
mu = -1.0 # chemical potential
fc = fermionchain.Spinful_Fermionic_Chain(n)
h = 0
for i in range(n-1): # hopping
    h = h + fc.Cdagup[i]*fc.Cup[i+1] + fc.Cdagdn[i]*fc.Cdn[i+1]
h = h + h.get_dagger()
for i in range(n): # chemical potential
    h = h + mu*(fc.Cdagup[i]*fc.Cup[i] + fc.Cdagdn[i]*fc.Cdn[i])
fc.set_hamiltonian(h)
fc.maxm = 40 # bond dimension
fc.nsweeps = 10 # number of DMRG sweeps

cs = [fc.vev(fc.Cdagup[0]*fc.Cup[j]) + fc.vev(fc.Cdagdn[0]*fc.Cdn[j]) for j in range(n)]

import matplotlib.pyplot as plt
plt.plot(range(n), np.real(cs), marker="o")
plt.xlabel("site j")
plt.ylabel("<Cdag_0 C_j>")
plt.show()
