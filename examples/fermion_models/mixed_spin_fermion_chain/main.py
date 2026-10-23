# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Kondo-like model in a mixed chain: a spin-1/2 moment between two
# spinful fermion sites, coupled by the exchange K to each fermion site,
# with hopping t between the fermion sites (across the spin site) and a
# field B. We compute the magnetization <Sz> of the moment against B.
import numpy as np
from dmrgpy import mixedchain

# fermion - spin - fermion, the labels "F" are spinful fermion sites
mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F"], itensor_version=3)
mc.maxm = 60 # bond dimension
mc.nsweeps = 20 # number of DMRG sweeps
K, t = 0.5, 0.4 # exchange and hopping

bs = np.linspace(0.05, 0.4, 6) # magnetic fields
szs = []
for b in bs:
    h = K*mc.SS(0, 1) + K*mc.SS(1, 2) # Kondo exchange
    hop = t*(mc.Cdagup[0]*mc.Cup[2] + mc.Cdagdn[0]*mc.Cdn[2])
    h = h + hop + hop.get_dagger() # hopping crossing the spin site
    h = h + b*(mc.Sz[0] + mc.Sz[1] + mc.Sz[2]) # magnetic field
    mc.set_hamiltonian(h)
    szs.append(mc.vev(mc.Sz[1]).real)

import matplotlib.pyplot as plt
plt.plot(bs, szs, marker="o")
plt.xlabel("magnetic field B")
plt.ylabel("<Sz> of the spin site")
plt.show()
