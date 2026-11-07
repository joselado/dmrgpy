# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Two local moments S=1/2 in a chain of spinful fermion sites, c S c S c.
# Each moment has an exchange K with its two neighbor fermion sites, the
# electrons hop across the moments and have a Hubbard U, at half filling.
# In a weak magnetic field B we compute the magnetization <Sz> of every
# site, with DMRG and with exact diagonalization (ED) of the same chain.
# The exchange is antiferromagnetic, so the moments and the electrons
# next to them point in opposite directions.
import numpy as np
from dmrgpy import mixedchain

mc = mixedchain.Mixed_Spin_Fermion_Chain(["F", "1/2", "F", "1/2", "F"])
mc.maxm = 60 # bond dimension
mc.nsweeps = 20 # number of DMRG sweeps
K, t, U, B = 0.6, 0.5, 2.0, 0.3 # exchange, hopping, Hubbard U and field

hop = 0
for i in [0, 2]: # hopping from fermion site i to i+2, across a moment
    hop = hop + t*(mc.Cdagup[i]*mc.Cup[i+2] + mc.Cdagdn[i]*mc.Cdn[i+2])
h = hop + hop.get_dagger()
for i in [1, 3]: h = h + K*mc.SS(i-1, i) + K*mc.SS(i, i+1) # exchange
for i in [0, 2, 4]: h = h + U*mc.Nup[i]*mc.Ndn[i] - U/2*mc.Ntot[i]
h = h - B*sum(mc.Sz) # magnetic field
mc.set_hamiltonian(h)

szdmrg = [mc.vev(mc.Sz[i], mode="DMRG").real for i in range(5)]
szed = [mc.vev(mc.Sz[i], mode="ED").real for i in range(5)]
print("Max difference DMRG-ED =", np.max(np.abs(np.array(szdmrg)-szed)))
assert np.max(np.abs(np.array(szdmrg)-szed)) < 1e-6

import matplotlib.pyplot as plt
plt.plot(range(5), szed, linewidth=4, label="ED")
plt.plot(range(5), szdmrg, "o--", label="DMRG")
plt.xticks(range(5), ["c", "S", "c", "S", "c"])
plt.xlabel("site")
plt.ylabel("<Sz>")
plt.legend()
plt.show()
