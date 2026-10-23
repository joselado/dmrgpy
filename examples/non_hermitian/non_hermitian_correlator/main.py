# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Interacting fermion chain with gain and loss, so H is non-Hermitian.
# We take the two states with the lowest Re(E), 0 and 1, and compute
# |<0|c_i^dag|1>|^2 + |<0|c_i|1>|^2 on every site i.
# The states come from ED and from DMRG with a non-Hermitian Arnoldi.
import numpy as np
from dmrgpy import fermionchain

ns = 8 # number of sites
fc = fermionchain.Fermionic_Chain(ns)
eta = 1.0 # strength of gain and loss
gl = [1j*eta, -1j*eta, -1j*eta, 1j*eta] # pattern of gain and loss

h = 0
for i in range(ns-1): # hopping and interaction
    h = h + fc.Cdag[i]*fc.C[i+1] + fc.Cdag[i+1]*fc.C[i]
    h = h + (fc.N[i]-0.5)*(fc.N[i+1]-0.5)
for i in range(ns): # gain and loss
    h = h + gl[i % 4]*fc.Cdag[i]*fc.C[i]
fc.set_hamiltonian(h)

(e0, wf0) = fc.get_excited_states(mode="ED", n=2)
(e1, wf1) = fc.get_excited_states(mode="DMRG", n=2, recursive=False)

def proj(wfs): # weight of the transition between the two states
    return [abs(wfs[0].dot(fc.Cdag[i]*wfs[1]))**2
            + abs(wfs[0].dot(fc.C[i]*wfs[1]))**2 for i in range(ns)]

import matplotlib.pyplot as plt
plt.plot(range(ns), proj(wf0), "o", markersize=12, label="ED")
plt.plot(range(ns), proj(wf1), "x", markersize=12, label="DMRG")
plt.xlabel("site i")
plt.ylabel("correlator")
plt.legend()
plt.show()
