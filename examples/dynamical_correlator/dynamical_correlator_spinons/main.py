# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Spinon chain: a Heisenberg chain written with spinon operators, which
# are fermions with one particle per site. We compute the spectral
# function A(w) = sum_n |<n|Cdag_up,0|GS>|^2 delta(w - E_n + E_0),
# which is the weight a tunneling electron or a spinon probe would see.
import numpy as np
from dmrgpy import fermionchain

n = 6 # number of sites
sc = fermionchain.Spinon_Chain(n)

# Heisenberg Hamiltonian, H = sum_i S_i . S_{i+1}
h = 0
for i in range(n-1):
    h = h + sc.Sx[i]*sc.Sx[i+1] + sc.Sy[i]*sc.Sy[i+1] + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
sc.maxm = 40 # bond dimension
sc.nsweeps = 10 # number of DMRG sweeps

es = np.linspace(-0.5, 6.0, 800) # energies
name = (sc.Cdagup[0], sc.Cup[0]) # the two operators of the correlator
(x, y) = sc.get_dynamical_correlator(mode="DMRG", name=name, es=es, delta=0.05)

import matplotlib.pyplot as plt
plt.plot(x, y.real, label="DMRG")
plt.xlabel("energy w")
plt.ylabel("A(w), first site")
plt.legend()
plt.show()
