# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Gap of the Hubbard chain at half filling, as a function of U,
# H = sum_{i,s} c^+_{i,s} c_{i+1,s} + h.c. + U sum_i (n_i - 1)^2.
# For U>0 a Mott gap opens, which transport or tunneling spectroscopy measures.
# The gap is E(N+1) + E(N-1) - 2 E(N), computed with DMRG.
import numpy as np
from dmrgpy import fermionchain

ns = 6 # number of spinful sites

def gap(U):
    fc = fermionchain.Fermionic_Chain(2*ns) # sites 2*i and 2*i+1 are spin up and down
    h = 0
    for i in range(ns-1):
        for s in range(2):
            h = h + fc.Cdag[2*i+s]*fc.C[2*(i+1)+s]
    h = h + h.get_dagger()
    for i in range(ns):
        d = fc.N[2*i] + fc.N[2*i+1]
        h = h + U*(d*d - 2*d) # Hubbard term
    fc.set_hamiltonian(h)
    fc.maxm = 40 # bond dimension
    fc.nsweeps = 10 # number of DMRG sweeps
    return fc.get_gap(mode="DMRG")

us = np.linspace(0.0, 10.0, 8)
gs = [gap(u) for u in us]

import matplotlib.pyplot as plt
plt.plot(us, gs, marker="o")
plt.xlabel("Hubbard U")
plt.ylabel("gap")
plt.show()
