# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Two spins coupled to a chain of four bosons with hopping t and
# Hubbard repulsion U, H_int = (Sx_0 + Sx_1)(b_2 + b_2^dagger).
# We compute the ground state and the occupation of each boson site,
# the average <N> and the probability P(n) of having n bosons.
import numpy as np
from dmrgpy import bosonchain

ns = 2 # number of spins
nb = 4 # number of bosons
bc = bosonchain.SpinBoson_Chain(["S=1/2"]*ns + ["B"]*nb)
t = 1.0 # hopping
U = 1.0 # Hubbard interaction
h = 0
for i in range(ns, ns+nb-1): # boson hopping
    h = h + t*bc.Adag[i]*bc.A[i+1]
for i in range(ns, ns+nb): # boson interaction
    h = h + U*(bc.N[i]-1.0)*(bc.N[i]-1.0)
h = h + bc.Sx[0]*bc.Sx[1] + bc.Sy[0]*bc.Sy[1] + bc.Sz[0]*bc.Sz[1] # spin-spin
h = h + (bc.Sx[0] + bc.Sx[1])*(bc.A[2] + bc.Adag[2]) # spin-boson
h = h + h.get_dagger() # make it Hermitian
bc.set_hamiltonian(h)
bc.maxm = 20 # bond dimension
bc.nsweeps = 10 # number of DMRG sweeps
print("Energy", bc.gs_energy())

sites = range(ns, ns+nb)
navg = [bc.vev(bc.N[i]).real for i in sites]
probs = [[bc.vev(D[i]).real for i in sites] for D in [bc.D0, bc.D1, bc.D2]]

import matplotlib.pyplot as plt
plt.plot(sites, navg, "o-", label="<N>")
for k in range(3):
    plt.plot(sites, probs[k], "o--", label="P(n=%d)" % k)
plt.xlabel("boson site")
plt.ylabel("occupation")
plt.legend()
plt.show()
