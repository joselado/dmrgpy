# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

# Bosonic chain with hopping t and onsite interaction U (N-1)^2, which
# favors one boson per site. We compute the average occupation <N_i> and
# the probability P(n) of finding n bosons on site i in the ground state.
# Occupations are what a quantum gas microscope measures.
import numpy as np
from dmrgpy import bosonchain

n = 4 # number of sites
t = 1.0 # hopping
U = 1.0 # onsite interaction
bc = bosonchain.Bosonic_Chain(n)
h = 0
for i in range(n-1): h = h + t*bc.Adag[i]*bc.A[i+1] # hopping
for i in range(n): h = h + U*(bc.N[i]-1.0)*(bc.N[i]-1.0) # interaction
h = h + h.get_dagger()
bc.set_hamiltonian(h)
print("Ground state energy =", bc.gs_energy())

sites = range(n)
occ = [bc.vev(bc.N[i]).real for i in sites] # average occupation
p0 = [bc.vev(bc.D0[i]).real for i in sites] # probability of n=0
p1 = [bc.vev(bc.D1[i]).real for i in sites] # probability of n=1
p2 = [bc.vev(bc.D2[i]).real for i in sites] # probability of n=2

import matplotlib.pyplot as plt
plt.plot(sites, occ, "o-", label="<N>")
plt.plot(sites, p0, "s-", label="P(n=0)")
plt.plot(sites, p1, "^-", label="P(n=1)")
plt.plot(sites, p2, "v-", label="P(n=2)")
plt.xlabel("site i")
plt.ylabel("occupation or probability")
plt.legend()
plt.show()
