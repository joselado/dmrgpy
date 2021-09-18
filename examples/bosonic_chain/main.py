# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import bosonchain
n = 4
bc = bosonchain.Bosonic_Chain(n) # create the bosonic chain
h = 0
t = 1.0 # hopping
U = 1.0 # Hubbard
mu = 0.0 # chemical potential
for i in range(n-1):  h = h +t*bc.Adag[i]*bc.A[i+1] # hopping
for i in range(n):  h = h + U*(bc.N[i]-1.0)*(bc.N[i]-1.0) # onsite interaction
for i in range(n):  h = h + mu*bc.N[i] # chemical potential
h = h + h.get_dagger()
bc.set_hamiltonian(h)
e = bc.gs_energy() # compute the ground state energy
print("Energy",e)
for i in range(n):
    print()
    print("Site #",i)
    print("Average occupation",np.round(bc.vev(bc.N[i]).real,2))
    print("No occupation",np.round(bc.vev(bc.D0[i]).real,2))
    print("Single occupation",np.round(bc.vev(bc.D1[i]).real,2))
    print("Double occupation",np.round(bc.vev(bc.D2[i]).real,2))
    print("Triple occupation",np.round(bc.vev(bc.D3[i]).real,2))
#print("Energy",e)









