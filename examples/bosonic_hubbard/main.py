# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import bosonchain
n = 3
bc = bosonchain.Bosonic_Chain(n) # create the chain

h = 0 # initialize
for i in range(n-1): # hopping
    h = h + bc.get_operator("Adag",i)*bc.get_operator("A",i+1)
    h = h + bc.get_operator("Adag",i+1)*bc.get_operator("A",i)

for i in range(n): # hopping
    den = bc.get_operator("Adag",i)*bc.get_operator("A",i)
    h = h + 0.2*den*den

bc.set_hamiltonian(h) # setup hamiltonian
print(bc.gs_energy(mode="ED"))

ds = [bc.vev(bc.get_operator("N",i),mode="ED").real for i in range(n)]
print(ds)









