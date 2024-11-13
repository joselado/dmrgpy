# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # generate a Heisenberg Hamiltonian
for i in range(n-1):
    h = h +sc.Sx[i]*sc.Sx[i+1]
    h = h +sc.Sy[i]*sc.Sy[i+1]
    h = h +sc.Sz[i]*sc.Sz[i+1]

# the dagger of an operator con be computed as
# B = A.get_dagger()

sc.set_hamiltonian(h) # set the Hamiltonian
wf0 = sc.get_gs(mode="DMRG") # compute ground state
wf0i = 1j*wf0 # imaginary unit times ground state


wf1 = wf0i.get_conjugate() # compute the conjugate

print("This should be 1")
print(wf0.dot(wf0))
print("This should be -i")
print(wf0.dot(wf1))

