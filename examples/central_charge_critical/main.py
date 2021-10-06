# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')
import numpy as np

from dmrgpy import spinchain
n = 100 # number of sites
spins = ["S=1/2" for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Chain(spins) # create the spin chain
h = 0 # initialize
for i in range(n-1): h = h + sc.Sz[i]*sc.Sz[i+1] # Ising coupling
for i in range(n): h = h + 0.5*sc.Sx[i] # transverse field
sc.set_hamiltonian(h) # set the Hamiltonian
sc.maxm = 200 # increase bond dimension for a critical system
wf = sc.get_gs() # compute ground state
print("Central charge",wf.get_central_charge()) # compute central charge
