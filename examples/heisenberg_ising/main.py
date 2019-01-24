# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
n = 40
# first one will be a Ising spin
spins = [3] + [3 for i in range(n)] # spin 1/2 chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def fj(i,j): 
    if j==(i+1):
        if i==0: return [[0.,0.,0.],[0.,0.,0.],[0.,0.,1.]] # ising
        else: return 1.0 # heisenberg
    return 0.0
sc.set_exchange(fj)
e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy


