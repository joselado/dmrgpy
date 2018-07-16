from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

ns = [10,100,1000]
spins = [2 for i in range(100)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
