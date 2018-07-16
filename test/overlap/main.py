from __future__ import print_function
import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

spins = [2 for i in range(10)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
e0 = sc.gs_energy(mode="DMRG") # compute the ground state energy
wf = sc.get_gs() # get the ground state as an MPS object
wf1 = wf*2.
wf2 = wf*3.
print(wf1.dot(wf2))
