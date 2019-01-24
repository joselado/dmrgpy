# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
n = 20
spins = [3 for i in range(n)] # spin 1/2 heisenberg chain
#spins = [np.random.randint(2,6) for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
#e = sc.gs_energy() # compute the ground state energy
#print(e)
#e = sc.entropy() # compute the ground state energy
#es = sc.get_excited(mode="DMRG",n=10)
sc.kpmmaxm = 10
es = sc.get_spismj(n=1000,mode="DMRG",i=0,j=0)
#es = sc.get_dos(n=10,mode="DMRG")
#print(es)
#print(np.round(es0-es1,4))
#print("Energy per site",e/n)
#e = sc.gs_energy(mode="aa") # compute the ground state energy
#print("Energy per site",e/n)


