import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 10
spins = [1 for i in range(n)] # spin 1/2 plus fermionic sites
#spins = [2,2]
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def ft(i,j):
    if i==j: return 1.0
    if abs(j-i)==1: return 1.0
    return 0.0
sc.set_hoppings(ft) # hoppings

#sc.set_fields(fb) # Add local magnetic fields

#sc.kpmmaxm = 20 # KPM maxm
import time
print(sc.gs_energy())
#exit()
(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=1,delta=0.02,
                                           name="cdc")


