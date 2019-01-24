# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 6
sc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
def ft(i,j):
    if abs(j-i)==1: return 1.0 #+ np.random.random()
    return 0.0


sc.set_hoppings(ft) # hoppings
#sc.set_fields(fb) # Add local magnetic fields
sc.maxm = 20 # maxm
sc.nsweeps = 15 # maxm


import time
e0 = sc.gs_energy() # compute ground state energy with DMRG
e1 = sc.gs_energy_free() # compute ground state energy for free electrons
print(sc.spinful)
print(e0,e1)
#exit()
pairs = [(0,i) for i in range(n)]
out = sc.get_correlator(pairs=pairs)
out2 = sc.get_correlator_free(pairs=pairs)
x = np.array(range(len(out2))) # positions
plt.plot(x,out,c="blue",label="DMRG")
plt.scatter(x,out2,marker="o",c="red",label="ED",s=40)
plt.legend()
#(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=1,delta=0.02,
#                                           name="cdc")
plt.show()


