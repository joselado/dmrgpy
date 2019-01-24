# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
sc = fermionchain.Fermionic_Hamiltonian(n) # create the chain
def ft(i,j):
#    if i==j: return 1.0
    if abs(j-i)==1: return 1.0 #+ np.random.random()
#    if i==j: return np.random.random()
#    if i==j: return 1.1
    return 0.0
sc.set_hoppings(ft) # hoppings
#sc.set_fields(fb) # Add local magnetic fields
#sc.kpmmaxm = 20 # KPM maxm
import time
e0 = sc.gs_energy() # compute ground state energy with DMRG
e1 = sc.gs_energy_free() # compute ground state energy for free electrons
print(e0,e1)
#exit()
pairs = [(0,i) for i in range(n)]
out = sc.correlator(pairs)
out2 = sc.correlator_free(pairs)
plt.plot(range(n),out,c="blue",label="DMRG")
plt.scatter(range(n),out2,marker="o",c="red",label="ED",s=40)
#(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=1,delta=0.02,
#                                           name="cdc")
plt.show()
