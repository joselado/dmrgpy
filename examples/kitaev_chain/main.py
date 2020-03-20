# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import matplotlib.pyplot as plt
from dmrgpy import fermionchain
n = 10
fc = fermionchain.Fermionic_Chain(n,spinful=False) # create the chain

def ft(i,j):
    if i==j: return 0.4
    if abs(j-i)==1: return 1.0 
    return 0.0

def fp(i,j):
    if abs(j-i)==1: 
        return i-j
    return 0.0


fc.set_hoppings(ft) # hoppings
#fc.set_pairings(fp) # pairings
fc.set_hubbard(lambda i,j: 4*abs(i-j)==1) # pairings
e0 = fc.gs_energy(mode="ED") # energy with exact diagonalization
e1 = fc.gs_energy(mode="DMRG") # energy with DMRG
print("Energy with ED",e0)
print("Energy with DMRG",e1)
#exit()

i = n//2
(x,y) = fc.get_dynamical_correlator(mode="ED",i=i,j=i,name="cdc")
(x2,y2) = fc.get_dynamical_correlator(mode="DMRG",i=i,j=i,name="cdc")

import matplotlib.pyplot as plt

plt.plot(x,y.real,c="blue",label="ED")
plt.plot(x2,y2.real,c="red",label="DMRG")
plt.legend()
plt.show()


