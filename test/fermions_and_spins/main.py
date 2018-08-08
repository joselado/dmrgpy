import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain
import spinchain
import random


sites = [1 for i in range(6)]
for i in random.sample(range(6), 3):
    sites[i] = 2

print(sites)

sc = spinchain.Spin_Hamiltonian(sites) # create the spin chain

def fh(i,j):
    if sites[i]==1 and sites[j]==1: 
        if i!=j: return 1.0
    return 0.0

def fj(i,j): return 0.3

sc.set_hoppings(fh) # Add hopping
sc.set_exchange(fj) # add exchange coupling

e0 = min([sc.gs_energy() for i in range(10)])
print(e0)
