import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 4
spins = [2] + [1 for i in range(n)] + [2] + [1] # spin 1/2 plus fermionic sites
#spins = [2,2]
sc = spinchain.Many_Body_Hamiltonian(spins) # create the spin chain
def ft(i,j):
    if i>0: # fermionic sites
      if abs(i-j)==2: return 1.0
    return 0.0
sc.set_hoppings(ft) # hoppings


def fj(i,j):
    if i==0 and j==i+1: return 1.0
    else: return 0.0

sc.set_exchange(fj) # set exchange couplings

def fb(i):
    if i==0: return [0.0,0.0,0.01]
    else: return [0.0,0.0,0.0]

sc.set_fields(fb) # Add local magnetic fields

#sc.kpmmaxm = 20 # KPM maxm
import time
print(sc.gs_energy())
#exit()
#(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=0,delta=0.02)


