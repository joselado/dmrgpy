# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
n = 39
U = 0.0
spins = [1 for i in range(n)] + [2] # fermionic sites plus 1/2
#spins = [2,2]
sc = spinchain.Many_Body_Hamiltonian(spins) # create the spin chain
def ft(i,j):
    if i<n and j<n: # fermionic sites
      if abs(i-j)==1: return 1.0
      if i==j: return -U
    return 0.0
sc.set_hoppings(ft) # hoppings
def fj(i,j):
    if i==(n//2) and j==n: return 0.5 # AF coupling between impurity and site
    else: return 0.0
sc.set_exchange(fj) # set exchange couplings
def fb(i):
    if i==n: return [0.0,0.0,0.0001] # in the spin
    else: return [0.0,0.0,0.0]
def fu(i,j):
    if i==j and i<5: return U/2. # small Hubbard
    return 0.0
sc.set_hubbard(fu)
sc.set_fields(fb) # Add local magnetic fields
#sc.kpmmaxm = 20 # KPM maxm
import time
print(sc.gs_energy())
#exit()
#(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=0,delta=0.02)


