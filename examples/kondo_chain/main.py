# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import fermionchain
import spinchain
n = 4 # number of sites (including fermions and spins)
sites = [] # list with the sites
for i in range(n):
    sites.append(1) # add fermion
for i in range(n):
    sites.append(2) # add spin
sc = spinchain.Spin_Hamiltonian(sites) # create the spin chain
def fh(i,j): # hopping
    if abs(i-j)==1: return 1.0
    return 0.0
def fj(i,j):
    if abs(i-j)==n: return 1.0 # coupling between fermion and spin
    if i>=n and j>=n and abs(i-j)==1: return 0.0 # coupling between spins
    return 0.0
sc.set_hoppings(fh) # Add hopping
sc.set_exchange(fj) # add exchange coupling
print("DMRG=",sc.gs_energy())
# Test 
m = np.matrix([[fh(i,j) for i in range(n)] for j in range(n)])
import scipy.linalg as lg
e = lg.eigvalsh(m)
print("ED=",2*np.sum(e[e<0.]))
