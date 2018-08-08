import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain
import spinchain




n = 5 # number of sites (including fermions and spins)
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
    if i%2==0 and j%2==1: # first one is a fermion, second one a spin
        if abs(i-j)==1: return 0.0
    return 0.0


sc.set_hoppings(fh) # Add hopping
sc.set_exchange(fj) # add exchange coupling


print("DMRG=",sc.gs_energy())


# Test 

m = np.matrix([[fh(i,j) for i in range(n)] for j in range(n)])
import scipy.linalg as lg
e = lg.eigvalsh(m)
print("ED=",2*np.sum(e[e<0.]))

