import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import fermionchain
import fermionchain




n = 4 # number of sites (including fermions and spins)
sites = [] # list with the sites
for i in range(n):
    sites.append(1) # add fermion
#for i in range(n):
#    sites.append(2) # add spin

sc = fermionchain.Fermionic_Hamiltonian(n) # create the spin chain

mr = np.random.random((n,n)) + 1j*np.random.random((n,n))
mr = np.matrix(mr)
mr = mr + mr.H

def fh(i,j): return mr[i,j] # hopping

def fj(i,j):
    if i%2==0 and j%2==1: # first one is a fermion, second one a spin
        if abs(i-j)==1: return 0.0
    return 0.0


sc.set_spinful_hoppings(mr) # Add hopping
#sc.set_exchange(fj) # add exchange coupling


print("DMRG=",sc.gs_energy())


# Test 

print("ED=",sc.gs_energy_free())

