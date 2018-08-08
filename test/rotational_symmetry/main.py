import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 5
spins = [1 for i in range(n)] +[2] # spin 1/2 plus fermionic sites
#spins = [2,2]
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def ft(i,j):
    if abs(i-j)==1: return 1.0
    return 0.0
sc.set_hoppings(ft) # hoppings


def fj(i,j):
    if i==0 and j==1: return 1.0
    else: return 0.0

sc.set_exchange(fj) # set exchange couplings

def fb(i):
    if i==n:
        d = np.random.random(3)
#        d[0] = 0.0
#        return [0.0,0.0,0.1]
        return d/np.sqrt(d.dot(d))
    else: return [0.0,0.0,0.0]



#sc.set_fields(fb) # Add local magnetic fields

#sc.kpmmaxm = 20 # KPM maxm
import time
print(sc.gs_energy())
#exit()
#(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=0,j=0,delta=0.02)


