import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
from dmrgpy import spinchain




ns = np.array(range(4,30,2))
es = []
n = 80
spins = [2 for i in range(n)] # S=1/2 chain
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
def fj(i,j): # function to define the exchange couplings
    if 0.1<abs(i-j)<1.1:
        if i==39 and j==40: return 0.2
        if i==40 and j==39: return 0.2
        return 1.0
    return 0.0
sc.set_exchange(fj)

# Pairs on sites on which to compute the correlator S_i S_j
pairs = [(i,i+1) for i in range(n-1)] # correlators, first neighbors
sc.gs_energy() # compute the ground state energy
cs = sc.correlator(pairs) # compute the correlator between these sites

# now plot the correlator
import matplotlib.pyplot as plt

plt.plot(range(len(cs)),cs,marker="o") # correlator using DMRG

plt.show()
