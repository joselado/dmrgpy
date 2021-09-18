# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,10,2))
es = []
n = 10
spins = [2 for i in range(n)]
sc = spinchain.Spin_Chain(spins) # create the chain
h = 0
for i in range(n-1): 
    h = h + sc.Sx[i]*sc.Sx[i+1]
    h = h + sc.Sy[i]*sc.Sy[i+1]
    h = h + sc.Sz[i]*sc.Sz[i+1]
sc.set_hamiltonian(h)
pairs = [sc.Sx[i]*sc.Sx[i+1] for i in range(n-1)] # correlators
#pairs = [(n//2,n//2+i) for i in range(1,n//2-1)] # correlators
#pairs = [(n//2,i+1) for i in range(n-1)] # correlators
sc.gs_energy() # compute the correlator between these sites
cs = [sc.vev(p) for p in pairs] # compute the correlator between these sites
import matplotlib.pyplot as plt
plt.plot(range(len(cs)),cs,marker="o") # correlator using DMRG
plt.show()











