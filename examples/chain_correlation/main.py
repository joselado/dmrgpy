import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain




ns = np.array(range(4,30,2))
es = []
n = 10
spins = [2 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
pairs = [(i,i+1) for i in range(n-1)] # correlators
#pairs = [(n//2,n//2+i) for i in range(1,n//2-1)] # correlators
#pairs = [(n//2,i+1) for i in range(n-1)] # correlators
sc.gs_energy() # compute the correlator between these sites
cs = sc.correlator(pairs) # compute the correlator between these sites

exit()  
import matplotlib.pyplot as plt

plt.plot(range(len(cs)),cs) # correlator using DMRG

plt.show()
