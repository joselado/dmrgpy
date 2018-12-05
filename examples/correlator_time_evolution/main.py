import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
from dmrgpy import spinchain




ns = np.array(range(4,30,2))
es = []
n = 2
spins = [2 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
def fj(i,j):
    if 0.9<abs(i-j)<1.1: return 1.0
    return 0.0
sc.set_exchange(fj) # set exchange couplings
sc.maxm = 10
sc.get_gs()
#(x,y) = sc.evolution(nt=300,dt=0.1)
es = np.linspace(-1.0,10.0,4000)
import time
print("Starting")
t0 = time.time()
(x1,y1) = sc.get_dynamical_correlator(es=es,use_kpm=True)
t1 = time.time()
sc.fit_td = True
(x2,y2) = sc.get_dynamical_correlator(es=es,use_kpm=False)
t2 = time.time()
print("Time in TD",t2-t1)
print("Time in KPM",t1-t0)

import matplotlib.pyplot as plt
plt.plot(x1,y1.real,label="KPM")
plt.plot(x2,np.abs(y2),label="TD")
plt.legend()
#plt.scatter(x,y.imag)
plt.show()


