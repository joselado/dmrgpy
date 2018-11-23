import os
import sys
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain




ns = np.array(range(4,30,2))
es = []
n = 5
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
(x1,y1) = sc.get_dynamical_correlator(es=es,use_kpm=True)
(x2,y2) = sc.get_dynamical_correlator(es=es,use_kpm=False)


import matplotlib.pyplot as plt
plt.plot(x1,y1.real)
plt.plot(x2,np.abs(y2))
#plt.scatter(x,y.imag)
plt.show()


