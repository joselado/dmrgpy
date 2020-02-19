# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

raise # this does not work yet

import numpy as np
from dmrgpy import spinchain
ns = np.array(range(4,30,2))
es = []
n = 6
spins = [2 for i in range(n)]
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
sc.clean()
sc = spinchain.Spin_Hamiltonian(spins) # create the chain
def fj(i,j):
    if 0.9<abs(i-j)<1.1: return 1.0
    return 0.0
sc.set_exchange(fj) # set exchange couplings
#sc.set_fields(lambda x: [0.2,0.2,0.2]) # set exchange couplings
sc.maxm = 10
sc.get_gs()
i = 0 ; j=5
(x,y) = sc.evolution(nt=1000,dt=0.03,mode="ED",name="ZZ",i=i,j=j)
(x1,y1) = sc.evolution(nt=500,dt=0.03,mode="DMRG",name="ZZ",i=i,j=j)
#(x2,y2) = sc.evolution(nt=500,dt=0.03,mode="DMRG",name="ZZ",
#        i=i,j=j,restart=False)
#y1 = np.concatenate([y1,y2])
#x1 = np.concatenate([x1,x2+np.max(x1)])
import matplotlib.pyplot as plt
plt.subplot(1,2,1)
plt.plot(x,y.real,c="red",label="ED")
plt.scatter(x1,y1.real,c="blue",label="DMRG")
plt.subplot(1,2,2)
plt.plot(x,y.imag,c="red",label="ED")
plt.scatter(x1,y1.imag,c="blue",label="DMRG")
plt.legend()
#plt.plot(x,y.imag)
#plt.scatter(x,y.imag)
plt.show()


