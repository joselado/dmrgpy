import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
from dmrgpy import spinchain

n = 4
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
#spins = [2,6]
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

def fj(i,j):
  if j==(i+1): return np.random.random((3,3))
  else: return 0.0
#sc.set_exchange(fj)
sc.set_fields(lambda x: np.random.random(3))

#sc.kpmmaxm = 20 # KPM maxm
import time


i = np.random.randint(n)
j = np.random.randint(n)

t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400))
t2 = time.time()
print("Time with DMRG",t2-t1)
(x3,y3) = sc.get_dynamical_correlator(n=300,mode="ED",i=i,j=j,delta=0.02)
t3 = time.time()

#print("Time in reduced energy window",t1-t0)
print("Time with ED",t3-t2)
#print("Time using exact KPM",t3-t2)

import matplotlib.pyplot as plt

plt.plot(x2,np.abs(y2),c="blue",label="DMRG")
plt.scatter(x3,np.abs(y3),c="green",label="ED")
plt.xlim([-0.5,4.5])
plt.legend()

plt.show()

#sc.clean()
