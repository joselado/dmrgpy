# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
import spinchain
n = 4
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
#spins = [2,6]
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def fj(i,j):
  if 0.9<abs(i-j)<1.1: 
      return 1.0
      return np.random.random((3,3))
  else: return 0.0
sc.set_exchange(fj)
sc.set_fields(lambda x: np.random.random(3))
#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
i,j = 0,0
t1 = time.time()
(es,xx) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400),name="XX")
(es,xy) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400),name="XY")
(es,yx) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400),name="YX")
(es,yy) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400),name="YY")
(es,pm) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=i,j=j,delta=0.02,
        es=np.linspace(-0.5,15.0,400),name="+-")
t2 = time.time()
print("Time with DMRG",t2-t1)
pm2 = xx+yy+1j*(xy-yx)
import matplotlib.pyplot as plt
plt.plot(es,pm.real,c="blue",label="Direct")
plt.scatter(es,pm2.real,c="red",label="Summing")
plt.xlim([-0.5,4.5])
plt.legend()
plt.show()
#sc.clean()


