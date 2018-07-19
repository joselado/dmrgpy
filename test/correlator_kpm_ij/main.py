import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 4
spins = [3 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

sc.set_fields(lambda x: [0.4,0.5,0.5])

def fc(i,j):
  if i==j+1: return 1.0
  else: return 0.0

sc.set_exchange(fc)

#sc.kpmmaxm = 10 # KPM max m
import time


t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",i=0,j=0,name="XZ")
t2 = time.time()
(x3,y3) = sc.get_dynamical_correlator(mode="fullKPM",i=0,j=0,name="XZ")
t3 = time.time()

print("Time with DMRG",t2-t1)
#print("Time in reduced energy window",t1-t0)
print("Time with ED",t3-t2)
#print("Time using exact KPM",t3-t2)

import matplotlib.pyplot as plt

plt.plot(x2,y2.imag,c="blue",label="DMRG")
plt.scatter(x3,y3.imag,c="green",label="ED")
plt.xlim([-0.5,4.5])
plt.legend()

plt.show()

