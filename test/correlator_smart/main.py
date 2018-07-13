import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 10
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

#sc.kpmmaxm = 10 # KPM max m
sc.kpmscale = 20.0
import time

sc.kpmmaxm = 10
t0 = time.time()
(x2,y2) = sc.get_spismj(n=1000,mode="DMRG",i=0,j=0,smart=False)
t1 = time.time()
(x1,y1) = sc.get_spismj(n=1000,mode="DMRG",i=0,j=0,smart=True)
t2 = time.time()

print("Time in full energy window",t2-t1)
print("Time in reduced energy window",t1-t0)

import matplotlib.pyplot as plt

plt.plot(x1,y1,c="red",label="Smart DMRG")
plt.scatter(x2,y2,c="blue",label="KPM DMRG")
#plt.xlim([-0.5,4.5])
plt.legend()

plt.show()

