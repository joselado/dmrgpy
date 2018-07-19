import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 6
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

sc.set_fields(lambda x: [0.,0.,1.])

#sc.kpmmaxm = 10 # KPM max m
sc.kpmscale = 10.0
import time


t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(n=600,mode="DMRG",i=2,j=2)
t2 = time.time()
(x3,y3) = sc.get_dynamical_correlator(n=300,mode="ED",i=2,j=2)
t3 = time.time()

print("Time with DMRG",t2-t1)
#print("Time in reduced energy window",t1-t0)
print("Time with ED",t3-t2)
#print("Time using exact KPM",t3-t2)

import matplotlib.pyplot as plt

plt.plot(x2,np.abs(y2),c="blue",label="DMRG")
plt.scatter(x3,np.abs(y3),c="green",label="ED")
plt.xlim([-0.5,4.5])
plt.legend()

plt.show()

sc.clean()
