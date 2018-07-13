import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 6
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

#sc.kpmmaxm = 10 # KPM max m
sc.kpmscale = 10.0
import time


t0 = time.time()
#(x1,y1) = sc.get_spismj(n=300,mode="DMRG",i=0,j=0,smart=True)
t1 = time.time()
(x2,y2) = sc.get_spismj(n=300,mode="DMRG",i=0,j=0,smart=False)
t2 = time.time()
(x3,y3) = sc.get_spismj(n=300,mode="full",i=0,j=0)
t3 = time.time()
#(x4,y4) = sc.get_spismj(n=300,mode="fullKPM",i=0,j=0)
t4 = time.time()

print("Time in full energy window",t2-t1)
#print("Time in reduced energy window",t1-t0)
print("Time using exact matrix inversion",t3-t2)
#print("Time using exact KPM",t3-t2)

import matplotlib.pyplot as plt

#plt.plot(x1,y1,c="red",label="Smart DMRG")
plt.plot(x2,y2,c="blue",label="KPM DMRG")
plt.scatter(x3,y3,c="green",label="Exact inversion")
#plt.scatter(x4,y4,c="black",label="KPM full")
plt.xlim([-0.5,4.5])
plt.legend()

plt.show()

