import sys
import os
import numpy as np
sys.path.append(os.environ["DMRGROOT"]) # root for dmrg
import spinchain

n = 20
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain

#sc.kpmmaxm = 10 # KPM max m
sc.kpmscale = 10.0
import time
import matplotlib.pyplot as plt

ns = [100,300,500,1000]

for n in ns:
  t1 = time.time()
  (x,y) = sc.get_spismj(n=n,mode="DMRG",i=0,j=0,smart=False)
  t2 = time.time()
  print("Time in",n,t2-t1)
  plt.plot(x,y,label=str(n))

plt.legend()
plt.xlim([-1,4])

plt.show()

