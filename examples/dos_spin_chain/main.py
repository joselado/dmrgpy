# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 7
# create a random spin chain
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = [3 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
def fj(i,j):
  if abs(i-j)==1: return 1.0
  else: return 0.0
sc.set_exchange(fj)
sc.get_gs()

#sc.kpmmaxm = 20 # KPM maxm
import time
t0 = time.time()
sc.kpmmaxm = 50
sc.kpm_scale = 2.0
(x1,y1) = sc.get_dos(mode="DMRG",delta=0.2)
t1 = time.time()
print("Time in DMRG",t1-t0)
(x,y) = sc.get_dos(mode="ED",delta=1.0)
t2 = time.time()
print("Time in ED",t2-t1)


# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x,y,c="blue",label="ED")
plt.plot(x1,y1,c="red",label="DMRG")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("DOS")
plt.show()



