# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../../src')

import numpy as np
from dmrgpy import spinchain
n = 6
# create a random spin chain
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Chain(spins) # create the spin chain
def fj(i,j):
  if abs(i-j)==1: return 1.0
  else: return 0.0
sc.set_exchange(fj)
#sc.get_gs()

#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
j = i

sc.get_gs() # compute ground state

t1 = time.time()
sc.kpm_accelerate = True # use acceleration in KPM mode
name = (sc.Sx[i],sc.Sx[j])
(x0,y0) = sc.get_dynamical_correlator(mode="DMRG",name=name)
t2 = time.time()
print("Time with acceleration",t2-t1)



sc.kpm_accelerate = False # use acceleration in KPM mode
(x1,y1) = sc.get_dynamical_correlator(mode="DMRG",name=name)
t3 = time.time()
print("Time without acceleration",t3-t2)





# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x0,y0.real,c="blue",label="With acceleration")
plt.scatter(x1,y1.real,c="red",label="Without acceleration")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()












