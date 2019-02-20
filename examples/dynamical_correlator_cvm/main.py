# Add the root path of the dmrgpy library
import os ; import sys ; sys.path.append(os.getcwd()+'/../../src')

import numpy as np
from dmrgpy import spinchain
n = 5
# create a random spin chain
spins = [np.random.randint(2,5) for i in range(n)] # spin 1/2 heisenberg chain
spins = [2 for i in range(n)] # spin 1/2 heisenberg chain


# create first neighbor exchange
sc = spinchain.Spin_Hamiltonian(spins) # create the spin chain
def fj(i,j):
  if abs(i-j)==1: return 1.0
  else: return 0.0
sc.set_exchange(fj)

sc.get_gs()

import time
i = np.random.randint(n)
j = np.random.randint(n)
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(submode="KPM",i=i,j=j,name="XX")
t2 = time.time()
print("Time with KPM",t2-t1)


(x3,y3) = sc.get_dynamical_correlator(submode="CVM",i=i,j=j,name="XX")
t3 = time.time()
print("Time with CVM",t3-t2)




# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,np.abs(y2),c="blue",label="DMRG")
plt.scatter(x3,np.abs(y3),c="green",label="ED")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



