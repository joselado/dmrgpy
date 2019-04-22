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
#sc.get_gs()

#sc.kpmmaxm = 20 # KPM maxm
import time
i = np.random.randint(n)
j = np.random.randint(n)
j = i
t1 = time.time()
(x2,y2) = sc.get_dynamical_correlator(mode="DMRG",i=i,j=j,name="ZZ")
t2 = time.time()
print("Time with DMRG",t2-t1)


(x3,y3) = sc.get_dynamical_correlator(mode="ED",submode="ED",
        i=i,j=j,name="ZZ")
t3 = time.time()
print("Time with ED",t3-t2)


(x4,y4) = sc.get_dynamical_correlator(mode="ED",submode="KPM",
        i=i,j=j,name="ZZ")
t4 = time.time()
print("Time with KPM-ED",t4-t3)


# plot the results
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.family'] = "Bitstream Vera Serif"
fig = plt.figure()
fig.subplots_adjust(0.2,0.2)
plt.plot(x2,np.abs(y2),c="blue",label="DMRG")
plt.scatter(x3,np.abs(y3),c="green",label="ED")
plt.scatter(x4,np.abs(y4),c="red",label="ED KPM")
plt.legend()
plt.xlabel("frequency [J]")
plt.ylabel("Dynamical correlator")
plt.xlim([-0.5,4.5])
plt.show()



